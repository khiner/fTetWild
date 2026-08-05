// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/mesh/mesh_reorder.cpp
// Copied rather than reimplemented: spatial sort; nth_element is not stable so the incoming
// arrangement decides ties.
//
// geogram parameterizes the sort over a comparator template, which is itself parameterized over
// the coordinate, the direction and a mesh class. Two things vary between the three sorts here:
// what coordinate an element compares by, and whether the order reverses along an axis, which is
// what separates a Hilbert curve from a Morton one. So those are the two parameters.

#include "geo_mesh_reorder.h"
#include "geo_parallel.h"

#include <algorithm>
#include <random>

namespace {

    using namespace floatTetWild::geo;

    using Iterator = vector<index_t>::iterator;

    // Partitions the sequence into two halves with the same number of elements, such that the
    // elements of the first are smaller than those of the second, and returns the middle.
    template <class CMP>
    inline Iterator reorder_split(Iterator begin, Iterator end, CMP cmp) {
        if(begin >= end) {
            return begin;
        }
        Iterator middle = begin + (end - begin) / 2;
        std::nth_element(begin, middle, end, cmp);
        return middle;
    }

    /************************************************************************/

    // The coordinate a vertex sorts by: its own.
    struct VertexCoord {
        const double* points;

        double operator() (index_t v, int coord) const {
            return points[3 * v + coord];
        }
    };

    // The coordinate a facet sorts by: the centre of its three vertices.
    struct FacetCoord {
        const Mesh& mesh;

        double operator() (index_t f, int coord) const {
            double result = 0.0;
            double s = 1.0 / 3.0;
            for(index_t lv = 0; lv < 3; ++lv) {
                result += s * mesh.point_ptr(mesh.facet_vertex(f, lv))[coord];
            }
            return result;
        }
    };

    // Compares two elements along coordinate COORD. UP picks direct or reverse order and is
    // ignored unless DIRECTED, which is true for the Hilbert order, that reverses along some
    // axes, and false for the Morton order, that does not.
    template <class COORD_OF, int COORD, bool UP, bool DIRECTED>
    struct Compare {
        const COORD_OF& coord_of;

        bool operator() (index_t i1, index_t i2) const {
            return (UP || !DIRECTED)
                ? coord_of(i1, COORD) < coord_of(i2, COORD)
                : coord_of(i1, COORD) > coord_of(i2, COORD);
        }
    };

    /************************************************************************/

    // Sorts elements in Hilbert or Morton order in 3d, after Delage and Devillers, "Spatial
    // Sorting", CGAL User and Reference Manual, 3.9 edition, 2011.
    template <class COORD_OF, bool DIRECTED>
    struct SpatialSort {

        template <int COORD, bool UP>
        static Compare<COORD_OF, COORD, UP, DIRECTED> cmp(const COORD_OF& coord_of) {
            return {coord_of};
        }

        // The recursion. COORDX is the first coordinate, 0, 1 or 2, and the second and third
        // are COORDX+1 and COORDX+2 modulo 3. UPX, UPY and UPZ are the direction along each.
        template <int COORDX, bool UPX, bool UPY, bool UPZ>
        static void sort(const COORD_OF& c, Iterator begin, Iterator end) {
            const int COORDY = (COORDX + 1) % 3, COORDZ = (COORDY + 1) % 3;
            if(end - begin <= 1) {
                return;
            }
            Iterator m0 = begin, m8 = end;
            Iterator m4 = reorder_split(m0, m8, cmp<COORDX,  UPX>(c));
            Iterator m2 = reorder_split(m0, m4, cmp<COORDY,  UPY>(c));
            Iterator m1 = reorder_split(m0, m2, cmp<COORDZ,  UPZ>(c));
            Iterator m3 = reorder_split(m2, m4, cmp<COORDZ, !UPZ>(c));
            Iterator m6 = reorder_split(m4, m8, cmp<COORDY, !UPY>(c));
            Iterator m5 = reorder_split(m4, m6, cmp<COORDZ,  UPZ>(c));
            Iterator m7 = reorder_split(m6, m8, cmp<COORDZ, !UPZ>(c));
            sort<COORDZ,  UPZ,  UPX,  UPY>(c, m0, m1);
            sort<COORDY,  UPY,  UPZ,  UPX>(c, m1, m2);
            sort<COORDY,  UPY,  UPZ,  UPX>(c, m2, m3);
            sort<COORDX,  UPX, !UPY, !UPZ>(c, m3, m4);
            sort<COORDX,  UPX, !UPY, !UPZ>(c, m4, m5);
            sort<COORDY, !UPY,  UPZ, !UPX>(c, m5, m6);
            sort<COORDY, !UPY,  UPZ, !UPX>(c, m6, m7);
            sort<COORDZ, !UPZ, !UPX,  UPY>(c, m7, m8);
        }

        // Sorts a sequence of element indices spatially, in place. The top of the recursion,
        // spread over threads. It is the body of sort() with COORDX 0 and every direction
        // direct, written out because a lambda cannot take a template parameter of the enclosing
        // function as a template argument on every compiler.
        static void run(const COORD_OF& c, Iterator b, Iterator e) {
            geo_debug_assert(e >= b);

            // If the sequence is small enough, sort it sequentially.
            if(e - b < 1024) {
                sort<0, false, false, false>(c, b, e);
                return;
            }

            // Parallel sorting (2 then 4 then 8 sorts in parallel)
            Iterator m0 = b, m8 = e;
            Iterator m1, m2, m3, m5, m6, m7;
            Iterator m4 = reorder_split(m0, m8, cmp<0, false>(c));

            parallel(
                [&]() { m2 = reorder_split(m0, m4, cmp<1, false>(c)); },
                [&]() { m6 = reorder_split(m4, m8, cmp<1, true >(c)); }
            );

            parallel(
                [&]() { m1 = reorder_split(m0, m2, cmp<2, false>(c)); },
                [&]() { m3 = reorder_split(m2, m4, cmp<2, true >(c)); },
                [&]() { m5 = reorder_split(m4, m6, cmp<2, false>(c)); },
                [&]() { m7 = reorder_split(m6, m8, cmp<2, true >(c)); }
            );

            parallel(
                [&]() { sort<2, false, false, false>(c, m0, m1); },
                [&]() { sort<1, false, false, false>(c, m1, m2); },
                [&]() { sort<1, false, false, false>(c, m2, m3); },
                [&]() { sort<0, false, true,  true >(c, m3, m4); },
                [&]() { sort<0, false, true,  true >(c, m4, m5); },
                [&]() { sort<1, true,  false, true >(c, m5, m6); },
                [&]() { sort<1, true,  false, true >(c, m6, m7); },
                [&]() { sort<2, true,  true,  false>(c, m7, m8); }
            );
        }
    };

    /************************************************************************/

    // The identity permutation over \p n elements, for a sort to rearrange.
    void trivial_indices(index_t n, vector<index_t>& sorted_indices) {
        sorted_indices.resize(n);
        for(index_t i = 0; i < n; ++i) {
            sorted_indices[i] = i;
        }
    }

    // Implementation of compute_BRIO_order(), over the [b,e) range of indices into \p vertices,
    // which holds three coordinates per vertex.
    void compute_BRIO_order_recursive(
        const double* vertices, Iterator b, Iterator e
    ) {
        geo_debug_assert(e > b);

        // Minimum size of interval to be sorted, and the splitting ratio between the current
        // interval and the rest.
        static const index_t threshold = 64;
        static const double ratio = 0.125;

        Iterator m = b;
        if(index_t(e - b) > threshold) {
            m = b + int32_t(double(e - b) * ratio);
            compute_BRIO_order_recursive(vertices, b, m);
        }

        SpatialSort<VertexCoord, true>::run(VertexCoord{vertices}, m, e);
    }
}

/****************************************************************************/

namespace floatTetWild {
namespace geo {


    void mesh_reorder(Mesh& M, vector<index_t>* facet_permutation) {

        // Step 1: reorder vertices
        {
            vector<index_t> sorted_indices;
            trivial_indices(M.nb_vertices(), sorted_indices);
            SpatialSort<VertexCoord, false>::run(
                VertexCoord{M.points.data()},
                sorted_indices.begin(), sorted_indices.end()
            );
            M.permute_vertices(sorted_indices);
        }

        // Step 2: reorder facets
        if(M.nb_facets() != 0) {
            vector<index_t> sorted_indices;
            trivial_indices(M.nb_facets(), sorted_indices);
            SpatialSort<FacetCoord, false>::run(
                FacetCoord{M},
                sorted_indices.begin(), sorted_indices.end()
            );
            if(facet_permutation != nullptr) {
                *facet_permutation = sorted_indices;
            }
            M.permute_facets(sorted_indices);
        }

    }


    void compute_BRIO_order(
        index_t nb_vertices, const double* vertices, vector<index_t>& sorted_indices
    ) {
        trivial_indices(nb_vertices, sorted_indices);

        // The fixed seed used to be applied to geogram's own copy by
        // cmake/patches/geogram_deterministic_shuffle.cmake. It is baked in here instead: geogram
        // seeds a fresh std::mt19937 from std::random_device on every call, so the insertion order
        // this picks would differ every run and every stage downstream of
        // FloatTetDelaunay::tetrahedralize would diverge with it. The order is still drawn from the
        // same distribution, just pinned to one draw.
        std::mt19937 urng(42);
        std::shuffle(sorted_indices.begin(), sorted_indices.end(), urng);
        compute_BRIO_order_recursive(
            vertices, sorted_indices.begin(), sorted_indices.end()
        );
    }
} }
