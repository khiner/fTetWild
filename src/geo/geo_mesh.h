// A stand-in for the slice of geogram/mesh/mesh.h that fTetWild and the vendored spatial sort,
// loaders and AABB use.
//
// Reimplemented rather than copied, unlike its neighbours here. A container addresses elements by
// index and has no comparisons of its own, so there is no tie-breaking to preserve; what it does
// have to match is the element order geogram's operations leave behind, which is why the permute
// functions in geo_mesh_reorder.cpp reproduce what geogram's permute_elements() leaves.
//
// geogram spreads this over MeshVertices, MeshFacetCorners and MeshFacets, and supports polygons,
// other dimensions and cells. Every facet built here is a triangle -- the loaders fan polygons as
// they read them -- every mesh is 3d, and fTetWild's only use of cells was Mesh::partition(),
// which has no live caller. So the whole thing is two flat arrays.

#pragma once

#include "geo_geometry.h"

namespace floatTetWild {
namespace geo {

    static const index_t NO_FACET = index_t(-1);

    // A triangle surface mesh: three coordinates per vertex, then three vertex indices per facet
    // in order around it.
    struct Mesh {
        std::vector<double> points;
        std::vector<index_t> corners;

        index_t nb_vertices() const {
            return index_t(points.size() / 3);
        }

        index_t nb_facets() const {
            return index_t(corners.size() / 3);
        }

        double* point_ptr(index_t v) {
            geo_debug_assert(v < nb_vertices());
            return points.data() + 3 * v;
        }

        const double* point_ptr(index_t v) const {
            geo_debug_assert(v < nb_vertices());
            return points.data() + 3 * v;
        }

        vec3& point(index_t v) {
            return *reinterpret_cast<vec3*>(point_ptr(v));
        }

        const vec3& point(index_t v) const {
            return *reinterpret_cast<const vec3*>(point_ptr(v));
        }

        index_t facet_vertex(index_t f, index_t lv) const {
            geo_debug_assert(f < nb_facets() && lv < 3);
            return corners[3 * f + lv];
        }

        void set_facet_vertex(index_t f, index_t lv, index_t v) {
            geo_debug_assert(f < nb_facets() && lv < 3);
            corners[3 * f + lv] = v;
        }

        const vec3& facet_point(index_t f, index_t lv) const {
            return point(facet_vertex(f, lv));
        }

        // Creates a contiguous chunk of vertices at the origin, and returns the first index.
        index_t create_vertices(index_t nb) {
            index_t first = nb_vertices();
            points.resize(size_t(first + nb) * 3, 0.0);
            return first;
        }

        // The same for triangles, with no vertices set yet.
        index_t create_triangles(index_t nb) {
            index_t first = nb_facets();
            corners.resize(size_t(first + nb) * 3, NO_INDEX);
            return first;
        }

        void clear() {
            points.clear();
            points.shrink_to_fit();
            corners.clear();
            corners.shrink_to_fit();
        }
    };
}
}
