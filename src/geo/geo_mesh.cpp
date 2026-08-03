#include "geo_mesh.h"
#include "geo_permutation.h"

namespace floatTetWild {
namespace geo {

    void MeshVertices::permute_elements(vector<index_t>& permutation) {
        // geogram keeps the coordinates in a vertex attribute, so its permute_elements() moves
        // them along with everything else bound to a vertex. Here they are a plain array, so
        // they are moved explicitly first.
        {
            std::vector<double> reordered(points_.size());
            for(index_t v = 0; v < permutation.size(); ++v) {
                const double* from = point_ptr(permutation[v]);
                std::copy(from, from + DIMENSION, reordered.begin() + v * DIMENSION);
            }
            points_.swap(reordered);
        }

        Permutation::invert(permutation);

        MeshFacetCorners& facet_corners = mesh_.facet_corners;
        for(index_t c = 0; c < facet_corners.nb(); ++c) {
            facet_corners.set_vertex(c, permutation[facet_corners.vertex(c)]);
        }
    }

    index_t MeshFacets::create_triangles(index_t nb_triangles) {
        index_t first_facet = nb();
        facet_corners_.create_sub_elements(nb_triangles * 3);
        nb_ += nb_triangles;
        return first_facet;
    }

    void MeshFacets::clear() {
        nb_ = 0;
        facet_corners_.clear();
    }

    void MeshFacets::permute_elements(vector<index_t>& permutation) {
        // Every facet is a triangle, so the corners move in place three at a time.
        Permutation::apply(
            facet_corners_.corner_vertex_.data(), permutation, index_t(sizeof(index_t) * 3)
        );
    }

    double bbox_diagonal(const Mesh& M) {
        double xyzmin[3];
        double xyzmax[3];
        for(index_t c = 0; c < 3; ++c) {
            xyzmin[c] = Numeric::max_float64();
            xyzmax[c] = Numeric::min_float64();
        }
        for(index_t v = 0; v < M.vertices.nb(); ++v) {
            const double* p = M.vertices.point_ptr(v);
            for(index_t c = 0; c < 3; ++c) {
                xyzmin[c] = std::min(xyzmin[c], p[c]);
                xyzmax[c] = std::max(xyzmax[c], p[c]);
            }
        }
        return ::sqrt(
            geo_sqr(xyzmax[0] - xyzmin[0]) +
            geo_sqr(xyzmax[1] - xyzmin[1]) +
            geo_sqr(xyzmax[2] - xyzmin[2])
        );
    }
}
}
