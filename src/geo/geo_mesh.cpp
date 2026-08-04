#include "geo_mesh.h"

namespace floatTetWild {
namespace geo {

    void Mesh::permute_vertices(const vector<index_t>& permutation) {
        std::vector<double> reordered(points.size());
        for(index_t v = 0; v < permutation.size(); ++v) {
            const double* from = point_ptr(permutation[v]);
            std::copy(from, from + 3, reordered.begin() + v * 3);
        }
        points.swap(reordered);

        // The vertex that was at permutation[v] is now at v, so a corner naming the old index
        // has to be sent the other way.
        vector<index_t> old_to_new(permutation.size());
        for(index_t v = 0; v < permutation.size(); ++v) {
            old_to_new[permutation[v]] = v;
        }
        for(index_t& corner : corners) {
            corner = old_to_new[corner];
        }
    }

    void Mesh::permute_facets(const vector<index_t>& permutation) {
        // Every facet is a triangle, so the corners move three at a time.
        std::vector<index_t> reordered(corners.size());
        for(index_t f = 0; f < permutation.size(); ++f) {
            for(index_t lv = 0; lv < 3; ++lv) {
                reordered[3 * f + lv] = corners[3 * permutation[f] + lv];
            }
        }
        corners.swap(reordered);
    }

    double bbox_diagonal(const Mesh& M) {
        double xyzmin[3];
        double xyzmax[3];
        for(index_t c = 0; c < 3; ++c) {
            xyzmin[c] = std::numeric_limits<double>::max();
            xyzmax[c] = std::numeric_limits<double>::lowest();
        }
        for(index_t v = 0; v < M.nb_vertices(); ++v) {
            const double* p = M.point_ptr(v);
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
