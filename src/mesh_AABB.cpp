// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see geo/LICENSE.geogram.
// Source: geogram/mesh/mesh_AABB.cpp
// Copied rather than reimplemented: the traversal order decides which of two equidistant facets a
// query lands on, and the mesher's envelope tests carry that facet forward as a hint.

#include <floattetwild/mesh_AABB.h>
#include <floattetwild/geo_mesh.h>
#include <floattetwild/geo_geometry.h>

namespace {

    using namespace floatTetWild::geo;

    // The smallest Box that encloses \p B1 and \p B2.
    inline void bbox_union(Box& target, const Box& B1, const Box& B2) {
        for(index_t c = 0; c < 3; ++c) {
            target.xyz_min[c] = std::min(B1.xyz_min[c], B2.xyz_min[c]);
            target.xyz_max[c] = std::max(B1.xyz_max[c], B2.xyz_max[c]);
        }
    }

    void get_facet_bbox(
        const Mesh& M, Box& B, index_t f
    ) {
        const double* p = M.point_ptr(M.facet_vertex(f, 0));
        for(index_t coord = 0; coord < 3; ++coord) {
            B.xyz_min[coord] = p[coord];
            B.xyz_max[coord] = p[coord];
        }
        for(index_t lv = 1; lv < 3; ++lv) {
            p = M.point_ptr(M.facet_vertex(f, lv));
            for(index_t coord = 0; coord < 3; ++coord) {
                B.xyz_min[coord] = std::min(B.xyz_min[coord], p[coord]);
                B.xyz_max[coord] = std::max(B.xyz_max[coord], p[coord]);
            }
        }
    }

    // The tree has no combinatorics of its own: the children of node n are 2n and 2n+1, and the
    // facets [b,e) under it are split down the middle. So the array of boxes has to be as long as
    // the largest node index the recursion below reaches.
    index_t max_node_index(index_t node_index, index_t b, index_t e) {
        geo_debug_assert(e > b);
        if(b + 1 == e) {
            return node_index;
        }
        index_t m = b + (e - b) / 2;
        index_t childl = 2 * node_index;
        index_t childr = 2 * node_index + 1;
        return std::max(
            max_node_index(childl, b, m),
            max_node_index(childr, m, e)
        );
    }

    void init_bboxes_recursive(
        const Mesh& M, vector<Box>& bboxes,
        index_t node_index,
        index_t b, index_t e
    ) {
        geo_debug_assert(node_index < bboxes.size());
        geo_debug_assert(b != e);
        if(b + 1 == e) {
            get_facet_bbox(M, bboxes[node_index], b);
            return;
        }
        index_t m = b + (e - b) / 2;
        index_t childl = 2 * node_index;
        index_t childr = 2 * node_index + 1;
        geo_debug_assert(childl < bboxes.size());
        geo_debug_assert(childr < bboxes.size());
        init_bboxes_recursive(M, bboxes, childl, b, m);
        init_bboxes_recursive(M, bboxes, childr, m, e);
        bbox_union(bboxes[node_index], bboxes[childl], bboxes[childr]);
    }

    // \pre p is inside B
    double inner_point_box_squared_distance(
        const vec3& p,
        const Box& B
    ) {
        double result = geo_sqr(p[0] - B.xyz_min[0]);
        result = std::min(result, geo_sqr(p[0] - B.xyz_max[0]));
        for(index_t c = 1; c < 3; ++c) {
            result = std::min(result, geo_sqr(p[c] - B.xyz_min[c]));
            result = std::min(result, geo_sqr(p[c] - B.xyz_max[c]));
        }
        return result;
    }

    // Negative when the point is inside the box.
    double point_box_signed_squared_distance(
        const vec3& p,
        const Box& B
    ) {
        bool inside = true;
        double result = 0.0;
        for(index_t c = 0; c < 3; c++) {
            if(p[c] < B.xyz_min[c]) {
                inside = false;
                result += geo_sqr(p[c] - B.xyz_min[c]);
            } else if(p[c] > B.xyz_max[c]) {
                inside = false;
                result += geo_sqr(p[c] - B.xyz_max[c]);
            }
        }
        if(inside) {
            result = -inner_point_box_squared_distance(p, B);
        }
        return result;
    }

    void get_point_facet_nearest_point(
        const Mesh& M,
        const vec3& p,
        index_t f,
        vec3& nearest_p,
        double& squared_dist
    ) {
        double lambda1, lambda2, lambda3;  // barycentric coords, not used.
        squared_dist = Geom::point_triangle_squared_distance(
            p, M.facet_point(f, 0), M.facet_point(f, 1), M.facet_point(f, 2),
            nearest_p, lambda1, lambda2, lambda3
        );
    }

    double point_box_center_squared_distance(
        const vec3& p, const Box& B
    ) {
        double result = 0.0;
        for(index_t c = 0; c < 3; ++c) {
            double d = p[c] - 0.5 * (B.xyz_min[c] + B.xyz_max[c]);
            result += geo_sqr(d);
        }
        return result;
    }

    // The recursive traversal behind both nearest facet queries. With \p Envelope set the walk
    // stops as soon as the facet found so far is within \p sq_epsilon, and never descends into a
    // box that is already outside it. Without it, \p sq_epsilon is unused and the walk finds the
    // nearest facet outright.
    template <bool Envelope>
    void nearest_recursive(
        const Mesh& M, const vector<Box>& bboxes,
        const vec3& p, double sq_epsilon,
        index_t& nearest_f, vec3& nearest_point, double& sq_dist,
        index_t n, index_t b, index_t e
    ) {
        geo_debug_assert(e > b);

        if(Envelope && sq_dist <= sq_epsilon) {
            return;
        }

        // If node is a leaf: compute point-facet distance
        // and replace current if nearer
        if(b + 1 == e) {
            vec3 cur_nearest_point;
            double cur_sq_dist;
            get_point_facet_nearest_point(M, p, b, cur_nearest_point, cur_sq_dist);
            if(cur_sq_dist < sq_dist) {
                nearest_f = b;
                nearest_point = cur_nearest_point;
                sq_dist = cur_sq_dist;
            }
            return;
        }
        index_t m = b + (e - b) / 2;
        index_t childl = 2 * n;
        index_t childr = 2 * n + 1;

        double dl = point_box_signed_squared_distance(p, bboxes[childl]);
        double dr = point_box_signed_squared_distance(p, bboxes[childr]);

        const auto descend = [&](double d, index_t child, index_t cb, index_t ce) {
            if(d < sq_dist && (!Envelope || d <= sq_epsilon)) {
                nearest_recursive<Envelope>(
                    M, bboxes, p, sq_epsilon,
                    nearest_f, nearest_point, sq_dist,
                    child, cb, ce
                );
            }
        };

        // Traverse the "nearest" child first, so that it has more chances
        // to prune the traversal of the other child.
        if(dl < dr) {
            descend(dl, childl, b, m);
            descend(dr, childr, m, e);
        } else {
            descend(dr, childr, m, e);
            descend(dl, childl, b, m);
        }
    }

}

/****************************************************************************/

namespace floatTetWild {

    MeshFacetsAABBWithEps::MeshFacetsAABBWithEps(const Mesh& M)
     : mesh_(M) {
        bboxes_.resize(
            max_node_index(
                1, 0, mesh_.nb_facets()
            ) + 1 // <-- this is because size == max_index + 1 !!!
        );
        init_bboxes_recursive(mesh_, bboxes_, 1, 0, mesh_.nb_facets());
    }

    void MeshFacetsAABBWithEps::get_nearest_facet_hint(
        const vec3& p,
        index_t& nearest_f, vec3& nearest_point, double& sq_dist
    ) const {
        index_t b = 0;
        index_t e = mesh_.nb_facets();
        index_t n = 1;
        while(e != b + 1) {
            index_t m = b + (e - b) / 2;
            index_t childl = 2 * n;
            index_t childr = 2 * n + 1;
            if(
                point_box_center_squared_distance(p, bboxes_[childl]) <
                point_box_center_squared_distance(p, bboxes_[childr])
            ) {
                e = m;
                n = childl;
            } else {
                b = m;
                n = childr;
            }
        }
        nearest_f = b;

        nearest_point = mesh_.facet_point(nearest_f, 0);
        sq_dist = Geom::distance2(p, nearest_point);
    }

    index_t MeshFacetsAABBWithEps::nearest_facet(
        const vec3& p, vec3& nearest_point, double& sq_dist
    ) const {
        index_t nearest_f;
        get_nearest_facet_hint(p, nearest_f, nearest_point, sq_dist);
        nearest_recursive<false>(
            mesh_, bboxes_, p, 0,
            nearest_f, nearest_point, sq_dist, 1, 0, mesh_.nb_facets()
        );
        return nearest_f;
    }

    void MeshFacetsAABBWithEps::facet_in_envelope_with_hint(
        const vec3& p, double sq_epsilon,
        index_t& nearest_f, vec3& nearest_point, double& sq_dist
    ) const {
        //  Measuring against the facet the caller already has is what makes the hint worth
        // carrying: for the samples along one triangle it usually still answers the query, and
        // then there is no descent at all. Returning here rather than letting the walk below
        // return at once is what keeps that case free of calls.
        if(nearest_f == NO_FACET) {
            get_nearest_facet_hint(p, nearest_f, nearest_point, sq_dist);
        } else {
            get_point_facet_nearest_point(mesh_, p, nearest_f, nearest_point, sq_dist);
            if(sq_dist <= sq_epsilon) {
                return;
            }
        }
        nearest_recursive<true>(
            mesh_, bboxes_, p, sq_epsilon,
            nearest_f, nearest_point, sq_dist, 1, 0, mesh_.nb_facets()
        );
    }

/****************************************************************************/

}
