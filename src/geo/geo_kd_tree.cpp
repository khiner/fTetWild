// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/points/kd_tree.cpp
// Copied rather than reimplemented: balanced kd-tree; nth_element splitting again, so ties follow the incoming order

#include "geo_kd_tree.h"
#include "geo_geometry.h"
#include "geo_parallel.h"

namespace floatTetWild {
namespace geo {

    void BalancedKdTree::set_points(
        index_t nb_points, const double* points
    ) {
        nb_points_ = nb_points;
        points_ = points;

        point_index_.resize(nb_points);
        for(index_t i = 0; i < nb_points; i++) {
            point_index_[i] = i;
        }

        // Compute the bounding box.
        for(coord_index_t c = 0; c < DIMENSION; ++c) {
            bbox_min_[c] =  std::numeric_limits<double>::max();
            bbox_max_[c] = -std::numeric_limits<double>::max();
        }
        for(index_t i = 0; i < nb_points; ++i) {
            const double* p = point_ptr(i);
            for(coord_index_t c = 0; c < DIMENSION; ++c) {
                bbox_min_[c] = std::min(bbox_min_[c], p[c]);
                bbox_max_[c] = std::max(bbox_max_[c], p[c]);
            }
        }

        build_tree();
    }

    double BalancedKdTree::nearest_sq_dist(const double* query_point) const {
        // Distance between the query point and the global bounding box, with that box copied to
        // stack locals. The traversal moves bbox_min and bbox_max as it descends, which is what
        // lets it measure the distance from the query point to the current node's box.
        double box_dist = 0.0;
        double bbox_min[DIMENSION];
        double bbox_max[DIMENSION];
        for(coord_index_t c = 0; c < DIMENSION; ++c) {
            bbox_min[c] = bbox_min_[c];
            bbox_max[c] = bbox_max_[c];
            if(query_point[c] < bbox_min_[c]) {
                box_dist += geo_sqr(bbox_min_[c] - query_point[c]);
            } else if(query_point[c] > bbox_max_[c]) {
                box_dist += geo_sqr(bbox_max_[c] - query_point[c]);
            }
        }

        // Root node is number 1. This is because "children at 2*n and 2*n+1" does not work with 0.
        double best_sq_dist = std::numeric_limits<double>::max();
        nearest_recursive(
            1, 0, nb_points(), bbox_min, bbox_max, box_dist, query_point, best_sq_dist
        );
        return best_sq_dist;
    }

    void BalancedKdTree::nearest_recursive(
        index_t node_index, index_t b, index_t e,
        double* bbox_min, double* bbox_max, double box_dist,
        const double* query_point, double& best_sq_dist
    ) const {
        geo_debug_assert(e > b);

        // Simple case (node is a leaf)
        if((e - b) <= MAX_LEAF_SIZE) {
            for(index_t i = b; i < e; ++i) {
                double sq_dist = Geom::distance2(
                    query_point, point_ptr(point_index_[i])
                );
                if(sq_dist < best_sq_dist) {
                    best_sq_dist = sq_dist;
                }
            }
            return;
        }

        const coord_index_t coord = splitting_coord_[node_index];
        const double val = splitting_val_[node_index];
        const index_t m = b + (e - b) / 2;
        const index_t left_node_index = 2 * node_index;
        const index_t right_node_index = 2 * node_index + 1;

        const double cut_diff = query_point[coord] - val;

        //  The subtree the query point falls in is visited first, with the cutting plane standing
        // in for the bound it replaces; \p bound is whichever of bbox_min and bbox_max that is.
        // The bound is restored on the way out, so the caller's box is unchanged.
        const auto descend = [&](index_t child, index_t cb, index_t ce, double* bound) {
            const double save = bound[coord];
            bound[coord] = val;
            nearest_recursive(
                child, cb, ce,
                bbox_min, bbox_max, box_dist, query_point, best_sq_dist
            );
            bound[coord] = save;
        };

        const bool near_is_left = cut_diff < 0.0;
        // The near subtree clamps the bound on the far side of the cut, and the far subtree the
        // near one.
        double* const near_bound = near_is_left ? bbox_max : bbox_min;
        double* const far_bound = near_is_left ? bbox_min : bbox_max;

        if(near_is_left) {
            descend(left_node_index, b, m, near_bound);
        } else {
            descend(right_node_index, m, e, near_bound);
        }

        // Update bbox distance (now measures the distance to the bbox of the far subtree)
        const double box_diff = near_is_left
            ? far_bound[coord] - query_point[coord]
            : query_point[coord] - far_bound[coord];
        if(box_diff > 0.0) {
            box_dist -= geo_sqr(box_diff);
        }
        box_dist += geo_sqr(cut_diff);

        // Traverse the far subtree, only if bbox distance is nearer than the nearest point so
        // far, else there is no chance that it contains anything nearer.
        if(box_dist <= best_sq_dist) {
            if(near_is_left) {
                descend(right_node_index, m, e, far_bound);
            } else {
                descend(left_node_index, b, m, far_bound);
            }
        }
    }

/****************************************************************************/

    void BalancedKdTree::build_tree() {
        index_t sz = max_node_index(1, 0, nb_points()) + 1;
        splitting_coord_.resize(sz);
        splitting_val_.resize(sz);

        // If there are more than 16*MAX_LEAF_SIZE (=256) points,
        // create the tree in parallel
        if(
            nb_points() < (16 * MAX_LEAF_SIZE) ||
            std::thread::hardware_concurrency() <= 1
        ) {
            create_kd_tree_recursive(1, 0, nb_points());
            return;
        }

        // Split points m0..m8 of the first four levels. parallel() joins before it returns, so
        // each level sees the one above it finished.
        index_t m0 = 0, m8 = nb_points();
        index_t m1, m2, m3, m4, m5, m6, m7;

        m4 = split_kd_node(1, m0, m8);

        parallel(
            [&]() { m2 = split_kd_node(2, m0, m4); },
            [&]() { m6 = split_kd_node(3, m4, m8); }
        );

        parallel(
            [&]() { m1 = split_kd_node(4, m0, m2); },
            [&]() { m3 = split_kd_node(5, m2, m4); },
            [&]() { m5 = split_kd_node(6, m4, m6); },
            [&]() { m7 = split_kd_node(7, m6, m8); }
        );

        parallel(
            [&]() { create_kd_tree_recursive(8 , m0, m1); },
            [&]() { create_kd_tree_recursive(9 , m1, m2); },
            [&]() { create_kd_tree_recursive(10, m2, m3); },
            [&]() { create_kd_tree_recursive(11, m3, m4); },
            [&]() { create_kd_tree_recursive(12, m4, m5); },
            [&]() { create_kd_tree_recursive(13, m5, m6); },
            [&]() { create_kd_tree_recursive(14, m6, m7); },
            [&]() { create_kd_tree_recursive(15, m7, m8); }
        );
    }

    index_t BalancedKdTree::split_kd_node(
        index_t node_index, index_t b, index_t e
    ) {

        geo_debug_assert(e > b);
        // Do not split leafs
        if(b + 1 == e) {
            return b;
        }

        coord_index_t splitting_coord = best_splitting_coord(b, e);
        index_t m = b + (e - b) / 2;
        geo_debug_assert(m < e);

        // sorts the indices in such a way that points's
        // coordinates splitting_coord in [b,m) are smaller
        // than m's and points in [m,e) are
        // greater or equal to m's
        std::nth_element(
            point_index_.begin() + std::ptrdiff_t(b),
            point_index_.begin() + std::ptrdiff_t(m),
            point_index_.begin() + std::ptrdiff_t(e),
            [this, splitting_coord](index_t i, index_t j) {
                return point_ptr(i)[splitting_coord] < point_ptr(j)[splitting_coord];
            }
        );

        // Initialize node's variables (splitting coord and
        // splitting value)
        splitting_coord_[node_index] = splitting_coord;
        splitting_val_[node_index] =
            point_ptr(point_index_[m])[splitting_coord];
        return m;
    }

    coord_index_t BalancedKdTree::best_splitting_coord(
        index_t b, index_t e
    ) {
        // Returns the coordinates that maximizes
        // point's spread. We should probably
        // use a tradeoff between spread and
        // bbox shape ratio, as done in ANN, but
        // this simple method seems to give good
        // results in our case.
        coord_index_t result = 0;
        double max_spread = spread(b, e, 0);
        for(coord_index_t c = 1; c < DIMENSION; ++c) {
            double coord_spread = spread(b, e, c);
            if(coord_spread > max_spread) {
                result = c;
                max_spread = coord_spread;
            }
        }
        return result;
    }

    /*************************************************************************/

} }
