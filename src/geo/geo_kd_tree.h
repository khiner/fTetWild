// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/points/kd_tree.h and geogram/points/nn_search.h
// Copied rather than reimplemented: balanced kd-tree; nth_element splitting again, so ties follow the incoming order
//
// geogram spread this over three classes: a NearestNeighborSearch interface, a KdTree base and a
// BalancedKdTree, because it also has an AdaptiveKdTree and a search over an ANN library. Only the
// balanced tree is vendored and it is used as a concrete type, so the three are one class here and
// nothing is virtual.
//
// geogram also returned the k nearest neighbours, through a sorted insertion buffer the caller
// supplied storage for. The mesher asks for one and reads only its distance, so the buffer is a
// single running minimum here. That is the same number: the traversal prunes on the distance to
// the furthest neighbour kept so far, which at k=1 is that minimum, and geogram's insert() kept
// the smaller of two equal distances just the same.
//
// geogram's points carry a dimension and a stride. Every tree built here is over 3d points packed
// three doubles to a point, so both are the constant below.

#pragma once

#include "geo_basic.h"
#include <algorithm>

namespace floatTetWild {
namespace geo {

    // Nearest neighbour search over a perfectly balanced kd-tree. No combinatorics is stored: the
    // two children of node n are 2n+1 and 2n+2. It works well for regular to moderately irregular
    // pointsets.
    class BalancedKdTree {
    public:
        // The dimension of the points, and the stride between two of them.
        static constexpr coord_index_t DIMENSION = 3;

        // \p points is an array of nb_points * DIMENSION doubles.
        void set_points(index_t nb_points, const double* points);

        // The squared distance from \p query_point, DIMENSION doubles, to the nearest point.
        double nearest_sq_dist(const double* query_point) const;

        index_t nb_points() const {
            return nb_points_;
        }

        const double* point_ptr(index_t i) const {
            geo_debug_assert(i < nb_points());
            return points_ + i * DIMENSION;
        }

    private:
        // Number of points stored in the leafs of the tree.
        static constexpr index_t MAX_LEAF_SIZE = 16;

        void build_tree();

        // Traverses the subtree under \p node_index, which holds the [b,e) point sequence,
        // lowering \p best_sq_dist as nearer points turn up. \p bbox_min and \p bbox_max are the
        // bounding box of that sequence, storage owned by the caller, modified here and restored
        // on exit, and \p box_dist is the squared distance from the query point to it, which is
        // what prunes subtrees that cannot hold a nearer point.
        void nearest_recursive(
            index_t node_index, index_t b, index_t e,
            double* bbox_min, double* bbox_max,
            double box_dist, const double* query_point,
            double& best_sq_dist
        ) const;

        // The largest node index the subtree rooted at \p node_id reaches, which is what the
        // per-node arrays have to be sized for.
        static index_t max_node_index(
            index_t node_id, index_t b, index_t e
        ) {
            if(e - b <= MAX_LEAF_SIZE) {
                return node_id;
            }
            index_t m = b + (e - b) / 2;
            return std::max(
                max_node_index(2 * node_id, b, m),
                max_node_index(2 * node_id + 1, m, e)
            );
        }

        // The coordinate the [b,e) sequence will be split along.
        coord_index_t best_splitting_coord(index_t b, index_t e);

        void create_kd_tree_recursive(
            index_t node_index, index_t b, index_t e
        ) {
            if(e - b <= MAX_LEAF_SIZE) {
                return;
            }
            index_t m = split_kd_node(node_index, b, e);
            create_kd_tree_recursive(2 * node_index, b, m);
            create_kd_tree_recursive(2 * node_index + 1, m, e);
        }

        // Computes and stores the splitting coordinate and value of \p node_index, which holds
        // the [b,e) sequence, and returns the m that splits it into [b,m) for the left child
        // (2*node_index) and [m,e) for the right (2*node_index+1).
        index_t split_kd_node(
            index_t node_index, index_t b, index_t e
        );

        index_t nb_points_ = 0;
        const double* points_ = nullptr;

        vector<index_t> point_index_;
        double bbox_min_[DIMENSION];
        double bbox_max_[DIMENSION];

        // One per node: the coordinate it splits along and the value it splits at.
        vector<coord_index_t> splitting_coord_;
        vector<double> splitting_val_;
    };

} }
