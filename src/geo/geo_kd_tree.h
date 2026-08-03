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

#pragma once

#include "geo_basic.h"
#include <algorithm>

/**
 * \file geogram/points/kd_tree.h
 * \brief Nearest neighbor search over a balanced kd-tree
 */

namespace floatTetWild {
namespace geo {

    /**
     * \brief Nearest neighbor search over a perfectly balanced Kd-tree.
     * \details No combinatorics is stored: the two children of node n are 2n+1 and 2n+2. For
     *  regular to moderately irregular pointsets it works well.
     */
    class BalancedKdTree {
    public:
        /**
         * \brief Creates a new BalancedKdTree.
         * \param[in] dim dimension of the points
         */
        BalancedKdTree(coord_index_t dim) :
            dimension_(dim), bbox_min_(dim), bbox_max_(dim) {
        }

        /**
         * \brief Sets the points and creates the search data structure.
         * \param[in] nb_points number of points
         * \param[in] points an array of nb_points * dimension() doubles
         */
        void set_points(index_t nb_points, const double* points);

        /**
         * \brief The squared distance from a query point to the nearest of the points.
         * \param[in] query_point as an array of dimension() doubles
         */
        double nearest_sq_dist(const double* query_point) const;

        coord_index_t dimension() const {
            return dimension_;
        }

        index_t nb_points() const {
            return nb_points_;
        }

        const double* point_ptr(index_t i) const {
            geo_debug_assert(i < nb_points());
            return points_ + i * dimension_;
        }

    private:
        /**
         * \brief Number of points stored in the leafs of the tree.
         */
        static constexpr index_t MAX_LEAF_SIZE = 16;

        /**
         * \brief Builds the tree.
         */
        void build_tree();

        /**
         * \brief The recursive function to implement KdTree traversal.
         * \details Traverses the subtree under the node_index node that corresponds to the
         *  [b,e) point sequence, lowering best_sq_dist as nearer points turn up.
         * \param[in] node_index index of the current node in the Kd tree
         * \param[in] b index of the first point in the subtree under node \p node_index
         * \param[in] e one position past the index of the last point in the subtree
         * \param[in,out] bbox_min , bbox_max the bounding box of the [b,e) sequence, allocated and
         *  managed by the caller. Modified by the function and restored on exit.
         * \param[in] box_dist squared distance between the query point and that bounding box. It
         *  is used to prune traversals that cannot hold a nearer point.
         * \param[in] query_point the query point
         * \param[in,out] best_sq_dist the nearest squared distance so far
         */
        void nearest_recursive(
            index_t node_index, index_t b, index_t e,
            double* bbox_min, double* bbox_max,
            double box_dist, const double* query_point,
            double& best_sq_dist
        ) const;

        /**
         * \brief Computes the minimum and maximum point coordinates
         *   along a coordinate.
         * \param[in] b first index of the point sequence
         * \param[in] e one position past the last index of the point sequence
         * \param[in] coord coordinate along which the extent is measured
         * \param[out] minval , maxval minimum and maximum
         */
        void get_minmax(
            index_t b, index_t e, coord_index_t coord,
            double& minval, double& maxval
        ) const {
            minval = Numeric::max_float64();
            maxval = Numeric::min_float64();
            for(index_t i = b; i < e; ++i) {
                double val = point_ptr(point_index_[i])[coord];
                minval = std::min(minval, val);
                maxval = std::max(maxval, val);
            }
        }

        /**
         * \brief Computes the extent of a point sequence
         *  along a given coordinate.
         * \param[in] b first index of the point sequence
         * \param[in] e one position past the last index of the point sequence
         * \param[in] coord coordinate along which the extent is measured
         * \return the extent of the sequence along the coordinate
         */
        double spread(index_t b, index_t e, coord_index_t coord) const {
            double minval,maxval;
            get_minmax(b,e,coord,minval,maxval);
            return maxval - minval;
        }

        /**
         * \brief Returns the maximum node index in subtree.
         * \param[in] node_id node index of the subtree
         * \param[in] b first index of the points sequence in the subtree
         * \param[in] e one position past the last index of the point
         *  sequence in the subtree
         */
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

        /**
         * \brief Computes the coordinate along which a point
         *   sequence will be split.
         * \param[in] b first index of the point sequence
         * \param[in] e one position past the last index of the point sequence
         */
        coord_index_t best_splitting_coord(index_t b, index_t e);

        /**
         * \brief Creates the subtree under a node.
         * \param[in] node_index index of the node that represents
         *  the subtree to create
         * \param[in] b first index of the point sequence in the subtree
         * \param[in] e one position past the last index of the point
         *  index in the subtree
         */
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

        /**
         * \brief Computes and stores the splitting coordinate
         *  and splitting value of the node node_index, that
         *  corresponds to the [b,e) points sequence.
         *
         * \return a node index m. The point sequences
         *  [b,m) and [m,e) correspond to the left
         *  child (2*node_index) and right child (2*node_index+1)
         *  of node_index.
         */
        index_t split_kd_node(
            index_t node_index, index_t b, index_t e
        );

        coord_index_t dimension_;
        index_t nb_points_ = 0;
        const double* points_ = nullptr;

        vector<index_t> point_index_;
        vector<double> bbox_min_;
        vector<double> bbox_max_;

        /**
         * \brief One per node, splitting coordinate.
         */
        vector<coord_index_t> splitting_coord_;

        /**
         * \brief One per node, splitting coordinate value.
         */
        vector<double> splitting_val_;
    };

    /*********************************************************************/
} }
