// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/points/kd_tree.cpp
// Copied rather than reimplemented: balanced kd-tree; nth_element splitting again, so ties follow the incoming order

#include "geo_kd_tree.h"
#include "geo_geometry_nd.h"
#include "geo_parallel.h"

namespace {

    using namespace floatTetWild::geo;

    /**
     * \brief Comparison functor used to
     *  sort the point indices. Used by
     *  BalancedKdTree.
     */
    class ComparePointCoord {
    public:
        /**
         * \brief Creates a new ComparePointCoord
         * \param[in] nb_points number of points
         * \param[in] points pointer to first point
         * \param[in] stride number of doubles between two
         *  consecutive points in array (=dimension if point
         *  array is compact).
         * \param[in] splitting_coord the coordinate to compare
         */
        ComparePointCoord(
            index_t nb_points,
            const double* points,
            index_t stride,
            coord_index_t splitting_coord
        ) :
            nb_points_(nb_points),
            points_(points),
            stride_(stride),
            splitting_coord_(splitting_coord) {
            geo_argused(nb_points_);
        }

        /**
         * \brief Compares to point indices (does the
         * indirection and coordinate lookup).
         * \param[in] i index of first point to compare
         * \param[in] j index of second point to compare
         * \return true if point \p i is before point \p j, false otherwise
         */
        bool operator() (index_t i, index_t j) const {
            geo_debug_assert(i < nb_points_);
            geo_debug_assert(j < nb_points_);
            return
                (points_ + i * stride_)[splitting_coord_] <
                (points_ + j * stride_)[splitting_coord_]
                ;
        }

    private:
        index_t nb_points_;
        const double* points_;
        index_t stride_;
        coord_index_t splitting_coord_;
    };
}

/****************************************************************************/

namespace floatTetWild {
namespace geo {

    KdTree::KdTree(coord_index_t dim) :
        NearestNeighborSearch(dim),
        bbox_min_(dim),
        bbox_max_(dim),
        root_(NO_INDEX) {
    }

    KdTree::~KdTree() {
    }

    bool KdTree::stride_supported() const {
        return true;
    }

    void KdTree::set_points(
        index_t nb_points, const double* points, index_t stride
    ) {
        nb_points_ = nb_points;
        points_ = points;
        stride_ = stride;

        point_index_.resize(nb_points);
        for(index_t i = 0; i < nb_points; i++) {
            point_index_[i] = i;
        }

        // Compute the bounding box.
        for(coord_index_t c = 0; c < dimension(); ++c) {
            bbox_min_[c] =  Numeric::max_float64();
            bbox_max_[c] = -Numeric::max_float64();
        }
        for(index_t i = 0; i < nb_points; ++i) {
            const double* p = point_ptr(i);
            for(coord_index_t c = 0; c < dimension(); ++c) {
                bbox_min_[c] = std::min(bbox_min_[c], p[c]);
                bbox_max_[c] = std::max(bbox_max_[c], p[c]);
            }
        }

        root_ = build_tree();
    }

    void KdTree::set_points(
        index_t nb_points, const double* points
    ) {
        set_points(nb_points, points, dimension());
    }


    void KdTree::get_nearest_neighbors(
        index_t nb_neighbors,
        const double* query_point,
        index_t* neighbors,
        double* neighbors_sq_dist
    ) const {

        geo_debug_assert(nb_neighbors <= nb_points());

        // Compute distance between query point and global bounding box
        // and copy global bounding box to local variables (bbox_min, bbox_max),
        // allocated on the stack. bbox_min and bbox_max are updated during the
        // traversal of the BalancedKdTree (see
        // get_nearest_neighbors_recursive()). They are necessary to
        // compute the distance between the query point and the
        // bbox of the current node.
        double box_dist = 0.0;
        double* bbox_min = (double*) (alloca(dimension() * sizeof(double)));
        double* bbox_max = (double*) (alloca(dimension() * sizeof(double)));
        init_bbox_and_bbox_dist_for_traversal(
            bbox_min, bbox_max, box_dist, query_point
        );
        NearestNeighbors NN(
            nb_neighbors,
            neighbors,
            neighbors_sq_dist,
            (index_t*)alloca(sizeof(index_t) * (nb_neighbors+1)),
            (double*)alloca(sizeof(double) * (nb_neighbors+1))
        );
        get_nearest_neighbors_recursive(
            root_, 0, nb_points(), bbox_min, bbox_max, box_dist, query_point, NN
        );
        NN.copy_to_user();
    }

    void KdTree::get_nearest_neighbors(
        index_t nb_neighbors,
        const double* query_point,
        index_t* neighbors,
        double* neighbors_sq_dist,
        KeepInitialValues KV
    ) const {
        geo_debug_assert(nb_neighbors <= nb_points());
        geo_argused(KV);
        // Compute distance between query point and global bounding box
        // and copy global bounding box to local variables (bbox_min, bbox_max),
        // allocated on the stack. bbox_min and bbox_max are updated during the
        // traversal of the BalancedKdTree
        // (see get_nearest_neighbors_recursive()). They
        // are necessary to compute the distance between the query point and the
        // bbox of the current node.
        double box_dist = 0.0;
        double* bbox_min = (double*) (alloca(dimension() * sizeof(double)));
        double* bbox_max = (double*) (alloca(dimension() * sizeof(double)));
        init_bbox_and_bbox_dist_for_traversal(
            bbox_min, bbox_max, box_dist, query_point
        );
        NearestNeighbors NN(
            nb_neighbors,
            neighbors,
            neighbors_sq_dist,
            (index_t*)alloca(sizeof(index_t) * (nb_neighbors+1)),
            (double*)alloca(sizeof(double) * (nb_neighbors+1))
        );
        NN.copy_from_user();
        get_nearest_neighbors_recursive(
            root_, 0, nb_points(), bbox_min, bbox_max, box_dist, query_point, NN
        );
        NN.copy_to_user();
    }

    void KdTree::get_nearest_neighbors(
        index_t nb_neighbors,
        index_t q_index,
        index_t* neighbors,
        double* neighbors_sq_dist
    ) const {
        // TODO: optimized version that uses the fact that
        // we know that query_point is in the search data
        // structure already.
        // (I tried something already, see in the Attic,
        //  but it did not give any significant speedup).
        get_nearest_neighbors(
            nb_neighbors, point_ptr(q_index),
            neighbors, neighbors_sq_dist
        );
    }

    void KdTree::get_nearest_neighbors_recursive(
        index_t node_index, index_t b, index_t e,
        double* bbox_min, double* bbox_max, double box_dist,
        const double* query_point, NearestNeighbors& NN
    ) const {
        geo_debug_assert(e > b);

        // Simple case (node is a leaf)
        if((e - b) <= MAX_LEAF_SIZE) {
            get_nearest_neighbors_leaf(node_index, b, e, query_point, NN);
            return;
        }

        // Get node attributes (virtual function call).

        index_t left_node_index;
        index_t right_node_index;
        coord_index_t coord;
        index_t m;
        double val;

        get_node(
            node_index, b, e,
            left_node_index, right_node_index,
            coord, m, val
        );

        double cut_diff = query_point[coord] - val;

        // If the query point is on the left side
        if(cut_diff < 0.0) {

            // Traverse left subtree
            {
                double bbox_max_save = bbox_max[coord];
                bbox_max[coord] = val;
                get_nearest_neighbors_recursive(
                    left_node_index, b, m,
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_max[coord] = bbox_max_save;
            }

            // Update bbox distance (now measures the
            // distance to the bbox of the right subtree)
            double box_diff = bbox_min[coord] - query_point[coord];
            if(box_diff > 0.0) {
                box_dist -= geo_sqr(box_diff);
            }
            box_dist += geo_sqr(cut_diff);

            // Traverse the right subtree, only if bbox
            // distance is nearer than furthest neighbor,
            // else there is no chance that the right
            // subtree contains points that will change
            // anything in the nearest neighbors NN.
            if(box_dist <= NN.furthest_neighbor_sq_dist()) {
                double bbox_min_save = bbox_min[coord];
                bbox_min[coord] = val;
                get_nearest_neighbors_recursive(
                    right_node_index, m, e,
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_min[coord] = bbox_min_save;
            }
        } else {
            // else the query point is on the right side
            // (then do the same with left and right subtree
            //  permutted).
            {
                double bbox_min_save = bbox_min[coord];
                bbox_min[coord] = val;
                get_nearest_neighbors_recursive(
                    right_node_index, m, e,
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_min[coord] = bbox_min_save;
            }

            // Update bbox distance (now measures the
            // distance to the bbox of the left subtree)
            double box_diff = query_point[coord] - bbox_max[coord];
            if(box_diff > 0.0) {
                box_dist -= geo_sqr(box_diff);
            }
            box_dist += geo_sqr(cut_diff);

            if(box_dist <= NN.furthest_neighbor_sq_dist()) {
                double bbox_max_save = bbox_max[coord];
                bbox_max[coord] = val;
                get_nearest_neighbors_recursive(
                    left_node_index, b, m,
                    bbox_min, bbox_max, box_dist, query_point, NN
                );
                bbox_max[coord] = bbox_max_save;
            }
        }
    }

    void KdTree::get_nearest_neighbors_leaf(
        index_t node_index, index_t b, index_t e,
        const double* query_point,
        NearestNeighbors& NN
    ) const {
        geo_argused(node_index);
        NN.nb_visited += (e-b);
        double R = NN.furthest_neighbor_sq_dist();
        index_t nb = e-b;
        const index_t* geo_restrict idx = &point_index_[b];

        // TODO: check generated ASM (I'd like to have AVX here).
        // We may need to dispatch according to dimension.

        index_t local_idx[MAX_LEAF_SIZE];
        double  local_sq_dist[MAX_LEAF_SIZE];

        // Cache indices and computed distances in local
        // array. I guess AVX likes that (to be checked).
        // Not sure, because access to p is indirect, maybe
        // I should copy the points to local memory before
        // computing the distances (or having another copy
        // of the points array that I pre-reorder so that
        // leaf's points are in a contiguous chunk of memory),
        // to be tested...
        for(index_t ii=0; ii<nb; ++ii) {
            index_t i = idx[ii];
            const double* geo_restrict p = point_ptr(i);
            double sq_dist = Geom::distance2(
                query_point, p, dimension()
            );
            local_idx[ii] = i;
            local_sq_dist[ii] = sq_dist;
        }

        // Now insert the points that are nearer to query
        // point than NN's bounding ball.
        for(index_t ii=0; ii<nb; ++ii) {
            double sq_dist = local_sq_dist[ii];
            if(sq_dist <= R) {
                NN.insert(local_idx[ii],sq_dist);
                R = NN.furthest_neighbor_sq_dist();
            }
        }
    }

    void KdTree::init_bbox_and_bbox_dist_for_traversal(
        double* bbox_min, double* bbox_max,
        double& box_dist, const double* query_point
    ) const {
        // Compute distance between query point and global bounding box
        // and copy global bounding box to local variables (bbox_min, bbox_max),
        // allocated on the stack. bbox_min and bbox_max are updated during the
        // traversal of the KdTree (see get_nearest_neighbors_recursive()). They
        // are necessary to compute the distance between the query point and the
        // bbox of the current node.
        box_dist = 0.0;
        for(coord_index_t c = 0; c < dimension(); ++c) {
            bbox_min[c] = bbox_min_[c];
            bbox_max[c] = bbox_max_[c];
            if(query_point[c] < bbox_min_[c]) {
                box_dist += geo_sqr(bbox_min_[c] - query_point[c]);
            } else if(query_point[c] > bbox_max_[c]) {
                box_dist += geo_sqr(bbox_max_[c] - query_point[c]);
            }
        }
    }

/****************************************************************************/

    BalancedKdTree::BalancedKdTree(coord_index_t dim) :
        KdTree(dim),
        m0_(max_index_t()),
        m1_(max_index_t()),
        m2_(max_index_t()),
        m3_(max_index_t()),
        m4_(max_index_t()),
        m5_(max_index_t()),
        m6_(max_index_t()),
        m7_(max_index_t()),
        m8_(max_index_t()) {
    }

    BalancedKdTree::~BalancedKdTree() {
    }

    index_t BalancedKdTree::build_tree() {
        index_t sz = max_node_index(1, 0, nb_points()) + 1;
        splitting_coord_.resize(sz);
        splitting_val_.resize(sz);

        // If there are more than 16*MAX_LEAF_SIZE (=256) points,
        // create the tree in parallel
        if(
            nb_points() >= (16 * MAX_LEAF_SIZE) &&
            Process::maximum_concurrent_threads() > 1
        ) {
            m0_ = 0;
            m8_ = nb_points();
            // Create the first level of the tree
            m4_ = split_kd_node(1, m0_, m8_);

            // Create the second level of the tree
            //  (using two threads)
            parallel(
                [this]() { m2_ = split_kd_node(2, m0_, m4_); },
                [this]() { m6_ = split_kd_node(3, m4_, m8_); }
            );

            // Create the third level of the tree
            //  (using four threads)
            parallel(
                [this]() { m1_ = split_kd_node(4, m0_, m2_); },
                [this]() { m3_ = split_kd_node(5, m2_, m4_); },
                [this]() { m5_ = split_kd_node(6, m4_, m6_); },
                [this]() { m7_ = split_kd_node(7, m6_, m8_); }
            );

            // Create the fourth level of the tree
            //  (using eight threads)
            parallel(
                [this]() { create_kd_tree_recursive(8 , m0_, m1_); },
                [this]() { create_kd_tree_recursive(9 , m1_, m2_); },
                [this]() { create_kd_tree_recursive(10, m2_, m3_); },
                [this]() { create_kd_tree_recursive(11, m3_, m4_); },
                [this]() { create_kd_tree_recursive(12, m4_, m5_); },
                [this]() { create_kd_tree_recursive(13, m5_, m6_); },
                [this]() { create_kd_tree_recursive(14, m6_, m7_); },
                [this]() { create_kd_tree_recursive(15, m7_, m8_); }
            );

        } else {
            create_kd_tree_recursive(1, 0, nb_points());
        }

        // Root node is number 1.
        // This is because "children at 2*n and 2*n+1" does not
        // work with 0 !!
        return 1;
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
            ComparePointCoord(
                nb_points_, points_, stride_, splitting_coord
            )
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
        for(coord_index_t c = 1; c < dimension(); ++c) {
            double coord_spread = spread(b, e, c);
            if(coord_spread > max_spread) {
                result = c;
                max_spread = coord_spread;
            }
        }
        return result;
    }

    void BalancedKdTree::get_node(
        index_t n, index_t b, index_t e,
        index_t& left_child, index_t& right_child,
        coord_index_t&  splitting_coord,
        index_t& m,
        double& splitting_val
    ) const {
        left_child = 2*n;
        right_child = 2*n+1;
        splitting_coord = splitting_coord_[n];
        m = b + (e - b) / 2;
        splitting_val = splitting_val_[n];
    }









    /*************************************************************************/

} }
