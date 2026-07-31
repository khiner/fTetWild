// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/points/nn_search.cpp
// Copied rather than reimplemented: nearest neighbour interface

#include "geo_nn_search.h"

/****************************************************************************/

namespace floatTetWild {
namespace geo {

    NearestNeighborSearch::NearestNeighborSearch(
        coord_index_t dimension
    ) :
        dimension_(dimension),
        nb_points_(0),
        stride_(0),
        points_(nullptr),
        exact_(true) {
    }

    void NearestNeighborSearch::get_nearest_neighbors(
        index_t nb_neighbors,
        const double* query_point,
        index_t* neighbors,
        double* neighbors_sq_dist,
        KeepInitialValues
    ) const {
        get_nearest_neighbors(
            nb_neighbors,
            query_point,
            neighbors,
            neighbors_sq_dist
        );
    }

    void NearestNeighborSearch::get_nearest_neighbors(
        index_t nb_neighbors,
        index_t query_point,
        index_t* neighbors,
        double* neighbors_sq_dist
    ) const {
        get_nearest_neighbors(
            nb_neighbors,
            point_ptr(query_point),
            neighbors,
            neighbors_sq_dist
        );
    }

    void NearestNeighborSearch::set_points(
        index_t nb_points, const double* points
    ) {
        nb_points_ = nb_points;
        points_ = points;
        stride_ = dimension_;
    }

    bool NearestNeighborSearch::stride_supported() const {
        return false;
    }

    void NearestNeighborSearch::set_points(
        index_t nb_points, const double* points, index_t stride
    ) {
        if(stride == index_t(dimension())) {
            set_points(nb_points, points);
            return;
        }
        geo_assert(stride_supported());
        nb_points_ = nb_points;
        points_ = points;
        stride_ = stride;
    }

    void NearestNeighborSearch::set_exact(bool x) {
        exact_ = x;
    }

    NearestNeighborSearch::~NearestNeighborSearch() {
    }

} }
