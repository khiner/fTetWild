// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/points/nn_search.h
// Copied rather than reimplemented: nearest neighbour interface

#pragma once

#include "geo_basic.h"

#include <string>

/**
 * \file geogram/points/nn_search.h
 * \brief Abstract interface for nearest neighbor searching
 */

namespace floatTetWild {
namespace geo {

    /**
     * \brief Abstract interface for nearest neighbor search algorithms.
     * \details
     * Given a point set in arbitrary dimension, creates a data
     * structure for efficient nearest neighbor queries.
     *
     */
    class NearestNeighborSearch {
    public:

        /**
         * \brief Sets the points and create the search data structure.
         * \param[in] nb_points number of points
         * \param[in] points an array of nb_points * dimension()
         */
        virtual void set_points(index_t nb_points, const double* points);

        /**
         * \brief Tests whether the stride variant of set_points() is supported
         * \return true if stride different from dimension can be used
         *  in set_points(), false otherwise
         */
        virtual bool stride_supported() const;

        /**
         * \brief Sets the points and create the search data structure.
         * \details This variant has a stride parameter. Call
         * stride_supported() before to check whether it is supported.
         *
         * \param[in] nb_points number of points
         * \param[in] points an array of nb_points * dimension()
         * \param[in] stride number of doubles between two consecutive
         *  points (stride=dimension() by default).
         */
        virtual void set_points(
            index_t nb_points, const double* points, index_t stride
        );

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  coordinates.
         * \param[in] nb_neighbors number of neighbors to be searched.
         *  Should be smaller or equal to nb_points() (else it triggers
         *  an assertion)
         * \param[in] query_point as an array of dimension() doubles
         * \param[out] neighbors array of nb_neighbors index_t
         * \param[out] neighbors_sq_dist array of nb_neighbors doubles
         */
        virtual void get_nearest_neighbors(
            index_t nb_neighbors,
            const double* query_point,
            index_t* neighbors,
            double* neighbors_sq_dist
        ) const = 0;


        /**
         * \brief A structure to discriminate between the two
         *  versions of get_nearest_neighbors()
         */
        struct KeepInitialValues {
        };

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  coordinates. Uses input neighbors and squared distance as
         *  an initialization.
         * \details Default implementation ignores the input values.
         *  Derived classes may have more efficient implementations.
         * \param[in] nb_neighbors number of neighbors to be searched.
         *  Should be smaller or equal to nb_points() (else it triggers
         *  an assertion)
         * \param[in] query_point as an array of dimension() doubles
         * \param[in,out] neighbors array of nb_neighbors index_t
         * \param[in,out] neighbors_sq_dist array of nb_neighbors doubles
         * \param[in] dummy a dummy parameter to discriminate between the
         *  two forms of get_nearest_neighbors()
         */
        virtual void get_nearest_neighbors(
            index_t nb_neighbors,
            const double* query_point,
            index_t* neighbors,
            double* neighbors_sq_dist,
            KeepInitialValues dummy
        ) const;

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  its index.
         * \details For some implementation, may be faster than
         *  nearest neighbor search by point coordinates.
         * \param[in] nb_neighbors number of neighbors to be searched.
         *  Should be smaller or equal to nb_points() (else it triggers
         *  an assertion)
         * \param[in] query_point as the index of one of the points that
         *  was inserted in this NearestNeighborSearch
         * \param[out] neighbors array of nb_neighbors index_t
         * \param[out] neighbors_sq_dist array of nb_neighbors doubles
         */
        virtual void get_nearest_neighbors(
            index_t nb_neighbors,
            index_t query_point,
            index_t* neighbors,
            double* neighbors_sq_dist
        ) const;

        /**
         * \brief Nearest neighbor search.
         * \param[in] query_point array of dimension() doubles
         * \return the index of the nearest neighbor from \p query_point
         */
        index_t get_nearest_neighbor(const double* query_point) const {
            index_t result;
            double sq_dist;
            get_nearest_neighbors(1, query_point, &result, &sq_dist);
	    geo_assert(result < nb_points());
            return result;
        }

        /**
         * \brief Gets the dimension of the points.
         * \return the dimension
         */
        coord_index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Gets the number of points.
         * \return the number of points
         */
        index_t nb_points() const {
            return nb_points_;
        }

        /**
         * \brief Gets a point by its index
         * \param[in] i index of the point
         * \return a const pointer to the coordinates of the point
         */
        const double* point_ptr(index_t i) const {
            geo_debug_assert(i < nb_points());
            return points_ + i * stride_;
        }

        /**
         * \brief Search can be exact or approximate. Approximate
         *  search may be faster.
         * \return true if nearest neighbor search is exact,
         *   false otherwise
         */
        bool exact() const {
            return exact_;
        }

        /**
         * \brief Search can be exact or approximate. Approximate
         *  search may be faster.
         * \param[in] x true if nearest neighbor search is exact,
         *   false otherwise. Default mode is exact.
         */
        virtual void set_exact(bool x);

    protected:
        /**
         * \brief Constructs a NearestNeighborSearch.
         * \param[in] dimension dimension of the points
         */
        NearestNeighborSearch(coord_index_t dimension);

        /**
         * \brief NearestNeighborSearch destructor
         * \details Public here, unlike in geogram. geogram made it protected because these are
         *  reference counted and are only ever destroyed through a NearestNeighborSearch_var,
         *  which is not the case any more.
         */
    public:
        virtual ~NearestNeighborSearch();
    protected:

    protected:
        coord_index_t dimension_;
        index_t nb_points_;
        index_t stride_;
        const double* points_;
        bool exact_;
    };
}
}
