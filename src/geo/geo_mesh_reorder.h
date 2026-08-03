// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/mesh/mesh_reorder.h
// Copied rather than reimplemented: spatial sort; nth_element is not stable so the incoming arrangement decides ties

#pragma once

#include "geo_mesh.h"

/**
 * \file geogram/mesh/mesh_reorder.h
 * \brief Reorders the elements in a mesh to improve data locality
 */

namespace floatTetWild {
namespace geo {


    /**
     * \brief Reorders all the elements of a mesh in Morton order, which improves data locality.
     * \details geogram also has a Hilbert ordering here, which has a continuous mapping between
     *  indices and space. No caller asked for it. The Hilbert curve itself is still what
     *  compute_BRIO_order() sorts along.
     * \param[in] M the mesh to reorder
     * \param[out] facet_permutation if non-null, receives the permutation the facets were
     *  reordered by: facet k is the one that was at facet_permutation[k]. Left alone when the
     *  mesh has no facets.
     */
    void mesh_reorder(Mesh& M, vector<index_t>* facet_permutation = nullptr);

    /**
     * \brief Computes the BRIO order for a set of 3D points.
     * \details It is used to accelerate incremental insertion
     *  in Delaunay triangulation. See the following reference:
     *  -Incremental constructions con brio. Nina Amenta, Sunghee Choi,
     *   Gunter Rote, Symposium on Computational Geometry conf. proc.,
     *   2003
     * \param[in] nb_vertices number of vertices to sort
     * \param[in] vertices pointer to the coordinates of the vertices
     * \param[out] sorted_indices a vector of element indices to
     *  be sorted spatially
     * \param[in] dimension number of vertices coordinates
     * \param[in] stride number of doubles between two consecutive vertices
     * \param[in] threshold minimum size of interval to be sorted
     * \param[in] ratio splitting ratio between current interval and
     *  the rest to be sorted
     * \param[out] levels if non-nullptr, indices that correspond to level l are
     *   in the range levels[l] (included) ... levels[l+1] (excluded)
     */
    void compute_BRIO_order(
        index_t nb_vertices, const double* vertices,
        vector<index_t>& sorted_indices,
        index_t dimension,
        index_t stride = 3,
        index_t threshold = 64,
        double ratio = 0.125,
        vector<index_t>* levels = nullptr
    );

}
}
