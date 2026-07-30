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
     * \brief Strategy for spatial sorting.
     */
    enum MeshOrder {
        /**
         * Hilbert ordering improves data locality and
         * has a continuous mapping between indices and space.
         */
        MESH_ORDER_HILBERT,
        /**
         * Morton ordering improves data locality and is
         * a bit simpler than Hilbert ordering.
         */
        MESH_ORDER_MORTON
    };

    /**
     * \brief Reorders all the elements of a mesh.
     * \details It is used for both improving data locality
     *  and for implementing mesh partitioning.
     * \param[in] M the mesh to reorder
     * \param[in] order the reordering scheme, one of MESH_ORDER_HILBERT,
     *  MESH_ORDER_MORTION
     */
    void mesh_reorder(
        Mesh& M, MeshOrder order = MESH_ORDER_HILBERT
    );

    /**
     * \brief Computes the Hilbert order for a set of 3D points.
     * \details
     *  This variant sorts a subsequence of the indices vector.
     *  The implementation is inspired by:
     *  - Christophe Delage and Olivier Devillers. Spatial Sorting.
     *   In CGAL User and Reference Manual. CGAL Editorial Board,
     *   3.9 edition, 2011
     * \param[in] total_nb_vertices total number of vertices
     *   in the \p vertices array, used to test indices in debug mode
     * \param[in] vertices pointer to the coordinates of the vertices
     * \param[in,out] sorted_indices a vector of vertex indices, sorted
     *  spatially on exit
     * \param[in] first index of the first element in \p sorted_indices
     *  to be sorted
     * \param[in] last one position past the index of the last element
     *  in \p sorted_indices to be sorted
     * \param[in] dimension number of vertices coordinates
     * \param[in] stride number of doubles between two consecutive vertices
     */

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
