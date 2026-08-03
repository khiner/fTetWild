/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */

#pragma once
/**
 * \file mesh_AABB.h
 * \brief Axis Aligned Bounding Box trees for accelerating
 *  geometric queries that operate on a Mesh.
 */

#include <floattetwild/geo_mesh.h>
#include <floattetwild/geo_geometry_nd.h>

namespace floatTetWild {

    /**
     * \brief Axis Aligned Bounding Box tree of mesh facets.
     * \details Used to quickly compute facet intersection and
     *  to locate the nearest facet from 3d query points.
     */
    class MeshFacetsAABBWithEps {
    public:
        /**
         * \brief Creates the Axis Aligned Bounding Boxes tree.
         * \param[in] M the input mesh, already triangulated and already Morton-ordered by
         *  mesh_reorder(): geogram did both here, and every caller now does them beforehand.
         * \pre M.facets.are_simplices()
         */
        MeshFacetsAABBWithEps(const geo::Mesh& M);

        /**
         * \brief Finds the nearest facet from an arbitrary 3d query point.
         * \param[in] p query point
         * \param[out] nearest_point nearest point on the surface
         * \param[out] sq_dist squared distance between p and the surface.
         * \return the index of the facet nearest to point p.
         */
        geo::index_t nearest_facet(
            const geo::vec3& p, geo::vec3& nearest_point, double& sq_dist
        ) const {
            geo::index_t nearest_facet;
            get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
            nearest_facet_recursive(
                p,
                nearest_facet, nearest_point, sq_dist,
                1, 0, mesh_.facets.nb()
            );
            return nearest_facet;
        }

        /*
         * Finds the nearest facet on the surface, but stops early if a
         * point within a given distance is found. Starting from a facet the
         * caller already has saves the search for a hint.
         */
        void facet_in_envelope_with_hint(
            const geo::vec3& p, double sq_epsilon,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist
        ) const {
            if(nearest_facet == geo::NO_FACET) {
                get_nearest_facet_hint(
                    p, nearest_facet, nearest_point, sq_dist
                );
            }
            facet_in_envelope_recursive(
                p, sq_epsilon,
                nearest_facet, nearest_point, sq_dist,
                1, 0, mesh_.facets.nb()
            );
        }

        /*
         * \overload With no facet to start from.
         */
        geo::index_t facet_in_envelope(
            const geo::vec3& p, double sq_epsilon, geo::vec3& nearest_point, double& sq_dist
        ) const {
            geo::index_t nearest_facet = geo::NO_FACET;
            facet_in_envelope_with_hint(p, sq_epsilon, nearest_facet, nearest_point, sq_dist);
            return nearest_facet;
        }

    protected:
        /**
         * \brief Computes a reasonable initialization for
         *  nearest facet search.
         *
         * \details A good initialization makes the algorithm faster,
         *  by allowing early pruning of subtrees that provably
         *  do not contain the nearest neighbor.
         *
         * \param[in] p query point
         * \param[out] nearest_facet a facet reasonably near p
         * \param[out] nearest_point a point in nearest_facet
         * \param[out] sq_dist squared distance between p and nearest_point
         */
        void get_nearest_facet_hint(
            const geo::vec3& p,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist
        ) const;

        /**
         * \brief The recursive function used by the implementation
         *  of nearest_facet().
         *
         * \details The first call may use get_nearest_facet_hint()
         * to initialize nearest_facet, nearest_point and sq_dist,
         * as done in nearest_facet().
         *
         * \param[in] p query point
         * \param[in,out] nearest_facet the nearest facet so far,
         * \param[in,out] nearest_point a point in nearest_facet
         * \param[in,out] sq_dist squared distance between p and nearest_point
         * \param[in] n index of the current node in the AABB tree
         * \param[in] b index of the first facet in the subtree under node \p n
         * \param[in] e one position past the index of the last facet in the
         *  subtree under node \p n
         */
        void nearest_facet_recursive(
            const geo::vec3& p,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist,
            geo::index_t n, geo::index_t b, geo::index_t e
        ) const;

        /*
         * Same as before, but stops early if a point within a given distance
         * is found.
         */
        void facet_in_envelope_recursive(
            const geo::vec3& p, double sq_epsilon,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist,
            geo::index_t n, geo::index_t b, geo::index_t e
        ) const;

        geo::vector<geo::Box> bboxes_;
        const geo::Mesh& mesh_;
    };

}

namespace floatTetWild {

inline void get_point_facet_nearest_point(
    const geo::Mesh& M,
    const geo::vec3& p,
    geo::index_t f,
    geo::vec3& nearest_p,
    double& squared_dist
) {
    using namespace floatTetWild::geo;
    geo_debug_assert(M.facets.nb_vertices(f) == 3);
    index_t c = M.facets.corners_begin(f);
    const vec3& p1 = M.vertices.point(M.facet_corners.vertex(c));
    ++c;
    const vec3& p2 = M.vertices.point(M.facet_corners.vertex(c));
    ++c;
    const vec3& p3 = M.vertices.point(M.facet_corners.vertex(c));
    double lambda1, lambda2, lambda3;  // barycentric coords, not used.
    squared_dist = Geom::point_triangle_squared_distance(
        p, p1, p2, p3, nearest_p, lambda1, lambda2, lambda3
    );
}
}
