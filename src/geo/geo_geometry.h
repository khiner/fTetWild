// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/basic/geometry.h
// Copied rather than reimplemented: vec2/3/4, Box and the Geom helpers, unchanged

#pragma once

#include "geo_vecg.h"

/**
 * \file geogram/basic/geometry.h
 * \brief Geometric functions in 2d and 3d
 */

namespace floatTetWild {
namespace geo {

    /**
     * \brief Represents points and vectors in 3d.
     * \details Syntax is (mostly) compatible with GLSL.
     */
    typedef vecng<3, Numeric::float64> vec3;

    /************************************************************************/

    /**
     * \brief Geometric functions and utilities.
     */
    namespace Geom {

        /**
         * \brief Computes the signed volume of a 3d tetrahedron
         * \param[in] p1 first vertex of the tetrahedron
         * \param[in] p2 second vertex of the tetrahedron
         * \param[in] p3 third vertex of the tetrahedron
         * \param[in] p4 fourth vertex of the tetrahedron
         * \return the signed volume of the tetrahedron
         *  (\p p1, \p p2, \p p3, \p p4)
         */
        inline double tetra_signed_volume(
            const vec3& p1, const vec3& p2,
            const vec3& p3, const vec3& p4
        ) {
            return dot(p2 - p1, cross(p3 - p1, p4 - p1)) / 6.0;
        }

        /**
         * \brief Computes the signed volume of a 3d tetrahedron
         * \param[in] p1 first vertex of the tetrahedron
         * \param[in] p2 second vertex of the tetrahedron
         * \param[in] p3 third vertex of the tetrahedron
         * \param[in] p4 fourth vertex of the tetrahedron
         * \return the signed volume of the tetrahedron
         *  (\p p1, \p p2, \p p3, \p p4)
         */
        inline double tetra_signed_volume(
            const double* p1, const double* p2,
            const double* p3, const double* p4
        ) {
            return tetra_signed_volume(
                *reinterpret_cast<const vec3*>(p1),
                *reinterpret_cast<const vec3*>(p2),
                *reinterpret_cast<const vec3*>(p3),
                *reinterpret_cast<const vec3*>(p4)
            );
        }

    }

    /*******************************************************************/

    /**
     * \brief Axis-aligned bounding box.
     */
    class Box {
    public:
        double xyz_min[3];
        double xyz_max[3];

        /**
         * \brief Tests whether a box contains a point.
         * \param[in] b the point
         * \return true if this box contains \p b, false otherwise
         */
        bool contains(const vec3& b) const {
            for(coord_index_t c = 0; c < 3; ++c) {
                if(b[c] < xyz_min[c]) {
                    return false;
                }
                if(b[c] > xyz_max[c]) {
                    return false;
                }
            }
            return true;
        }
    };

    /**
     * \brief Computes the smallest Box that encloses two Boxes.
     * \param[out] target the smallest axis-aligned box
     *  that encloses \p B1 and \p B2
     * \param[in] B1 first box
     * \param[in] B2 second box
     */
    inline void bbox_union(Box& target, const Box& B1, const Box& B2) {
        for(coord_index_t c = 0; c < 3; ++c) {
            target.xyz_min[c] = std::min(B1.xyz_min[c], B2.xyz_min[c]);
            target.xyz_max[c] = std::max(B1.xyz_max[c], B2.xyz_max[c]);
        }
    }

    /******************************************************************/

} }

