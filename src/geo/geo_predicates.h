// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/numerics/predicates.h
// Copied rather than reimplemented: the symbolic-perturbation tie-breaking is what decides
// Delaunay's output on the cospherical configurations a regular background grid is full of.
//
// Trimmed to the predicates fTetWild reaches. The declarations are unchanged.

#pragma once

#include "geo_geometry.h"

namespace floatTetWild {
namespace geo {

    /**
     * \brief PCK (Predicate Construction Kit) implements a set of
     *  geometric predicates. PCK uses arithmetic filters (Meyer and Pion),
     *  expansion arithmetics (Shewchuk) and simulation of simplicity
     *  (Edelsbrunner).
     */
    namespace PCK {

        /**
         * \brief Tests whether a point is in the circumscribed sphere of
         *  four other points.
         * \details If the tetrahedron \p p0, \p p1, \p p2, \p p3 is oriented
         *  negatively, then the result is inversed. Symbolic perturbation
         *  is applied whenever the five points are cospherical.
         * \param[in] p0 , p1 , p2 , p3 vertices of the tetrahedron
         * \param[in] p4 point to be tested
         * \retval POSITIVE whenever \p p4 is in the circumscribed sphere
         *  of \p p0, \p p1, \p p2, \p p3
         * \retval NEGATIVE whenever \p p4 is outside the circumscribed sphere
         *  of \p p0, \p p1, \p p2, \p p3
         */
        Sign in_sphere_3d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4
        );

        /**
         * \brief Computes the orientation predicate in 3d.
         * \details Computes the sign of the signed volume of
         *  the tetrahedron p0, p1, p2, p3.
         * \param[in] p0 , p1 , p2 , p3 vertices of the tetrahedron
         * \retval POSITIVE if the tetrahedron is oriented positively
         * \retval ZERO if the tetrahedron is flat
         * \retval NEGATIVE if the tetrahedron is oriented negatively
         */
        Sign orient_3d(
            const double* p0, const double* p1,
            const double* p2, const double* p3
        );

        /**
         * \brief Computes the (approximate) orientation predicate in 3d.
         * \details Computes the sign of the (approximate) signed volume of
         *  the tetrahedron p0, p1, p2, p3.
         * \param[in] p0 first vertex of the tetrahedron
         * \param[in] p1 second vertex of the tetrahedron
         * \param[in] p2 third vertex of the tetrahedron
         * \param[in] p3 fourth vertex of the tetrahedron
         * \retval POSITIVE if the tetrahedron is oriented positively
         * \retval ZERO if the tetrahedron is flat
         * \retval NEGATIVE if the tetrahedron is oriented negatively
         */
        inline Sign orient_3d_inexact(
            const double* p0, const double* p1,
            const double* p2, const double* p3
        ) {
            double a11 = p1[0] - p0[0] ;
            double a12 = p1[1] - p0[1] ;
            double a13 = p1[2] - p0[2] ;

            double a21 = p2[0] - p0[0] ;
            double a22 = p2[1] - p0[1] ;
            double a23 = p2[2] - p0[2] ;

            double a31 = p3[0] - p0[0] ;
            double a32 = p3[1] - p0[1] ;
            double a33 = p3[2] - p0[2] ;

            double Delta = det3x3(
                a11,a12,a13,
                a21,a22,a23,
                a31,a32,a33
            );

            return geo_sgn(Delta);
        }

        /**
         * \brief Tests whether two 3d points are identical.
         * \retval true if \p p1 and \p p2 have exactly the same coordinates
         * \retval false otherwise
         */
        bool points_are_identical_3d(
            const double* p1,
            const double* p2
        );

        /**
         * \brief Tests whether three 3d points are colinear.
         * \retval true if \p p1, \p p2 and \p p3 are colinear
         * \retval false otherwise
         */
        bool points_are_colinear_3d(
            const double* p1,
            const double* p2,
            const double* p3
        );

    }
}
}
