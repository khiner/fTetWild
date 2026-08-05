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

    // The Predicate Construction Kit: geometric predicates built from arithmetic filters (Meyer
    // and Pion), expansion arithmetic (Shewchuk) and simulation of simplicity (Edelsbrunner).
    namespace PCK {

        // POSITIVE when \p p4 is inside the sphere circumscribed to the tetrahedron
        // \p p0, \p p1, \p p2, \p p3, NEGATIVE when it is outside, inverted when the
        // tetrahedron is oriented negatively. Symbolic perturbation breaks ties whenever the five
        // points are cospherical.
        Sign in_sphere_3d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4
        );

        // The sign of the signed volume of the tetrahedron \p p0, \p p1, \p p2, \p p3, so ZERO
        // when it is flat.
        Sign orient_3d(
            const double* p0, const double* p1,
            const double* p2, const double* p3
        );

        // The same, from the plain floating point determinant.
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

        // Whether \p p1 and \p p2 have exactly the same coordinates.
        bool points_are_identical_3d(
            const double* p1,
            const double* p2
        );

        bool points_are_colinear_3d(
            const double* p1,
            const double* p2,
            const double* p3
        );

    }
}
}
