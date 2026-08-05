// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/numerics/predicates.h
// Copied rather than reimplemented: the symbolic-perturbation tie-breaking is what decides
// Delaunay's output on the cospherical configurations a regular background grid is full of.
//
// Trimmed to the predicates fTetWild reaches. The declarations are unchanged.

#pragma once

#include "geo_basic.h"

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

        bool points_are_colinear_3d(
            const double* p1,
            const double* p2,
            const double* p3
        );

    }
}
}
