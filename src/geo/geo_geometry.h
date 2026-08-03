// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/basic/geometry.h
// Copied rather than reimplemented: vec3 and Box, unchanged. The Geom helpers that lived here as
// well, tetra_signed_volume in both its shapes, had no caller and are gone; what is left of Geom
// is in geo_geometry_nd.h.

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

    /*******************************************************************/

    /**
     * \brief Axis-aligned bounding box.
     */
    class Box {
    public:
        double xyz_min[3];
        double xyz_max[3];
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
