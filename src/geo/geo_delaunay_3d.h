// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/delaunay/delaunay_3d.h
// Copied rather than reimplemented: create_first_tetrahedron scans raw indices and SOS breaks ties by point address, so insertion order decides the mesh

#pragma once

#include "geo_basic.h"

#include <vector>

namespace floatTetWild {
namespace geo {

    // The Delaunay tetrahedralization of the given points, three doubles each, as four vertex
    // indices per tetrahedron. Empty when the points are all coplanar or all colinear. Degenerate
    // input aside, the result covers the convex hull of the points.
    std::vector<index_t> delaunay_3d(index_t nb_vertices, const double* vertices);

} }
