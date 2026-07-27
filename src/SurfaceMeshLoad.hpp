// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/geo_mesh.h>

#include <string>

namespace floatTetWild {

// Read a surface mesh in one of the formats -i advertises, chosen by extension.
//
// This replaces GEO::mesh_load. geogram's readers for these formats merge nothing and sort
// nothing: they append each vertex as it appears in the file and each facet in file order, and
// mesh_load's only post-pass fills facet adjacency, which nothing downstream of here reads.
// Reproducing the file order is therefore all that is required for the mesh handed to
// mesh_reorder to be identical, and the corpus checks that for .stl. Duplicate vertices are
// expected here and are resolved later by remove_duplicates() in Simplification.cpp.
//
// Polygons with more than three vertices are fanned into triangles as they are read, so the
// result is always simplicial. geogram left them as polygons and fTetWild then called
// GEO::mesh_repair(TRIANGULATE|QUIET), which also dropped degenerate facets, reoriented facets
// anti-Moebius, split non-manifold vertices and oriented normals. Those passes are gone: no
// input the corpus or the tests exercise is polygonal, and none of them ran for triangle input
// in the first place.
bool load_surface_mesh(const std::string& path, geo::Mesh& mesh);
}
