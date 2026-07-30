// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/Types.hpp>

#include <vector>

namespace floatTetWild {

// Run the meshing pipeline on an input surface: preprocessing, tetrahedralization, cutting,
// optimization, and tracked surface orientation. Reads and writes no files.
//
// Settings come from mesh.params, which the caller fills in before calling, and which this
// completes via Parameters::init. input_vertices and input_faces are modified in place by the
// preprocessing step, and input_tags alongside them.
//
// Interior extraction is deliberately left to the caller, which is why this stops after
// orientation. Follow it with one of filter_outside, filter_outside_floodfill or
// boolean_operation, so that the mesh can be observed in between. See main.cpp.
//
// Returns 0 on success.
int tetrahedralization(AABBWrapper&           tree,
                       std::vector<Vector3>&  input_vertices,
                       std::vector<Vector3i>& input_faces,
                       std::vector<int>&      input_tags,
                       Mesh&                  mesh,
                       bool                   skip_simplify = false);

}
