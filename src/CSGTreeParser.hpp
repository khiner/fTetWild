// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include "Types.hpp"

#include "geo/geo_mesh.h"

#include <istream>
#include <memory>
#include <string>
#include <vector>

namespace floatTetWild {

// A csg tree is a binary tree of boolean operations over surface meshes, read from a file like
//
//   {"operation": "difference",
//    "left": {"operation": "union", "left": "a.stl", "right": "b.stl"},
//    "right": "c.stl"}
//
// where a child is either the name of a mesh file or another node. assign_mesh_ids() then gives
// every distinct name an index into the mesh list, which is the form the rest of the code
// consumes.
struct CSGTree
{
    struct Child
    {
        std::string              mesh;     // the mesh file name, when this child is not a subtree
        int                      id = -1;  // its index into the mesh list, after assign_mesh_ids()
        std::shared_ptr<CSGTree> subtree;  // set when this child is another node instead
    };

    std::string operation;  // "union", "intersection" or "difference"
    Child       left, right;
};

namespace CSGTreeParser {

// Reads a tree, leaving the caller to open the file. Returns false with a description in
// error rather than throwing.
bool parse(std::istream& in, CSGTree& tree, std::string& error);

// Numbers every mesh in the tree and returns the names, so that meshes[child.id] is the name the
// child was parsed with. Ids run dense from 0 in the order the names are first seen walking the
// tree left to right, and a name used twice keeps the id it was given first.
std::vector<std::string> assign_mesh_ids(CSGTree& tree);

bool keep_tet(const CSGTree& csg_tree, int t_id, const std::vector<MatrixXd>& w);

void merge(const std::vector<std::vector<Vector3>>&  Vs,
           const std::vector<std::vector<Vector3i>>& Fs,
           std::vector<Vector3>&                     V,
           std::vector<Vector3i>&                    F,
           geo::Mesh&                                sf_mesh,
           std::vector<int>&                         tags);

}  // namespace CSGTreeParser
}  // namespace floatTetWild
