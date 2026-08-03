// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Types.hpp>

#include <floattetwild/geo_mesh.h>

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
// where a child is either the name of a mesh file or another node. get_meshes() collects the names
// and returns the same tree with every name replaced by an index into that list, which is the form
// the rest of the code consumes.
struct CSGTree
{
    struct Child
    {
        enum Kind
        {
            Mesh,    // a mesh file name, as it appears in the file
            Id,      // an index into the mesh list, after get_meshes()
            Subtree  // another node
        };

        Kind                     kind = Mesh;
        std::string              mesh;
        int                      id = -1;
        std::shared_ptr<CSGTree> subtree;
    };

    std::string operation;  // "union", "intersection" or "difference"
    Child       left, right;
};

namespace CSGTreeParser {

// Reads a tree, leaving the caller to open the file. Returns false with a description in
// error rather than throwing.
bool parse(std::istream& in, CSGTree& tree, std::string& error);

// The mesh names in the tree, in the order they are first seen, and the same tree with each name
// replaced by its index into that list. A name used twice keeps the index it was given first.
void get_meshes(const CSGTree&            csg_tree,
                std::vector<std::string>& meshes,
                CSGTree&                  csg_tree_with_ids);

// The largest mesh index in the tree, or -1 when it holds none.
int get_max_id(const CSGTree& csg_tree_with_ids);

bool keep_tet(const CSGTree& csg_tree_with_ids, const int t_id, const std::vector<MatrixXd>& w);

void merge(const std::vector<std::vector<Vector3>>&  Vs,
           const std::vector<std::vector<Vector3i>>& Fs,
           std::vector<Vector3>&                     V,
           std::vector<Vector3i>&                    F,
           geo::Mesh&                                sf_mesh,
           std::vector<int>&                         tags);

}  // namespace CSGTreeParser
}  // namespace floatTetWild
