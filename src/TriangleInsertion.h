// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_TRIANGLEINSERTION_H
#define FLOATTETWILD_TRIANGLEINSERTION_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/CutMesh.h>

namespace floatTetWild {
    void insert_triangles(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                          const std::vector<int> &input_tags, Mesh &mesh,
                          std::vector<bool> &is_face_inserted, AABBWrapper &tree, bool is_again);

    // The edges of the input surface that the mesher has to preserve: its open boundary, plus the
    // creases between faces that are not coplanar. b_edge_infos and is_on_cut_edges describe the
    // subset that still needs inserting and are left empty when every face is already in.
    void find_boundary_edges(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                             const std::vector<bool> &is_face_inserted, const std::vector<bool> &old_is_face_inserted,
                             std::vector<std::pair<std::array<int, 2>, std::vector<int>>> &b_edge_infos,
                             std::vector<bool> &is_on_cut_edges,
                             std::vector<std::array<int, 2>> &b_edges);
}

#endif //FLOATTETWILD_TRIANGLEINSERTION_H
