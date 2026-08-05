// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LocalOperations.h"

namespace floatTetWild {
namespace {

// Parked in MeshTet::scalar to mark a tet whose quality the pass still owes an update, and cleared
// again once it has. Nothing else reads scalar between the two.
constexpr Scalar TET_MODIFIED = 100;

bool split_an_edge(Mesh& mesh, int v1_id, int v2_id, bool is_repush,
        std::vector<std::array<int, 2>>& new_edges, std::vector<bool>& is_splittable, const AABBWrapper& tree) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    MeshVertex new_v;
    new_v.pos = (tet_vertices[v1_id].pos + tet_vertices[v2_id].pos) / 2;
    const int v_id = get_new_vertex_slot(mesh, new_v);

    std::vector<int> old_t_ids;
    set_intersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, old_t_ids);
    for (int t_id: old_t_ids) {
        if(!is_splittable[t_id]){
            tet_vertices[v_id].is_removed = true;
            return false;
        }
        for (int j = 0; j < 4; j++) {
            if (tets[t_id][j] == v1_id || tets[t_id][j] == v2_id) {
                if(is_inverted(mesh, t_id, j, new_v.pos)) {
                    for (int t_id1: old_t_ids)
                        is_splittable[t_id1] = false;
                    tet_vertices[v_id].is_removed = true;
                    return false;
                }
            }
        }
    }

    tet_vertices[v_id].sizing_scalar = (tet_vertices[v1_id].sizing_scalar + tet_vertices[v2_id].sizing_scalar) / 2;
    tet_vertices[v_id].is_on_bbox = is_bbox_edge(mesh, v1_id, v2_id, old_t_ids);
    tet_vertices[v_id].is_on_surface = is_surface_edge(mesh, v1_id, v2_id, old_t_ids);
    tet_vertices[v_id].is_on_boundary = is_boundary_edge(mesh, v1_id, v2_id, tree);
    if(!mesh.is_input_all_inserted && tet_vertices[v_id].is_on_boundary) {
        tet_vertices[v_id].is_on_cut = (tet_vertices[v1_id].is_on_cut && tet_vertices[v2_id].is_on_cut);
    }

    std::vector<int> new_t_ids;
    split_tets_at(mesh, v_id, v1_id, v2_id, old_t_ids, new_t_ids);

    // The face opposite the end each half gave away is new, so its marks go, and both halves are
    // marked for a quality update.
    const auto finish = [&](int t_id, int dropped_v_id) {
        clear_face_marks(tets[t_id], tets[t_id].find(dropped_v_id));
        tets[t_id].scalar = TET_MODIFIED;
    };
    for (int i = 0; i < old_t_ids.size(); i++) {
        finish(old_t_ids[i], v2_id);
        finish(new_t_ids[i], v1_id);
    }

    if(mesh.tets.size()!=is_splittable.size())
        is_splittable.resize(mesh.tets.size(), true);

    if (!is_repush)
        return true;

    std::vector<int> n_v_ids = {v1_id, v2_id};
    n_v_ids.reserve(old_t_ids.size()*3);
    for(int t_id: old_t_ids){
        for(int j=0;j<4;j++)
            if(tets[t_id][j]!=v_id)
                n_v_ids.push_back(tets[t_id][j]);
    }
    vector_unique(n_v_ids);

    for (int n_v_id: n_v_ids) {
        new_edges.push_back({{v_id, n_v_id}});
    }

    return true;
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::edge_splitting(Mesh& mesh, const AABBWrapper& tree) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    mesh.reset_t_empty_start();
    mesh.reset_v_empty_start();

    LongestFirstQueue es_queue;

    // An edge is worth splitting once it is longer than the threshold scaled by the sizing field
    // averaged over its two ends.
    const auto push_if_long = [&](const std::array<int, 2>& e) {
        Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
        Scalar sizing_scalar = avg_sizing_scalar(mesh, e[0], e[1]);
        if (l_2 > mesh.params.split_threshold_2 * sizing_scalar * sizing_scalar)
            es_queue.push({e, l_2});
    };

    for (auto& e : get_all_edges(mesh))
        push_if_long(e);

    // Every split adds a vertex and turns each incident tet into two, so reserve for the queue as
    // it stands and leave the empty slots already in the mesh to be reused.
    const int v_slots = int(tet_vertices.size()) - mesh.get_v_num();
    const int t_slots = int(tets.size()) - mesh.get_t_num();
    if (v_slots < es_queue.size() * 2)
        tet_vertices.reserve(tet_vertices.size() + es_queue.size() * 2 - v_slots);
    if (t_slots < es_queue.size() * 6 * 2)
        tets.reserve(tets.size() + es_queue.size() * 6 * 2 - t_slots + 1);

    std::vector<bool> is_splittable(mesh.tets.size(), true);
    bool is_repush = true;
    while (!es_queue.empty()) {
        std::array<int, 2> v_ids = es_queue.top().v_ids;
        es_queue.pop();

        if(tet_vertices[v_ids[0]].is_freezed && tet_vertices[v_ids[1]].is_freezed)
            continue;

        std::vector<std::array<int, 2>> new_edges;
        if (!split_an_edge(mesh, v_ids[0], v_ids[1], is_repush, new_edges, is_splittable, tree))
            is_repush = false;

        for (auto &e:new_edges)
            push_if_long(e);
    }

    for(auto& t: tets){
        if(t.is_removed)
            continue;
        if(t.scalar == TET_MODIFIED) {
            t.quality = get_quality(mesh, t);
            t.scalar = 0;
        }
    }
}
