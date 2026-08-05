// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "CutMesh.h"
#include "LocalOperations.h"
#include "Predicates.hpp"
#include "intersections.h"
#include "Logger.hpp"

void floatTetWild::CutMesh::construct(const std::vector<int>& cut_t_ids) {
    collect_tet_vertices(mesh, cut_t_ids, v_ids);

    for (int i = 0; i < v_ids.size(); i++)
        map_v_ids[v_ids[i]] = i;

    tets.resize(cut_t_ids.size());
    for (int i = 0; i < cut_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++)
            tets[i][j] = map_v_ids[mesh.tets[cut_t_ids[i]][j]];
    }
}

floatTetWild::Scalar floatTetWild::CutMesh::get_signed_plane_dist(const Vector3 &p, bool &snaps) const {
    const int ori = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2], p);
    if (ori == Predicates::ORI_ZERO) {
        snaps = false;
        return 0;
    }

    Scalar dist = get_to_plane_dist(p);
    if ((ori == Predicates::ORI_POSITIVE && dist > 0) || (ori == Predicates::ORI_NEGATIVE && dist < 0))
        dist = -dist;
    snaps = std::fabs(dist) < mesh.params.eps_coplanar;
    return dist;
}

bool floatTetWild::CutMesh::snap_to_plane() {
    bool snapped = false;
    to_plane_dists.resize(map_v_ids.size());
    is_snapped.resize(map_v_ids.size());
    for (auto &v:map_v_ids) {
        bool is_v_snapped = false;
        to_plane_dists[v.second] = get_signed_plane_dist(mesh.tet_vertices[v.first].pos, is_v_snapped);
        if (is_v_snapped) {
            is_snapped[v.second] = true;
            snapped = true;
        }
    }

    revert_totally_snapped_tets();

    return snapped;
}

void floatTetWild::CutMesh::expand_new(std::vector<int> &cut_t_ids) {
    const int t = get_t(p_vs[0], p_vs[1], p_vs[2]);
    const std::array<Vector2, 3> tri_2d = {{to_2d(p_vs[0], t), to_2d(p_vs[1], t), to_2d(p_vs[2], t)}};

    std::vector<bool> is_in_cutmesh(mesh.tets.size(), false);
    for (int t_id:cut_t_ids)
        is_in_cutmesh[t_id] = true;

    std::vector<bool> is_interior(v_ids.size(), false);

    // The local id of a mesh vertex the cut has reached, taking it on with its plane distance the
    // first time it is asked for.
    const auto local_id = [&](int gv_id) {
        const auto it = map_v_ids.find(gv_id);
        if (it != map_v_ids.end())
            return it->second;
        v_ids.push_back(gv_id);
        is_interior.push_back(false);
        bool is_v_snapped = false;
        to_plane_dists.push_back(get_signed_plane_dist(mesh.tet_vertices[gv_id].pos, is_v_snapped));
        is_snapped.push_back(is_v_snapped);
        is_projected.push_back(false);
        return map_v_ids[gv_id] = int(v_ids.size()) - 1;
    };

    while (true) {
        std::vector<bool> is_visited(mesh.tets.size(), false);
        for (int t_id:cut_t_ids)
            is_visited[t_id] = true;

        int old_cut_t_ids = cut_t_ids.size();
        for (const auto &m: map_v_ids) {
            int gv_id = m.first;
            int lv_id = m.second;

            if (is_interior[lv_id])
                continue;
            if (!is_snapped[lv_id])
                continue;

            bool is_in = true;
            for (int gt_id: mesh.tet_vertices[gv_id].conn_tets) {
                if (is_in_cutmesh[gt_id])
                    continue;
                is_in = false;

                if (is_visited[gt_id])
                    continue;
                is_visited[gt_id] = true;

                int cnt = 0;
                for (int j = 0; j < 4; j++)
                    cnt += map_v_ids.count(mesh.tets[gt_id][j]);
                if (cnt < 3)
                    continue;
                int oris[4];
                for (int j = 0; j < 4; j++)
                    oris[j] = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2],
                                                    tet_pos(mesh, gt_id, j));
                int cnt_pos, cnt_neg;
                count_orientations(oris, 4, cnt_pos, cnt_neg);
                if (cnt_neg == 0 || cnt_pos == 0)  // does not straddle the plane
                    continue;

                bool is_overlapped = false;
                std::array<Vector2, 4> tet_2d;
                for (int j = 0; j < 4; j++)
                    tet_2d[j] = to_2d(tet_pos(mesh, gt_id, j), p_n, p_vs[0], t);
                for(int j=0;j<4;j++) {
                    if (is_tri_tri_cutted_2d({{tet_2d[(j + 1) % 4], tet_2d[(j + 2) % 4], tet_2d[(j + 3) % 4]}},
                                             tri_2d)) {
                        is_overlapped = true;
                        break;
                    }
                }
                if(!is_overlapped)
                    continue;

                cut_t_ids.push_back(gt_id);
                is_in_cutmesh[gt_id] = true;

                tets.emplace_back();
                for (int j = 0; j < 4; j++)
                    tets.back()[j] = local_id(mesh.tets[gt_id][j]);
            }
            if (is_in)
                is_interior[lv_id] = true;
        }
        if (cut_t_ids.size() == old_cut_t_ids)
            break;
    }
    revert_totally_snapped_tets();
}

void floatTetWild::CutMesh::project_to_plane(int input_vertices_size) {
    is_projected.resize(v_ids.size(), false);

    for (int i = 0; i < is_snapped.size(); i++) {
        if (!is_snapped[i] || is_projected[i])
            continue;
        if (v_ids[i] < input_vertices_size)
            continue;
        Scalar dist = get_to_plane_dist(mesh.tet_vertices[v_ids[i]].pos);
        Vector3 proj_p = mesh.tet_vertices[v_ids[i]].pos - p_n * dist;
        bool is_snappable = true;
        for (int t_id: mesh.tet_vertices[v_ids[i]].conn_tets) {
            int j = mesh.tets[t_id].find(v_ids[i]);
            if (is_inverted(mesh, t_id, j, proj_p)) {
                is_snappable = false;
                break;
            }
        }
        if (is_snappable) {
            mesh.tet_vertices[v_ids[i]].pos = proj_p;
            is_projected[i] = true;
            to_plane_dists[i] = get_to_plane_dist(proj_p);
        }
    }
}

// A tet with all four corners on the plane would be cut into nothing, so the corner furthest from
// the plane gives up its snap.
void floatTetWild::CutMesh::revert_totally_snapped_tets() {
    for (const auto &t : tets) {
        if (!is_v_on_plane(t[0]) || !is_v_on_plane(t[1]) || !is_v_on_plane(t[2]) ||
            !is_v_on_plane(t[3]))
            continue;

        auto tmp_t = t;
        std::sort(tmp_t.begin(), tmp_t.end(), [&](int a, int b) {
            return fabs(to_plane_dists[a]) > fabs(to_plane_dists[b]);
        });
        for (int j = 0; j < 3; j++) {
            if (is_snapped[tmp_t[j]]) {
                is_snapped[tmp_t[j]] = false;
                break;
            }
        }
    }
}

bool floatTetWild::CutMesh::get_intersecting_edges_and_points(std::vector<Vector3> &points,
                                                              std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
                                                              std::vector<int>& subdivide_t_ids) {
    std::vector<std::array<int, 2>> edges;
    collect_tet_edges(tets, edges);

    std::vector<int> e_v_ids;
    for (int i = 0; i < edges.size(); i++) {
        const auto &e = edges[i];
        if (is_v_on_plane(e[0]) || is_v_on_plane(e[1]))
            continue;
        if (to_plane_dists[e[0]] * to_plane_dists[e[1]] >= 0)
            continue;

        int v1_id = v_ids[e[0]];
        int v2_id = v_ids[e[1]];
        Vector3 p;
        bool is_result = seg_plane_intersection(mesh.tet_vertices[v1_id].pos, mesh.tet_vertices[v2_id].pos,
                                                p_vs[0], p_n, p);
        if (!is_result) {
            logger().error("seg_plane_intersection no result!");
            return false;
        }

        points.push_back(p);
        map_edge_to_intersecting_point[sorted_edge(v1_id, v2_id)] = points.size() - 1;

        e_v_ids.push_back(v1_id);
        e_v_ids.push_back(v2_id);
    }
    vector_unique(e_v_ids);
    for (int v_id: e_v_ids)
        subdivide_t_ids.insert(subdivide_t_ids.end(), mesh.tet_vertices[v_id].conn_tets.begin(),
                               mesh.tet_vertices[v_id].conn_tets.end());
    vector_unique(subdivide_t_ids);

    return true;
}

