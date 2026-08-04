// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/EdgeCollapsing.h>
#include <floattetwild/LocalOperations.h>

namespace floatTetWild {
namespace {
bool is_edge_freezed(Mesh& mesh, int v1_id, int v2_id)
{
    return mesh.tet_vertices[v1_id].is_freezed || mesh.tet_vertices[v2_id].is_freezed;
}

bool is_collapsable_bbox(Mesh& mesh, int v1_id, int v2_id)
{
    if (!mesh.tet_vertices[v1_id].is_on_bbox)
        return true;
    else if (!mesh.tet_vertices[v2_id].is_on_bbox)
        return false;

    std::vector<int> bbox_fs2;
    for (int t_id : mesh.tet_vertices[v2_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v2_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX)
                bbox_fs2.push_back(mesh.tets[t_id].is_bbox_fs[j]);
        }
    }
    vector_unique(bbox_fs2);

    for (int t_id : mesh.tet_vertices[v1_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX) {
                if (std::find(bbox_fs2.begin(), bbox_fs2.end(), mesh.tets[t_id].is_bbox_fs[j]) ==
                    bbox_fs2.end())
                    return false;
            }
        }
    }

    return true;
}

bool is_collapsable_length(Mesh& mesh, int v1_id, int v2_id, Scalar l_2)
{
    Scalar sizing_scalar = avg_sizing_scalar(mesh, v1_id, v2_id);
    return l_2 <= mesh.params.collapse_threshold_2 * sizing_scalar * sizing_scalar;
}

bool is_collapsable_boundary(Mesh&              mesh,
                             int                v1_id,
                             int                v2_id,
                             const AABBWrapper& tree)
{
    return !mesh.tet_vertices[v1_id].is_on_boundary || is_boundary_edge(mesh, v1_id, v2_id, tree);
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::edge_collapsing(Mesh& mesh, const AABBWrapper& tree)
{
    auto& tets         = mesh.tets;
    auto& tet_vertices = mesh.tet_vertices;

    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);

    ShortestFirstQueue ec_queue;
    // A collapse moves the first end onto the second, so an edge goes in once per direction.
    const auto push_both_ends = [&](const std::array<int, 2>& e, Scalar l_2) {
        ec_queue.push({e, l_2});
        ec_queue.push({{{e[1], e[0]}}, l_2});
    };
    for (auto& e : edges) {
        Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
        if (is_collapsable_length(mesh, e[0], e[1], l_2) &&
            is_collapsable_boundary(mesh, e[0], e[1], tree))
            push_both_ends(e, l_2);
    }
    edges.clear();

    int                             ts = 0;
    std::vector<std::array<int, 2>> inf_es;
    std::vector<int>                inf_e_tss;
    std::vector<int>                tet_tss;
    tet_tss.assign(tets.size(), 0);

    for (;;) {
        int suc_counter = 0;
        while (!ec_queue.empty()) {
            std::array<int, 2> v_ids      = ec_queue.top().v_ids;
            Scalar             old_weight = ec_queue.top().weight;
            ec_queue.pop();

            pop_duplicates(ec_queue, v_ids);

            if (is_edge_freezed(mesh, v_ids[0], v_ids[1]))
                continue;

            if (!is_valid_edge(mesh, v_ids[0], v_ids[1]))
                continue;

            if (!is_collapsable_boundary(mesh, v_ids[0], v_ids[1], tree))
                continue;

            Scalar weight = get_edge_length_2(mesh, v_ids[0], v_ids[1]);
            if (weight != old_weight || !is_collapsable_length(mesh, v_ids[0], v_ids[1], weight))
                continue;

            if (!is_collapsable_bbox(mesh, v_ids[0], v_ids[1]))
                continue;

            std::vector<std::array<int, 2>> new_edges;
            if (collapse_an_edge(mesh, v_ids[0], v_ids[1], tree, new_edges, ts, tet_tss)) {
                suc_counter++;

                for (auto& e : new_edges) {
                    if (is_edge_freezed(mesh, e[0], e[1]))
                        continue;

                    Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
                    if (is_collapsable_length(mesh, e[0], e[1], l_2))
                        push_both_ends(e, l_2);
                }
            }
            else {
                inf_es.push_back(v_ids);
                inf_e_tss.push_back(ts);
            }
        }

        if (suc_counter == 0)
            break;

        std::vector<std::array<int, 2>> tmp_inf_es;
        const unsigned int              inf_es_size = inf_es.size();
        tmp_inf_es.reserve(inf_es_size / 4.0 + 1);
        for (unsigned int i = 0; i < inf_es_size; i++) {
            const std::array<int, 2>& e = inf_es[i];
            if (is_edge_freezed(mesh, e[0], e[1]))
                continue;
            if (!is_valid_edge(mesh, e[0], e[1]))
                continue;

            Scalar weight = get_edge_length_2(mesh, e[0], e[1]);
            if (!is_collapsable_length(mesh, e[0], e[1], weight))
                continue;

            if (!is_collapsable_bbox(mesh, e[0], e[1]))
                continue;

            // Worth another try only if something around it has moved since it last failed.
            bool is_recal = false;
            for (int t_id : tet_vertices[e[0]].conn_tets) {
                if (tet_tss[t_id] > inf_e_tss[i]) {
                    is_recal = true;
                    break;
                }
            }
            if (is_recal)
                ec_queue.push({e, weight});
            else
                tmp_inf_es.push_back(e);
        }
        vector_unique(tmp_inf_es);
        inf_es = tmp_inf_es;

        ts++;
        inf_e_tss = std::vector<int>(inf_es.size(), ts);
    }
}

// Moves v1 onto v2, or leaves the mesh untouched and returns false when the collapse would
// invert a tet, worsen the worst quality around v1, or take v1's surface or boundary out of its
// envelope.
bool floatTetWild::collapse_an_edge(Mesh&                            mesh,
                                    int                              v1_id,
                                    int                              v2_id,
                                    const AABBWrapper&               tree,
                                    std::vector<std::array<int, 2>>& new_edges,
                                    int                              ts,
                                    std::vector<int>&                tet_tss,
                                    bool                             is_check_quality,
                                    bool                             is_update_tss)
{
    auto& tet_vertices = mesh.tet_vertices;
    auto& tets         = mesh.tets;

    if (tet_vertices[v1_id].is_on_surface && is_isolate_surface_point(mesh, v1_id)) {
        tet_vertices[v1_id].is_on_surface  = false;
        tet_vertices[v1_id].is_on_boundary = false;
    }
    if (tet_vertices[v1_id].is_on_boundary &&
        is_point_out_boundary_envelope(
          mesh,
          tet_vertices[v2_id].pos,
          tree))  // todo: you should check/unmark is_on_boundary around here
        return false;
    if (tet_vertices[v1_id].is_on_surface &&
        tree.is_out_sf_envelope(tet_vertices[v2_id].pos, mesh.params.eps_2))
        return false;

    std::vector<int> n12_t_ids;
    set_intersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, n12_t_ids);
    if (n12_t_ids.empty())
        return false;
    std::vector<int> n1_t_ids;  // v1.conn_tets - n12_t_ids
    std::sort(tet_vertices[v1_id].conn_tets.begin(), tet_vertices[v1_id].conn_tets.end());
    std::sort(n12_t_ids.begin(), n12_t_ids.end());
    std::set_difference(tet_vertices[v1_id].conn_tets.begin(),
                        tet_vertices[v1_id].conn_tets.end(),
                        n12_t_ids.begin(),
                        n12_t_ids.end(),
                        std::back_inserter(n1_t_ids));

    std::vector<int> js_n1_t_ids;
    for (int t_id : n1_t_ids) {
        int j = tets[t_id].find(v1_id);
        js_n1_t_ids.push_back(j);
        assert(j < 4);
        if (is_inverted(mesh, t_id, j, tet_vertices[v2_id].pos))
            return false;
    }

    Scalar old_max_quality = 0;
    if (mesh.is_coarsening)
        old_max_quality = mesh.params.stop_energy;
    else if (is_check_quality)
        old_max_quality = get_max_quality(mesh, tet_vertices[v1_id].conn_tets);
    std::vector<Scalar> new_qs;
    new_qs.reserve(tet_vertices[v1_id].conn_tets.size());
    for (size_t i = 0; i < n1_t_ids.size(); i++) {
        const int t_id = n1_t_ids[i];
        const int j    = js_n1_t_ids[i];
        Scalar new_q = get_quality(tet_vertices[v2_id].pos,
                                   tet_pos(mesh, t_id, (j + 1) % 4),
                                   tet_pos(mesh, t_id, (j + 2) % 4),
                                   tet_pos(mesh, t_id, (j + 3) % 4));
        if (is_check_quality && new_q > old_max_quality)
            return false;
        new_qs.push_back(new_q);
    }

    Scalar l = get_edge_length_2(mesh, v1_id, v2_id);
    if (l > 0) {
        if (tet_vertices[v1_id].is_on_boundary) {
            if (is_out_boundary_envelope(mesh, v1_id, tet_vertices[v2_id].pos, tree))
                return false;
        }
        if (tet_vertices[v1_id].is_on_surface) {
            if (is_out_envelope(mesh, v1_id, tet_vertices[v2_id].pos, tree))
                return false;
        }
    }

    tet_vertices[v1_id].is_removed = true;
    tet_vertices[v2_id].is_on_bbox =
      tet_vertices[v1_id].is_on_bbox || tet_vertices[v2_id].is_on_bbox;
    tet_vertices[v2_id].is_on_surface =
      tet_vertices[v1_id].is_on_surface || tet_vertices[v2_id].is_on_surface;
    tet_vertices[v2_id].is_on_boundary =
      tet_vertices[v1_id].is_on_boundary || tet_vertices[v2_id].is_on_boundary;

    std::vector<int> n1_v_ids;
    collect_tet_vertices(mesh, n1_t_ids, n1_v_ids);

    // The marks on the face of a tet that faces one end of the edge.
    struct FaceMarks
    {
        int j    = -1;
        int sf   = NOT_SURFACE;
        int tag  = NO_SURFACE_TAG;
        int bbox = NOT_BBOX;
    };

    for (int t_id : n12_t_ids) {
        const MeshTet& t = tets[t_id];

        // facing[0] is the face across from v1 and facing[1] the one across from v2. The collapse
        // merges the pair onto the two faces that connect the two ends.
        std::array<FaceMarks, 2> facing;
        for (int j = 0; j < 4; j++) {
            const int i = t[j] == v1_id ? 0 : (t[j] == v2_id ? 1 : -1);
            if (i < 0)
                continue;
            facing[i] = {j, t.is_surface_fs[j], t.surface_tags[j], t.is_bbox_fs[j]};
        }

        std::array<int, 2> sf_connecting_v12  = {{NOT_SURFACE, NOT_SURFACE}};
        std::array<int, 2> tag_connecting_v12 = {{NO_SURFACE_TAG, NO_SURFACE_TAG}};
        if (facing[1].sf != NOT_SURFACE && facing[0].sf != NOT_SURFACE) {
            const int new_tag     = facing[0].sf - facing[1].sf;
            sf_connecting_v12[0]  = (new_tag > 0) - (new_tag < 0);
            tag_connecting_v12[0] = facing[0].tag;
        }
        else if (facing[1].sf != NOT_SURFACE) {
            sf_connecting_v12[0]  = -facing[1].sf;
            tag_connecting_v12[0] = facing[1].tag;
        }
        else if (facing[0].sf != NOT_SURFACE) {
            sf_connecting_v12[0]  = facing[0].sf;
            tag_connecting_v12[0] = facing[0].tag;
        }
        if (sf_connecting_v12[0] != NOT_SURFACE) {
            sf_connecting_v12[1]  = -sf_connecting_v12[0];
            tag_connecting_v12[1] = tag_connecting_v12[0];
        }

        const int bbox_connecting_v12 =
          facing[1].bbox != NOT_BBOX ? facing[1].bbox : facing[0].bbox;

        for (int i = 0; i < 2; i++) {
            // Face i is the one facing the other end, so it is the one the collapse keeps, and
            // the marks go onto its far side.
            int opp_t_id, j;
            if (!get_opp_face(mesh, t_id, facing[(i + 1) % 2].j, opp_t_id, j))
                continue;

            tets[opp_t_id].is_surface_fs[j] = sf_connecting_v12[i];
            tets[opp_t_id].surface_tags[j]  = tag_connecting_v12[i];
            tets[opp_t_id].is_bbox_fs[j]    = bbox_connecting_v12;
        }
    }

    ts++;
    for (size_t i = 0; i < n1_t_ids.size(); i++) {
        const int t_id     = n1_t_ids[i];
        tets[t_id][js_n1_t_ids[i]] = v2_id;
        tets[t_id].quality = new_qs[i];
        tet_vertices[v2_id].conn_tets.push_back(t_id);
        if (is_update_tss)
            tet_tss[t_id] = ts;
    }
    for (int t_id : n12_t_ids) {
        tets[t_id].is_removed = true;
        for (int j = 0; j < 4; j++) {
            if (tets[t_id][j] != v1_id)
                vector_erase(tet_vertices[tets[t_id][j]].conn_tets, t_id);
        }
    }

    tet_vertices[v1_id].conn_tets.clear();

    for (int v_id : n1_v_ids) {
        if (v_id != v1_id)
            new_edges.push_back({{v2_id, v_id}});
    }

    return true;
}
