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
void edge_collapsing_aux(Mesh&                            mesh,
                         const AABBWrapper&               tree,
                         std::vector<std::array<int, 2>>& edges)
{
    auto& tets         = mesh.tets;
    auto& tet_vertices = mesh.tet_vertices;

    int suc_counter = 0;

    ////init
    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_s> ec_queue;
    // A collapse moves the first end onto the second, so an edge goes in once per direction.
    const auto push_both_ends = [&](const std::array<int, 2>& e, Scalar l_2) {
        ec_queue.push(ElementInQueue(e, l_2));
        ec_queue.push(ElementInQueue({{e[1], e[0]}}, l_2));
    };
    for (auto& e : edges) {
        Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
        if (is_collapsable_length(mesh, e[0], e[1], l_2) &&
            is_collapsable_boundary(mesh, e[0], e[1], tree))
            push_both_ends(e, l_2);
    }
    edges.clear();

    ////collapse
    int                             ts = 0;
    std::vector<std::array<int, 2>> inf_es;
    std::vector<int>                inf_e_tss;
    std::vector<int>                tet_tss;
    tet_tss.assign(tets.size(), 0);

    do {
        suc_counter = 0;
        while (!ec_queue.empty()) {
            std::array<int, 2> v_ids      = ec_queue.top().v_ids;
            Scalar             old_weight = ec_queue.top().weight;
            ec_queue.pop();

            while (!ec_queue.empty()) {
                if (ec_queue.top().v_ids == v_ids)
                    ec_queue.pop();
                else
                    break;
            }

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

        ////postprocess
        std::vector<std::array<int, 2>> tmp_inf_es;
        const unsigned int              inf_es_size = inf_es.size();
        tmp_inf_es.reserve(inf_es_size / 4.0 + 1);
        for (unsigned int i = 0; i < inf_es_size; i++) {
            if (is_edge_freezed(mesh, inf_es[i][0], inf_es[i][1]))
                continue;
            if (!is_valid_edge(mesh, inf_es[i][0], inf_es[i][1]))
                continue;

            Scalar weight = get_edge_length_2(mesh, inf_es[i][0], inf_es[i][1]);
            if (!is_collapsable_length(mesh, inf_es[i][0], inf_es[i][1], weight))
                continue;

            if (!is_collapsable_bbox(mesh, inf_es[i][0], inf_es[i][1]))
                continue;

            bool is_recal = false;
            for (int t_id : tet_vertices[inf_es[i][0]].conn_tets) {
                if (tet_tss[t_id] > inf_e_tss[i]) {
                    is_recal = true;
                    break;
                }
            }
            if (is_recal)
                ec_queue.push(ElementInQueue(inf_es[i], weight));
            else
                tmp_inf_es.push_back(inf_es[i]);
        }
        std::sort(tmp_inf_es.begin(), tmp_inf_es.end());
        tmp_inf_es.erase(std::unique(tmp_inf_es.begin(), tmp_inf_es.end()),
                         tmp_inf_es.end());  // it's better
        inf_es = tmp_inf_es;

        ts++;
        inf_e_tss = std::vector<int>(inf_es.size(), ts);
    } while (suc_counter > 0);
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::edge_collapsing(Mesh& mesh, const AABBWrapper& tree)
{
    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);
    edge_collapsing_aux(mesh, tree, edges);
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

    ////check vertices
    // check isolate surface points
    if (tet_vertices[v1_id].is_on_surface && is_isolate_surface_point(mesh, v1_id)) {
        tet_vertices[v1_id].is_on_surface  = false;
        tet_vertices[v1_id].is_on_boundary = false;
    }
    // check boundary/surface
    if (tet_vertices[v1_id].is_on_boundary &&
        is_point_out_boundary_envelope(
          mesh,
          tet_vertices[v2_id].pos,
          tree))  // todo: you should check/unmark is_on_boundary around here
        return false;
    if (tet_vertices[v1_id].is_on_surface &&
        is_point_out_envelope(mesh, tet_vertices[v2_id].pos, tree))
        return false;

    ////check tets
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

    // inversion
    std::vector<int> js_n1_t_ids;
    for (int t_id : n1_t_ids) {
        int j = tets[t_id].find(v1_id);
        js_n1_t_ids.push_back(j);
        assert(j < 4);
        if (is_inverted(mesh, t_id, j, tet_vertices[v2_id].pos))
            return false;
    }

    // quality
    Scalar old_max_quality = 0;
    if (mesh.is_coarsening)
        old_max_quality = mesh.params.stop_energy;
    else if (is_check_quality)
        old_max_quality = get_max_quality(mesh, tet_vertices[v1_id].conn_tets);
    std::vector<Scalar> new_qs;
    new_qs.reserve(tet_vertices[v1_id].conn_tets.size());
    int ii = 0;
    for (int t_id : n1_t_ids) {
        int    j     = js_n1_t_ids[ii++];
        Scalar new_q = get_quality(tet_vertices[v2_id],
                                   tet_vertices[tets[t_id][mod4(j + 1)]],
                                   tet_vertices[tets[t_id][mod4(j + 2)]],
                                   tet_vertices[tets[t_id][mod4(j + 3)]]);
        if (is_check_quality && new_q > old_max_quality)
            return false;
        new_qs.push_back(new_q);
    }

    // envelope
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

    ////real update
    // vertex
    tet_vertices[v1_id].is_removed = true;
    tet_vertices[v2_id].is_on_bbox =
      tet_vertices[v1_id].is_on_bbox || tet_vertices[v2_id].is_on_bbox;
    tet_vertices[v2_id].is_on_surface =
      tet_vertices[v1_id].is_on_surface || tet_vertices[v2_id].is_on_surface;
    tet_vertices[v2_id].is_on_boundary =
      tet_vertices[v1_id].is_on_boundary || tet_vertices[v2_id].is_on_boundary;
    if (tet_vertices[v1_id].on_boundary_e_id >= 0)
        tet_vertices[v2_id].on_boundary_e_id = tet_vertices[v1_id].on_boundary_e_id;

    // tets
    // update quality
    int q_i = 0;
    for (int t_id : n1_t_ids) {
        tets[t_id].quality = new_qs[q_i++];
    }

    // n_v_id for repush
    //    std::set_difference(n1_v_ids.begin(), n1_v_ids.end(), n12_v_ids.begin(), n12_v_ids.end(),
    std::vector<int> n1_v_ids;
    n1_v_ids.reserve(n1_t_ids.size() * 4);
    for (int t_id : n1_t_ids) {
        for (int j = 0; j < 4; j++)
            n1_v_ids.push_back(tets[t_id][j]);
    }
    vector_unique(n1_v_ids);

    // update tags

    for (int t_id : n12_t_ids) {
        // The marks on the two faces of this tet that face v1 and v2, indexed 0 for v1 and 1 for
        // v2. The collapse merges them onto the two faces that connect the two ends.
        std::array<int, 2> j12         = {{-1, -1}};
        std::array<int, 2> sf_facing   = {{NOT_SURFACE, NOT_SURFACE}};
        std::array<int, 2> tag_facing  = {{NO_SURFACE_TAG, NO_SURFACE_TAG}};
        std::array<int, 2> bbox_facing = {{NOT_BBOX, NOT_BBOX}};
        for (int j = 0; j < 4; j++) {
            const int i = tets[t_id][j] == v1_id ? 0 : (tets[t_id][j] == v2_id ? 1 : -1);
            if (i < 0)
                continue;
            sf_facing[i]   = tets[t_id].is_surface_fs[j];
            tag_facing[i]  = tets[t_id].surface_tags[j];
            bbox_facing[i] = tets[t_id].is_bbox_fs[j];
            j12[i]         = j;
        }

        std::array<int, 2> sf_connecting_v12  = {{NOT_SURFACE, NOT_SURFACE}};
        std::array<int, 2> tag_connecting_v12 = {{NO_SURFACE_TAG, NO_SURFACE_TAG}};
        if (sf_facing[1] != NOT_SURFACE && sf_facing[0] != NOT_SURFACE) {
            const int new_tag     = sf_facing[0] - sf_facing[1];
            sf_connecting_v12[0]  = (new_tag > 0) - (new_tag < 0);
            tag_connecting_v12[0] = tag_facing[0];
        }
        else if (sf_facing[1] != NOT_SURFACE) {
            sf_connecting_v12[0]  = -sf_facing[1];
            tag_connecting_v12[0] = tag_facing[1];
        }
        else if (sf_facing[0] != NOT_SURFACE) {
            sf_connecting_v12[0]  = sf_facing[0];
            tag_connecting_v12[0] = tag_facing[0];
        }
        if (sf_connecting_v12[0] != NOT_SURFACE) {
            sf_connecting_v12[1]  = -sf_connecting_v12[0];
            tag_connecting_v12[1] = tag_connecting_v12[0];
        }

        int bbox_connecting_v12 = NOT_BBOX;
        if (bbox_facing[1] != NOT_BBOX)
            bbox_connecting_v12 = bbox_facing[1];
        else if (bbox_facing[0] != NOT_BBOX)
            bbox_connecting_v12 = bbox_facing[0];

        for (int i = 0; i < 2; i++) {
            // Face i is the one facing the other end, so it is the one the collapse keeps.
            const int          jj = j12[(i + 1) % 2];
            std::array<int, 3> f  = {{tets[t_id][mod4(jj + 1)],
                                      tets[t_id][mod4(jj + 2)],
                                      tets[t_id][mod4(jj + 3)]}};

            std::vector<int> pair;
            set_intersection(tet_vertices[f[0]].conn_tets,
                             tet_vertices[f[1]].conn_tets,
                             tet_vertices[f[2]].conn_tets,
                             pair);
            if (pair.size() <= 1)
                continue;

            const int opp_t_id = pair[0] == t_id ? pair[1] : pair[0];
            const int j = get_local_f_id(opp_t_id, f[0], f[1], f[2], mesh);
            tets[opp_t_id].is_surface_fs[j] = sf_connecting_v12[i];
            tets[opp_t_id].surface_tags[j]  = tag_connecting_v12[i];
            tets[opp_t_id].is_bbox_fs[j]    = bbox_connecting_v12;
        }
    }

    // update connectivity
    ts++;
    ii = 0;
    for (int t_id : n1_t_ids) {
        int j         = js_n1_t_ids[ii++];
        tets[t_id][j] = v2_id;
        tet_vertices[v2_id].conn_tets.push_back(t_id);
        if (is_update_tss)
            tet_tss[t_id] = ts;  // update timestamp
    }
    for (int t_id : n12_t_ids) {
        tets[t_id].is_removed = true;
        for (int j = 0; j < 4; j++) {
            if (tets[t_id][j] != v1_id)
                vector_erase(tet_vertices[tets[t_id][j]].conn_tets, t_id);
        }
    }

    tet_vertices[v1_id].conn_tets.clear();

    ////re-push
    for (int v_id : n1_v_ids) {
        if (v_id != v1_id)
            new_edges.push_back({{v2_id, v_id}});
    }

    return true;
}

bool floatTetWild::is_edge_freezed(Mesh& mesh, int v1_id, int v2_id)
{
    return mesh.tet_vertices[v1_id].is_freezed || mesh.tet_vertices[v2_id].is_freezed;
}

bool floatTetWild::is_collapsable_bbox(Mesh& mesh, int v1_id, int v2_id)
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

bool floatTetWild::is_collapsable_length(Mesh& mesh, int v1_id, int v2_id, Scalar l_2)
{
    Scalar sizing_scalar = avg_sizing_scalar(mesh, v1_id, v2_id);
    return l_2 <= mesh.params.collapse_threshold_2 * sizing_scalar * sizing_scalar;
}

bool floatTetWild::is_collapsable_boundary(Mesh&              mesh,
                                           int                v1_id,
                                           int                v2_id,
                                           const AABBWrapper& tree)
{
    return !mesh.tet_vertices[v1_id].is_on_boundary || is_boundary_edge(mesh, v1_id, v2_id, tree);
}
