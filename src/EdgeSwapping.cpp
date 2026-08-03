// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/EdgeSwapping.h>
#include <floattetwild/LocalOperations.h>

namespace floatTetWild {
namespace {

// The surface, tag and bbox marks of the faces the swap keeps, keyed by their sorted corner triple
// so that the rebuilt tets can pick their own faces' marks back up.
struct FaceMarks
{
    std::vector<std::array<int, 3>> fs;
    std::vector<char>               is_sf, tags, is_bx;

    // Records the faces of t that contain v1_id or v2_id: those are the ones a swap preserves.
    void record(const MeshTet& t, int v1_id, int v2_id)
    {
        for (int j = 0; j < 4; j++) {
            if (t[j] != v1_id && t[j] != v2_id)
                continue;
            fs.push_back(sorted_face(t, j));
            is_sf.push_back(t.is_surface_fs[j]);
            tags.push_back(t.surface_tags[j]);
            is_bx.push_back(t.is_bbox_fs[j]);
        }
    }

    // Copies the recorded marks onto face j of t, or clears them when that face was not recorded.
    void apply(MeshTet& t, int j) const
    {
        const auto it = std::find(fs.begin(), fs.end(), sorted_face(t, j));
        if (it == fs.end()) {
            clear_face_marks(t, j);
            return;
        }
        const size_t k     = it - fs.begin();
        t.is_surface_fs[j] = is_sf[k];
        t.surface_tags[j]  = tags[k];
        t.is_bbox_fs[j]    = is_bx[k];
    }

    // Puts the marks back on every face of t. kept(v) says whether the face opposite corner v
    // survived the swap. A face that did not is new, so it is cleared instead.
    template<typename Kept>
    void reapply(MeshTet& t, Kept&& kept) const
    {
        for (int j = 0; j < 4; j++) {
            if (kept(t[j]))
                apply(t, j);
            else
                clear_face_marks(t, j);
        }
    }
};

// The two corners of t that are not on the edge, in local order.
std::array<int, 2> corners_off_edge(const MeshTet& t, int v1_id, int v2_id)
{
    std::array<int, 2> off;
    int                cnt = 0;
    for (int j = 0; j < 4; j++) {
        if (t[j] != v1_id && t[j] != v2_id)
            off[cnt++] = t[j];
    }
    return off;
}

// Walks the ring of tets around the edge, so that n12_v_ids[i] and n12_v_ids[i+1] are the two
// corners of n12_t_ids[i] that are not on the edge. The 4-4 and 5-6 swaps both need this.
void order_edge_ring(const Mesh&             mesh,
                     int                     v1_id,
                     int                     v2_id,
                     const std::vector<int>& old_t_ids,
                     std::vector<int>&       n12_v_ids,
                     std::vector<int>&       n12_t_ids)
{
    const int N = old_t_ids.size();

    // {the two corners off the edge, the tet}
    std::vector<std::array<int, 3>> n12_es;
    n12_es.reserve(N);
    for (int t_id : old_t_ids) {
        const std::array<int, 2> off = corners_off_edge(mesh.tets[t_id], v1_id, v2_id);
        n12_es.push_back({{off[0], off[1], t_id}});
    }

    n12_v_ids.push_back(n12_es[0][0]);
    n12_v_ids.push_back(n12_es[0][1]);
    n12_t_ids.push_back(n12_es[0][2]);
    std::vector<bool> is_visited(N, false);
    is_visited[0] = true;
    for (int i = 0; i < N - 2; i++) {
        for (int j = 0; j < N; j++) {
            if (is_visited[j])
                continue;
            if (n12_es[j][0] == n12_v_ids.back()) {
                is_visited[j] = true;
                n12_v_ids.push_back(n12_es[j][1]);
            }
            // else if, not a second if: the first arm moves back(), and the tet would
            // then match on its other corner too.
            else if (n12_es[j][1] == n12_v_ids.back()) {
                is_visited[j] = true;
                n12_v_ids.push_back(n12_es[j][0]);
            }
            if (is_visited[j]) {
                n12_t_ids.push_back(n12_es[j][2]);
                break;
            }
        }
    }
    n12_t_ids.push_back(
      n12_es[std::find(is_visited.begin(), is_visited.end(), false) - is_visited.begin()][2]);
}

bool remove_an_edge_32(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges){
    auto& tet_vertices = mesh.tet_vertices;
    auto& tets = mesh.tets;

    const std::array<int, 2> v_ids = corners_off_edge(tets[old_t_ids[0]], v1_id, v2_id);
    // The tet that carries v_ids[0] comes first, so the two rebuilt tets take the ends in a
    // fixed order.
    const bool in_order = tets[old_t_ids[1]].find(v_ids[0]) >= 0;
    const std::array<int, 2> t_ids = in_order ? std::array<int, 2>{{old_t_ids[1], old_t_ids[2]}}
                                              : std::array<int, 2>{{old_t_ids[2], old_t_ids[1]}};
    std::array<MeshTet, 2> new_tets = {{tets[t_ids[0]], tets[t_ids[1]]}};
    new_tets[0][new_tets[0].find(v1_id)] = v_ids[1];
    new_tets[1][new_tets[1].find(v2_id)] = v_ids[0];

    for(auto& t:new_tets) {
        if (is_inverted(mesh, t))
            return false;
    }
    std::array<Scalar, 2> new_qs;
    Scalar old_max_quality = get_max_quality(mesh, old_t_ids);
    for (int i = 0; i < 2; i++) {
        new_qs[i] = get_quality(mesh, new_tets[i]);
        if (new_qs[i] >= old_max_quality)
            return false;
    }

    FaceMarks marks;
    for (int t_id : old_t_ids)
        marks.record(tets[t_id], v1_id, v2_id);

    tets[old_t_ids[0]].is_removed = true;
    tets[t_ids[0]] = new_tets[0];
    tets[t_ids[1]] = new_tets[1];

    for (int i = 0; i < 2; i++)
        tets[t_ids[i]].quality = new_qs[i];

    marks.reapply(tets[t_ids[0]], [&](int v_id) { return v_id != v2_id; });
    marks.reapply(tets[t_ids[1]], [&](int v_id) { return v_id != v1_id; });

    vector_erase(tet_vertices[v_ids[0]].conn_tets, old_t_ids[0]);
    vector_erase(tet_vertices[v_ids[1]].conn_tets, old_t_ids[0]);

    tet_vertices[v_ids[0]].conn_tets.push_back(t_ids[1]);
    tet_vertices[v_ids[1]].conn_tets.push_back(t_ids[0]);

    vector_erase(tet_vertices[v1_id].conn_tets, old_t_ids[0]);
    vector_erase(tet_vertices[v2_id].conn_tets, old_t_ids[0]);

    vector_erase(tet_vertices[v1_id].conn_tets, t_ids[0]);
    vector_erase(tet_vertices[v2_id].conn_tets, t_ids[1]);

    collect_tet_edges(new_tets, new_edges);

    return true;
}

bool remove_an_edge_44(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    order_edge_ring(mesh, v1_id, v2_id, old_t_ids, n12_v_ids, n12_t_ids);

    bool is_valid = false;
    std::vector<Vector4i> new_tets;
    std::array<int, 2> v_ids;
    std::vector<Scalar> new_qs;
    Scalar old_max_quality = get_max_quality(mesh, old_t_ids);
    Scalar new_max_quality = 0;
    // The two candidate diagonals. Either loop below stopping short means the diagonal is out,
    // which is what the size checks after them say.
    for (int i = 0; i < 2; i++) {
        const std::array<int, 2> tmp_v_ids = {{n12_v_ids[i], n12_v_ids[2 + i]}};

        std::vector<Vector4i> tmp_new_tets;
        for (int t_id : old_t_ids) {
            // A tet that already has the diagonal's first end gives up its v2 corner to the
            // second end; one that does not gives up its v1 corner to the first.
            MeshTet t = tets[t_id];
            if (t.find(tmp_v_ids[0]) >= 0)
                t[t.find(v2_id)] = tmp_v_ids[1];
            else
                t[t.find(v1_id)] = tmp_v_ids[0];
            if (is_inverted(mesh, t))
                break;
            tmp_new_tets.push_back(t.indices);
        }
        if (tmp_new_tets.size() != old_t_ids.size())
            continue;

        std::vector<Scalar> tmp_new_qs;
        for (auto &t: tmp_new_tets) {
            Scalar q = get_quality(mesh, t);
            if (q >= old_max_quality)
                break;
            if (q > new_max_quality)
                new_max_quality = q;
            tmp_new_qs.push_back(q);
        }
        if (tmp_new_qs.size() != tmp_new_tets.size())
            continue;

        is_valid = true;
        old_max_quality = new_max_quality;
        new_tets = tmp_new_tets;
        new_qs = tmp_new_qs;
        v_ids = tmp_v_ids;
    }
    if (!is_valid)
        return false;

    FaceMarks marks;
    for (int t_id : old_t_ids)
        marks.record(tets[t_id], v1_id, v2_id);

    for (int j = 0; j < new_tets.size(); j++) {
        // Which end of the edge this tet gave up, read off the old tet, which is still in place
        // here, so it is the same side it was rebuilt on above.
        const int gave_up = tets[old_t_ids[j]].find(v_ids[0]) >= 0 ? 1 : 0;
        vector_erase(tet_vertices[gave_up == 0 ? v1_id : v2_id].conn_tets, old_t_ids[j]);
        tet_vertices[v_ids[gave_up]].conn_tets.push_back(old_t_ids[j]);
        tets[old_t_ids[j]].indices = new_tets[j];
        tets[old_t_ids[j]].quality = new_qs[j];
    }

    for (int t_id : old_t_ids) {  // now holding the rebuilt tets
        marks.reapply(tets[t_id],
                      [&](int v_id) { return v_id == v_ids[0] || v_id == v_ids[1]; });
    }

    collect_tet_edges(new_tets, new_edges);

    return true;
}

bool remove_an_edge_56(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    order_edge_ring(mesh, v1_id, v2_id, old_t_ids, n12_v_ids, n12_t_ids);

    Scalar old_max_quality = get_max_quality(mesh, old_t_ids);
    Scalar new_max_quality = 0;

    // The pair a tet of the ring turns into: one copy with v1 moved onto nv_id and one with v2
    // moved onto it. False, leaving out untouched, when either copy would be inverted.
    const auto split_pair = [&](const MeshTet& t, int nv_id, std::array<Vector4i, 2>& out) {
        for (int k = 0; k < 2; k++) {
            MeshTet nt = t;
            nt[nt.find(k == 0 ? v1_id : v2_id)] = nv_id;
            if (is_inverted(mesh, nt))
                return false;
            out[k] = nt.indices;
        }
        return true;
    };

    // Slots 0..4 hold the pair each ring position splits into, and 5..9 the pair for the tet
    // across from it, which is only filled for the position finally selected.
    std::array<std::array<Scalar, 2>, 10> tet_qs = {};
    std::array<std::array<Vector4i, 2>, 10> new_tets = {};
    std::vector<bool> is_v_valid(5, true);
    for (int i = 0; i < n12_v_ids.size(); i++) {
        if (!is_v_valid[(i + 1) % 5] && !is_v_valid[(i - 1 + 5) % 5])
            continue;

        std::array<Vector4i, 2> new_ts;
        if (!split_pair(tets[n12_t_ids[i]], n12_v_ids[(i - 1 + 5) % 5], new_ts)) {
            is_v_valid[(i + 1) % 5] = false;
            is_v_valid[(i - 1 + 5) % 5] = false;
            continue;
        }

        new_tets[i] = new_ts;
        for (int k = 0; k < 2; k++)
            tet_qs[i][k] = get_quality(mesh, new_ts[k]);
    }
    if (std::count(is_v_valid.begin(), is_v_valid.end(), true) == 0)
        return false;

    int selected_id = -1;
    for (int i = 0; i < is_v_valid.size(); i++) {
        if (!is_v_valid[i])
            continue;

        std::array<Vector4i, 2> new_ts;
        if (!split_pair(tets[n12_t_ids[(i + 2) % 5]], n12_v_ids[i], new_ts))
            continue;

        // The swap rebuilds six tets: this pair and the pairs at the two neighbouring positions.
        // new_max_quality carries over between iterations, as it always has.
        std::array<Scalar, 2> qs;
        for (int k = 0; k < 2; k++)
            qs[k] = get_quality(mesh, new_ts[k]);
        for (int k = 0; k < 2; k++) {
            for (Scalar q : {qs[k], tet_qs[(i + 1) % 5][k], tet_qs[(i - 1 + 5) % 5][k]}) {
                if (q > new_max_quality)
                    new_max_quality = q;
            }
        }
        if (new_max_quality >= old_max_quality)
            continue;

        old_max_quality = new_max_quality;
        selected_id = i;
        tet_qs[i + 5] = qs;
        new_tets[i + 5] = new_ts;
    }
    if (selected_id < 0)
        return false;

    FaceMarks marks;
    for (int t_id : old_t_ids)
        marks.record(tets[t_id], v1_id, v2_id);

    // Every slot below, the freshly claimed one included, is overwritten outright.
    std::vector<int> new_t_ids = old_t_ids;
    get_new_tet_slots(mesh, 1, new_t_ids);

    for (int i = 0; i < 2; i++) {
        tets[new_t_ids[i]] = new_tets[(selected_id + 1) % 5][i];
        tets[new_t_ids[i + 2]] = new_tets[(selected_id - 1 + 5) % 5][i];
        tets[new_t_ids[i + 4]] = new_tets[selected_id + 5][i];

        tets[new_t_ids[i]].quality = tet_qs[(selected_id + 1) % 5][i];
        tets[new_t_ids[i + 2]].quality = tet_qs[(selected_id - 1 + 5) % 5][i];
        tets[new_t_ids[i + 4]].quality = tet_qs[selected_id + 5][i];
    }

    // A face opposite one of the four vertices the swap touched is new.
    const int nv1 = n12_v_ids[(selected_id + 1) % 5];
    const int nv2 = n12_v_ids[(selected_id - 1 + 5) % 5];
    for (int t_id : new_t_ids) {
        marks.reapply(tets[t_id], [&](int v_id) {
            return v_id != v1_id && v_id != v2_id && v_id != nv1 && v_id != nv2;
        });
    }

    for (int i = 0; i < n12_v_ids.size(); i++) {
        vector_erase(tet_vertices[n12_v_ids[i]].conn_tets, n12_t_ids[i]);
        vector_erase(tet_vertices[n12_v_ids[i]].conn_tets, n12_t_ids[(i - 1 + 5) % 5]);
    }
    for (int i = 0; i < n12_t_ids.size(); i++) {
        vector_erase(tet_vertices[v1_id].conn_tets, n12_t_ids[i]);
        vector_erase(tet_vertices[v2_id].conn_tets, n12_t_ids[i]);
    }

    for (int t_id : new_t_ids) {
        for (int j = 0; j < 4; j++)
            tet_vertices[tets[t_id][j]].conn_tets.push_back(t_id);
    }

    collect_tet_edges(mesh, new_t_ids, new_edges);

    return true;
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::edge_swapping(Mesh& mesh) {
    auto &tet_vertices = mesh.tet_vertices;

    mesh.reset_t_empty_start();
    mesh.reset_v_empty_start();

    auto is_swappable = [&](int v1_id, int v2_id, const std::vector<int> &n12_t_ids) {
        if (n12_t_ids.size() < 3 || n12_t_ids.size() > 5)
            return false;
        if (!is_valid_edge(mesh, v1_id, v2_id))
            return false;
        if (is_surface_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        if (is_bbox_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        return true;
    };

    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);

    LongestFirstQueue es_queue;
    for (auto &e:edges) {
        es_queue.push({e, get_edge_length_2(mesh, e[0], e[1])});
    }
    edges.clear();

    while (!es_queue.empty()) {
        std::array<int, 2> v_ids = es_queue.top().v_ids;
        es_queue.pop();

        if(tet_vertices[v_ids[0]].is_freezed && tet_vertices[v_ids[1]].is_freezed)
            continue;

        std::vector<int> n12_t_ids;
        set_intersection(tet_vertices[v_ids[0]].conn_tets, tet_vertices[v_ids[1]].conn_tets, n12_t_ids);
        if (!is_swappable(v_ids[0], v_ids[1], n12_t_ids))
            continue;

        pop_duplicates(es_queue, v_ids);

        // is_swappable has already pinned the ring to one of these three sizes.
        std::vector<std::array<int, 2>> new_edges;
        if (n12_t_ids.size() == 3)
            remove_an_edge_32(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges);
        else if (n12_t_ids.size() == 4)
            remove_an_edge_44(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges);
        else
            remove_an_edge_56(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges);

        for (auto &e:new_edges) {
            es_queue.push({e, get_edge_length_2(mesh, e[0], e[1])});
        }
    }
}
