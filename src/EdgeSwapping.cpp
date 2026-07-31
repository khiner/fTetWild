// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/EdgeSwapping.h>
#include <floattetwild/LocalOperations.h>
#include <unordered_map>

namespace floatTetWild {
namespace {

// The surface, tag and bbox marks of the faces the swap keeps, keyed by their sorted corner triple
// so that the rebuilt tets can pick their own faces' marks back up.
struct FaceMarks
{
    std::vector<std::array<int, 3>> fs;
    std::vector<char>               is_sf, tags, is_bx;

    static std::array<int, 3> corners(const MeshTet& t, int j)
    {
        std::array<int, 3> f = {{t[mod4(j + 1)], t[mod4(j + 2)], t[mod4(j + 3)]}};
        std::sort(f.begin(), f.end());
        return f;
    }

    // Records the faces of t that contain v1_id or v2_id: those are the ones a swap preserves.
    void record(const MeshTet& t, int v1_id, int v2_id)
    {
        for (int j = 0; j < 4; j++) {
            if (t[j] != v1_id && t[j] != v2_id)
                continue;
            fs.push_back(corners(t, j));
            is_sf.push_back(t.is_surface_fs[j]);
            tags.push_back(t.surface_tags[j]);
            is_bx.push_back(t.is_bbox_fs[j]);
        }
    }

    static void clear(MeshTet& t, int j)
    {
        t.is_surface_fs[j] = NOT_SURFACE;
        t.surface_tags[j]  = NO_SURFACE_TAG;
        t.is_bbox_fs[j]    = NOT_BBOX;
    }

    // Copies the recorded marks onto face j of t, or clears them when that face was not recorded.
    void apply(MeshTet& t, int j) const
    {
        const auto it = std::find(fs.begin(), fs.end(), corners(t, j));
        if (it == fs.end()) {
            clear(t, j);
            return;
        }
        const size_t k     = it - fs.begin();
        t.is_surface_fs[j] = is_sf[k];
        t.surface_tags[j]  = tags[k];
        t.is_bbox_fs[j]    = is_bx[k];
    }
};

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
        std::array<int, 3> e;
        int                cnt = 0;
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id][j] != v2_id)
                e[cnt++] = mesh.tets[t_id][j];
        }
        e[cnt] = t_id;
        n12_es.push_back(e);
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
            else if (n12_es[j][1] == n12_v_ids.back()) {  // else if!!!!!!!!!!
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

}  // namespace
}  // namespace floatTetWild

void floatTetWild::edge_swapping(Mesh& mesh) {
    auto &tet_vertices = mesh.tet_vertices;

    mesh.reset_t_empty_start();
    mesh.reset_v_empty_start();

    auto is_swappable = [&](int v1_id, int v2_id, const std::vector<int> &n12_t_ids) {
        if (n12_t_ids.size() < 3 || n12_t_ids.size() > 5)
            return false;
        if (!is_valid_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        if (is_surface_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        if (is_bbox_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        return true;
    };

    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);

    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_l> es_queue;
    for (auto &e:edges) {
        es_queue.push(ElementInQueue(e, get_edge_length_2(mesh, e[0], e[1])));
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

        while (!es_queue.empty()) {
            if (es_queue.top().v_ids == v_ids)
                es_queue.pop();
            else
                break;
        }

        std::vector<std::array<int, 2>> new_edges;
        if (n12_t_ids.size() == 3)
            remove_an_edge_32(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges);
        if (n12_t_ids.size() == 4)
            remove_an_edge_44(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges);
        if (n12_t_ids.size() == 5)
            remove_an_edge_56(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges);

        for (auto &e:new_edges) {
            es_queue.push(ElementInQueue(e, get_edge_length_2(mesh, e[0], e[1])));
        }
    }
}

bool floatTetWild::remove_an_edge_32(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges){
    if(old_t_ids.size()!=3)
        return false;

    auto& tet_vertices = mesh.tet_vertices;
    auto& tets = mesh.tets;

    ////construct
    std::array<int, 2> v_ids;
    std::vector<MeshTet> new_tets;
    std::array<int, 2> t_ids;
    int cnt = 0;
    for (int i = 0; i < 4; i++) {
        if (tets[old_t_ids[0]][i] != v1_id && tets[old_t_ids[0]][i] != v2_id) {
            v_ids[cnt++] = tets[old_t_ids[0]][i];
        }
    }
    int found_j = tets[old_t_ids[1]].find(v_ids[0]);
    if(found_j>=0){
        new_tets.push_back(tets[old_t_ids[1]]);
        new_tets.push_back(tets[old_t_ids[2]]);
        t_ids = {{old_t_ids[1], old_t_ids[2]}};
    } else {
        new_tets.push_back(tets[old_t_ids[2]]);
        new_tets.push_back(tets[old_t_ids[1]]);
        t_ids = {{old_t_ids[2], old_t_ids[1]}};
    }
    found_j = new_tets[0].find(v1_id);
    new_tets[0][found_j] = v_ids[1];
    found_j = new_tets[1].find(v2_id);
    new_tets[1][found_j] = v_ids[0];

    ////check
    for(auto& t:new_tets) {
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]))
            return false;
    }
    std::vector<Scalar> new_qs;
    Scalar old_max_quality = 0;
    for(int t_id: old_t_ids) {
        if (tets[t_id].quality > old_max_quality)
            old_max_quality = tets[t_id].quality;
    }
    for(auto& t:new_tets) {
        Scalar q = get_quality(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]);
        if (q >= old_max_quality)//or use > ???
            return false;
        new_qs.push_back(q);
    }

    ////real update
    FaceMarks marks;
    for (int t_id : old_t_ids)
        marks.record(tets[t_id], v1_id, v2_id);

    tets[old_t_ids[0]].is_removed = true;
    tets[t_ids[0]] = new_tets[0];//v2
    tets[t_ids[1]] = new_tets[1];//v1

    for (int i = 0; i < 2; i++)
        tets[t_ids[i]].quality = new_qs[i];

    for(int i=0;i<4;i++) {
        // The face opposite the vertex the swap moved is new, so it carries no marks.
        if (tets[t_ids[0]][i] != v2_id)
            marks.apply(tets[t_ids[0]], i);
        else
            FaceMarks::clear(tets[t_ids[0]], i);

        if (tets[t_ids[1]][i] != v1_id)
            marks.apply(tets[t_ids[1]], i);
        else
            FaceMarks::clear(tets[t_ids[1]], i);
    }

    vector_erase(tet_vertices[v_ids[0]].conn_tets, old_t_ids[0]);
    vector_erase(tet_vertices[v_ids[1]].conn_tets, old_t_ids[0]);

    tet_vertices[v_ids[0]].conn_tets.push_back(t_ids[1]);
    tet_vertices[v_ids[1]].conn_tets.push_back(t_ids[0]);

    vector_erase(tet_vertices[v1_id].conn_tets, old_t_ids[0]);
    vector_erase(tet_vertices[v2_id].conn_tets, old_t_ids[0]);

    vector_erase(tet_vertices[v1_id].conn_tets, t_ids[0]);
    vector_erase(tet_vertices[v2_id].conn_tets, t_ids[1]);

    ////re-push

    new_edges.reserve(new_tets.size() * 6);
    for (const auto &t : new_tets)
        push_tet_edges(t, new_edges);
    vector_unique(new_edges);

    return true;
}

bool floatTetWild::remove_an_edge_44(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges) {

    const int N = 4;
    if (old_t_ids.size() != N)
        return false;

    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    ////construct
    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    order_edge_ring(mesh, v1_id, v2_id, old_t_ids, n12_v_ids, n12_t_ids);

    ////check
    bool is_valid = false;
    std::vector<Vector4i> new_tets;
    new_tets.reserve(4);
    std::vector<int> tags;
    std::array<int, 2> v_ids;
    std::vector<Scalar> new_qs;
    Scalar old_max_quality = 0;
    Scalar new_max_quality = 0;
    for (int t_id: old_t_ids) {
        if (tets[t_id].quality > old_max_quality)
            old_max_quality = tets[t_id].quality;
    }
    for (int i = 0; i < 2; i++) {
        std::vector<Vector4i> tmp_new_tets;
        std::vector<int> tmp_tags;
        std::array<int, 2> tmp_v_ids;
        tmp_v_ids = {{n12_v_ids[0 + i], n12_v_ids[2 + i]}};
        bool is_break = false;
        for (int j = 0; j < old_t_ids.size(); j++) {
            auto t = tets[old_t_ids[j]];
            int ii = t.find(tmp_v_ids[0]);
            if (ii >= 0) {
                int jt = t.find(v2_id);
                t[jt] = tmp_v_ids[1];
                tmp_tags.push_back(1);
            } else {
                int jt = t.find(v1_id);
                t[jt] = tmp_v_ids[0];
                tmp_tags.push_back(0);
            }
            if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]])) {
                is_break = true;
                break;
            }
            tmp_new_tets.push_back(t.indices);
        }
        if (is_break)
            continue;

        std::vector<Scalar> tmp_new_qs;
        for (auto &t: tmp_new_tets) {
            Scalar q = get_quality(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]);
            if (q >= old_max_quality) {
                is_break = true;
                break;
            }
            if (q > new_max_quality)
                new_max_quality = q;
            tmp_new_qs.push_back(q);
        }
        if (is_break)
            continue;

        is_valid = true;
        old_max_quality = new_max_quality;
        new_tets = tmp_new_tets;
        tags = tmp_tags;
        new_qs = tmp_new_qs;
        v_ids = tmp_v_ids;
    }
    if (!is_valid)
        return false;

    ////real update
    FaceMarks marks;
    for (int t_id : old_t_ids)
        marks.record(tets[t_id], v1_id, v2_id);

    for (int j = 0; j < new_tets.size(); j++) {
        if (tags[j] == 0) {
            vector_erase(tet_vertices[v1_id].conn_tets, old_t_ids[j]);
            tet_vertices[v_ids[0]].conn_tets.push_back(old_t_ids[j]);
        } else {
            vector_erase(tet_vertices[v2_id].conn_tets, old_t_ids[j]);
            tet_vertices[v_ids[1]].conn_tets.push_back(old_t_ids[j]);
        }
        tets[old_t_ids[j]].indices = new_tets[j];
        tets[old_t_ids[j]].quality = new_qs[j];
    }

    for (int t_id : old_t_ids) {//old_t_ids contains new tets
        for (int j = 0; j < 4; j++) {
            if (tets[t_id][j] == v_ids[0] || tets[t_id][j] == v_ids[1])
                marks.apply(tets[t_id], j);
            else
                FaceMarks::clear(tets[t_id], j);
        }
    }

    ////re-push
    new_edges.reserve(new_tets.size() * 6);
    for (const auto &t : new_tets)
        push_tet_edges(t, new_edges);
    vector_unique(new_edges);

    return true;
}

bool floatTetWild::remove_an_edge_56(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges) {
    if (old_t_ids.size() != 5)
        return false;

    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    ////construct
    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    order_edge_ring(mesh, v1_id, v2_id, old_t_ids, n12_v_ids, n12_t_ids);

    ////check
    Scalar old_max_quality = 0;
    Scalar new_max_quality = 0;
    for (int t_id: old_t_ids) {
        if (tets[t_id].quality > old_max_quality)
            old_max_quality = tets[t_id].quality;
    }

    std::unordered_map<int, std::array<Scalar, 2>> tet_qs;
    std::unordered_map<int, std::array<Vector4i, 2>> new_tets;
    std::vector<bool> is_v_valid(5, true);
    for (int i = 0; i < n12_v_ids.size(); i++) {
        if (!is_v_valid[(i + 1) % 5] && !is_v_valid[(i - 1 + 5) % 5])
            continue;

        std::vector<Vector4i> new_ts;
        new_ts.reserve(6);
        auto t = tets[n12_t_ids[i]];
        int it = t.find(v1_id);
        t[it] = n12_v_ids[(i - 1 + 5) % 5];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]])) {
            is_v_valid[(i + 1) % 5] = false;
            is_v_valid[(i - 1 + 5) % 5] = false;
            continue;
        }
        new_ts.push_back(t.indices);

        t = tets[n12_t_ids[i]];
        it = t.find(v2_id);
        t[it] = n12_v_ids[(i - 1 + 5) % 5];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]])) {
            is_v_valid[(i + 1) % 5] = false;
            is_v_valid[(i - 1 + 5) % 5] = false;
            continue;
        }
        new_ts.push_back(t.indices);
        new_tets[i] = std::array<Vector4i, 2>({{new_ts[0], new_ts[1]}});

        std::vector<Scalar> qs;
        for (auto &nt: new_ts) {
            qs.emplace_back(
                    get_quality(tet_vertices[nt[0]], tet_vertices[nt[1]], tet_vertices[nt[2]], tet_vertices[nt[3]]));
        }
        tet_qs[i] = std::array<Scalar, 2>({{qs[0], qs[1]}});
    }
    if (std::count(is_v_valid.begin(), is_v_valid.end(), true) == 0)
        return false;

    int selected_id = -1;
    for (int i = 0; i < is_v_valid.size(); i++) {
        if (!is_v_valid[i])
            continue;

        std::vector<Vector4i> new_ts;
        new_ts.reserve(6);
        auto t = tets[n12_t_ids[(i + 2) % 5]];
        int it = t.find(v1_id);
        t[it] = n12_v_ids[i];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]))
            continue;
        new_ts.push_back(t.indices);
        t = tets[n12_t_ids[(i + 2) % 5]];
        it = t.find(v2_id);
        t[it] = n12_v_ids[i];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]))
            continue;
        new_ts.push_back(t.indices);

        std::vector<Scalar> qs;
        for (auto &nt:new_ts) {
            qs.push_back(get_quality(tet_vertices[nt[0]], tet_vertices[nt[1]], tet_vertices[nt[2]], tet_vertices[nt[3]]));
        }
        for (int j = 0; j < 2; j++) {
            qs.push_back(tet_qs[(i + 1) % 5][j]);
            qs.push_back(tet_qs[(i - 1 + 5) % 5][j]);
        }
        if (qs.size() != 6) {
            assert("qs.size() != 6");
        }
        for (auto &q:qs) {
            if (q > new_max_quality)
                new_max_quality = q;
        }
        if (new_max_quality >= old_max_quality)
            continue;

        old_max_quality = new_max_quality;
        selected_id = i;
        tet_qs[i + 5] = std::array<Scalar, 2>({{qs[0], qs[1]}});
        new_tets[i + 5] = std::array<Vector4i, 2>({{new_ts[0], new_ts[1]}});
    }
    if (selected_id < 0)
        return false;

    ////real update
    //update on surface -- 1
    FaceMarks marks;
    for (int t_id : old_t_ids)
        marks.record(tets[t_id], v1_id, v2_id);

    std::vector<int> new_t_ids = old_t_ids;
    get_new_tet_slots(mesh, 1, new_t_ids);
    tets[new_t_ids.back()].reset();

    for (int i = 0; i < 2; i++) {
        tets[new_t_ids[i]] = new_tets[(selected_id + 1) % 5][i];
        tets[new_t_ids[i + 2]] = new_tets[(selected_id - 1 + 5) % 5][i];
        tets[new_t_ids[i + 4]] = new_tets[selected_id + 5][i];

        tets[new_t_ids[i]].quality = tet_qs[(selected_id + 1) % 5][i];
        tets[new_t_ids[i + 2]].quality = tet_qs[(selected_id - 1 + 5) % 5][i];
        tets[new_t_ids[i + 4]].quality = tet_qs[selected_id + 5][i];
    }

    //update on_surface -- 2
    for (int t_id : new_t_ids) {
        for (int j = 0; j < 4; j++) {
            // A face opposite one of the four vertices the swap touched is new.
            if (tets[t_id][j] != v1_id && tets[t_id][j] != v2_id
                && tets[t_id][j] != n12_v_ids[(selected_id + 1) % 5]
                && tets[t_id][j] != n12_v_ids[(selected_id - 1 + 5) % 5])
                marks.apply(tets[t_id], j);
            else
                FaceMarks::clear(tets[t_id], j);
        }
    }

    //update conn_tets
    for (int i = 0; i < n12_v_ids.size(); i++) {
        vector_erase(tet_vertices[n12_v_ids[i]].conn_tets, n12_t_ids[i]);
        vector_erase(tet_vertices[n12_v_ids[i]].conn_tets, n12_t_ids[(i - 1 + 5) % 5]);
    }
    for (int i = 0; i < n12_t_ids.size(); i++) {
        vector_erase(tet_vertices[v1_id].conn_tets, n12_t_ids[i]);
        vector_erase(tet_vertices[v2_id].conn_tets, n12_t_ids[i]);
    }

    //add
    for (int t_id : new_t_ids) {
        for (int j = 0; j < 4; j++)
            tet_vertices[tets[t_id][j]].conn_tets.push_back(t_id);
    }

    ////re-push
    new_edges.reserve(new_t_ids.size() * 6);
    for (int t_id : new_t_ids)
        push_tet_edges(tets[t_id], new_edges);
    vector_unique(new_edges);

    return true;
}