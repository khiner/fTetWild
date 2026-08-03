// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/LocalOperations.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Logger.hpp>

#include <floattetwild/ParallelFor.hpp>
#include <floattetwild/MeshCleanup.hpp>

#include <numeric>
#include <random>

namespace floatTetWild {
namespace {

// A Fisher-Yates shuffle over one process-wide generator, seeded so that a run is repeatable.
// Only one_ring_edge_set() shuffles, and the sequence it draws decides which collapses a round
// makes, so the generator has to stay exactly this one.
void shuffle_ids(std::vector<int>& vec)
{
    static std::mt19937 gen(42);
    for (int i = vec.size() - 1; i > 0; --i) {
        const int index = (gen() - std::mt19937::min()) /
                          double(std::mt19937::max() - std::mt19937::min()) * (i + 1);
        std::swap(vec[i], vec[index]);
    }
}

// A random subset of the edges whose incident faces are pairwise disjoint, so that collapsing
// all of them at once is the same as collapsing them one at a time.
void one_ring_edge_set(const std::vector<std::array<int, 2>>&      edges,
                       const std::vector<char>&                    v_is_removed,
                       const std::vector<char>&                    f_is_removed,
                       const std::vector<std::unordered_set<int>>& conn_fs,
                       std::vector<int>&                           safe_set)
{
    std::vector<int> indices(edges.size());
    std::iota(std::begin(indices), std::end(indices), 0);
    shuffle_ids(indices);

    std::vector<bool> unsafe_face(f_is_removed.size(), false);
    safe_set.clear();
    for (const int e_id : indices) {
        const auto e = edges[e_id];
        if (v_is_removed[e[0]] || v_is_removed[e[1]])
            continue;

        bool ok = true;
        for (const int v : e) {
            for (const int f : conn_fs[v]) {
                if (!f_is_removed[f] && unsafe_face[f]) {
                    ok = false;
                    break;
                }
            }
            if (!ok)
                break;
        }
        if (!ok)
            continue;

        safe_set.push_back(e_id);

        for (const int v : e) {
            for (const int f : conn_fs[v]) {
                if (!f_is_removed[f])
                    unsafe_face[f] = true;
            }
        }
    }
}

Scalar get_angle_cos(const Vector3& p,
                     const Vector3& p1,
                     const Vector3& p2)
{
    Vector3 v1  = p1 - p;
    Vector3 v2  = p2 - p;
    Scalar  res = v1.dot(v2) / (v1.norm() * v2.norm());
    if (res > 1)
        return 1;
    if (res < -1)
        return -1;
    return res;
}

bool is_out_envelope(const std::array<Vector3, 3>& vs,
                     const AABBWrapper&            tree,
                     const Parameters&             params)
{
    geo::index_t prev_facet = geo::NO_FACET;
    return sample_triangle_and_check_is_out(
      vs, params.dd_simplification, params.eps_2_simplification, tree, prev_facet);
}

bool remove_duplicates(std::vector<Vector3>&  input_vertices,
                       std::vector<Vector3i>& input_faces,
                       std::vector<int>&      input_tags,
                       const Parameters&      params)
{
    MatrixXs V_tmp = to_matrix(input_vertices), V_in;
    MatrixXi F_tmp = to_matrix(input_faces), F_in;

    MatrixXi IV, _;
    remove_duplicate_vertices(
      V_tmp, F_tmp, SCALAR_ZERO * params.bbox_diag_length, V_in, IV, _, F_in);
    for (int i = 0; i < F_in.rows(); i++) {
        int j_min = 0;
        for (int j = 1; j < 3; j++) {
            if (F_in(i, j) < F_in(i, j_min))
                j_min = j;
        }
        if (j_min == 0)
            continue;
        int v0_id = F_in(i, j_min);
        int v1_id = F_in(i, (j_min + 1) % 3);
        int v2_id = F_in(i, (j_min + 2) % 3);
        F_in.row(i) << v0_id, v1_id, v2_id;
    }
    F_tmp.resize(0, 0);
    MatrixXi IF;
    unique_rows(F_in, F_tmp, IF, _);
    F_in = F_tmp;
    // The faces were reordered and deduplicated, so their tags follow them.
    std::vector<int> face_tags(IF.rows());
    for (int i = 0; i < IF.rows(); i++)
        face_tags[i] = input_tags[IF(i)];
    if (V_in.rows() == 0 || F_in.rows() == 0) {
        input_tags = face_tags;
        return false;
    }

    logger().info("remove duplicates: ");
    logger().info("#v: {} -> {}", input_vertices.size(), V_in.rows());
    logger().info("#f: {} -> {}", input_faces.size(), F_in.rows());

    input_vertices.resize(V_in.rows());
    input_faces.clear();
    input_faces.reserve(F_in.rows());
    input_tags.clear();
    for (int i = 0; i < V_in.rows(); i++)
        input_vertices[i] = V_in.row(i);
    for (int i = 0; i < F_in.rows(); i++) {
        if (F_in(i, 0) == F_in(i, 1) || F_in(i, 0) == F_in(i, 2) || F_in(i, 2) == F_in(i, 1))
            continue;
        if (i > 0 && (F_in(i, 0) == F_in(i - 1, 0) && F_in(i, 1) == F_in(i - 1, 2) &&
                      F_in(i, 2) == F_in(i - 1, 1)))
            continue;
        const Vector3 p0 = V_in.row(F_in(i, 0));
        const Vector3 p1 = V_in.row(F_in(i, 1));
        const Vector3 p2 = V_in.row(F_in(i, 2));
        if ((p1 - p0).cross(p2 - p0).norm() / 2 <= SCALAR_ZERO * params.bbox_diag_length)
            continue;
        input_faces.push_back(F_in.row(i));
        input_tags.push_back(face_tags[i]);
    }

    return true;
}

void collapsing(std::vector<Vector3>&                 input_vertices,
                std::vector<Vector3i>&                input_faces,
                const AABBWrapper&                    tree,
                const Parameters&                     params,
                std::vector<char>&                    v_is_removed,
                std::vector<char>&                    f_is_removed,
                std::vector<std::unordered_set<int>>& conn_fs)
{
    std::vector<std::array<int, 2>> edges;

    const auto build_edges = [&]() {
        edges = parallel_collect<std::array<int, 2>>(
          0, input_faces.size(), [&](size_t f_id, std::vector<std::array<int, 2>>& out) {
              if (f_is_removed[f_id])
                  return;

              push_tri_edges(input_faces[f_id], out);
          });

        vector_unique(edges);
    };

    auto is_onering_clean = [&](int v_id) {
        std::vector<int> v_ids;
        v_ids.reserve(conn_fs[v_id].size() * 2);
        for (const auto& f_id : conn_fs[v_id]) {
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] != v_id)
                    v_ids.push_back(input_faces[f_id][j]);
            }
        }
        std::sort(v_ids.begin(), v_ids.end());

        if (v_ids.size() % 2 != 0)
            return false;
        for (int i = 0; i < v_ids.size(); i += 2) {
            if (v_ids[i] != v_ids[i + 1])
                return false;
        }

        return true;
    };

    // Moves v1 and v2 onto their midpoint projected back to the surface, or returns false and
    // leaves the mesh alone when that would flip a face or leave the envelope.
    auto remove_an_edge = [&](int v1_id, int v2_id, const std::vector<int>& n12_f_ids) {
        if (!is_onering_clean(v1_id) || !is_onering_clean(v2_id))
            return false;

        std::vector<int> new_f_ids;
        for (int f_id : conn_fs[v1_id]) {
            if (f_id != n12_f_ids[0] && f_id != n12_f_ids[1])
                new_f_ids.push_back(f_id);
        }
        for (int f_id : conn_fs[v2_id]) {
            if (f_id != n12_f_ids[0] && f_id != n12_f_ids[1])
                new_f_ids.push_back(f_id);
        }
        vector_unique(new_f_ids);

        Vector3 p = (input_vertices[v1_id] + input_vertices[v2_id]) / 2;
        tree.project_to_sf(p);

        // The corner of the face that sits on the edge, and so moves onto p. -1 when the face has
        // none, which the two checks below then skip.
        const auto moved_corner = [&](int f_id) {
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] == v1_id || input_faces[f_id][j] == v2_id)
                    return j;
            }
            return -1;
        };

        for (int f_id : new_f_ids) {
            const int j = moved_corner(f_id);
            if (j < 0)
                continue;
            const Vector3i& f = input_faces[f_id];
            const Vector3 old_nv = tri_normal(
              input_vertices[f[0]], input_vertices[f[1]], input_vertices[f[2]]);
            const Vector3 new_nv =
              tri_normal(p, input_vertices[f[(j + 1) % 3]], input_vertices[f[(j + 2) % 3]]);
            if (old_nv.dot(new_nv) <= 0)
                return false;
            if (new_nv.norm() / 2 <= SCALAR_ZERO_2)
                return false;
        }

        for (int f_id : new_f_ids) {
            const int j = moved_corner(f_id);
            if (j < 0)
                continue;
            const std::array<Vector3, 3> tri = {
              {p,
               input_vertices[input_faces[f_id][(j + 1) % 3]],
               input_vertices[input_faces[f_id][(j + 2) % 3]]}};
            if (is_out_envelope(tri, tree, params))
                return false;
        }

        v_is_removed[v1_id]   = true;
        input_vertices[v2_id] = p;
        // conn_fs is not pruned here. The cleanup pass after the parallel loop below drops the
        // faces marked removed, which is what keeps this safe to run on several edges at once.
        for (int f_id : n12_f_ids)
            f_is_removed[f_id] = true;
        for (int f_id : conn_fs[v1_id]) {
            if (f_is_removed[f_id])
                continue;
            conn_fs[v2_id].insert(f_id);
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] == v1_id)
                    input_faces[f_id][j] = v2_id;
            }
        }

        // No queue to push new edges onto. The loop below rebuilds the edge list and a fresh safe
        // set on every round instead.
        return true;
    };

    std::atomic<int> cnt(0);
    int              cnt_suc = 0;

    const int stopping = input_vertices.size() / 10000.;

    std::vector<int> safe_set;
    do {
        build_edges();
        one_ring_edge_set(edges, v_is_removed, f_is_removed, conn_fs, safe_set);
        cnt = 0;

        // safe_set holds edges with disjoint neighbourhoods, so these collapses do not interact.
        parallel_for(size_t(0), safe_set.size(), [&](size_t i) {
            std::array<int, 2>& v_ids = edges[safe_set[i]];

            std::vector<int> n12_f_ids;
            set_intersection(conn_fs[v_ids[0]], conn_fs[v_ids[1]], n12_f_ids);

            if (n12_f_ids.size() != 2)
                return;

            if (remove_an_edge(v_ids[0], v_ids[1], n12_f_ids))
                cnt++;
        });

        parallel_for(size_t(0), conn_fs.size(), [&](size_t i) {
            if (v_is_removed[i])
                return;
            std::vector<int> r_f_ids;
            for (int f_id : conn_fs[i]) {
                if (f_is_removed[f_id])
                    r_f_ids.push_back(f_id);
            }
            for (int f_id : r_f_ids)
                conn_fs[i].erase(f_id);
        });

        cnt_suc += cnt;
    } while (cnt > stopping);

    logger().debug("{}  faces are collapsed!!", cnt_suc);
}

void swapping(std::vector<Vector3>&                 input_vertices,
              std::vector<Vector3i>&                input_faces,
              const AABBWrapper&                    tree,
              const Parameters&                     params,
              const std::vector<char>&              f_is_removed,
              std::vector<std::unordered_set<int>>& conn_fs)
{
    std::vector<std::array<int, 2>> edges;
    edges.reserve(input_faces.size() * 6);
    for (int i = 0; i < input_faces.size(); i++) {
        if (f_is_removed[i])
            continue;
        push_tri_edges(input_faces[i], edges);
    }
    vector_unique(edges);

    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_l> sm_queue;
    for (auto& e : edges) {
        Scalar weight = (input_vertices[e[0]] - input_vertices[e[1]]).squaredNorm();
        sm_queue.push(ElementInQueue(e, weight));
        sm_queue.push(ElementInQueue(std::array<int, 2>({{e[1], e[0]}}), weight));
    }

    int cnt = 0;
    while (!sm_queue.empty()) {
        int v1_id = sm_queue.top().v_ids[0];
        int v2_id = sm_queue.top().v_ids[1];
        sm_queue.pop();

        std::vector<int> n12_f_ids;
        set_intersection(conn_fs[v1_id], conn_fs[v2_id], n12_f_ids);
        if (n12_f_ids.size() != 2)
            continue;

        // The first corner of each face that is off the edge.
        const auto corner_off_edge = [&](int f_id) {
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] != v1_id && input_faces[f_id][j] != v2_id)
                    return input_faces[f_id][j];
            }
            return -1;
        };
        const std::array<int, 2> n_v_ids = {
          {corner_off_edge(n12_f_ids[0]), corner_off_edge(n12_f_ids[1])}};

        Scalar cos_a0 =
          get_angle_cos(input_vertices[n_v_ids[0]], input_vertices[v1_id], input_vertices[v2_id]);
        Scalar cos_a1 =
          get_angle_cos(input_vertices[n_v_ids[1]], input_vertices[v1_id], input_vertices[v2_id]);
        std::array<Vector3, 2> old_nvs;
        for (int f = 0; f < 2; f++) {
            const Vector3i& fv = input_faces[n12_f_ids[f]];
            old_nvs[f] = tri_normal(input_vertices[fv[0]], input_vertices[fv[1]],
                                    input_vertices[fv[2]]).normalized();
        }
        if (cos_a0 > -0.999) {                          // maybe it's for avoiding numerical issue
            if (old_nvs[0].dot(old_nvs[1]) < 1 - 1e-6)  // not coplanar
                continue;
        }

        // check inversion: the two faces must not already face opposite ways
        if (old_nvs[0].dot(old_nvs[1]) < 0)
            continue;

        Scalar cos_a0_new = get_angle_cos(
          input_vertices[v1_id], input_vertices[n_v_ids[0]], input_vertices[n_v_ids[1]]);
        Scalar cos_a1_new = get_angle_cos(
          input_vertices[v2_id], input_vertices[n_v_ids[0]], input_vertices[n_v_ids[1]]);
        if (std::min(cos_a0_new, cos_a1_new) <= std::min(cos_a0, cos_a1))
            continue;

        if (is_out_envelope(
              {{input_vertices[v1_id], input_vertices[n_v_ids[0]], input_vertices[n_v_ids[1]]}},
              tree,
              params) ||
            is_out_envelope(
              {{input_vertices[v2_id], input_vertices[n_v_ids[0]], input_vertices[n_v_ids[1]]}},
              tree,
              params)) {
            continue;
        }

        for (int j = 0; j < 3; j++) {
            if (input_faces[n12_f_ids[0]][j] == v2_id)
                input_faces[n12_f_ids[0]][j] = n_v_ids[1];
            if (input_faces[n12_f_ids[1]][j] == v1_id)
                input_faces[n12_f_ids[1]][j] = n_v_ids[0];
        }
        conn_fs[v1_id].erase(n12_f_ids[1]);
        conn_fs[v2_id].erase(n12_f_ids[0]);
        conn_fs[n_v_ids[0]].insert(n12_f_ids[1]);
        conn_fs[n_v_ids[1]].insert(n12_f_ids[0]);
        cnt++;
    }

    logger().debug("{}  faces are swapped!!", cnt);
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::simplify(std::vector<Vector3>&  input_vertices,
                            std::vector<Vector3i>& input_faces,
                            std::vector<int>&      input_tags,
                            const AABBWrapper&     tree,
                            const Parameters&      params,
                            bool                   skip_simplify)
{
    remove_duplicates(input_vertices, input_faces, input_tags, params);
    if (skip_simplify)
        return;

    std::vector<char>                    v_is_removed(input_vertices.size(), false);
    std::vector<char>                    f_is_removed(input_faces.size(), false);
    std::vector<std::unordered_set<int>> conn_fs(input_vertices.size());
    for (int i = 0; i < input_faces.size(); i++) {
        for (int j = 0; j < 3; j++)
            conn_fs[input_faces[i][j]].insert(i);
    }

    collapsing(input_vertices, input_faces, tree, params, v_is_removed, f_is_removed, conn_fs);
    swapping(input_vertices, input_faces, tree, params, f_is_removed, conn_fs);

    std::vector<int>     map_v_ids(input_vertices.size(), -1);
    std::vector<Vector3> new_input_vertices;
    new_input_vertices.reserve(input_vertices.size());
    for (int i = 0; i < input_vertices.size(); i++) {
        if (v_is_removed[i] || conn_fs[i].empty())
            continue;
        map_v_ids[i] = new_input_vertices.size();
        new_input_vertices.push_back(input_vertices[i]);
    }
    input_vertices = new_input_vertices;

    std::vector<Vector3i> new_input_faces;
    std::vector<int>      new_input_tags;
    new_input_faces.reserve(input_faces.size());
    new_input_tags.reserve(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++) {
        if (f_is_removed[i])
            continue;
        Vector3i f;
        for (int j = 0; j < 3; j++)
            f[j] = map_v_ids[input_faces[i][j]];
        new_input_faces.push_back(f);
        new_input_tags.push_back(input_tags[i]);
    }
    input_faces = new_input_faces;
    input_tags  = new_input_tags;

    remove_duplicates(input_vertices, input_faces, input_tags, params);

    logger().info("#v = {}", input_vertices.size());
    logger().info("#f = {}", input_faces.size());
}
