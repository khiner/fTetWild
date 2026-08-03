// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <algorithm>
#include <random>

#include <floattetwild/TriangleInsertion.h>

#include <floattetwild/geo_geometry_nd.h>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/intersections.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/auto_table.hpp>

#include <floattetwild/ParallelFor.hpp>

#include <bitset>
#include <floattetwild/EdgeCollapsing.h>

// The faces the current insertion laid onto the cut plane. subdivide_tets fills it and
// simplify_subdivision_result consumes it, both under one insert_one_triangle, which is serial.
static std::vector<std::array<int, 3>> covered_tet_fs;

namespace floatTetWild {
namespace {

// The live tets whose bounding box overlaps [min, max] and that keep() accepts, in index order.
// The queue order below decides the result, and parallel_collect returns index order, which is
// the order the loop would push in if it ran on one thread.
template<typename Keep>
std::vector<int> tets_overlapping_bbox(const Mesh&    mesh,
                                       const Vector3& min,
                                       const Vector3& max,
                                       Keep&&         keep)
{
    return parallel_collect<int>(0, mesh.tets.size(), [&](size_t t_id, std::vector<int>& out) {
        if (mesh.tets[t_id].is_removed || !keep(t_id))
            return;

        Vector3 min_t, max_t;
        get_bbox({tet_pos(mesh, t_id, 0), tet_pos(mesh, t_id, 1), tet_pos(mesh, t_id, 2),
                  tet_pos(mesh, t_id, 3)},
                 min_t,
                 max_t);
        if (is_bbox_intersected(min, max, min_t, max_t))
            out.push_back(t_id);
    });
}

// Walks outward from seeds, visiting each tet once. visit(t_id, is_cut_vs) does the work for one
// tet and marks the corners the cut reached there; the walk carries on through every tet around
// those, and stops where nothing was reached.
template<typename Visit>
void flood_from_cut(const Mesh& mesh, const std::vector<int>& seeds, Visit&& visit)
{
    std::vector<bool> is_visited(mesh.tets.size(), false);
    std::queue<int>   queue;
    const auto        push = [&](int t_id) {
        if (!is_visited[t_id]) {
            is_visited[t_id] = true;
            queue.push(t_id);
        }
    };

    for (int t_id : seeds)
        push(t_id);

    while (!queue.empty()) {
        const int t_id = queue.front();
        queue.pop();

        std::array<bool, 4> is_cut_vs = {{false, false, false, false}};
        visit(t_id, is_cut_vs);

        for (int j = 0; j < 4; j++) {
            if (!is_cut_vs[j])
                continue;
            for (int n_t_id : mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets)
                push(n_t_id);
        }
    }
}

// The three corners of input face f_id.
std::array<Vector3, 3> face_points(const std::vector<Vector3>&  input_vertices,
                                   const std::vector<Vector3i>& input_faces,
                                   int                          f_id)
{
    const Vector3i& f = input_faces[f_id];
    return {{input_vertices[f[0]], input_vertices[f[1]], input_vertices[f[2]]}};
}

// The unit normal of input face f_id.
Vector3 input_face_normal(const std::vector<Vector3>&  input_vertices,
                          const std::vector<Vector3i>& input_faces,
                          int                          f_id)
{
    const std::array<Vector3, 3> p = face_points(input_vertices, input_faces, f_id);
    return tri_normal(p[0], p[1], p[2]).normalized();
}

// The unit normal of the plane through the three points, wound from the first corner. This is not
// tri_normal()'s spelling and does not round to the same value, so the two are kept apart.
Vector3 plane_normal(const std::array<Vector3, 3>& vs)
{
    Vector3 n = (vs[1] - vs[0]).cross(vs[2] - vs[0]);
    n.normalize();
    return n;
}

void match_surface_fs(const Mesh&                                   mesh,
                      const std::vector<Vector3i>&                  input_faces,
                      std::vector<bool>&                            is_face_inserted,
                      std::vector<std::array<std::vector<int>, 4>>& track_surface_fs)
{
    auto comp = [](const std::array<int, 4>& a, const std::array<int, 4>& b) {
        return std::tuple<int, int, int>(a[0], a[1], a[2]) <
               std::tuple<int, int, int>(b[0], b[1], b[2]);
    };

    std::vector<std::array<int, 4>> input_fs(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++) {
        input_fs[i] = {{input_faces[i][0], input_faces[i][1], input_faces[i][2], i}};
        std::sort(input_fs[i].begin(), input_fs[i].begin() + 3);
    }
    std::sort(input_fs.begin(), input_fs.end(), comp);

    for (int i = 0; i < mesh.tets.size(); i++) {
        auto& t = mesh.tets[i];
        for (int j = 0; j < 4; j++) {
            const std::array<int, 3> f = sorted_face(t, j);
            auto bounds = std::equal_range(
              input_fs.begin(), input_fs.end(), std::array<int, 4>({{f[0], f[1], f[2], -1}}), comp);
            for (auto it = bounds.first; it != bounds.second; ++it) {
                int f_id               = (*it)[3];
                is_face_inserted[f_id] = true;
                track_surface_fs[i][j].push_back(f_id);
            }
        }
    }
}

void sort_input_faces(const std::vector<Vector3i>& input_faces,
                      const Mesh&                  mesh,
                      std::vector<int>&            sorted_f_ids)
{
    sorted_f_ids.resize(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++)
        sorted_f_ids[i] = i;

    if (mesh.params.not_sort_input)
        return;
    std::mt19937 g(42);  // fixed seed, so the same input gives the same insertion order
    std::shuffle(sorted_f_ids.begin(), sorted_f_ids.end(), g);
}

void push_new_tets(Mesh&                                         mesh,
                   std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
                   std::vector<Vector3>&                         points,
                   std::vector<MeshTet>&                         new_tets,
                   std::vector<std::array<std::vector<int>, 4>>& new_track_surface_fs,
                   std::vector<int>&                             modified_t_ids)
{
    const int old_v_size = mesh.tet_vertices.size();
    mesh.tet_vertices.resize(mesh.tet_vertices.size() + points.size());
    for (int i = 0; i < points.size(); i++)
        mesh.tet_vertices[old_v_size + i].pos = points[i];

    // The first new_tets replace the subdivided ones in place. The rest are appended, so their ids
    // start at the current end of mesh.tets. subdivide_tets already set their tags.
    for (int i = 0; i < modified_t_ids.size(); i++) {
        const int t_id = modified_t_ids[i];
        for (int j = 0; j < 4; j++)
            vector_erase(mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets, t_id);

        mesh.tets[t_id]        = new_tets[i];
        track_surface_fs[t_id] = new_track_surface_fs[i];
        for (int j = 0; j < 4; j++)
            mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets.push_back(t_id);
    }
    for (int i = modified_t_ids.size(); i < new_tets.size(); i++) {
        const int t_id = mesh.tets.size() + i - modified_t_ids.size();
        for (int j = 0; j < 4; j++)
            mesh.tet_vertices[new_tets[i][j]].conn_tets.push_back(t_id);
    }

    mesh.tets.insert(mesh.tets.end(), new_tets.begin() + modified_t_ids.size(), new_tets.end());
    track_surface_fs.insert(track_surface_fs.end(),
                            new_track_surface_fs.begin() + modified_t_ids.size(),
                            new_track_surface_fs.end());
    modified_t_ids.clear();
}

void find_cutting_tets(int                           f_id,
                       const std::vector<Vector3>&   input_vertices,
                       const std::vector<Vector3i>&  input_faces,
                       const std::array<Vector3, 3>& vs,
                       Mesh&                         mesh,
                       std::vector<int>&             cut_t_ids,
                       bool                          is_again)
{
    std::vector<int> seeds;
    if (!is_again) {
        for (int j = 0; j < 3; j++) {
            seeds.insert(seeds.end(),
                         mesh.tet_vertices[input_faces[f_id][j]].conn_tets.begin(),
                         mesh.tet_vertices[input_faces[f_id][j]].conn_tets.end());
        }
        vector_unique(seeds);
    }
    else {
        const std::array<Vector3, 3> fps = face_points(input_vertices, input_faces, f_id);
        Vector3                      min_f, max_f;
        get_bbox({fps[0], fps[1], fps[2]}, min_f, max_f);
        seeds = tets_overlapping_bbox(mesh, min_f, max_f, [](size_t) { return true; });
    }

    flood_from_cut(mesh, seeds, [&](int t_id, std::array<bool, 4>& is_cut_vs) {
        if (is_again) {
            bool is_cut = false;
            for (int j = 0; j < 3; j++) {
                if (is_point_inside_tet(input_vertices[input_faces[f_id][j]],
                                        tet_pos(mesh, t_id, 0),
                                        tet_pos(mesh, t_id, 1),
                                        tet_pos(mesh, t_id, 2),
                                        tet_pos(mesh, t_id, 3))) {
                    is_cut = true;
                    break;
                }
            }
            if (is_cut) {
                cut_t_ids.push_back(t_id);
                is_cut_vs = {{true, true, true, true}};
                return;
            }
        }

        std::array<int, 4> oris;
        for (int j = 0; j < 4; j++)
            oris[j] = Predicates::orient_3d(vs[0], vs[1], vs[2], tet_pos(mesh, t_id, j));

        bool is_cutted = false;
        for (int j = 0; j < 4; j++) {
            const int face_oris[3] = {
              oris[(j + 1) % 4], oris[(j + 2) % 4], oris[(j + 3) % 4]};
            int cnt_pos, cnt_neg;
            count_orientations(face_oris, 3, cnt_pos, cnt_neg);
            const int cnt_on = 3 - cnt_pos - cnt_neg;

            const Vector3& tp1 = tet_pos(mesh, t_id, (j + 1) % 4);
            const Vector3& tp2 = tet_pos(mesh, t_id, (j + 2) % 4);
            const Vector3& tp3 = tet_pos(mesh, t_id, (j + 3) % 4);
            const auto cuts = [&](int hint) {
                return is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, hint) == hint;
            };

            // The face lies in the plane, or straddles it, so the whole face is on the cut.
            if (cnt_on == 3 || (cnt_pos > 0 && cnt_neg > 0)) {
                if (cuts(cnt_on == 3 ? CUT_COPLANAR : CUT_FACE)) {
                    is_cutted = true;
                    for (int k = 1; k < 4; k++)
                        is_cut_vs[(j + k) % 4] = true;
                }
            }
            // Exactly one of the face's edges is in the plane, so only its two ends are.
            else if (cnt_on == 2) {
                for (int e = CUT_EDGE_0; e <= CUT_EDGE_2; e++) {
                    const int a = (j + 1 + e) % 4;
                    const int b = (j + 1 + (e + 1) % 3) % 4;
                    if (oris[a] != Predicates::ORI_ZERO || oris[b] != Predicates::ORI_ZERO)
                        continue;
                    if (cuts(e)) {
                        is_cut_vs[a] = true;
                        is_cut_vs[b] = true;
                    }
                    break;
                }
            }
        }
        if (is_cutted)
            cut_t_ids.push_back(t_id);
    });
}

bool subdivide_tets(
  int                                           insert_f_id,
  Mesh&                                         mesh,
  CutMesh&                                      cut_mesh,
  std::vector<Vector3>&                         points,
  std::map<std::array<int, 2>, int>&            map_edge_to_intersecting_point,
  std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
  std::vector<int>&                             subdivide_t_ids,
  std::vector<bool>&                            is_mark_surface,
  std::vector<MeshTet>&                         new_tets,
  std::vector<std::array<std::vector<int>, 4>>& new_track_surface_fs,
  std::vector<int>&                             modified_t_ids)
{
    static const std::array<std::array<int, 2>, 6> t_es = {
      {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
    static const std::array<std::array<int, 3>, 4> t_f_es = {
      {{{1, 5, 4}}, {{5, 3, 2}}, {{3, 0, 4}}, {{0, 1, 2}}}};
    static const std::array<std::array<int, 3>, 4> t_f_vs = {
      {{{3, 1, 2}}, {{0, 2, 3}}, {{1, 3, 0}}, {{2, 0, 1}}}};

    covered_tet_fs.clear();
    for (int I = 0; I < subdivide_t_ids.size(); I++) {
        int  t_id       = subdivide_t_ids[I];
        bool is_mark_sf = is_mark_surface[I];

        // Per tet edge, {the local vertex id the cut point takes, its index in points}, or
        // {-1, -1} where the edge is not cut. The local ids carry on from the four corners.
        std::bitset<6>                     config_bit;
        std::array<std::pair<int, int>, 6> on_edge_p_ids;
        on_edge_p_ids.fill(std::make_pair(-1, -1));
        int cnt = 4;
        for (int i = 0; i < t_es.size(); i++) {
            const auto& le = t_es[i];
            const auto  it = map_edge_to_intersecting_point.find(
              sorted_edge(mesh.tets[t_id][le[0]], mesh.tets[t_id][le[1]]));
            if (it == map_edge_to_intersecting_point.end())
                continue;
            on_edge_p_ids[i] = std::make_pair(cnt++, it->second);
            config_bit.set(i);
        }
        int config_id = config_bit.to_ulong();
        if (config_id == 0) {  // no intersection
            if (is_mark_sf) {
                for (int j = 0; j < 4; j++) {
                    int cnt_on = 0;
                    for (int k = 0; k < 3; k++) {
                        if (cut_mesh.is_v_on_plane(
                              cut_mesh.map_v_ids[mesh.tets[t_id][(j + k + 1) % 4]])) {
                            cnt_on++;
                        }
                    }

                    if (cnt_on == 3) {
                        new_tets.push_back(mesh.tets[t_id]);
                        new_track_surface_fs.push_back(track_surface_fs[t_id]);
                        (new_track_surface_fs.back())[j].push_back(insert_f_id);
                        modified_t_ids.push_back(t_id);

                        covered_tet_fs.push_back({{mesh.tets[t_id][(j + 1) % 4],
                                                   mesh.tets[t_id][(j + 2) % 4],
                                                   mesh.tets[t_id][(j + 3) % 4]}});

                        break;
                    }
                }
            }
            continue;
        }

        const auto& configs = CutTable::get_tet_confs(config_id);
        if (configs.empty())
            continue;

        std::vector<Vector2i> my_diags;
        for (int j = 0; j < 4; j++) {
            std::vector<int> le_ids;
            for (int k = 0; k < 3; k++) {
                if (on_edge_p_ids[t_f_es[j][k]].first < 0)
                    continue;
                le_ids.push_back(k);
            }
            if (le_ids.size() != 2)  // no ambiguity
                continue;

            // The diagonal runs from whichever of the two cut edges carries the larger point id.
            const int cut = on_edge_p_ids[t_f_es[j][le_ids[0]]].second >
                                on_edge_p_ids[t_f_es[j][le_ids[1]]].second
                              ? le_ids[0]
                              : le_ids[1];
            my_diags.emplace_back();
            auto& diag = my_diags.back();
            diag << on_edge_p_ids[t_f_es[j][cut]].first, t_f_vs[j][cut];
            if (diag[0] > diag[1])
                std::swap(diag[0], diag[1]);
        }
        std::sort(my_diags.begin(), my_diags.end(), [](const Vector2i& a, const Vector2i& b) {
            return std::make_tuple(a[0], a[1]) < std::make_tuple(b[0], b[1]);
        });

        std::map<int, int> map_lv_to_v_id;
        const int          v_size  = mesh.tet_vertices.size();
        const int          vp_size = mesh.tet_vertices.size() + points.size();
        for (int i = 0; i < 4; i++)
            map_lv_to_v_id[i] = mesh.tets[t_id][i];
        cnt = 0;
        for (int i = 0; i < t_es.size(); i++) {
            if (config_bit[i] == 0)
                continue;
            map_lv_to_v_id[4 + cnt] = v_size + on_edge_p_ids[i].second;
            cnt++;
        }

        // A local vertex id is either one the mesh already has or one of the intersection points
        // appended after them.
        const auto pos_of = [&](int lv_id) -> const Vector3& {
            const int v_id = map_lv_to_v_id[lv_id];
            return v_id < v_size ? mesh.tet_vertices[v_id].pos : points[v_id - v_size];
        };

        auto get_centroid = [&](const std::vector<Vector4i>& config, int lv_id, Vector3& c) {
            std::vector<int> n_ids;
            for (const auto& tet : config) {
                std::vector<int> tmp;
                for (int j = 0; j < 4; j++) {
                    if (tet[j] != lv_id)
                        tmp.push_back(tet[j]);
                }
                if (tmp.size() == 4)
                    continue;
                n_ids.insert(n_ids.end(), tmp.begin(), tmp.end());
            }
            vector_unique(n_ids);
            c << 0, 0, 0;
            for (int n_id : n_ids)
                c += pos_of(n_id);
            c /= n_ids.size();
        };

        auto check_config = [&](int                                   diag_config_id,
                                std::vector<std::pair<int, Vector3>>& centroids) {
            const std::vector<Vector4i>& config = CutTable::get_tet_conf(config_id, diag_config_id);
            Scalar                       min_q  = -666;
            int                          cnt    = 0;
            std::map<int, int>           map_lv_to_c;
            for (const auto& tet : config) {
                std::array<Vector3, 4> vs;
                for (int j = 0; j < 4; j++) {
                    if (map_lv_to_v_id.find(tet[j]) == map_lv_to_v_id.end()) {
                        if (map_lv_to_c.find(tet[j]) == map_lv_to_c.end()) {
                            get_centroid(config, tet[j], vs[j]);
                            map_lv_to_c[tet[j]] = centroids.size();
                            centroids.push_back(std::make_pair(tet[j], vs[j]));
                        }
                        vs[j] = centroids[map_lv_to_c[tet[j]]].second;
                    }
                    else
                        vs[j] = pos_of(tet[j]);
                }

                Scalar volume = Predicates::orient_3d_volume(vs[0], vs[1], vs[2], vs[3]);

                if (cnt++ == 0 || volume < min_q)
                    min_q = volume;
            }

            return min_q;
        };

        int                                  diag_config_id = 0;
        std::vector<std::pair<int, Vector3>> centroids;
        if (!my_diags.empty()) {
            const auto&                         all_diags = CutTable::get_diag_confs(config_id);
            std::vector<std::pair<int, Scalar>> min_qualities(all_diags.size());
            std::vector<std::vector<std::pair<int, Vector3>>> all_centroids(all_diags.size());
            for (int i = 0; i < all_diags.size(); i++) {
                if (my_diags != all_diags[i]) {
                    min_qualities[i] = std::make_pair(i, -1);
                    continue;
                }

                std::vector<std::pair<int, Vector3>> tmp_centroids;
                Scalar                               min_q = check_config(i, tmp_centroids);
                min_qualities[i] = std::make_pair(i, min_q);
                all_centroids[i] = tmp_centroids;
            }
            std::sort(min_qualities.begin(),
                      min_qualities.end(),
                      [](const std::pair<int, Scalar>& a, const std::pair<int, Scalar>& b) {
                          return a.second < b.second;
                      });

            if (min_qualities.back().second < SCALAR_ZERO_3)
                return false;

            diag_config_id = min_qualities.back().first;
            centroids      = all_centroids[diag_config_id];
        }
        else if (check_config(diag_config_id, centroids) < SCALAR_ZERO_3)
            return false;

        for (int i = 0; i < centroids.size(); i++) {
            map_lv_to_v_id[centroids[i].first] = vp_size + i;
            points.push_back(centroids[i].second);
        }

        const auto& config            = CutTable::get_tet_conf(config_id, diag_config_id);
        const auto& new_is_surface_fs = CutTable::get_surface_conf(config_id, diag_config_id);
        const auto& new_local_f_ids   = CutTable::get_face_id_conf(config_id, diag_config_id);
        for (int i = 0; i < config.size(); i++) {
            const auto& lt = config[i];
            new_tets.push_back(MeshTet(map_lv_to_v_id[lt[0]],
                                       map_lv_to_v_id[lt[1]],
                                       map_lv_to_v_id[lt[2]],
                                       map_lv_to_v_id[lt[3]]));
            new_track_surface_fs.emplace_back();
            MeshTet&                         nt  = new_tets.back();
            std::array<std::vector<int>, 4>& nfs = new_track_surface_fs.back();

            for (int j = 0; j < 4; j++) {
                if (new_is_surface_fs[i][j] && is_mark_sf) {
                    nfs[j].push_back(insert_f_id);
                    covered_tet_fs.push_back(
                      {{nt[(j + 1) % 4], nt[(j + 3) % 4], nt[(j + 2) % 4]}});
                }

                // A face of the new tet that lies on one of the old tet's faces inherits it.
                const int old_local_f_id = new_local_f_ids[i][j];
                if (old_local_f_id < 0)
                    continue;
                const std::vector<int>& old_fs = track_surface_fs[t_id][old_local_f_id];
                nfs[j].insert(nfs[j].end(), old_fs.begin(), old_fs.end());
                nt.is_bbox_fs[j]    = mesh.tets[t_id].is_bbox_fs[old_local_f_id];
                nt.is_surface_fs[j] = mesh.tets[t_id].is_surface_fs[old_local_f_id];
                nt.surface_tags[j]  = mesh.tets[t_id].surface_tags[old_local_f_id];
            }
        }
        modified_t_ids.push_back(t_id);
    }

    return true;
}

void simplify_subdivision_result(
  int                                           insert_f_id,
  int                                           input_v_size,
  Mesh&                                         mesh,
  AABBWrapper&                                  tree,
  std::vector<std::array<std::vector<int>, 4>>& track_surface_fs)
{
    if (covered_tet_fs.empty())
        return;

    for (auto& f : covered_tet_fs)
        std::sort(f.begin(), f.end());
    vector_unique(covered_tet_fs);

    // The edges of the covered faces, each tagged with the face it came from.
    std::vector<std::array<int, 3>> edges;
    for (int i = 0; i < covered_tet_fs.size(); i++) {
        const auto& f = covered_tet_fs[i];
        for (int j = 0; j < 3; j++) {
            const std::array<int, 2> e = sorted_edge(f[j], f[(j + 1) % 3]);
            edges.push_back({{e[0], e[1], i}});
        }
    }
    std::sort(
      edges.begin(), edges.end(), [](const std::array<int, 3>& a, const std::array<int, 3>& b) {
          return std::make_tuple(a[0], a[1]) < std::make_tuple(b[0], b[1]);
      });
    std::unordered_set<int> freezed_v_ids;
    bool                    is_duplicated = false;
    for (int i = 0; i < edges.size() - 1; i++) {
        if (edges[i][0] < input_v_size)
            freezed_v_ids.insert(edges[i][0]);
        if (edges[i][1] < input_v_size)
            freezed_v_ids.insert(edges[i][1]);

        if (edges[i][0] == edges[i + 1][0] && edges[i][1] == edges[i + 1][1]) {
            is_duplicated = true;
            edges.erase(edges.begin() + i);
            i--;
        }
        else {
            if (!is_duplicated) {
                freezed_v_ids.insert(edges[i][0]);
                freezed_v_ids.insert(edges[i][1]);
            }
            is_duplicated = false;
        }
    }
    if (!is_duplicated) {
        freezed_v_ids.insert(edges.back()[0]);
        freezed_v_ids.insert(edges.back()[1]);
    }

    ShortestFirstQueue ec_queue;
    // Collapsing an edge moves its first end onto its second, so an edge goes in once per end that
    // is free to move.
    const auto push_free_ends = [&](int v0_id, int v1_id) {
        Scalar l_2 = get_edge_length_2(mesh, v0_id, v1_id);
        if (freezed_v_ids.find(v0_id) == freezed_v_ids.end())
            ec_queue.push({{{v0_id, v1_id}}, l_2});
        if (freezed_v_ids.find(v1_id) == freezed_v_ids.end())
            ec_queue.push({{{v1_id, v0_id}}, l_2});
    };

    for (const auto& e : edges)
        push_free_ends(e[0], e[1]);
    if (ec_queue.empty())
        return;
    std::unordered_set<int> all_v_ids;
    for (const auto& e : edges) {
        all_v_ids.insert(e[0]);
        all_v_ids.insert(e[1]);
    }
    // A collapse here may only touch tets that this insertion owns outright: no face tracked to
    // another input face, and no surface or bbox face other than the ones facing v_id.
    const auto only_this_insertion = [&](int v_id) {
        for (int t_id : mesh.tet_vertices[v_id].conn_tets) {
            for (int j = 0; j < 4; j++) {
                if ((!track_surface_fs[t_id][j].empty() &&
                     !vector_contains(track_surface_fs[t_id][j], insert_f_id)) ||
                    (mesh.tets[t_id][j] != v_id &&
                     (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE ||
                      mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX)))
                    return false;
            }
        }
        return true;
    };

    std::vector<int> tet_tss;  // the timestamps collapse_an_edge keeps, unused here
    while (!ec_queue.empty()) {
        std::array<int, 2> v_ids      = ec_queue.top().v_ids;
        Scalar             old_weight = ec_queue.top().weight;
        ec_queue.pop();

        pop_duplicates(ec_queue, v_ids);

        if (!is_valid_edge(mesh, v_ids[0], v_ids[1]))
            continue;
        if (freezed_v_ids.find(v_ids[0]) != freezed_v_ids.end())
            continue;
        Scalar weight = get_edge_length_2(mesh, v_ids[0], v_ids[1]);
        if (weight != old_weight)
            continue;

        if (!only_this_insertion(v_ids[0]))
            continue;

        std::vector<std::array<int, 2>> new_edges;
        for (int t_id : mesh.tet_vertices[v_ids[0]].conn_tets)
            mesh.tets[t_id].quality = get_quality(mesh, t_id);

        if (collapse_an_edge(
              mesh, v_ids[0], v_ids[1], tree, new_edges, 0, tet_tss, true, false)) {
            for (const auto& e : new_edges) {
                if (all_v_ids.count(e[0]) && all_v_ids.count(e[1]))
                    push_free_ends(e[0], e[1]);
            }
        }
    }
}

bool is_uninserted_face_covered(int                          uninserted_f_id,
                                const std::vector<Vector3>&  input_vertices,
                                const std::vector<Vector3i>& input_faces,
                                const std::vector<int>&      cut_t_ids,
                                Mesh&                        mesh)
{
    std::vector<geo::vec3> ps;
    sample_triangle(face_points(input_vertices, input_faces, uninserted_f_id), ps, mesh.params.dd);

    std::vector<int> n_t_ids;
    for (int t_id : cut_t_ids) {
        for (int j = 0; j < 4; j++)
            n_t_ids.insert(n_t_ids.end(),
                           mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets.begin(),
                           mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets.end());
    }
    vector_unique(n_t_ids);

    std::vector<std::array<int, 3>> faces;
    for (int t_id : n_t_ids) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE)
                faces.push_back(sorted_face(mesh.tets[t_id], j));
        }
    }
    vector_unique(faces);

    for (auto& p : ps) {
        bool is_valid = false;
        for (auto& f : faces) {
            double dis_2 =
              geo::Geom::point_triangle_squared_distance(p,
                                                         to_geo_p(mesh.tet_vertices[f[0]].pos),
                                                         to_geo_p(mesh.tet_vertices[f[1]].pos),
                                                         to_geo_p(mesh.tet_vertices[f[2]].pos));
            if (dis_2 < mesh.params.eps_2) {
                is_valid = true;
                break;
            }
        }
        if (!is_valid)
            return false;
    }

    return true;
}

bool insert_one_triangle(
  int                                           insert_f_id,
  const std::vector<Vector3>&                   input_vertices,
  const std::vector<Vector3i>&                  input_faces,
  Mesh&                                         mesh,
  std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
  AABBWrapper&                                  tree,
  bool                                          is_again)
{
    const std::array<Vector3, 3> vs = face_points(input_vertices, input_faces, insert_f_id);
    const Vector3                n  = plane_normal(vs);

    std::vector<int> cut_t_ids;
    find_cutting_tets(insert_f_id, input_vertices, input_faces, vs, mesh, cut_t_ids, is_again);

    if (cut_t_ids.empty())
        return false;

    // A step that could not cut. On a re-run the face may already be covered by the surface
    // around the tets that were reached, in which case there is nothing left to insert.
    const auto failed = [&](const char* what) {
        if (is_again &&
            is_uninserted_face_covered(insert_f_id, input_vertices, input_faces, cut_t_ids, mesh))
            return true;
        logger().warn(what);
        return false;
    };

    CutMesh cut_mesh(mesh, n, vs);
    cut_mesh.construct(cut_t_ids);

    if (cut_mesh.snap_to_plane()) {
        cut_mesh.project_to_plane(input_vertices.size());
        cut_mesh.expand_new(cut_t_ids);
        cut_mesh.project_to_plane(input_vertices.size());
    }

    std::vector<Vector3>              points;
    std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
    std::vector<int>                  subdivide_t_ids;
    if (!cut_mesh.get_intersecting_edges_and_points(
          points, map_edge_to_intersecting_point, subdivide_t_ids))
        return failed("FAIL get_intersecting_edges_and_points");
    vector_unique(cut_t_ids);
    std::vector<int> tmp;
    std::set_difference(subdivide_t_ids.begin(),
                        subdivide_t_ids.end(),
                        cut_t_ids.begin(),
                        cut_t_ids.end(),
                        std::back_inserter(tmp));
    std::vector<bool> is_mark_surface(cut_t_ids.size(), true);
    cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
    is_mark_surface.resize(is_mark_surface.size() + tmp.size(), false);

    std::vector<MeshTet>                         new_tets;
    std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
    std::vector<int>                             modified_t_ids;
    if (!subdivide_tets(insert_f_id,
                        mesh,
                        cut_mesh,
                        points,
                        map_edge_to_intersecting_point,
                        track_surface_fs,
                        cut_t_ids,
                        is_mark_surface,
                        new_tets,
                        new_track_surface_fs,
                        modified_t_ids))
        return failed("FAIL subdivide_tets");

    push_new_tets(
      mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids);

    simplify_subdivision_result(insert_f_id, input_vertices.size(), mesh, tree, track_surface_fs);

    return true;
}

void pair_track_surface_fs(
  Mesh&                                         mesh,
  std::vector<std::array<std::vector<int>, 4>>& track_surface_fs)
{
    std::vector<std::array<bool, 4>> is_visited(track_surface_fs.size(),
                                                {{false, false, false, false}});
    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
            if (is_visited[t_id][j])
                continue;
            is_visited[t_id][j] = true;
            if (track_surface_fs[t_id][j].empty())
                continue;
            int opp_t_id, k;
            if (!get_opp_face(mesh, t_id, j, opp_t_id, k))
                continue;
            is_visited[opp_t_id][k] = true;
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
            std::vector<int> f_ids;
            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k]) {
                std::set_union(track_surface_fs[t_id][j].begin(),
                               track_surface_fs[t_id][j].end(),
                               track_surface_fs[opp_t_id][k].begin(),
                               track_surface_fs[opp_t_id][k].end(),
                               std::back_inserter(f_ids));
                track_surface_fs[t_id][j]     = f_ids;
                track_surface_fs[opp_t_id][k] = f_ids;
            }
        }
    }
}

bool insert_boundary_edges_get_intersecting_edges_and_points(
  const std::vector<std::vector<std::pair<int, int>>>& covered_fs_infos,
  const std::vector<Vector3>&                          input_vertices,
  const std::vector<Vector3i>&                         input_faces,
  const std::array<int, 2>&                            e,
  const std::vector<int>&                              n_f_ids,
  std::vector<std::array<std::vector<int>, 4>>&        track_surface_fs,
  Mesh&                                                mesh,
  std::vector<Vector3>&                                points,
  std::map<std::array<int, 2>, int>&                   map_edge_to_intersecting_point,
  std::vector<int>&                                    snapped_v_ids,
  std::vector<std::array<int, 3>>&                     cut_fs,
  bool                                                 is_again)
{

    // The edge is cut in the plane of the first of the faces it belongs to.
    const std::array<Vector3, 3> fps = face_points(input_vertices, input_faces, n_f_ids.front());
    const int                    t   = get_t(fps[0], fps[1], fps[2]);
    std::array<Vector2, 2>       evs_2d = {
      {to_2d(input_vertices[e[0]], t), to_2d(input_vertices[e[1]], t)}};
    const Vector3& pp = fps[0];
    const Vector3  n  = plane_normal(fps);

    /// find seed t_ids: the tets already covering one of the faces, plus the tets around the edge
    std::vector<int> t_ids;
    for (int f_id : n_f_ids) {
        for (const auto& info : covered_fs_infos[f_id])
            t_ids.push_back(info.first);
    }
    if (!is_again) {
        for (int v_id : e)
            t_ids.insert(t_ids.end(),
                         mesh.tet_vertices[v_id].conn_tets.begin(),
                         mesh.tet_vertices[v_id].conn_tets.end());
    }
    else {
        Vector3 min_e, max_e;
        get_bbox({input_vertices[e[0]], input_vertices[e[0]], input_vertices[e[1]]}, min_e, max_e);

        const std::vector<int> bbox_t_ids =
          tets_overlapping_bbox(mesh, min_e, max_e, [&](size_t t_id) {
              return !track_surface_fs[t_id][0].empty() || !track_surface_fs[t_id][1].empty() ||
                     !track_surface_fs[t_id][2].empty() || !track_surface_fs[t_id][3].empty();
          });
        t_ids.insert(t_ids.end(), bbox_t_ids.begin(), bbox_t_ids.end());
    }
    vector_unique(t_ids);

    std::vector<int> v_oris(mesh.tet_vertices.size(), Predicates::ORI_UNKNOWN);
    flood_from_cut(mesh, t_ids, [&](int t_id, std::array<bool, 4>& is_cut_vs) {
        for (int j = 0; j < 4; j++) {
            if (track_surface_fs[t_id][j].empty())
                continue;
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            std::vector<int> tmp;
            std::set_intersection(track_surface_fs[t_id][j].begin(),
                                  track_surface_fs[t_id][j].end(),
                                  n_f_ids.begin(),
                                  n_f_ids.end(),
                                  std::back_inserter(tmp));
            if (tmp.empty())
                continue;

            std::array<int, 3>     f_v_ids = face_corners(mesh.tets[t_id], j);
            std::array<Vector2, 3> fvs_2d  = {{to_2d(mesh.tet_vertices[f_v_ids[0]].pos, n, pp, t),
                                               to_2d(mesh.tet_vertices[f_v_ids[1]].pos, n, pp, t),
                                               to_2d(mesh.tet_vertices[f_v_ids[2]].pos, n, pp, t)}};
            int                    cnt_pos = 0;
            int                    cnt_neg = 0;
            int                    cnt_on  = 0;
            for (int k = 0; k < 3; k++) {
                int& ori = v_oris[f_v_ids[k]];
                if (ori == Predicates::ORI_UNKNOWN)
                    ori = Predicates::orient_2d(evs_2d[0], evs_2d[1], fvs_2d[k]);
                if (ori == Predicates::ORI_ZERO) {
                    cnt_on++;
                }
                else {
                    Scalar dis_2 = p_seg_squared_dist_3d(mesh.tet_vertices[f_v_ids[k]].pos,
                                                         input_vertices[e[0]],
                                                         input_vertices[e[1]]);
                    if (dis_2 < mesh.params.eps_2_coplanar) {
                        ori = Predicates::ORI_ZERO;
                        cnt_on++;
                        continue;
                    }
                    if (ori == Predicates::ORI_POSITIVE)
                        cnt_pos++;
                    else
                        cnt_neg++;
                }
            }
            if (cnt_on >= 2) {
                cut_fs.push_back(sorted_face(mesh.tets[t_id], j));
                continue;
            }
            if (cnt_neg == 0 || cnt_pos == 0)
                continue;

            bool is_intersected = false;
            for (auto& p : evs_2d) {
                if (is_p_inside_tri_2d(p, fvs_2d)) {
                    is_intersected = true;
                    break;
                }
            }
            if (!is_intersected) {
                for (int k = 0; k < 3; k++) {
                    if (!is_crossing(v_oris[f_v_ids[k]], v_oris[f_v_ids[(k + 1) % 3]]))
                        continue;
                    std::array<int, 2> tri_e = sorted_edge(f_v_ids[k], f_v_ids[(k + 1) % 3]);
                    if (map_edge_to_intersecting_point.find(tri_e) !=
                        map_edge_to_intersecting_point.end()) {
                        is_intersected = true;
                        break;
                    }
                    double t2 = -1;
                    if (seg_seg_intersection_2d(evs_2d, {{fvs_2d[k], fvs_2d[(k + 1) % 3]}}, t2)) {
                        Vector3 p = (1 - t2) * mesh.tet_vertices[f_v_ids[k]].pos +
                                    t2 * mesh.tet_vertices[f_v_ids[(k + 1) % 3]].pos;
                        double dis1 = (p - mesh.tet_vertices[f_v_ids[k]].pos).squaredNorm();
                        double dis2 =
                          (p - mesh.tet_vertices[f_v_ids[(k + 1) % 3]].pos).squaredNorm();
                        if (dis1 < mesh.params.eps_2_coplanar) {
                            v_oris[f_v_ids[k]] = Predicates::ORI_ZERO;
                            is_intersected     = true;
                            break;
                        }
                        if (dis2 < mesh.params.eps_2_coplanar) {
                            v_oris[f_v_ids[k]] = Predicates::ORI_ZERO;
                            is_intersected     = true;
                            break;
                        }
                        points.push_back(p);
                        map_edge_to_intersecting_point[tri_e] = points.size() - 1;
                        is_intersected                        = true;
                        break;
                    }
                }
                if (!is_intersected)  /// no need to return false here
                    continue;
            }

            std::sort(f_v_ids.begin(), f_v_ids.end());
            cut_fs.push_back(f_v_ids);
            is_cut_vs[(j + 1) % 4] = true;
            is_cut_vs[(j + 2) % 4] = true;
            is_cut_vs[(j + 3) % 4] = true;
        }
    });
    vector_unique(cut_fs);

    std::vector<std::array<int, 2>> tet_edges;
    for (const auto& f : cut_fs)
        push_tri_edges(f, tet_edges);
    vector_unique(tet_edges);
    for (const auto& tet_e : tet_edges) {
        bool is_snapped = false;
        for (int j = 0; j < 2; j++) {
            if (v_oris[tet_e[j]] == Predicates::ORI_ZERO) {
                snapped_v_ids.push_back(tet_e[j]);
                is_snapped = true;
            }
        }
        if (is_snapped)
            continue;
        if (map_edge_to_intersecting_point.find(tet_e) != map_edge_to_intersecting_point.end())
            continue;
        std::array<Vector2, 2> tri_evs_2d = {{to_2d(mesh.tet_vertices[tet_e[0]].pos, n, pp, t),
                                              to_2d(mesh.tet_vertices[tet_e[1]].pos, n, pp, t)}};
        Scalar                 t_seg      = -1;
        if (seg_line_intersection_2d(tri_evs_2d, evs_2d, t_seg)) {
            points.push_back((1 - t_seg) * mesh.tet_vertices[tet_e[0]].pos +
                             t_seg * mesh.tet_vertices[tet_e[1]].pos);
            map_edge_to_intersecting_point[tet_e] = points.size() - 1;
        }
    }
    vector_unique(snapped_v_ids);

    return true;
}

void insert_boundary_edges(
  const std::vector<Vector3>&                                   input_vertices,
  const std::vector<Vector3i>&                                  input_faces,
  std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos,
  std::vector<bool>&                                            is_on_cut_edges,
  std::vector<std::array<std::vector<int>, 4>>&                 track_surface_fs,
  Mesh&                                                         mesh,
  AABBWrapper&                                                  tree,
  std::vector<bool>&                                            is_face_inserted,
  bool                                                          is_again,
  std::vector<std::array<int, 3>>&                              known_surface_fs,
  std::vector<std::array<int, 3>>&                              known_not_surface_fs)
{
    auto mark_known_surface_fs = [&](const std::array<int, 3>& f, int tag) {
        std::vector<int> n_t_ids;
        get_face_tets(mesh, f[0], f[1], f[2], n_t_ids);
        if (n_t_ids.size() != 2)
            return;

        for (int t_id : n_t_ids) {
            int j                            = get_local_f_id(t_id, f[0], f[1], f[2], mesh);
            mesh.tets[t_id].is_surface_fs[j] = tag;
        }
    };

    auto record_boundary_info = [&](const std::vector<Vector3>& points,
                                    const std::vector<int>&     snapped_v_ids,
                                    const std::array<int, 2>&   e,
                                    bool                        is_on_cut) {
        const int tet_vertices_size = mesh.tet_vertices.size();
        for (int i = points.size(); i > 0; i--) {
            mesh.tet_vertices[tet_vertices_size - i].is_on_boundary = true;
            mesh.tet_vertices[tet_vertices_size - i].is_on_cut      = is_on_cut;
        }

        for (int v_id : snapped_v_ids) {
            Scalar dis_2 = p_seg_squared_dist_3d(
              mesh.tet_vertices[v_id].pos, input_vertices[e[0]], input_vertices[e[1]]);
            if (dis_2 <= mesh.params.eps_2) {
                mesh.tet_vertices[v_id].is_on_boundary = true;
                mesh.tet_vertices[v_id].is_on_cut      = is_on_cut;
            }
        }
    };

    std::vector<std::vector<std::pair<int, int>>> covered_fs_infos(input_faces.size());
    for (int i = 0; i < track_surface_fs.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            for (int f_id : track_surface_fs[i][j])
                covered_fs_infos[f_id].push_back(std::make_pair(i, j));
        }
    }
    int cnt = 0;
    for (int I = 0; I < b_edge_infos.size(); I++) {
        const auto& e         = b_edge_infos[I].first;
        auto&       n_f_ids   = b_edge_infos[I].second;  /// it is sorted
        bool        is_on_cut = is_on_cut_edges[I];
        if (!is_again) {
            mesh.tet_vertices[e[0]].is_on_boundary = true;
            mesh.tet_vertices[e[1]].is_on_boundary = true;
            mesh.tet_vertices[e[0]].is_on_cut      = is_on_cut;
            mesh.tet_vertices[e[1]].is_on_cut      = is_on_cut;
        }

        for (int i = 0; i < n_f_ids.size(); i++) {
            if (!is_face_inserted[n_f_ids[i]]) {
                n_f_ids.erase(n_f_ids.begin() + i);
                i--;
                break;
            }
        }
        if (n_f_ids.empty()) {
            logger().warn("FAIL n_f_ids.empty()");
            continue;
        }

        std::vector<Vector3>              points;
        std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
        std::vector<int>                  snapped_v_ids;
        std::vector<std::array<int, 3>>   cut_fs;
        if (!insert_boundary_edges_get_intersecting_edges_and_points(covered_fs_infos,
                                                                     input_vertices,
                                                                     input_faces,
                                                                     e,
                                                                     n_f_ids,
                                                                     track_surface_fs,
                                                                     mesh,
                                                                     points,
                                                                     map_edge_to_intersecting_point,
                                                                     snapped_v_ids,
                                                                     cut_fs,
                                                                     is_again)) {
            for (int f_id : n_f_ids)
                is_face_inserted[f_id] = false;

            logger().warn("FAIL insert_boundary_edges_get_intersecting_edges_and_points");
            continue;
        }
        if (points.empty()) {  // everything snapped, nothing to subdivide
            record_boundary_info(points, snapped_v_ids, e, is_on_cut);
            cnt++;
            continue;
        }

        std::vector<int> cut_t_ids;
        for (const auto& m : map_edge_to_intersecting_point) {
            const auto&      tet_e = m.first;
            std::vector<int> tmp;
            set_intersection(
              mesh.tet_vertices[tet_e[0]].conn_tets, mesh.tet_vertices[tet_e[1]].conn_tets, tmp);
            cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
        }
        vector_unique(cut_t_ids);
        std::vector<bool> is_mark_surface(cut_t_ids.size(), false);
        CutMesh           empty_cut_mesh(mesh, Vector3(0, 0, 0), std::array<Vector3, 3>());
        std::vector<MeshTet>                         new_tets;
        std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
        std::vector<int>                             modified_t_ids;
        if (!subdivide_tets(-1,
                            mesh,
                            empty_cut_mesh,
                            points,
                            map_edge_to_intersecting_point,
                            track_surface_fs,
                            cut_t_ids,
                            is_mark_surface,
                            new_tets,
                            new_track_surface_fs,
                            modified_t_ids)) {
            bool is_inside_envelope = true;
            for (auto& f : cut_fs) {
                geo::index_t prev_facet = geo::NO_FACET;
                if (sample_triangle_and_check_is_out({{mesh.tet_vertices[f[0]].pos,
                                                       mesh.tet_vertices[f[1]].pos,
                                                       mesh.tet_vertices[f[2]].pos}},
                                                     mesh.params.dd,
                                                     mesh.params.eps_2,
                                                     tree,
                                                     prev_facet)) {
                    is_inside_envelope = false;
                    break;
                }
            }
            if (!is_inside_envelope) {
                for (int f_id : n_f_ids)
                    is_face_inserted[f_id] = false;
                for (auto& f : cut_fs) {
                    mark_known_surface_fs(f, KNOWN_NOT_SURFACE);
                }
                known_not_surface_fs.insert(
                  known_not_surface_fs.end(), cut_fs.begin(), cut_fs.end());
                logger().warn("FAIL subdivide_tets");
            }
            else {
                for (auto& f : cut_fs)
                    mark_known_surface_fs(f, KNOWN_SURFACE);
                known_surface_fs.insert(known_surface_fs.end(), cut_fs.begin(), cut_fs.end());
                logger().warn("SEMI-FAIL subdivide_tets");
            }

            continue;
        }

        push_new_tets(
          mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids);

        // points went on the end of mesh.tet_vertices, which is what this reads back.
        record_boundary_info(points, snapped_v_ids, e, is_on_cut);
        cnt++;
    }

    logger().info("uninsert boundary #e = {}/{}", b_edge_infos.size() - cnt, b_edge_infos.size());
}

void mark_surface_fs(const std::vector<Vector3>&                   input_vertices,
                     const std::vector<Vector3i>&                  input_faces,
                     const std::vector<int>&                       input_tags,
                     std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
                     const std::vector<std::array<int, 3>>&        known_surface_fs,
                     const std::vector<std::array<int, 3>>& known_not_surface_fs,
                     std::vector<std::array<int, 2>>&       b_edges,
                     Mesh&                                  mesh,
                     AABBWrapper&                           tree)
{
    std::vector<std::array<bool, 4>> is_visited(track_surface_fs.size(),
                                                {{false, false, false, false}});
    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        MeshTet& t = mesh.tets[t_id];
        for (int j = 0; j < 4; j++) {

            int ff_id    = -1;
            int opp_t_id = -1;
            int k        = -1;
            if (t.is_surface_fs[j] == KNOWN_SURFACE || t.is_surface_fs[j] == KNOWN_NOT_SURFACE) {
                if (!get_opp_face(mesh, t_id, j, opp_t_id, k)) {
                    t.is_surface_fs[j] = NOT_SURFACE;
                    continue;
                }
                is_visited[t_id][j]     = true;
                is_visited[opp_t_id][k] = true;
                if (t.is_surface_fs[j] == KNOWN_NOT_SURFACE ||
                    track_surface_fs[t_id][j].empty()) {
                    t.is_surface_fs[j]                   = NOT_SURFACE;
                    mesh.tets[opp_t_id].is_surface_fs[k] = NOT_SURFACE;
                    continue;
                }
                ff_id = track_surface_fs[t_id][j].front();
            }
            else {
                if (t.is_surface_fs[j] != NOT_SURFACE || is_visited[t_id][j])
                    continue;
                is_visited[t_id][j] = true;
                if (track_surface_fs[t_id][j].empty())
                    continue;

                auto& tp1_3d = mesh.tet_vertices[t[(j + 1) % 4]].pos;
                auto& tp2_3d = mesh.tet_vertices[t[(j + 2) % 4]].pos;
                auto& tp3_3d = mesh.tet_vertices[t[(j + 3) % 4]].pos;

                double eps_2 = (mesh.params.eps + mesh.params.eps_simplification) / 2;
                double dd    = (mesh.params.dd + mesh.params.dd_simplification) / 2;
                eps_2 *= eps_2;
                geo::index_t prev_facet = geo::NO_FACET;
                if (sample_triangle_and_check_is_out(
                      {{tp1_3d, tp2_3d, tp3_3d}}, dd, eps_2, tree, prev_facet))
                    continue;
                ff_id = track_surface_fs[t_id][j].front();

                if (!get_opp_face(mesh, t_id, j, opp_t_id, k)) {
                    t.is_surface_fs[j] = NOT_SURFACE;
                    continue;
                }
                is_visited[opp_t_id][k] = true;
            }

            MeshTet& opp = mesh.tets[opp_t_id];
            // The face carries opposite signs on its two sides.
            const auto set_sides = [&](int side, int opp_side) {
                t.is_surface_fs[j]   = side;
                opp.is_surface_fs[k] = opp_side;
            };

            t.surface_tags[j]   = input_tags[ff_id];
            opp.surface_tags[k] = input_tags[ff_id];

            auto& fv1   = input_vertices[input_faces[ff_id][0]];
            auto& fv2   = input_vertices[input_faces[ff_id][1]];
            auto& fv3   = input_vertices[input_faces[ff_id][2]];
            int   ori   = Predicates::orient_3d(fv1, fv2, fv3, mesh.tet_vertices[t[j]].pos);
            int opp_ori = Predicates::orient_3d(fv1, fv2, fv3, mesh.tet_vertices[opp[k]].pos);
            // The two corners straddle the face, so each side takes the sign it is on.
            if (ori != Predicates::ORI_ZERO && opp_ori != Predicates::ORI_ZERO &&
                ori != opp_ori) {
                set_sides(ori, opp_ori);
                continue;
            }
            // One corner is in the face's plane, so it takes the other's sign flipped.
            if (ori == Predicates::ORI_ZERO && opp_ori != Predicates::ORI_ZERO) {
                set_sides(-opp_ori, opp_ori);
                continue;
            }
            if (opp_ori == Predicates::ORI_ZERO && ori != Predicates::ORI_ZERO) {
                set_sides(ori, -ori);
                continue;
            }
            // Both corners are on the same side, or both are in the plane.
            Vector3 n = (fv2 - fv1).cross(fv3 - fv1);
            n.normalize();
            if (ori == Predicates::ORI_ZERO) {
                Vector3 nt = get_face_normal(mesh, t_id, j);
                nt.normalize();
                if (n.dot(nt) > 0)
                    set_sides(-1, 1);
                else
                    set_sides(1, -1);
            }
            else {
                // The nearer corner keeps its side and the further one is flipped.
                const Scalar dist     = n.dot(mesh.tet_vertices[t[j]].pos - fv1);
                const Scalar opp_dist = n.dot(mesh.tet_vertices[opp[k]].pos - fv1);
                if (dist < opp_dist)
                    set_sides(-ori, opp_ori);
                else
                    set_sides(ori, -opp_ori);
            }
        }
    }

    if (known_surface_fs.empty() && known_not_surface_fs.empty())
        return;

    std::vector<std::array<int, 2>> tmp_edges;
    for (const auto& f : known_surface_fs)
        push_tri_edges(f, tmp_edges);
    for (const auto& f : known_not_surface_fs)
        push_tri_edges(f, tmp_edges);
    vector_unique(tmp_edges);

    for (auto& e : tmp_edges) {
        int              cnt = 0;
        std::vector<int> n_t_ids;
        set_intersection(
          mesh.tet_vertices[e[0]].conn_tets, mesh.tet_vertices[e[1]].conn_tets, n_t_ids);
        for (int t_id : n_t_ids) {
            for (int j = 0; j < 4; j++) {
                if (mesh.tets[t_id][j] != e[0] && mesh.tets[t_id][j] != e[1]) {
                    if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE) {
                        cnt++;
                        if (cnt > 2)
                            break;
                    }
                }
                if (cnt > 2)
                    break;
            }
        }
        if (cnt == 2) {
            b_edges.push_back(e);
            mesh.tet_vertices[e[0]].is_on_boundary = true;
        }
    }
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::find_boundary_edges(
  const std::vector<Vector3>&                                   input_vertices,
  const std::vector<Vector3i>&                                  input_faces,
  const std::vector<bool>&                                      is_face_inserted,
  const std::vector<bool>&                                      old_is_face_inserted,
  std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos,
  std::vector<bool>&                                            is_on_cut_edges,
  std::vector<std::array<int, 2>>&                              b_edges)
{
    std::vector<std::array<int, 2>> edges;
    std::vector<std::vector<int>>   conn_tris(input_vertices.size());
    std::vector<std::vector<int>>   uninserted_conn_tris(input_vertices.size());
    for (int i = 0; i < input_faces.size(); i++) {
        if (!is_face_inserted[i]) {
            for (int j = 0; j < 3; j++)
                uninserted_conn_tris[input_faces[i][j]].push_back(i);
            continue;
        }
        const auto& f = input_faces[i];
        push_tri_edges(f, edges);
        for (int j = 0; j < 3; j++)
            conn_tris[f[j]].push_back(i);
    }
    vector_unique(edges);

    // An edge that has to be preserved is recorded with the faces it belongs to, and marked as
    // being on a cut when any of those faces has not been inserted yet.
    const auto record = [&](const std::array<int, 2>& e,
                            const std::vector<int>&   n12_f_ids,
                            bool                      needs_preserve,
                            const std::vector<int>&   uninserted_n12_f_ids) {
        b_edges.push_back(e);
        if (!needs_preserve)
            return;
        b_edge_infos.push_back(std::make_pair(e, n12_f_ids));
        is_on_cut_edges.push_back(!uninserted_n12_f_ids.empty());
    };

    for (const auto& e : edges) {
        std::vector<int> n12_f_ids;
        std::set_intersection(conn_tris[e[0]].begin(),
                              conn_tris[e[0]].end(),
                              conn_tris[e[1]].begin(),
                              conn_tris[e[1]].end(),
                              std::back_inserter(n12_f_ids));
        std::vector<int> uninserted_n12_f_ids;
        std::set_intersection(uninserted_conn_tris[e[0]].begin(),
                              uninserted_conn_tris[e[0]].end(),
                              uninserted_conn_tris[e[1]].begin(),
                              uninserted_conn_tris[e[1]].end(),
                              std::back_inserter(uninserted_n12_f_ids));

        bool needs_preserve = false;
        for (int f_id : n12_f_ids) {
            if (!old_is_face_inserted[f_id]) {
                needs_preserve = true;
                break;
            }
        }

        if (n12_f_ids.size() == 1) {  // open boundary
            record(e, n12_f_ids, needs_preserve, uninserted_n12_f_ids);
        }
        else {
            int f_id = n12_f_ids[0];
            int j    = 0;
            for (; j < 3; j++) {
                if ((input_faces[f_id][j] == e[0] && input_faces[f_id][(j + 1) % 3] == e[1]) ||
                    (input_faces[f_id][j] == e[1] && input_faces[f_id][(j + 1) % 3] == e[0]))
                    break;
            }

            const Vector3 n = input_face_normal(input_vertices, input_faces, f_id);
            const int     t = get_t(input_vertices[input_faces[f_id][0]],
                                input_vertices[input_faces[f_id][1]],
                                input_vertices[input_faces[f_id][2]]);

            // An edge between two faces that are not parallel does not need preserving.
            bool is_fine = false;
            for (int k = 0; k < n12_f_ids.size(); k++) {
                if (n12_f_ids[k] == f_id)
                    continue;
                const Vector3 n1 = input_face_normal(input_vertices, input_faces, n12_f_ids[k]);
                if (abs(n1.dot(n)) < 1 - SCALAR_ZERO) {
                    is_fine = true;
                    break;
                }
            }
            if (is_fine)
                continue;

            // Which side of the edge each neighbouring face's third corner falls on. Two of them
            // on different sides means the edge does not need preserving.
            const Vector2 e0_2d = to_2d(input_vertices[input_faces[f_id][j]], t);
            const Vector2 e1_2d = to_2d(input_vertices[input_faces[f_id][(j + 1) % 3]], t);
            is_fine             = false;
            int ori             = Predicates::ORI_UNKNOWN;
            for (int k = 0; k < n12_f_ids.size(); k++) {
                for (int r = 0; r < 3; r++) {
                    if (input_faces[n12_f_ids[k]][r] == input_faces[f_id][j] ||
                        input_faces[n12_f_ids[k]][r] == input_faces[f_id][(j + 1) % 3])
                        continue;
                    const int new_ori = Predicates::orient_2d(
                      e0_2d, e1_2d, to_2d(input_vertices[input_faces[n12_f_ids[k]][r]], t));
                    if (k == 0)
                        ori = new_ori;
                    else if (new_ori != ori)
                        is_fine = true;
                    break;
                }
                if (is_fine)
                    break;
            }
            if (is_fine)
                continue;

            record(e, n12_f_ids, needs_preserve, uninserted_n12_f_ids);
        }
    }
}

void floatTetWild::insert_triangles(const std::vector<Vector3>&  input_vertices,
                                    const std::vector<Vector3i>& input_faces,
                                    const std::vector<int>&      input_tags,
                                    Mesh&                        mesh,
                                    std::vector<bool>&           is_face_inserted,
                                    AABBWrapper&                 tree,
                                    bool                         is_again)
{
    std::vector<bool> old_is_face_inserted = is_face_inserted;

    logger().info("triangle insertion start, #f = {}, #v = {}, #t = {}",
                  input_faces.size(),
                  mesh.tet_vertices.size(),
                  mesh.tets.size());
    std::vector<std::array<std::vector<int>, 4>> track_surface_fs(mesh.tets.size());
    if (!is_again) {
        match_surface_fs(mesh, input_faces, is_face_inserted, track_surface_fs);
    }
    int cnt_matched = std::count(is_face_inserted.begin(), is_face_inserted.end(), true);
    logger().info(
      "matched #f = {}, uninserted #f = {}", cnt_matched, is_face_inserted.size() - cnt_matched);

    std::vector<int> sorted_f_ids;
    sort_input_faces(input_faces, mesh, sorted_f_ids);

    for (int f_id : sorted_f_ids) {
        if (is_face_inserted[f_id])
            continue;
        if (insert_one_triangle(f_id,
                                input_vertices,
                                input_faces,
                                mesh,
                                track_surface_fs,
                                tree,
                                is_again))
            is_face_inserted[f_id] = true;
    }
    logger().info(
      "insert_one_triangle * n done, #v = {}, #t = {}", mesh.tet_vertices.size(), mesh.tets.size());
    logger().info("uninserted #f = {}/{}",
                  std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
                  is_face_inserted.size() - cnt_matched);

    pair_track_surface_fs(mesh, track_surface_fs);
    logger().info("pair_track_surface_fs done");

    std::vector<std::array<int, 2>>                              b_edges1;
    std::vector<std::pair<std::array<int, 2>, std::vector<int>>> b_edge_infos;
    std::vector<bool>                                            is_on_cut_edges;
    find_boundary_edges(input_vertices,
                        input_faces,
                        is_face_inserted,
                        old_is_face_inserted,
                        b_edge_infos,
                        is_on_cut_edges,
                        b_edges1);
    logger().info("find_boundary_edges done");
    std::vector<std::array<int, 3>> known_surface_fs;
    std::vector<std::array<int, 3>> known_not_surface_fs;
    insert_boundary_edges(input_vertices,
                          input_faces,
                          b_edge_infos,
                          is_on_cut_edges,
                          track_surface_fs,
                          mesh,
                          tree,
                          is_face_inserted,
                          is_again,
                          known_surface_fs,
                          known_not_surface_fs);
    logger().info("uninserted #f = {}/{}",
                  std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
                  is_face_inserted.size() - cnt_matched);

    std::vector<std::array<int, 2>> b_edges2;
    mark_surface_fs(input_vertices,
                    input_faces,
                    input_tags,
                    track_surface_fs,
                    known_surface_fs,
                    known_not_surface_fs,
                    b_edges2,
                    mesh,
                    tree);
    logger().info("mark_surface_fs done");

    tree.init_tmp_b_mesh_and_tree(input_vertices, b_edges1, mesh, b_edges2);
    parallel_for(size_t(0), mesh.tets.size(), [&](size_t i) {
        auto& t = mesh.tets[i];
        if (!t.is_removed) {
            t.quality = get_quality(mesh, t);
        }
    });

    if (std::count(is_face_inserted.begin(), is_face_inserted.end(), false) == 0)
        mesh.is_input_all_inserted = true;
    logger().info("#b_edge1 = {}, #b_edges2 = {}", b_edges1.size(), b_edges2.size());
}
