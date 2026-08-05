// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/MeshImprovement.h>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/EdgeSplitting.h>
#include <floattetwild/EdgeCollapsing.h>
#include <floattetwild/EdgeSwapping.h>
#include <floattetwild/VertexSmoothing.h>
#include <floattetwild/CSGTreeParser.hpp>

#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/Statistics.h>

#include <floattetwild/Logger.hpp>

#include <floattetwild/Timer.h>
#include <floattetwild/ParallelFor.hpp>
#include <floattetwild/geo_kd_tree.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/MeshCleanup.hpp>
#include <floattetwild/FastWindingNumber.h>

namespace floatTetWild {
namespace {

// Defined at the end of this namespace, next to the rest of the pass driving.
void run_passes(Mesh &mesh, AABBWrapper &tree, const std::array<int, 4> &ops);

// One row per live tet, holding its barycentre. This is where the winding number gets evaluated.
MatrixXd tet_barycenters(const Mesh &mesh) {
    MatrixXd C(mesh.get_t_num(), 3);
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += tet_pos(mesh, i, j).cast<double>();
        C.row(index) /= 4.0;
        index++;
    }
    return C;
}

// The winding number of every query point in C against the surface (V, F), by the approximation
// in "Fast Winding Numbers for Soups and Clouds" [Barill et al. 2018]. Flipping the faces flips
// the sign convention, which is the fallback when the input turned out to be inside out.
//
// libigl also exposed the two halves of this, building the hierarchy once and querying it many
// times. Nothing here reuses a hierarchy, so the two are one call.
MatrixXd winding_numbers(const MatrixXs &V, const MatrixXi &F, const MatrixXd &C,
                         bool invert_faces) {
    using FastWindingNumber::Vec3f;

    std::vector<Vec3f> points(V.rows());
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < 3; j++)
            points[i][j] = V(i, j);
    }
    std::vector<int> faces(F.size());
    for (int f = 0; f < F.rows(); f++) {
        faces[3 * f] = F(f, 0);
        faces[3 * f + 1] = F(f, invert_faces ? 2 : 1);
        faces[3 * f + 2] = F(f, invert_faces ? 1 : 2);
    }
    std::vector<Vec3f> queries(C.rows());
    for (int i = 0; i < C.rows(); i++) {
        for (int j = 0; j < 3; j++)
            queries[i][j] = C(i, j);
    }

    const std::vector<float> angles = FastWindingNumber::solid_angles(points, faces, queries);

    // igl::PI, and the 4 pi a solid angle is divided by to become a winding number.
    constexpr double Pi = 3.1415926535897932384626433832795;

    MatrixXd W(C.rows(), 1);
    for (int i = 0; i < C.rows(); i++)
        W(i) = angles[i] / (4.0 * Pi);
    return W;
}

// The same, against the mesh's own tracked surface carrying tag c_id.
MatrixXd tracked_surface_wn(Mesh &mesh, const MatrixXd &C, int c_id) {
    MatrixXs V;
    MatrixXi F;
    get_tracked_surface(mesh, V, F, c_id);
    return winding_numbers(V, F, C, false);
}

// Which live tets have their barycentre inside the given surface, one entry per tet. A tet that
// was already removed is false. When that left nothing inside, the faces are flipped and the whole
// thing is redone, which is what an inside-out input needs.
std::vector<bool> tets_inside_surface(const Mesh &mesh, const MatrixXs &V, const MatrixXi &F) {
    const MatrixXd C = tet_barycenters(mesh);
    for (int pass = 0; pass < 2; pass++) {
        const MatrixXd W = winding_numbers(V, F, C, pass == 1);

        std::vector<bool> is_inside(mesh.tets.size(), false);
        int index = 0;
        int n_tets = 0;
        for (size_t t_id = 0; t_id < mesh.tets.size(); t_id++) {
            if (mesh.tets[t_id].is_removed)
                continue;
            if (W(index++) > 0.5) {
                is_inside[t_id] = true;
                n_tets++;
            }
        }
        if (n_tets > 0)
            return is_inside;

        if (pass == 0)
            logger().debug("Empty mesh trying to reverse the faces");
        else
            logger().error("Empty mesh, problem with inverted faces");
    }
    return std::vector<bool>(mesh.tets.size(), false);
}

// The faces carried by exactly one tet, as {v0, v1, v2, t_id, j}: the corners sorted, then the tet
// the face belongs to and its local index there. The result is sorted by corner triple.
// skip picks which tets take part, is_outside for the passes that run after mark_outside and
// is_removed for the rest.
std::vector<std::array<int, 5>> boundary_faces(const Mesh &mesh,
                                               bool MeshTet::*skip) {
    std::vector<std::array<int, 5>> faces;
    faces.reserve(mesh.tets.size() * 4);
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        const auto &t = mesh.tets[i];
        if (t.*skip)
            continue;
        for (int j = 0; j < 4; j++) {
            const std::array<int, 3> f = sorted_face(t, j);
            faces.push_back({{f[0], f[1], f[2], int(i), j}});
        }
    }
    const auto corners = [](const std::array<int, 5> &f) {
        return std::make_tuple(f[0], f[1], f[2]);
    };
    std::sort(faces.begin(), faces.end(),
              [&](const std::array<int, 5> &a, const std::array<int, 5> &b) {
                  return corners(a) < corners(b);
              });

    std::vector<std::array<int, 5>> b_faces;
    bool is_boundary = true;
    for (size_t i = 0; i + 1 < faces.size(); i++) {
        if (corners(faces[i]) == corners(faces[i + 1]))
            is_boundary = false;
        else {
            if (is_boundary)
                b_faces.push_back(faces[i]);
            is_boundary = true;
        }
    }
    if (!faces.empty() && is_boundary)
        b_faces.push_back(faces.back());
    return b_faces;
}

// Splits ids into connected groups. neighbours(id, reached) calls reached(n_id) for everything
// adjacent to id, and anything outside ids is ignored. More than one group means the mesh is
// non-manifold at whatever the ids were collected around.
template <typename Neighbours>
std::vector<std::vector<int>> group_connected(const std::vector<int> &ids,
                                              Neighbours &&neighbours) {
    std::map<int, bool> is_visited;
    for (int id : ids)
        is_visited[id] = false;

    std::vector<std::vector<int>> groups;
    for (int id : ids) {
        if (is_visited[id])
            continue;
        is_visited[id] = true;

        groups.emplace_back();
        std::queue<int> queue;
        queue.push(id);
        while (!queue.empty()) {
            const int cur_id = queue.front();
            queue.pop();
            groups.back().push_back(cur_id);

            neighbours(cur_id, [&](int n_id) {
                const auto it = is_visited.find(n_id);
                if (it != is_visited.end() && !it->second) {
                    queue.push(it->first);
                    it->second = true;
                }
            });
        }
    }
    return groups;
}

// group_connected over tets, where two are connected when they share a face that crosses(t_id, j)
// accepts.
template <typename Crosses>
std::vector<std::vector<int>> group_connected_tets(const Mesh &mesh, const std::vector<int> &t_ids,
                                                   Crosses &&crosses) {
    return group_connected(t_ids, [&](int t_id, auto &&reached) {
        for (int j = 0; j < 4; j++) {
            if (crosses(t_id, j))
                reached(get_opp_t_id(mesh, t_id, j));
        }
    });
}

// The vertices the given faces use, each once, ascending.
template <typename Face>
std::vector<int> face_vertices(const std::vector<Face> &faces) {
    std::vector<int> v_ids;
    v_ids.reserve(faces.size() * 3);
    for (const auto &f : faces) {
        for (int j = 0; j < 3; j++)
            v_ids.push_back(f[j]);
    }
    vector_unique(v_ids);
    return v_ids;
}

// Drops the sizing field back to uniform.
void reset_sizing_field(Mesh &mesh) {
    for (auto &v : mesh.tet_vertices) {
        if (!v.is_removed)
            v.sizing_scalar = 1;
    }
}

// Where each entry ends up once the removed ones before end_id are dropped. Everything from
// end_id on is kept, which is what makes the compaction partial.
template <typename T>
std::vector<int> compaction_map(const std::vector<T> &v, int end_id) {
    std::vector<int> map_ids(v.size(), -1);
    int cnt = 0;
    for (int i = 0; i < end_id; i++) {
        if (!v[i].is_removed)
            map_ids[i] = cnt++;
    }
    for (int i = end_id; i < int(v.size()); i++)
        map_ids[i] = cnt++;
    return map_ids;
}

template <typename T>
void erase_removed(std::vector<T> &v, int end_id) {
    v.erase(std::remove_if(v.begin(), v.begin() + end_id,
                           [](const T &e) { return e.is_removed; }),
            v.begin() + end_id);
}

// A vertex is gone once every tet around it is.
void drop_isolated_vertices(Mesh &mesh) {
    for (auto &v : mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        bool is_remove = true;
        for (int t_id : v.conn_tets) {
            if (!mesh.tets[t_id].is_removed) {
                is_remove = false;
                break;
            }
        }
        v.is_removed = is_remove;
    }
}

// Steps of the optimization loop and the manifold fix-ups. Everything above and below is
// reached only from this file; MeshImprovement.h carries the entry points.

// Compacts the leading 70% of both lists, which is enough to keep a big mesh's free slots from
// piling up without paying to renumber the whole thing.
void cleanup_empty_slots(Mesh &mesh) {
    if (mesh.tets.size() < 9e5)
        return;
    const double percentage = 0.7;
    const int v_end_id = mesh.tet_vertices.size() * percentage;
    const int t_end_id = mesh.tets.size() * percentage;
    const std::vector<int> map_v_ids = compaction_map(mesh.tet_vertices, v_end_id);
    const std::vector<int> map_t_ids = compaction_map(mesh.tets, t_end_id);

    erase_removed(mesh.tet_vertices, v_end_id);
    erase_removed(mesh.tets, t_end_id);

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        for (auto &t_id: v.conn_tets)
            t_id = map_t_ids[t_id];
    }
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            t[j] = map_v_ids[t[j]];
    }
}

bool update_scaling_field(Mesh &mesh, Scalar max_energy) {

    bool is_hit_min_edge_length = false;

    Scalar radius0 = mesh.params.ideal_edge_length * 1.8;//increasing the radius would increase the #v in output

    // Not static: it reads mesh.params, and tetrahedralization() can be called more than once in
    // a process with different parameters.
    const Scalar stop_filter_energy = mesh.params.stop_energy * 0.8;

    const Scalar filter_energy = std::min(std::max(max_energy / 100, stop_filter_energy), Scalar(100));

    const Scalar recover = 1.5;
    std::vector<Scalar> scale_multipliers(mesh.tet_vertices.size(), recover);
    Scalar refine_scale = 0.5;
    Scalar min_refine_scale = mesh.params.min_edge_len_rel / mesh.params.ideal_edge_length_rel;

    const int N = -int(std::log2(min_refine_scale) - 1);
    std::vector<std::vector<int>> v_ids(N, std::vector<int>());
    for (int i = 0; i < mesh.tet_vertices.size(); i++) {
        auto &v = mesh.tet_vertices[i];
        if (v.is_removed)
            continue;

        bool is_refine = false;
        for (int t_id: v.conn_tets) {
            if (mesh.tets[t_id].quality > filter_energy) {
                is_refine = true;
                break;
            }
        }
        if (!is_refine)
            continue;

        int n = -int(std::log2(v.sizing_scalar) - 0.5);
        if (n >= N)
            n = N - 1;
        v_ids[n].push_back(i);
    }

    for (int n = 0; n < N; n++) {
        if (v_ids[n].size() == 0)
            continue;

        Scalar radius = radius0 / std::pow(2, n);

        std::unordered_set<int> is_visited;
        std::queue<int> v_queue;

        std::vector<double> pts;//the kd-tree needs double []
        pts.reserve(v_ids[n].size() * 3);
        for (int v_id : v_ids[n]) {
            const Vector3 &p = mesh.tet_vertices[v_id].pos;
            pts.insert(pts.end(), {p[0], p[1], p[2]});

            v_queue.push(v_id);
            is_visited.insert(v_id);
            scale_multipliers[v_id] = refine_scale;
        }
        geo::BalancedKdTree nnsearch;
        nnsearch.set_points(int(v_ids[n].size()), pts.data());

        while (!v_queue.empty()) {
            int v_id = v_queue.front();
            v_queue.pop();

            for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
                for (int j = 0; j < 4; j++) {
                    const int n_v_id = mesh.tets[t_id][j];
                    if (is_visited.find(n_v_id) != is_visited.end())
                        continue;
                    const Vector3 &np = mesh.tet_vertices[n_v_id].pos;

                    const double p[3] = {np[0], np[1], np[2]};
                    Scalar dis = sqrt(nnsearch.nearest_sq_dist(p));

                    if (dis < radius) {
                        v_queue.push(n_v_id);
                        Scalar new_ss = (dis / radius) * (1 - refine_scale) + refine_scale;
                        if (new_ss < scale_multipliers[n_v_id])
                            scale_multipliers[n_v_id] = new_ss;
                    }
                    is_visited.insert(n_v_id);
                }
            }
        }
    }

    for (int i=0;i< mesh.tet_vertices.size();i++) {
        auto& v = mesh.tet_vertices[i];
        if (v.is_removed)
            continue;
        Scalar new_scale = v.sizing_scalar * scale_multipliers[i];
        if (new_scale > 1)
            v.sizing_scalar = 1;
        else if (new_scale < min_refine_scale) {
            is_hit_min_edge_length = true;
            v.sizing_scalar = min_refine_scale;
        } else
            v.sizing_scalar = new_scale;
    }

    return is_hit_min_edge_length;
}

int get_max_p(const Mesh &mesh)
{
    const Scalar scaling = 1.0 / (mesh.params.bbox_max - mesh.params.bbox_min).maxCoeff();
    const Scalar B = 3;
    const int p_ref = 1;

    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);
    Scalar h_ref = 0;
    for(const auto &e : edges){
        const Vector3 edge = (mesh.tet_vertices[e[0]].pos - mesh.tet_vertices[e[1]].pos)*scaling;
        h_ref += edge.norm();
    }
    h_ref /= edges.size();

    const Scalar rho_ref = sqrt(6.)/12.*h_ref;
    const Scalar sigma_ref = rho_ref / h_ref;

    int max_p = 1;

    for(const auto &t : mesh.tets)
    {
        if(t.is_removed)
            continue;

        const auto &v0 = mesh.tet_vertices[t[0]].pos;
        const auto &v1 = mesh.tet_vertices[t[1]].pos;
        const auto &v2 = mesh.tet_vertices[t[2]].pos;
        const auto &v3 = mesh.tet_vertices[t[3]].pos;

        std::array<Vector3, 6> e;
        e[0] = (v0 - v1) * scaling;
        e[1] = (v1 - v2) * scaling;
        e[2] = (v2 - v0) * scaling;

        e[3] = (v0 - v3) * scaling;
        e[4] = (v1 - v3) * scaling;
        e[5] = (v2 - v3) * scaling;

        const Scalar S = (e[0].cross(e[1]).norm() + e[0].cross(e[4]).norm() + e[4].cross(e[1]).norm() + e[2].cross(e[5]).norm()) / 2;
        const Scalar V = std::abs(e[3].dot(e[2].cross(-e[0])))/6;
        const Scalar rho = 3 * V / S;
        Scalar h = e[0].norm();
        for (int i = 1; i < 6; i++)
            h = std::max(h, e[i].norm());

        const Scalar sigma = rho / h;

        const Scalar ptmp = (std::log(B*std::pow(h_ref, p_ref + 1)*sigma*sigma/sigma_ref/sigma_ref) - std::log(h))/std::log(h);
        const int p = (int)std::round(ptmp);
        max_p = std::max(p, max_p);
    }

    return max_p;
}

void apply_coarsening(Mesh& mesh, AABBWrapper& tree) {
    mesh.is_coarsening = true;

    reset_sizing_field(mesh);

    int tets_size = mesh.get_t_num();
    int stop_size = tets_size * 0.001;
    for (int i = 0; i < 20; i++) {
        run_passes(mesh, tree, {{0, 1, 1, 0}});
        int new_size = mesh.get_t_num();
        if (abs(new_size - tets_size) < stop_size)
            break;
        tets_size = new_size;
    }

    mesh.is_coarsening = false;
}

void mark_outside(Mesh& mesh){
    MatrixXs V;
    MatrixXi F;
    get_tracked_surface(mesh, V, F);
    const std::vector<bool> is_inside = tets_inside_surface(mesh, V, F);
    for (size_t t_id = 0; t_id < mesh.tets.size(); ++t_id)
        mesh.tets[t_id].is_outside = !is_inside[t_id];
}

void untangle(Mesh &mesh) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;
    static const Scalar zero_area = 1e2 * SCALAR_ZERO_2;
    static const std::vector<std::array<int, 4>> face_pairs = {{{0, 1, 2, 3}},
                                                               {{0, 2, 1, 3}},
                                                               {{0, 3, 1, 2}}};

    for (int t_id = 0; t_id < tets.size(); t_id++) {
        auto &t = tets[t_id];
        if (t.is_removed)
            continue;
        if (t.quality < 1e10)
            continue;

        const auto is_surface = [&](int j) { return t.is_surface_fs[j] != NOT_SURFACE; };

        // Drop the surface tag from face j of this tet and from the same face on its neighbour.
        const auto unmark_surface = [&](int j) {
            t.is_surface_fs[j] = NOT_SURFACE;
            int opp_t_id, k;
            if (get_opp_face(mesh, t_id, j, opp_t_id, k))
                tets[opp_t_id].is_surface_fs[k] = NOT_SURFACE;
        };

        int cnt_on_surface = 0;
        bool has_degenerate_face = false;
        std::array<double, 4> areas;
        double max_area = 0;
        int max_j = -1;
        for (int j = 0; j < 4; j++) {
            if (is_surface(j))
                cnt_on_surface++;
            areas[j] = get_area(tet_vertices[t[(j + 1) % 4]].pos,
                                tet_vertices[t[(j + 2) % 4]].pos,
                                tet_vertices[t[(j + 3) % 4]].pos);
            if (areas[j] < zero_area)
                has_degenerate_face = true;
            if (areas[j] > max_area) {
                max_area = areas[j];
                max_j = j;
            }
        }
        if (cnt_on_surface == 0)
            continue;

        // The tet has flattened onto its largest face, so nothing else about it is surface.
        if (has_degenerate_face) {
            if (is_surface(max_j) && max_area > zero_area) {
                for (int j = 0; j < 4; j++) {
                    if (j != max_j)
                        unmark_surface(j);
                }
            }
            continue;
        }
        if (cnt_on_surface < 2)
            continue;

        // The three smaller faces exactly cover the largest one, same thing.
        if (std::abs(max_area - areas[(max_j + 1) % 4] - areas[(max_j + 2) % 4] -
                     areas[(max_j + 3) % 4]) < zero_area) {
            if (!is_surface(max_j))
                continue;
            for (int j = 0; j < 4; j++) {
                if (j != max_j && is_surface(j))
                    unmark_surface(j);
            }
            continue;
        }

        // Otherwise the tet may have flattened onto a pair of opposite edges. fp[0] and fp[1] are
        // the two corners on one side of the edge through fp[2] and fp[3]; the two sides folding
        // together is what makes one of the two face pairs redundant.
        for (const auto &fp: face_pairs) {
            std::array<Vector3, 2> ns;
            auto &p1 = tet_vertices[t[fp[2]]].pos;
            auto &p2 = tet_vertices[t[fp[3]]].pos;
            Vector3 v = (p2 - p1).normalized();
            for (int j = 0; j < 2; j++) {
                auto &p = tet_vertices[t[fp[j]]].pos;
                Vector3 q = p1 + ((p - p1).dot(v)) * v;
                ns[j] = p - q;
            }
            if (ns[0].dot(ns[1]) > 0)
                continue;

            if (std::abs(areas[fp[0]] + areas[fp[1]] - areas[fp[2]] - areas[fp[3]]) > zero_area)
                continue;

            std::array<int, 2> js = {{-1, -1}};
            if (is_surface(fp[0]) && is_surface(fp[1]))
                js = {{fp[2], fp[3]}};
            else if (is_surface(fp[2]) && is_surface(fp[3]))
                js = {{fp[0], fp[1]}};

            for (int j: js) {
                if (j >= 0 && is_surface(j))
                    unmark_surface(j);
            }
            break;
        }
    }
}

void manifold_edges(Mesh& mesh) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    auto split = [&](int v1_id, int v2_id, std::vector<int> &old_t_ids, std::vector<int> &new_t_ids) {
        MeshVertex new_v;
        new_v.pos = (tet_vertices[v1_id].pos + tet_vertices[v2_id].pos) / 2;
        const int v_id = get_new_vertex_slot(mesh, new_v);

        set_intersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, old_t_ids);
        for (int t_id: old_t_ids) {
            for (int j = 0; j < 4; j++) {
                if (tets[t_id][j] == v1_id || tets[t_id][j] == v2_id) {
                    if (is_inverted(mesh, t_id, j, new_v.pos)) {
                        tet_vertices[v_id].is_removed = true;
                        return -1;
                    }
                }
            }
        }

        split_tets_at(mesh, v_id, v1_id, v2_id, old_t_ids, new_t_ids);
        for (int i = 0; i < old_t_ids.size(); i++) {
            tets[new_t_ids[i]].quality = get_quality(mesh, new_t_ids[i]);
            tets[old_t_ids[i]].quality = get_quality(mesh, old_t_ids[i]);
        }

        return v_id;
    };

    const std::vector<std::array<int, 5>> b_faces = boundary_faces(mesh, &MeshTet::is_removed);

    std::vector<std::array<int, 2>> b_edges;
    for (const auto &f : b_faces)
        push_tri_edges(f, b_edges);
    vector_unique(b_edges);

    std::queue<std::array<int, 2>> edge_queue;
    for (auto &e: b_edges)
        edge_queue.push(e);

    while (!edge_queue.empty()) {
        auto e = edge_queue.front();
        edge_queue.pop();

        std::vector<int> n_t_ids;
        set_intersection(tet_vertices[e[0]].conn_tets, tet_vertices[e[1]].conn_tets, n_t_ids);
        if (n_t_ids.empty())
            continue;

        // The two faces of a tet that contain the edge are the ones a walk around it may cross.
        const std::vector<std::vector<int>> tet_groups =
            group_connected_tets(mesh, n_t_ids, [&](int t0_id, int j) {
                return tets[t0_id][j] != e[0] && tets[t0_id][j] != e[1];
            });
        if (tet_groups.size() < 2)
            continue;

        logger().debug("find non-manifold edge {} {}", e[0], e[1]);

        std::vector<int> new_t_ids;
        std::vector<int> old_t_ids;
        int v_id = split(e[0], e[1], old_t_ids, new_t_ids);
        if (v_id < 0)
            continue;
        std::map<int, int> old_t_ids_map;
        for (int i = 0; i < old_t_ids.size(); i++) {
            old_t_ids_map[old_t_ids[i]] = i;
        }

        for (int i = 0; i < tet_groups.size(); i++) {
            int dup_v_id = v_id;
            if(i > 0) {
                tet_vertices.push_back(tet_vertices[v_id]);
                dup_v_id = tet_vertices.size() - 1;
            }
            tet_vertices[dup_v_id].conn_tets.clear();
            for (int old_t_id: tet_groups[i]) {
                int new_t_id = new_t_ids[old_t_ids_map[old_t_id]];
                int j = tets[old_t_id].find(v_id);
                tets[old_t_id][j] = dup_v_id;
                j = tets[new_t_id].find(v_id);
                tets[new_t_id][j] = dup_v_id;
                tet_vertices[dup_v_id].conn_tets.push_back(old_t_id);
                tet_vertices[dup_v_id].conn_tets.push_back(new_t_id);
            }
        }

        old_t_ids.insert(old_t_ids.end(), new_t_ids.begin(), new_t_ids.end());
        std::vector<std::array<int, 2>> new_edges;
        collect_tet_edges(mesh, old_t_ids, new_edges);
        for (auto &ne : new_edges)
            edge_queue.push(ne);
    }
}

void manifold_vertices(Mesh& mesh){
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    const std::vector<std::array<int, 5>> b_faces = boundary_faces(mesh, &MeshTet::is_removed);

    for (int b_v_id : face_vertices(b_faces)) {
        std::vector<int> n_t_ids;
        for (int t_id: tet_vertices[b_v_id].conn_tets) {
            if (!tets[t_id].is_removed)
                n_t_ids.push_back(t_id);
        }

        // The three faces of a tet that contain the vertex are the ones a walk around it may cross.
        const std::vector<std::vector<int>> tet_groups =
            group_connected_tets(mesh, n_t_ids, [&](int t0_id, int j) {
                return tets[t0_id][j] != b_v_id;
            });
        if (tet_groups.size() < 2)
            continue;

        logger().debug("find non-manifold vertex {}", b_v_id);

        for (int i = 1; i < tet_groups.size(); i++) {
            tet_vertices.push_back(tet_vertices[b_v_id]);
            for (int t_id: tet_groups[i]) {
                int j = tets[t_id].find(b_v_id);
                tets[t_id][j] = tet_vertices.size() - 1;
            }
        }
    }
}

// Re-derives the per-vertex surface and bbox marks from the per-face ones, after an operation
// that could have moved a face off either.
void init(Mesh &mesh) {
    logger().info("initializing...");

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.is_on_surface = false;
        v.is_on_bbox = false;
    }
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            const bool on_surface = t.is_surface_fs[j] <= 0;
            const bool on_bbox = t.is_bbox_fs[j] != NOT_BBOX;
            if (!on_surface && !on_bbox)
                continue;
            for (int k = 1; k < 4; k++) {
                auto &v = mesh.tet_vertices[t[(j + k) % 4]];
                v.is_on_surface = v.is_on_surface || on_surface;
                v.is_on_bbox = v.is_on_bbox || on_bbox;
            }
        }
    }

    //todo: after all faces are inserted, mark true is_on_boundary
}

void run_passes(Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops){
    // Smoothing is the only one that does not untangle first.
    struct Pass {
        const char *name;
        int id;
        bool untangle_first;
        void (*run)(Mesh &, AABBWrapper &);
    };
    static const std::array<Pass, 4> passes = {{
        {"edge splitting", StateInfo::splitting_id, true,
         [](Mesh &m, AABBWrapper &t) { edge_splitting(m, t); }},
        {"edge collapsing", StateInfo::collapsing_id, true,
         [](Mesh &m, AABBWrapper &t) { edge_collapsing(m, t); }},
        {"edge swapping", StateInfo::swapping_id, true,
         [](Mesh &m, AABBWrapper &) { edge_swapping(m); }},
        {"vertex smoothing", StateInfo::smoothing_id, false,
         [](Mesh &m, AABBWrapper &t) { vertex_smoothing(m, t); }},
    }};

    Timer igl_timer;
    for (int p = 0; p < 4; p++) {
        for (int i = 0; i < ops[p]; i++) {
            igl_timer.start();
            logger().info(passes[p].name);
            if (passes[p].untangle_first)
                untangle(mesh);
            passes[p].run(mesh, tree);
            const double time = igl_timer.getElapsedTimeInSec();
            Scalar max_energy, avg_energy;
            mesh.get_max_avg_energy(max_energy, avg_energy);
            stats().push_back({passes[p].id, time, mesh.get_v_num(), mesh.get_t_num(),
                               max_energy, avg_energy});
        }
    }
}

// One round: the local operation passes ops asks for, then, if reinsert_triangles, another go at
// the input faces that are still missing.
void operation(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
        Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops, bool reinsert_triangles) {
    run_passes(mesh, tree, ops);

    if (!reinsert_triangles || mesh.is_input_all_inserted)
        return;

    Timer igl_timer;
    igl_timer.start();
    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, true);
    init(mesh);
    stats().push_back({StateInfo::cutting_id, igl_timer.getElapsedTimeInSec(),
                       mesh.get_v_num(), mesh.get_t_num(),
                       mesh.get_max_energy(), mesh.get_avg_energy(),
                       int(std::count(is_face_inserted.begin(), is_face_inserted.end(), false))});

    // A closed input that is now fully inserted has no cut left to preserve. Otherwise the
    // boundary marks only survive where they still sit inside the boundary envelope.
    const bool drop_all = mesh.is_input_all_inserted && mesh.is_closed;
    for (auto &v : mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        if (!drop_all) {
            if (!v.is_on_boundary)
                continue;
            geo::index_t prev_facet;
            if (!tree.is_out_tmp_b_envelope(v.pos, mesh.params.eps_2, prev_facet))
                continue;
        }
        v.is_on_boundary = false;
        v.is_on_cut = false;
    }
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::optimization(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
        Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops) {
    init(mesh);

    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{0, 1, 0, 0}}, false);
    cleanup_empty_slots(mesh);
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{0, 0, 0, 0}}, true);

    const int M = 5;
    const int N = 5;

    int it_after_al_inserted = 0;
    bool is_just_after_update = false;
    bool is_hit_min_edge_length = false;
    std::vector<std::array<Scalar, 2>> quality_queue;
    int cnt_increase_epsilon = mesh.params.stage - 1;
    const auto enlarge_envelope = [&]() {
        mesh.params.eps += mesh.params.eps_delta;
        mesh.params.eps_2 = mesh.params.eps * mesh.params.eps;
        cnt_increase_epsilon--;
        logger().info("enlarge envelope, eps = {}", mesh.params.eps);
    };
    for (int it = 0; it < mesh.params.max_its; it++) {
        if (mesh.is_input_all_inserted)
            it_after_al_inserted++;

        Scalar max_energy, avg_energy;
        mesh.get_max_avg_energy(max_energy, avg_energy);
        if (max_energy <= mesh.params.stop_energy && mesh.is_input_all_inserted)
            break;

        if (mesh.params.stop_p > 0) {
            int p = get_max_p(mesh);
            if (p <= mesh.params.stop_p && mesh.is_input_all_inserted)
                break;
        }

        logger().info("pass {}", it);
        // Every third pass also tries to insert what is still missing.
        operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, ops,
                  it % 3 == 2);

        if (it > mesh.params.max_its / 4 && max_energy > 1e3) {
            if (cnt_increase_epsilon > 0 && cnt_increase_epsilon == mesh.params.stage - 1)
                enlarge_envelope();
        }

        Scalar new_max_energy, new_avg_energy;
        mesh.get_max_avg_energy(new_max_energy, new_avg_energy);
        if (!is_just_after_update) {
            if (max_energy - new_max_energy < 5e-1 && (avg_energy - new_avg_energy) / avg_energy < 0.1) {
                is_hit_min_edge_length = update_scaling_field(mesh, new_max_energy) || is_hit_min_edge_length;
                is_just_after_update = true;
                if (cnt_increase_epsilon > 0)
                    enlarge_envelope();
            }
        } else
            is_just_after_update = false;

        quality_queue.push_back(std::array<Scalar, 2>({{new_max_energy, new_avg_energy}}));
        if (is_hit_min_edge_length && mesh.is_input_all_inserted && it_after_al_inserted > M && it > M + N) {
            if (quality_queue[it][0] - quality_queue[it - N][0] >= SCALAR_ZERO
                && quality_queue[it][1] - quality_queue[it - N][1] >= SCALAR_ZERO)
                break;
        }
    }

    logger().info("postprocessing");
    reset_sizing_field(mesh);
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{0, 1, 0, 0}}, false);

    if(mesh.params.coarsen){
        apply_coarsening(mesh, tree);
    }
}

void floatTetWild::get_tracked_surface(Mesh& mesh, MatrixXs &V_sf, MatrixXi &F_sf, int c_id) {
    auto &tet_vertices = mesh.tet_vertices;

    const auto is_tracked = [c_id](const MeshTet &t, int j) {
        return t.is_surface_fs[j] <= 0 && t.surface_tags[j] == c_id;
    };

    int cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (is_tracked(t, j))
                cnt++;
        }
    }

    V_sf.resize(cnt * 3, 3);
    F_sf.resize(cnt, 3);
    cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (is_tracked(t, j)) {
                const std::array<int, 3> f = face_corners(t, j);
                for (int k = 0; k < 3; k++)
                    V_sf.row(cnt * 3 + k) = tet_vertices[f[k]].pos;
                // Wind the face away from the corner it is opposite to.
                const bool flip =
                  Predicates::orient_3d(tet_vertices[f[0]].pos,
                                        tet_vertices[f[1]].pos,
                                        tet_vertices[f[2]].pos,
                                        tet_vertices[t[j]].pos) == Predicates::ORI_POSITIVE;
                F_sf.row(cnt) = Vector3i(cnt * 3, cnt * 3 + (flip ? 2 : 1), cnt * 3 + (flip ? 1 : 2));
                cnt++;
            }
        }
    }

    MatrixXd V;
    MatrixXi F;
    MatrixXi _1, _2;
    remove_duplicate_vertices(V_sf, F_sf, -1, V, _1, _2, F);
    V_sf = V;
    F_sf.resize(0, 3);
    bfs_orient(F, F_sf, _1);
}

void floatTetWild::correct_tracked_surface_orientation(Mesh &mesh, AABBWrapper& tree){
    std::vector<std::array<bool, 4>> is_visited(mesh.tets.size(), {{false, false, false, false}});
    for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
        auto &t = mesh.tets[t_id];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] == NOT_SURFACE || is_visited[t_id][j])
                continue;
            is_visited[t_id][j] = true;
            int opp_t_id, k;
            if (!get_opp_face(mesh, t_id, j, opp_t_id, k)) {
                t.is_surface_fs[j] = NOT_SURFACE;
                continue;
            }
            is_visited[opp_t_id][k] = true;
            Vector3 c = (tet_pos(mesh, t_id, (j + 1) % 4) + tet_pos(mesh, t_id, (j + 2) % 4) +
                         tet_pos(mesh, t_id, (j + 3) % 4)) / 3;
            int f_id = tree.get_nearest_face_sf(c);
            const auto &fv1 = tree.sf_mesh.point(tree.sf_mesh.facet_vertex(f_id, 0));
            const auto &fv2 = tree.sf_mesh.point(tree.sf_mesh.facet_vertex(f_id, 1));
            const auto &fv3 = tree.sf_mesh.point(tree.sf_mesh.facet_vertex(f_id, 2));
            auto nf = geo::cross((fv2 - fv1), (fv3 - fv1));
            const Vector3 n(nf[0], nf[1], nf[2]);
            // The face carries opposite signs on its two sides.
            const char side = get_face_normal(mesh, t_id, j).dot(n) > 0 ? 1 : -1;
            t.is_surface_fs[j] = side;
            mesh.tets[opp_t_id].is_surface_fs[k] = -side;
        }
    }
}

void floatTetWild::boolean_operation(Mesh&                                     mesh,
                                    const CSGTree&                            csg_tree_with_ids,
                                    const std::vector<std::vector<Vector3>>&  Vs,
                                    const std::vector<std::vector<Vector3i>>& Fs)
{
    const MatrixXd C = tet_barycenters(mesh);

    const int max_id = CSGTreeParser::get_max_id(csg_tree_with_ids);
    std::vector<MatrixXd> w(max_id + 1);

    // An empty operand list means the operands are the tracked surfaces of the mesh itself.
    for (int i = 0; i <= max_id; ++i) {
        w[i] = Vs.empty() ? tracked_surface_wn(mesh, C, i)
                          : winding_numbers(to_matrix(Vs[i]), to_matrix(Fs[i]), C, false);
    }

    // A tet survives if the tree says so, and is tagged with the highest-numbered operand it is
    // inside, which is what separate_components writes out.
    int cnt = 0;
    for (auto &t : mesh.tets) {
        if(t.is_removed)
            continue;

        t.is_removed = !CSGTreeParser::keep_tet(csg_tree_with_ids, cnt, w);
        int tid = 0;
        for (int id = 0; id <= max_id; ++id) {
            if (w[id][cnt] > 0.5)
                tid = std::max(id + 1, tid);
        }
        t.scalar = tid;
        cnt++;
    }
}

void floatTetWild::boolean_operation(Mesh& mesh, int op){
    enum { OP_UNION = 0, OP_INTERSECTION = 1, OP_DIFFERENCE = 2 };

    const MatrixXd C  = tet_barycenters(mesh);
    const MatrixXd w1 = tracked_surface_wn(mesh, C, 1);
    const MatrixXd w2 = tracked_surface_wn(mesh, C, 2);

    int cnt = 0;
    for (auto &t : mesh.tets) {
        if(t.is_removed)
            continue;
        const bool in1 = w1(cnt) > 0.5;
        const bool in2 = w2(cnt) > 0.5;
        // Anything else the command line let through keeps every tet, as it always has.
        t.is_removed = !(op == OP_UNION        ? in1 || in2
                       : op == OP_INTERSECTION ? in1 && in2
                       : op == OP_DIFFERENCE   ? in1 && !in2
                                               : true);
        cnt++;
    }
}

void floatTetWild::filter_outside(Mesh& mesh, const MatrixXs& V, const MatrixXi& F) {
    const std::vector<bool> is_inside = tets_inside_surface(mesh, V, F);
    for (size_t t_id = 0; t_id < mesh.tets.size(); ++t_id) {
        if (!is_inside[t_id])
            mesh.tets[t_id].is_removed = true;
    }

    drop_isolated_vertices(mesh);
}

void floatTetWild::filter_outside(Mesh& mesh, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces){
    const MatrixXd W = winding_numbers(
      to_matrix(input_vertices), to_matrix(input_faces), tet_barycenters(mesh), false);

    // Note that index only advances for tets that survive, so a removed one shifts every later
    // lookup. That is what this path has always done.
    int index = 0;
    for (auto &t : mesh.tets) {
        if (t.is_removed)
            continue;
        if (W(index) <= 0.5)
            t.is_removed = true;
        else
            index++;
    }

    drop_isolated_vertices(mesh);
}

void floatTetWild::filter_outside_floodfill(Mesh& mesh) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    std::queue<int> t_queue;
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (tets[i].is_bbox_fs[j] != NOT_BBOX) {
                t_queue.push(i);
                tets[i].is_removed = true;
                break;
            }
        }
    }

    while (!t_queue.empty()) {
        int t_id = t_queue.front();
        t_queue.pop();

        for (int j = 0; j < 4; j++) {
            if (tets[t_id].is_bbox_fs[j] != NOT_BBOX || tets[t_id].is_surface_fs[j] != NOT_SURFACE)
                continue;
            int n_t_id = get_opp_t_id(mesh, t_id, j);
            if (n_t_id < 0 || tets[n_t_id].is_removed)
                continue;
            tets[n_t_id].is_removed = true;
            t_queue.push(n_t_id);
        }
        for (int j = 0; j < 4; j++) {
            vector_erase(tet_vertices[tets[t_id][j]].conn_tets, t_id);
        }
    }

    // The walk above unlinked every tet it removed, so a vertex with none left goes too.
    for (auto &v : tet_vertices) {
        if (v.conn_tets.empty())
            v.is_removed = true;
    }
}

void floatTetWild::smooth_open_boundary(Mesh& mesh, const AABBWrapper& tree) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    mark_outside(mesh);
    const std::vector<std::array<int, 5>> faces = boundary_faces(mesh, &MeshTet::is_outside);
    if (faces.empty())
        return;

    // Only vertices off the tracked surface move, so only they collect their incident faces.
    std::vector<std::vector<int>> conn_b_fs(tet_vertices.size());
    for (int i = 0; i < faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            if (!tet_vertices[faces[i][j]].is_on_surface)
                conn_b_fs[faces[i][j]].push_back(i);
        }
    }

    const int IT = 8;
    for (int it = 0; it < IT; it++) {
        // Move each free vertex towards the centroid of its boundary neighbours, halving the
        // step until it keeps every incident tet valid and under the energy bound.
        for (int v_id = 0; v_id < tet_vertices.size(); v_id++) {
            if (conn_b_fs[v_id].empty())
                continue;

            tet_vertices[v_id].is_freezed = true;

            std::vector<int> n_v_ids;
            for (auto &f_id: conn_b_fs[v_id]) {
                for (int j = 0; j < 3; j++) {
                    if (faces[f_id][j] != v_id)
                        n_v_ids.push_back(faces[f_id][j]);
                }
            }
            vector_unique(n_v_ids);

            Vector3 c(0, 0, 0);
            for (int n_v_id: n_v_ids) {
                c += tet_vertices[n_v_id].pos;
            }
            c /= n_v_ids.size();

            double dis = (c - tet_vertices[v_id].pos).norm();
            Vector3 v = (c - tet_vertices[v_id].pos).normalized();
            static const int N = 7;
            Vector3 p;
            for (int n = 0; n < N; n++) {
                p = tet_vertices[v_id].pos + dis / pow(2, n) * v;
                bool is_valid = true;
                for (int t_id: tet_vertices[v_id].conn_tets) {
                    int j = tets[t_id].find(v_id);
                    if (is_inverted(mesh, t_id, j, p)) {
                        is_valid = false;
                        break;
                    }
                    if (get_quality(p, tet_pos(mesh, t_id, (j + 1) % 4),
                                    tet_pos(mesh, t_id, (j + 2) % 4),
                                    tet_pos(mesh, t_id, (j + 3) % 4)) > mesh.params.stop_energy) {
                        is_valid = false;
                        break;
                    }
                }
                if (!is_valid)
                    continue;

                tet_vertices[v_id].pos = p;
                for (int t_id: tet_vertices[v_id].conn_tets)
                    tets[t_id].quality = get_quality(mesh, t_id);
                break;
            }
        }

        for(auto& v: tet_vertices){
            if(v.is_on_surface)
                v.is_freezed = true;
        }
        edge_collapsing(mesh, tree);
        vertex_smoothing(mesh, tree);

        for (auto& v: tet_vertices)
            v.is_freezed = false;
    }
}

void floatTetWild::get_surface(Mesh& mesh, MatrixXd& V, MatrixXi& F) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    const std::vector<std::array<int, 5>> faces = boundary_faces(mesh, &MeshTet::is_removed);
    if (faces.empty())
        return;

    // Wind each face away from the tet it belongs to, so the surface faces outwards.
    std::vector<std::array<int, 3>> b_faces;
    b_faces.reserve(faces.size());
    for (const auto &f : faces) {
        b_faces.push_back({{f[0], f[1], f[2]}});
        if (!is_inverted(tet_vertices[tets[f[3]][f[4]]],
                         tet_vertices[f[0]],
                         tet_vertices[f[1]],
                         tet_vertices[f[2]]))
            std::swap(b_faces.back()[1], b_faces.back()[2]);
    }

    const std::vector<int> b_v_ids = face_vertices(b_faces);

    V.resize(b_v_ids.size(), 3);
    F.resize(b_faces.size(), 3);
    std::map<int, int> map_v_ids;
    for (int i = 0; i < b_v_ids.size(); i++) {
        map_v_ids[b_v_ids[i]] = i;
        V.row(i) = tet_vertices[b_v_ids[i]].pos;
    }
    for (int i = 0; i < b_faces.size(); i++) {
        F.row(i) = Vector3i(map_v_ids[b_faces[i][0]], map_v_ids[b_faces[i][1]],
                            map_v_ids[b_faces[i][2]]);
    }
}

void floatTetWild::manifold_surface(Mesh& mesh, MatrixXd& V, MatrixXi& F) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    for (auto &v: tet_vertices) {
        if (v.is_removed)
            continue;

        v.conn_tets.erase(std::remove_if(v.conn_tets.begin(), v.conn_tets.end(),
                                         [&](int t_id) { return tets[t_id].is_removed; }),
                          v.conn_tets.end());
        if (v.conn_tets.empty())
            v.is_removed = true;
    }

    manifold_edges(mesh);
    manifold_vertices(mesh);

    get_surface(mesh, V, F);
    std::vector<std::vector<int>> conn_f4v(V.rows());
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++)
            conn_f4v[F(i, j)].push_back(i);
    }
    int V_size = V.rows();
    for (int v_id = 0; v_id < V_size; v_id++) {
        if (conn_f4v[v_id].empty())
            continue;
        // Two faces around the vertex are connected when they share an edge off it that no third
        // face uses. More than one group means the surface pinches here.
        const std::vector<std::vector<int>> f_groups =
            group_connected(conn_f4v[v_id], [&](int f_id, auto &&reached) {
                for (int j = 0; j < 3; j++) {
                    if (F(f_id, j) == v_id)
                        continue;
                    std::vector<int> tmp;
                    set_intersection(conn_f4v[F(f_id, (j + 1) % 3)],
                                     conn_f4v[F(f_id, (j + 2) % 3)], tmp);
                    if (tmp.size() == 2)
                        reached(tmp[0] == f_id ? tmp[1] : tmp[0]);
                }
            });
        if (f_groups.size() < 2)
            continue;

        V.conservativeResize(V.rows() + 1, V.cols());
        V.row(V.rows() - 1) = V.row(v_id);
        for (int f_id: f_groups[0]) {
            for (int j = 0; j < 3; j++) {
                if (F(f_id, j) == v_id) {
                    F(f_id, j) = V.rows() - 1;
                    break;
                }
            }
        }
        conn_f4v.push_back(f_groups[0]);
    }
}
