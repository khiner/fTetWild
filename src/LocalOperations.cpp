// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>

#include <floattetwild/ParallelFor.hpp>
#include <floattetwild/geo_multi_precision.h>

namespace floatTetWild {
bool        use_old_energy       = false;
}  // namespace floatTetWild

using floatTetWild::Scalar;

int floatTetWild::get_opp_t_id(const Mesh& mesh, int t_id, int j)
{
    std::vector<int> pair;
    get_face_tets(mesh,
                  mesh.tets[t_id][(j + 1) % 4],
                  mesh.tets[t_id][(j + 2) % 4],
                  mesh.tets[t_id][(j + 3) % 4],
                  pair);
    if (pair.size() == 2)
        return pair[0] == t_id ? pair[1] : pair[0];
    return OPP_T_ID_BOUNDARY;
}

void floatTetWild::get_all_edges(const Mesh& mesh, std::vector<std::array<int, 2>>& edges)
{
    edges = parallel_collect<std::array<int, 2>>(
      0, mesh.tets.size(), [&](size_t i, std::vector<std::array<int, 2>>& out) {
          if (!mesh.tets[i].is_removed)
              push_tet_edges(mesh.tets[i], out);
      });
    vector_unique(edges);
}

Scalar floatTetWild::get_edge_length(const Mesh& mesh, int v1_id, int v2_id)
{
    return (mesh.tet_vertices[v1_id].pos - mesh.tet_vertices[v2_id].pos).norm();
}

Scalar floatTetWild::get_edge_length_2(const Mesh& mesh, int v1_id, int v2_id)
{
    return (mesh.tet_vertices[v1_id].pos - mesh.tet_vertices[v2_id].pos).squaredNorm();
}

namespace floatTetWild {
namespace {
// True when some tet around the edge carries a tagged face that does not contain the edge. Both
// callers differ only in which per-face tag they look at and what "untagged" means for it.
bool has_tagged_face_off_edge(const Mesh&                       mesh,
                              int                               v1_id,
                              int                               v2_id,
                              const std::vector<int>&           n12_t_ids,
                              std::array<char, 4> MeshTet::*    tags,
                              int                               untagged)
{
    for (int t_id : n12_t_ids) {
        const auto& t = mesh.tets[t_id];
        for (int j = 0; j < 4; j++) {
            if (t[j] != v1_id && t[j] != v2_id && (t.*tags)[j] != untagged)
                return true;
        }
    }
    return false;
}

// Appends the interior samples of the segment from ps[from_id] to ps.back(), spaced by dd.
void sample_segment(std::vector<geo::vec3>& ps, int from_id, Scalar length, Scalar dd)
{
    const int to_id = ps.size() - 1;
    const int N     = length / dd + 1;
    for (Scalar j = 1; j < N - 1; j++)
        ps.push_back(ps[from_id] * (j / N) + ps[to_id] * (1 - j / N));
}
}  // namespace
}  // namespace floatTetWild

bool floatTetWild::is_bbox_edge(const Mesh&             mesh,
                                int                     v1_id,
                                int                     v2_id,
                                const std::vector<int>& n12_t_ids)
{
    if (!mesh.tet_vertices[v1_id].is_on_bbox || !mesh.tet_vertices[v2_id].is_on_bbox)
        return false;
    return has_tagged_face_off_edge(
      mesh, v1_id, v2_id, n12_t_ids, &MeshTet::is_bbox_fs, NOT_BBOX);
}

bool floatTetWild::is_surface_edge(const Mesh&             mesh,
                                   int                     v1_id,
                                   int                     v2_id,
                                   const std::vector<int>& n12_t_ids)
{
    if (!mesh.tet_vertices[v1_id].is_on_surface || !mesh.tet_vertices[v2_id].is_on_surface)
        return false;
    return has_tagged_face_off_edge(
      mesh, v1_id, v2_id, n12_t_ids, &MeshTet::is_surface_fs, NOT_SURFACE);
}

bool floatTetWild::is_boundary_edge(const Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree)
{
    if (!mesh.tet_vertices[v1_id].is_on_boundary || !mesh.tet_vertices[v2_id].is_on_boundary)
        return false;

    std::vector<geo::vec3> ps;
    ps.push_back(to_geo_p(mesh.tet_vertices[v1_id].pos));
    ps.push_back(to_geo_p(mesh.tet_vertices[v2_id].pos));
    sample_segment(ps, 0, get_edge_length(mesh, v1_id, v2_id), mesh.params.dd);

    if (!mesh.is_input_all_inserted)
        return !tree.is_out_tmp_b_envelope(ps, mesh.params.eps_2);
    return !tree.is_out_b_envelope(ps, mesh.params.eps_2);
}

bool floatTetWild::is_valid_edge(const Mesh& mesh, int v1_id, int v2_id)
{
    return !mesh.tet_vertices[v1_id].is_removed && !mesh.tet_vertices[v2_id].is_removed;
}

bool floatTetWild::is_valid_edge(const Mesh&             mesh,
                                 int                     v1_id,
                                 int                     v2_id,
                                 const std::vector<int>& n12_t_ids)
{
    return !n12_t_ids.empty() && is_valid_edge(mesh, v1_id, v2_id);
}

bool floatTetWild::is_isolate_surface_point(const Mesh& mesh, int v_id)
{
    for (int t_id : mesh.tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE)
                return false;
        }
    }

    return true;
}

bool floatTetWild::is_point_out_envelope(const Mesh&        mesh,
                                         const Vector3&     p,
                                         const AABBWrapper& tree)
{
    geo::index_t prev_facet;
    return tree.is_out_sf_envelope(p, mesh.params.eps_2, prev_facet);
}

bool floatTetWild::is_point_out_boundary_envelope(const Mesh&        mesh,
                                                  const Vector3&     p,
                                                  const AABBWrapper& tree)
{
    if (mesh.is_input_all_inserted)
        return false;

    geo::index_t prev_facet;
    return tree.is_out_tmp_b_envelope(p, mesh.params.eps_2, prev_facet);
}

Scalar floatTetWild::get_quality(const Mesh& mesh, const MeshTet& t)
{
    return get_quality(mesh.tet_vertices[t[0]].pos,
                       mesh.tet_vertices[t[1]].pos,
                       mesh.tet_vertices[t[2]].pos,
                       mesh.tet_vertices[t[3]].pos);
}

Scalar floatTetWild::get_quality(const Mesh& mesh, int t_id)
{
    return get_quality(mesh, mesh.tets[t_id]);
}

Scalar floatTetWild::get_quality(const MeshVertex& v0,
                                 const MeshVertex& v1,
                                 const MeshVertex& v2,
                                 const MeshVertex& v3)
{
    return get_quality(v0.pos, v1.pos, v2.pos, v3.pos);
}

Scalar floatTetWild::get_quality(const Vector3& v0,
                                 const Vector3& v1,
                                 const Vector3& v2,
                                 const Vector3& v3)
{
    const std::array<Scalar, 12> T = {
      {v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]}};
    return AMIPS_energy(T);
}

Scalar floatTetWild::get_max_quality(const Mesh& mesh, const std::vector<int>& t_ids)
{
    Scalar max_quality = 0;
    for (int t_id : t_ids) {
        if (mesh.tets[t_id].quality > max_quality)
            max_quality = mesh.tets[t_id].quality;
    }
    return max_quality;
}

bool floatTetWild::is_inverted(const Mesh& mesh, int t_id)
{
    return is_inverted(mesh, mesh.tets[t_id]);
}

bool floatTetWild::is_inverted(const Mesh& mesh, const MeshTet& t)
{
    return is_inverted(mesh.tet_vertices[t[0]].pos,
                       mesh.tet_vertices[t[1]].pos,
                       mesh.tet_vertices[t[2]].pos,
                       mesh.tet_vertices[t[3]].pos);
}

// The tet with its j-th corner moved to new_p.
bool floatTetWild::is_inverted(const Mesh& mesh, int t_id, int j, const Vector3& new_p)
{
    const Vector3* ps[4];
    for (int k = 0; k < 4; k++)
        ps[k] = k == j ? &new_p : &mesh.tet_vertices[mesh.tets[t_id][k]].pos;
    return is_inverted(*ps[0], *ps[1], *ps[2], *ps[3]);
}

bool floatTetWild::is_inverted(const MeshVertex& v0,
                               const MeshVertex& v1,
                               const MeshVertex& v2,
                               const MeshVertex& v3)
{
    return is_inverted(v0.pos, v1.pos, v2.pos, v3.pos);
}

bool floatTetWild::is_inverted(const Vector3& v0,
                               const Vector3& v1,
                               const Vector3& v2,
                               const Vector3& v3)
{
    return Predicates::orient_3d(v0, v1, v2, v3) != Predicates::ORI_POSITIVE;
}

bool floatTetWild::is_degenerate(const Vector3& v0,
                                 const Vector3& v1,
                                 const Vector3& v2,
                                 const Vector3& v3)
{
    return Predicates::orient_3d(v0, v1, v2, v3) == Predicates::ORI_ZERO;
}

bool floatTetWild::is_out_boundary_envelope(const Mesh&        mesh,
                                            int                v_id,
                                            const Vector3&     new_pos,
                                            const AABBWrapper& tree)
{
    if (mesh.is_input_all_inserted)
        return false;
    if (!mesh.tet_vertices[v_id].is_on_cut)
        return false;

    geo::index_t prev_facet;
    if (tree.is_out_tmp_b_envelope(new_pos, mesh.params.eps_2 / 100, prev_facet))
        return true;

    std::vector<int> tmp_b_v_ids;
    for (int t_id : mesh.tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] <= 0) {
                for (int k = 0; k < 3; k++) {
                    int b_v_id = mesh.tets[t_id][(j + 1 + k) % 4];
                    if (b_v_id != v_id && mesh.tet_vertices[b_v_id].is_on_boundary)
                        tmp_b_v_ids.push_back(b_v_id);
                }
            }
        }
    }
    vector_unique(tmp_b_v_ids);

    std::vector<int> b_v_ids;
    b_v_ids.reserve(tmp_b_v_ids.size());
    for (int b_v_id : tmp_b_v_ids) {
        if (is_boundary_edge(mesh, v_id, b_v_id, tree))
            b_v_ids.push_back(b_v_id);
    }
    if (b_v_ids.empty())
        return false;

    std::vector<geo::vec3> ps;
    ps.push_back(to_geo_p(new_pos));
    for (int b_v_id : b_v_ids) {
        ps.push_back(to_geo_p(mesh.tet_vertices[b_v_id].pos));
        sample_segment(ps, 0, get_edge_length(mesh, v_id, b_v_id), mesh.params.dd);
    }

    return tree.is_out_tmp_b_envelope(ps, mesh.params.eps_2 / 100, prev_facet);
}

bool floatTetWild::is_out_envelope(Mesh&              mesh,
                                   int                v_id,
                                   const Vector3&     new_pos,
                                   const AABBWrapper& tree)
{
    geo::index_t prev_facet;
    if (tree.is_out_sf_envelope(new_pos, mesh.params.eps_2, prev_facet))
        return true;

    std::vector<geo::vec3> ps;
    for (int t_id : mesh.tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] <= 0) {
                std::array<Vector3, 3> vs;
                for (int k = 0; k < 3; k++) {
                    if (mesh.tets[t_id][mod4(j + 1 + k)] == v_id)
                        vs[k] = new_pos;
                    else
                        vs[k] = mesh.tet_vertices[mesh.tets[t_id][mod4(j + 1 + k)]].pos;
                }
                if (sample_triangle_and_check_is_out(
                      vs, mesh.params.dd, mesh.params.eps_2, tree, prev_facet))
                    return true;
            }
        }
    }

    return false;
}

namespace floatTetWild {
namespace {
// Walks the sample points of a triangle spaced by sampling_dist and hands each one to visit(),
// stopping early and returning true as soon as visit() does. The two callers below differ only in
// what they do with a point: collect it, or test it against the envelope.
// The samples form a triangular lattice laid out from the longest edge, followed by that edge's own
// subdivision and the other two edges'.
template<typename Visit>
bool walk_triangle_samples(const std::array<Vector3, 3>& vs, Scalar sampling_dist, Visit&& visit)
{
    const Scalar sqrt3_2 = std::sqrt(3) / 2;

    std::array<Scalar, 3> ls;
    for (int i = 0; i < 3; i++)
        ls[i] = (vs[i] - vs[mod3(i + 1)]).squaredNorm();
    // minmax_element, not max_element: on a tie it picks the last of the equal lengths, and which
    // corner the lattice starts from decides the samples.
    const int max_i = std::minmax_element(ls.begin(), ls.end()).second - ls.begin();

    Scalar N = sqrt(ls[max_i]) / sampling_dist;
    if (N <= 1) {
        for (int i = 0; i < 3; i++) {
            if (visit(to_geo_p(vs[i])))
                return true;
        }
        return false;
    }
    if (N == int(N))
        N -= 1;

    const geo::vec3 v0 = to_geo_p(vs[max_i]);
    const geo::vec3 v1 = to_geo_p(vs[mod3(max_i + 1)]);
    const geo::vec3 v2 = to_geo_p(vs[mod3(max_i + 2)]);

    const geo::vec3 n_v0v1 = geo::normalize(v1 - v0);
    for (int n = 0; n <= N; n++) {
        if (visit(v0 + n_v0v1 * sampling_dist * n))
            return true;
    }
    if (visit(v1))
        return true;

    const Scalar h = geo::distance(geo::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
    const int    M = h / (sqrt3_2 * sampling_dist);
    if (M < 1)
        return visit(v2);

    const geo::vec3 n_v0v2 = geo::normalize(v2 - v0);
    const geo::vec3 n_v1v2 = geo::normalize(v2 - v1);
    const Scalar    sin_v0 = geo::length(geo::cross((v2 - v0), (v1 - v0))) /
                          (geo::distance(v0, v2) * geo::distance(v0, v1));
    const Scalar tan_v0 =
      geo::length(geo::cross((v2 - v0), (v1 - v0))) / geo::dot((v2 - v0), (v1 - v0));
    const Scalar sin_v1 = geo::length(geo::cross((v2 - v1), (v0 - v1))) /
                          (geo::distance(v1, v2) * geo::distance(v0, v1));

    for (int m = 1; m <= M; m++) {
        int n  = sqrt3_2 / tan_v0 * m + 0.5;
        int n1 = sqrt3_2 / tan_v0 * m;
        if (m % 2 == 0 && n == n1) {
            n += 1;
        }
        geo::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
        geo::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
        if (geo::distance(v0_m, v1_m) <= sampling_dist)
            break;

        Scalar    delta_d = ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * sampling_dist;
        geo::vec3 v       = v0_m + delta_d * n_v0v1;
        int       N1      = geo::distance(v, v1_m) / sampling_dist;
        for (int i = 0; i <= N1; i++) {
            if (visit(v + i * n_v0v1 * sampling_dist))
                return true;
        }
    }
    if (visit(v2))
        return true;

    // The two shorter edges, walked from v1 towards v2 and from v2 towards v0.
    const std::array<geo::vec3, 2> from = {{v1, v2}};
    const std::array<geo::vec3, 2> dir  = {{n_v1v2, geo::normalize(v0 - v2)}};
    for (int e = 0; e < 2; e++) {
        N = sqrt(ls[mod3(max_i + 1 + e)]) / sampling_dist;
        if (N <= 1)
            continue;
        if (N == int(N))
            N -= 1;
        for (int n = 1; n <= N; n++) {
            if (visit(from[e] + dir[e] * sampling_dist * n))
                return true;
        }
    }

    return false;
}
}  // namespace
}  // namespace floatTetWild

void floatTetWild::sample_triangle(const std::array<Vector3, 3>& vs,
                                   std::vector<geo::vec3>&       ps,
                                   Scalar                        sampling_dist)
{
    walk_triangle_samples(vs, sampling_dist, [&](const geo::vec3& p) {
        ps.push_back(p);
        return false;
    });
}

bool floatTetWild::sample_triangle_and_check_is_out(const std::array<Vector3, 3>& vs,
                                                    Scalar                        sampling_dist,
                                                    Scalar                        eps_2,
                                                    const AABBWrapper&            tree,
                                                    geo::index_t&                 prev_facet)
{
    geo::vec3 nearest_point;
    double    sq_dist = std::numeric_limits<double>::max();
    return walk_triangle_samples(vs, sampling_dist, [&](const geo::vec3& p) {
        return tree.is_out_sf_envelope(p, eps_2, prev_facet, sq_dist, nearest_point);
    });
}

void floatTetWild::get_new_tet_slots(Mesh& mesh, int n, std::vector<int>& new_conn_tets)
{
    int cnt = 0;
    for (int i = mesh.t_empty_start; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed) {
            new_conn_tets.push_back(i);
            cnt++;
            if (cnt == n) {
                mesh.t_empty_start = i + 1;
                break;
            }
        }
    }
    if (cnt < n) {
        for (int i = 0; i < n - cnt; i++)
            new_conn_tets.push_back(mesh.tets.size() + i);
        mesh.tets.resize(mesh.tets.size() + n - cnt);
        mesh.t_empty_start = mesh.tets.size();
    }
}

// Reuses the first free vertex slot at or after v_empty_start, appending when there is none.
int floatTetWild::get_new_vertex_slot(Mesh& mesh, const MeshVertex& v)
{
    while (mesh.v_empty_start < mesh.tet_vertices.size() &&
           !mesh.tet_vertices[mesh.v_empty_start].is_removed)
        mesh.v_empty_start++;

    const int v_id = mesh.v_empty_start;
    if (v_id < mesh.tet_vertices.size())
        mesh.tet_vertices[v_id] = v;
    else
        mesh.tet_vertices.push_back(v);
    return v_id;
}

// After a split has duplicated each tet of old_t_ids into new_t_ids, hand the connectivity over:
// the new vertex joins both halves, v1_id keeps only the new copies, and the corners that are on
// neither end of the split edge gain them.
void floatTetWild::relink_split_tets(Mesh&                   mesh,
                                     int                     v_id,
                                     int                     v1_id,
                                     int                     v2_id,
                                     const std::vector<int>& old_t_ids,
                                     const std::vector<int>& new_t_ids)
{
    auto& tet_vertices = mesh.tet_vertices;
    auto& tets         = mesh.tets;

    auto& conn = tet_vertices[v_id].conn_tets;
    conn.insert(conn.end(), old_t_ids.begin(), old_t_ids.end());
    conn.insert(conn.end(), new_t_ids.begin(), new_t_ids.end());
    for (int i = 0; i < old_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] != v_id && tets[old_t_ids[i]][j] != v2_id)
                tet_vertices[tets[old_t_ids[i]][j]].conn_tets.push_back(new_t_ids[i]);
        }
        auto& conn1 = tet_vertices[v1_id].conn_tets;
        conn1.erase(std::find(conn1.begin(), conn1.end(), old_t_ids[i]));
        conn1.push_back(new_t_ids[i]);
    }
}

// Iterating the smaller set and probing the larger one is what makes this linear in the smaller.
void floatTetWild::set_intersection(const std::unordered_set<int>& s1,
                                    const std::unordered_set<int>& s2,
                                    std::vector<int>&              v)
{
    const std::unordered_set<int>& small = s1.size() <= s2.size() ? s1 : s2;
    const std::unordered_set<int>& large = s1.size() <= s2.size() ? s2 : s1;
    v.clear();
    v.reserve(small.size());
    for (int x : small) {
        if (large.count(x))
            v.push_back(x);
    }
}

// The conn_tets lists these run over are unsorted, so both take sorted copies.
void floatTetWild::set_intersection(const std::vector<int>& s11,
                                    const std::vector<int>& s22,
                                    std::vector<int>&       v)
{
    std::vector<int> s1 = s11;
    std::vector<int> s2 = s22;
    std::sort(s1.begin(), s1.end());
    std::sort(s2.begin(), s2.end());
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(v));
}

void floatTetWild::get_face_tets(const Mesh&       mesh,
                                 int               v1_id,
                                 int               v2_id,
                                 int               v3_id,
                                 std::vector<int>& t_ids)
{
    std::vector<int> s3 = mesh.tet_vertices[v3_id].conn_tets;
    std::sort(s3.begin(), s3.end());
    set_intersection(mesh.tet_vertices[v1_id].conn_tets, mesh.tet_vertices[v2_id].conn_tets, t_ids);
    auto it =
      std::set_intersection(t_ids.begin(), t_ids.end(), s3.begin(), s3.end(), t_ids.begin());
    t_ids.resize(it - t_ids.begin());
}

namespace {
// The AMIPS energy of a tet is
//     P / cbrt(16 * det^2)
// where P is the sum of the six squared edge lengths and det is the determinant of the three edge
// vectors leaving the first vertex, so six times the signed volume. Both are polynomials in the
// twelve coordinates with integer coefficients.
// That is the same value the rational version computed. It evaluated 27/16 * Q^3 / det^2 under a
// cube root, and 3Q is P, so 27/16 * Q^3 / det^2 is P^3 / (16 det^2). The cube root then collapses:
// 16 det^2 is positive wherever det is nonzero and cbrt is odd, so cbrt(P^3 / (16 det^2)) is
// P / cbrt(16 det^2) whatever the sign of P. Nothing has to be cubed.
// det is why this path exists. It is an orient3d determinant and loses itself to cancellation in
// double precision exactly when the tet is near-degenerate, which is when the caller falls through
// to here. Expansions carry it and P exactly, and the result is rounded only at the end.
Scalar AMIPS_energy_exact(const std::array<Scalar, 12>& T)
{
    using floatTetWild::geo::expansion;

    const double* p0 = T.data();
    const double* p1 = T.data() + 3;
    const double* p2 = T.data() + 6;
    const double* p3 = T.data() + 9;

    // These macros allocate in this frame with alloca, so every expansion has to be built here and
    // none of it can move into a helper that returns one.
    const expansion& a11 = expansion_diff(p1[0], p0[0]);
    const expansion& a12 = expansion_diff(p1[1], p0[1]);
    const expansion& a13 = expansion_diff(p1[2], p0[2]);
    const expansion& a21 = expansion_diff(p2[0], p0[0]);
    const expansion& a22 = expansion_diff(p2[1], p0[1]);
    const expansion& a23 = expansion_diff(p2[2], p0[2]);
    const expansion& a31 = expansion_diff(p3[0], p0[0]);
    const expansion& a32 = expansion_diff(p3[1], p0[1]);
    const expansion& a33 = expansion_diff(p3[2], p0[2]);

    // Capacity 192, under the 1024 that new_expansion_on_stack() allows.
    const expansion& det = expansion_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33);
    if (det.sign() == floatTetWild::geo::ZERO)
        return std::numeric_limits<double>::infinity();

    const expansion& d01  = expansion_sq_dist(p0, p1, 3);
    const expansion& d02  = expansion_sq_dist(p0, p2, 3);
    const expansion& d03  = expansion_sq_dist(p0, p3, 3);
    const expansion& d12  = expansion_sq_dist(p1, p2, 3);
    const expansion& d13  = expansion_sq_dist(p1, p3, 3);
    const expansion& d23  = expansion_sq_dist(p2, p3, 3);
    const expansion& sum4 = expansion_sum4(d01, d02, d03, d12);
    const expansion& P    = expansion_sum3(sum4, d13, d23);

    const double d = det.estimate();
    return P.estimate() / std::cbrt(16.0 * d * d);
}
}  // namespace

Scalar floatTetWild::AMIPS_energy(const std::array<Scalar, 12>& T)
{
    const Scalar res = AMIPS_energy_aux(T);
    // !(res > 1e8), not res <= 1e8: a NaN falls through unchanged, as it always has.
    if (use_old_energy || !(res > 1e8))
        return res;

    // The cheap form has lost the determinant to cancellation, so redo it exactly.
    if (is_degenerate(Vector3(T[0], T[1], T[2]),
                      Vector3(T[3], T[4], T[5]),
                      Vector3(T[6], T[7], T[8]),
                      Vector3(T[9], T[10], T[11])))
        return std::numeric_limits<double>::infinity();
    return AMIPS_energy_exact(T);
}

Scalar floatTetWild::AMIPS_energy_aux(const std::array<Scalar, 12>& T)
{
    Scalar helper_1 = T[2];
    Scalar helper_2 = T[11];
    Scalar helper_3 = T[0];
    Scalar helper_4 = T[3];
    Scalar helper_5 = T[9];
    Scalar helper_6 =
      0.577350269189626 * helper_3 - 1.15470053837925 * helper_4 + 0.577350269189626 * helper_5;
    Scalar helper_7  = T[1];
    Scalar helper_8  = T[4];
    Scalar helper_9  = T[7];
    Scalar helper_10 = T[10];
    Scalar helper_11 = 0.408248290463863 * helper_10 + 0.408248290463863 * helper_7 +
                       0.408248290463863 * helper_8 - 1.22474487139159 * helper_9;
    Scalar helper_12 =
      0.577350269189626 * helper_10 + 0.577350269189626 * helper_7 - 1.15470053837925 * helper_8;
    Scalar helper_13 = T[6];
    Scalar helper_14 = -1.22474487139159 * helper_13 + 0.408248290463863 * helper_3 +
                       0.408248290463863 * helper_4 + 0.408248290463863 * helper_5;
    Scalar helper_15 = T[5];
    Scalar helper_16 = T[8];
    Scalar helper_17 = 0.408248290463863 * helper_1 + 0.408248290463863 * helper_15 -
                       1.22474487139159 * helper_16 + 0.408248290463863 * helper_2;
    Scalar helper_18 =
      0.577350269189626 * helper_1 - 1.15470053837925 * helper_15 + 0.577350269189626 * helper_2;
    Scalar helper_19 = 0.5 * helper_13 + 0.5 * helper_4;
    Scalar helper_20 = 0.5 * helper_8 + 0.5 * helper_9;
    Scalar helper_21 = 0.5 * helper_15 + 0.5 * helper_16;
    Scalar helper_22 = (helper_1 - helper_2) * (helper_11 * helper_6 - helper_12 * helper_14) -
                       (-helper_10 + helper_7) * (-helper_14 * helper_18 + helper_17 * helper_6) +
                       (helper_3 - helper_5) * (-helper_11 * helper_18 + helper_12 * helper_17);
    Scalar res =
      -(helper_1 * (-1.5 * helper_1 + 0.5 * helper_2 + helper_21) +
        helper_10 * (-1.5 * helper_10 + helper_20 + 0.5 * helper_7) +
        helper_13 * (-1.5 * helper_13 + 0.5 * helper_3 + 0.5 * helper_4 + 0.5 * helper_5) +
        helper_15 * (0.5 * helper_1 - 1.5 * helper_15 + 0.5 * helper_16 + 0.5 * helper_2) +
        helper_16 * (0.5 * helper_1 + 0.5 * helper_15 - 1.5 * helper_16 + 0.5 * helper_2) +
        helper_2 * (0.5 * helper_1 - 1.5 * helper_2 + helper_21) +
        helper_3 * (helper_19 - 1.5 * helper_3 + 0.5 * helper_5) +
        helper_4 * (0.5 * helper_13 + 0.5 * helper_3 - 1.5 * helper_4 + 0.5 * helper_5) +
        helper_5 * (helper_19 + 0.5 * helper_3 - 1.5 * helper_5) +
        helper_7 * (0.5 * helper_10 + helper_20 - 1.5 * helper_7) +
        helper_8 * (0.5 * helper_10 + 0.5 * helper_7 - 1.5 * helper_8 + 0.5 * helper_9) +
        helper_9 * (0.5 * helper_10 + 0.5 * helper_7 + 0.5 * helper_8 - 1.5 * helper_9)) /
      std::cbrt(helper_22 * helper_22);
    return res;
}

void floatTetWild::AMIPS_jacobian(const std::array<Scalar, 12>& T, Vector3& result_0)
{
    Scalar helper_1 = T[1];
    Scalar helper_2 = T[10];
    Scalar helper_3 = helper_1 - helper_2;
    Scalar helper_4 = T[0];
    Scalar helper_5 = T[3];
    Scalar helper_6 = T[9];
    Scalar helper_7 =
      0.577350269189626 * helper_4 - 1.15470053837925 * helper_5 + 0.577350269189626 * helper_6;
    Scalar helper_8  = T[2];
    Scalar helper_9  = 0.408248290463863 * helper_8;
    Scalar helper_10 = T[5];
    Scalar helper_11 = 0.408248290463863 * helper_10;
    Scalar helper_12 = T[8];
    Scalar helper_13 = 1.22474487139159 * helper_12;
    Scalar helper_14 = T[11];
    Scalar helper_15 = 0.408248290463863 * helper_14;
    Scalar helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
    Scalar helper_17 = 0.577350269189626 * helper_8;
    Scalar helper_18 = 1.15470053837925 * helper_10;
    Scalar helper_19 = 0.577350269189626 * helper_14;
    Scalar helper_20 = helper_17 - helper_18 + helper_19;
    Scalar helper_21 = T[6];
    Scalar helper_22 = -1.22474487139159 * helper_21 + 0.408248290463863 * helper_4 +
                       0.408248290463863 * helper_5 + 0.408248290463863 * helper_6;
    Scalar helper_23 = helper_16 * helper_7 - helper_20 * helper_22;
    Scalar helper_24 = -helper_14 + helper_8;
    Scalar helper_25 = 0.408248290463863 * helper_1;
    Scalar helper_26 = T[4];
    Scalar helper_27 = 0.408248290463863 * helper_26;
    Scalar helper_28 = T[7];
    Scalar helper_29 = 1.22474487139159 * helper_28;
    Scalar helper_30 = 0.408248290463863 * helper_2;
    Scalar helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
    Scalar helper_32 = helper_31 * helper_7;
    Scalar helper_33 = 0.577350269189626 * helper_1;
    Scalar helper_34 = 1.15470053837925 * helper_26;
    Scalar helper_35 = 0.577350269189626 * helper_2;
    Scalar helper_36 = helper_33 - helper_34 + helper_35;
    Scalar helper_37 = helper_22 * helper_36;
    Scalar helper_38 = helper_4 - helper_6;
    Scalar helper_39 = helper_23 * helper_3 - helper_24 * (helper_32 - helper_37) -
                       helper_38 * (helper_16 * helper_36 - helper_20 * helper_31);
    Scalar helper_40 = pow(pow(helper_39, 2), -0.333333333333333);
    Scalar helper_41 = 0.707106781186548 * helper_10 - 0.707106781186548 * helper_12;
    Scalar helper_42 = 0.707106781186548 * helper_26 - 0.707106781186548 * helper_28;
    Scalar helper_43 = 0.5 * helper_21 + 0.5 * helper_5;
    Scalar helper_44 = 0.5 * helper_26 + 0.5 * helper_28;
    Scalar helper_45 = 0.5 * helper_10 + 0.5 * helper_12;
    Scalar helper_46 =
      0.666666666666667 *
      (helper_1 * (-1.5 * helper_1 + 0.5 * helper_2 + helper_44) +
       helper_10 * (-1.5 * helper_10 + 0.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8) +
       helper_12 * (0.5 * helper_10 - 1.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8) +
       helper_14 * (-1.5 * helper_14 + helper_45 + 0.5 * helper_8) +
       helper_2 * (0.5 * helper_1 - 1.5 * helper_2 + helper_44) +
       helper_21 * (-1.5 * helper_21 + 0.5 * helper_4 + 0.5 * helper_5 + 0.5 * helper_6) +
       helper_26 * (0.5 * helper_1 + 0.5 * helper_2 - 1.5 * helper_26 + 0.5 * helper_28) +
       helper_28 * (0.5 * helper_1 + 0.5 * helper_2 + 0.5 * helper_26 - 1.5 * helper_28) +
       helper_4 * (-1.5 * helper_4 + helper_43 + 0.5 * helper_6) +
       helper_5 * (0.5 * helper_21 + 0.5 * helper_4 - 1.5 * helper_5 + 0.5 * helper_6) +
       helper_6 * (0.5 * helper_4 + helper_43 - 1.5 * helper_6) +
       helper_8 * (0.5 * helper_14 + helper_45 - 1.5 * helper_8)) /
      helper_39;
    Scalar helper_47 = -0.707106781186548 * helper_21 + 0.707106781186548 * helper_5;
    result_0[0] =
      -helper_40 *
      (1.0 * helper_21 - 3.0 * helper_4 +
       helper_46 *
         (helper_41 * (-helper_1 + helper_2) - helper_42 * (helper_14 - helper_8) -
          (-helper_17 + helper_18 - helper_19) * (-helper_25 - helper_27 + helper_29 - helper_30) +
          (-helper_33 + helper_34 - helper_35) * (-helper_11 + helper_13 - helper_15 - helper_9)) +
       1.0 * helper_5 + 1.0 * helper_6);
    result_0[1] =
      helper_40 * (3.0 * helper_1 - 1.0 * helper_2 - 1.0 * helper_26 - 1.0 * helper_28 +
                   helper_46 * (helper_23 + helper_24 * helper_47 - helper_38 * helper_41));
    result_0[2] =
      helper_40 *
      (-1.0 * helper_10 - 1.0 * helper_12 - 1.0 * helper_14 +
       helper_46 * (-helper_3 * helper_47 - helper_32 + helper_37 + helper_38 * helper_42) +
       3.0 * helper_8);
}

void floatTetWild::AMIPS_hessian(const std::array<Scalar, 12>& T, Matrix3& result_0)
{
    Scalar helper_1  = T[2];
    Scalar helper_2  = T[11];
    Scalar helper_3  = helper_1 - helper_2;
    Scalar helper_4  = T[0];
    Scalar helper_5  = 0.577350269189626 * helper_4;
    Scalar helper_6  = T[3];
    Scalar helper_7  = 1.15470053837925 * helper_6;
    Scalar helper_8  = T[9];
    Scalar helper_9  = 0.577350269189626 * helper_8;
    Scalar helper_10 = helper_5 - helper_7 + helper_9;
    Scalar helper_11 = T[1];
    Scalar helper_12 = 0.408248290463863 * helper_11;
    Scalar helper_13 = T[4];
    Scalar helper_14 = 0.408248290463863 * helper_13;
    Scalar helper_15 = T[7];
    Scalar helper_16 = 1.22474487139159 * helper_15;
    Scalar helper_17 = T[10];
    Scalar helper_18 = 0.408248290463863 * helper_17;
    Scalar helper_19 = helper_12 + helper_14 - helper_16 + helper_18;
    Scalar helper_20 = helper_10 * helper_19;
    Scalar helper_21 = 0.577350269189626 * helper_11;
    Scalar helper_22 = 1.15470053837925 * helper_13;
    Scalar helper_23 = 0.577350269189626 * helper_17;
    Scalar helper_24 = helper_21 - helper_22 + helper_23;
    Scalar helper_25 = 0.408248290463863 * helper_4;
    Scalar helper_26 = 0.408248290463863 * helper_6;
    Scalar helper_27 = T[6];
    Scalar helper_28 = 1.22474487139159 * helper_27;
    Scalar helper_29 = 0.408248290463863 * helper_8;
    Scalar helper_30 = helper_25 + helper_26 - helper_28 + helper_29;
    Scalar helper_31 = helper_24 * helper_30;
    Scalar helper_32 = helper_3 * (helper_20 - helper_31);
    Scalar helper_33 = helper_4 - helper_8;
    Scalar helper_34 = 0.408248290463863 * helper_1;
    Scalar helper_35 = T[5];
    Scalar helper_36 = 0.408248290463863 * helper_35;
    Scalar helper_37 = T[8];
    Scalar helper_38 = 1.22474487139159 * helper_37;
    Scalar helper_39 = 0.408248290463863 * helper_2;
    Scalar helper_40 = helper_34 + helper_36 - helper_38 + helper_39;
    Scalar helper_41 = helper_24 * helper_40;
    Scalar helper_42 = 0.577350269189626 * helper_1;
    Scalar helper_43 = 1.15470053837925 * helper_35;
    Scalar helper_44 = 0.577350269189626 * helper_2;
    Scalar helper_45 = helper_42 - helper_43 + helper_44;
    Scalar helper_46 = helper_19 * helper_45;
    Scalar helper_47 = helper_41 - helper_46;
    Scalar helper_48 = helper_33 * helper_47;
    Scalar helper_49 = helper_11 - helper_17;
    Scalar helper_50 = helper_10 * helper_40;
    Scalar helper_51 = helper_30 * helper_45;
    Scalar helper_52 = helper_50 - helper_51;
    Scalar helper_53 = helper_49 * helper_52;
    Scalar helper_54 = helper_32 + helper_48 - helper_53;
    Scalar helper_55 = pow(helper_54, 2);
    Scalar helper_56 = pow(helper_55, -0.333333333333333);
    Scalar helper_57 = 1.0 * helper_27 - 3.0 * helper_4 + 1.0 * helper_6 + 1.0 * helper_8;
    Scalar helper_58 = 0.707106781186548 * helper_13;
    Scalar helper_59 = 0.707106781186548 * helper_15;
    Scalar helper_60 = helper_58 - helper_59;
    Scalar helper_61 = helper_3 * helper_60;
    Scalar helper_62 = 0.707106781186548 * helper_35 - 0.707106781186548 * helper_37;
    Scalar helper_63 = helper_49 * helper_62;
    Scalar helper_64 = helper_47 + helper_61 - helper_63;
    Scalar helper_65 = 1.33333333333333 / helper_54;
    Scalar helper_66 = 1.0 / helper_55;
    Scalar helper_67 = 0.5 * helper_27 + 0.5 * helper_6;
    Scalar helper_68 = -1.5 * helper_4 + helper_67 + 0.5 * helper_8;
    Scalar helper_69 = 0.5 * helper_4 + helper_67 - 1.5 * helper_8;
    Scalar helper_70 = -1.5 * helper_27 + 0.5 * helper_4 + 0.5 * helper_6 + 0.5 * helper_8;
    Scalar helper_71 = 0.5 * helper_27 + 0.5 * helper_4 - 1.5 * helper_6 + 0.5 * helper_8;
    Scalar helper_72 = 0.5 * helper_13 + 0.5 * helper_15;
    Scalar helper_73 = -1.5 * helper_11 + 0.5 * helper_17 + helper_72;
    Scalar helper_74 = 0.5 * helper_11 - 1.5 * helper_17 + helper_72;
    Scalar helper_75 = 0.5 * helper_11 + 0.5 * helper_13 - 1.5 * helper_15 + 0.5 * helper_17;
    Scalar helper_76 = 0.5 * helper_11 - 1.5 * helper_13 + 0.5 * helper_15 + 0.5 * helper_17;
    Scalar helper_77 = 0.5 * helper_35 + 0.5 * helper_37;
    Scalar helper_78 = -1.5 * helper_1 + 0.5 * helper_2 + helper_77;
    Scalar helper_79 = 0.5 * helper_1 - 1.5 * helper_2 + helper_77;
    Scalar helper_80 = 0.5 * helper_1 + 0.5 * helper_2 + 0.5 * helper_35 - 1.5 * helper_37;
    Scalar helper_81 = 0.5 * helper_1 + 0.5 * helper_2 - 1.5 * helper_35 + 0.5 * helper_37;
    Scalar helper_82 = helper_1 * helper_78 + helper_11 * helper_73 + helper_13 * helper_76 +
                       helper_15 * helper_75 + helper_17 * helper_74 + helper_2 * helper_79 +
                       helper_27 * helper_70 + helper_35 * helper_81 + helper_37 * helper_80 +
                       helper_4 * helper_68 + helper_6 * helper_71 + helper_69 * helper_8;
    Scalar helper_83 = 0.444444444444444 * helper_66 * helper_82;
    Scalar helper_84 = helper_66 * helper_82;
    Scalar helper_85 = -helper_32 - helper_48 + helper_53;
    Scalar helper_86 = 1.0 / helper_85;
    Scalar helper_87 = helper_86 * pow(pow(helper_85, 2), -0.333333333333333);
    Scalar helper_88 = 0.707106781186548 * helper_6;
    Scalar helper_89 = 0.707106781186548 * helper_27;
    Scalar helper_90 = helper_88 - helper_89;
    Scalar helper_91 =
      0.666666666666667 * helper_10 * helper_40 + 0.666666666666667 * helper_3 * helper_90 -
      0.666666666666667 * helper_30 * helper_45 - 0.666666666666667 * helper_33 * helper_62;
    Scalar helper_92 = -3.0 * helper_11 + 1.0 * helper_13 + 1.0 * helper_15 + 1.0 * helper_17;
    Scalar helper_93 = -helper_11 + helper_17;
    Scalar helper_94 = -helper_1 + helper_2;
    Scalar helper_95 = -helper_21 + helper_22 - helper_23;
    Scalar helper_96 = -helper_34 - helper_36 + helper_38 - helper_39;
    Scalar helper_97 = -helper_42 + helper_43 - helper_44;
    Scalar helper_98 = -helper_12 - helper_14 + helper_16 - helper_18;
    Scalar helper_99 =
      -0.666666666666667 * helper_60 * helper_94 + 0.666666666666667 * helper_62 * helper_93 +
      0.666666666666667 * helper_95 * helper_96 - 0.666666666666667 * helper_97 * helper_98;
    Scalar helper_100 = helper_3 * helper_90;
    Scalar helper_101 = helper_33 * helper_62;
    Scalar helper_102 = helper_100 - helper_101 + helper_52;
    Scalar helper_103 = -helper_60 * helper_94 + helper_62 * helper_93 + helper_95 * helper_96 -
                        helper_97 * helper_98;
    Scalar helper_104 = 0.444444444444444 * helper_102 * helper_103 * helper_82 * helper_86 +
                        helper_57 * helper_91 - helper_92 * helper_99;
    Scalar helper_105 =
      1.85037170770859e-17 * helper_1 * helper_78 + 1.85037170770859e-17 * helper_11 * helper_73 +
      1.85037170770859e-17 * helper_13 * helper_76 + 1.85037170770859e-17 * helper_15 * helper_75 +
      1.85037170770859e-17 * helper_17 * helper_74 + 1.85037170770859e-17 * helper_2 * helper_79 +
      1.85037170770859e-17 * helper_27 * helper_70 + 1.85037170770859e-17 * helper_35 * helper_81 +
      1.85037170770859e-17 * helper_37 * helper_80 + 1.85037170770859e-17 * helper_4 * helper_68 +
      1.85037170770859e-17 * helper_6 * helper_71 + 1.85037170770859e-17 * helper_69 * helper_8;
    Scalar helper_106 = helper_64 * helper_82 * helper_86;
    Scalar helper_107 =
      -0.666666666666667 * helper_10 * helper_19 + 0.666666666666667 * helper_24 * helper_30 +
      0.666666666666667 * helper_33 * helper_60 - 0.666666666666667 * helper_49 * helper_90;
    Scalar helper_108 = -3.0 * helper_1 + 1.0 * helper_2 + 1.0 * helper_35 + 1.0 * helper_37;
    Scalar helper_109 = -helper_20 + helper_31 + helper_33 * helper_60 - helper_49 * helper_90;
    Scalar helper_110 = 0.444444444444444 * helper_109 * helper_82 * helper_86;
    Scalar helper_111 = helper_103 * helper_110 + helper_107 * helper_57 - helper_108 * helper_99;
    Scalar helper_112 = -helper_4 + helper_8;
    Scalar helper_113 = -helper_88 + helper_89;
    Scalar helper_114 = -helper_5 + helper_7 - helper_9;
    Scalar helper_115 = -helper_25 - helper_26 + helper_28 - helper_29;
    Scalar helper_116 = helper_82 * helper_86 *
                        (helper_112 * helper_62 + helper_113 * helper_94 + helper_114 * helper_96 -
                         helper_115 * helper_97);
    Scalar helper_117 = -helper_100 + helper_101 - helper_50 + helper_51;
    Scalar helper_118 = -helper_102 * helper_110 + helper_107 * helper_92 + helper_108 * helper_91;
    Scalar helper_119 = helper_82 * helper_86 *
                        (helper_112 * (-helper_58 + helper_59) - helper_113 * helper_93 -
                         helper_114 * helper_98 + helper_115 * helper_95);
    result_0(0, 0) =
      helper_56 * (helper_57 * helper_64 * helper_65 - pow(helper_64, 2) * helper_83 +
                   0.666666666666667 * helper_64 * helper_84 *
                     (-helper_41 + helper_46 - helper_61 + helper_63) +
                   3.0);
    result_0(0, 1) = helper_87 * (helper_104 - helper_105 * helper_35 + helper_106 * helper_91);
    result_0(0, 2) = helper_87 * (helper_106 * helper_107 + helper_111);
    result_0(1, 0) = helper_87 * (helper_104 + helper_116 * helper_99);
    result_0(1, 1) =
      helper_56 * (-pow(helper_117, 2) * helper_83 + helper_117 * helper_65 * helper_92 +
                   helper_117 * helper_84 * helper_91 + 3.0);
    result_0(1, 2) = helper_87 * (-helper_105 * helper_6 - helper_107 * helper_116 + helper_118);
    result_0(2, 0) = helper_87 * (-helper_105 * helper_13 + helper_111 + helper_119 * helper_99);
    result_0(2, 1) = helper_87 * (helper_118 - helper_119 * helper_91);
    result_0(2, 2) = helper_56 * (-helper_108 * helper_109 * helper_65 -
                                  1.11111111111111 * pow(helper_109, 2) * helper_84 + 3.0);
}
