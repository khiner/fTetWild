// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LocalOperations.h"
#include "VertexSmoothing.h"

#include "ParallelFor.hpp"

namespace floatTetWild {
namespace {

// The vertices grouped into sets that can be smoothed at the same time, plus the ones that have
// to go one at a time: the removed ones, and the colours with fewer than two members. Greedy
// graph colouring over the one-ring adjacency, so that same-coloured vertices are never
// neighbours. -1 marks a removed vertex.
void one_ring_vertex_sets(const Mesh&                    mesh,
                          std::vector<std::vector<int>>& concurrent_sets,
                          std::vector<int>&              serial_set)
{
    const auto& tet_vertices = mesh.tet_vertices;

    std::vector<int> colors(tet_vertices.size(), -1);
    colors[0] = 0;

    std::vector<bool> available(tet_vertices.size(), true);

    std::vector<int> ring;

    for (int i = 1; i < tet_vertices.size(); ++i) {
        const auto& v = tet_vertices[i];
        if (v.is_removed)
            continue;

        ring.clear();
        collect_tet_vertices(mesh, v.conn_tets, ring);

        for (const auto n : ring) {
            if (colors[n] != -1)
                available[colors[n]] = false;
        }

        colors[i] = std::find(available.begin(), available.end(), true) - available.begin();

        for (const auto n : ring) {
            if (colors[n] != -1)
                available[colors[n]] = true;
        }
    }

    int max_c = -1;
    for (const auto c : colors)
        max_c = std::max(max_c, c);

    concurrent_sets.clear();
    concurrent_sets.resize(max_c + 1);
    serial_set.clear();

    for (size_t i = 0; i < colors.size(); ++i) {
        const int col = colors[i];
        if (col < 0)
            serial_set.push_back(i);
        else
            concurrent_sets[col].push_back(i);
    }

    for (int i = concurrent_sets.size() - 1; i >= 0; --i) {
        if (concurrent_sets[i].size() < 2) {
            serial_set.insert(
              serial_set.end(), concurrent_sets[i].begin(), concurrent_sets[i].end());
            concurrent_sets.erase(concurrent_sets.begin() + i);
        }
    }
}

bool project_and_check(Mesh&                mesh,
                       int                  v_id,
                       Vector3&             p,
                       const AABBWrapper&   tree,
                       bool                 is_sf,
                       std::vector<Scalar>& new_qs)
{
    if (is_sf)
        tree.project_to_sf(p);
    else {
        if (mesh.is_input_all_inserted)
            tree.project_to_b(p);
        else
            tree.project_to_tmp_b(p);
    }

    const Scalar max_q = get_max_quality(mesh, mesh.tet_vertices[v_id].conn_tets);
    for (int t_id : mesh.tet_vertices[v_id].conn_tets) {
        auto& t = mesh.tets[t_id];
        int   j = t.find(v_id);
        if (is_inverted(mesh, t_id, j, p))
            return false;
        Scalar new_q = get_quality(p,
                                   mesh.tet_vertices[t[(j + 1) % 4]].pos,
                                   mesh.tet_vertices[t[(j + 2) % 4]].pos,
                                   mesh.tet_vertices[t[(j + 3) % 4]].pos);
        if (new_q > max_q)
            return false;
        new_qs.push_back(new_q);
    }

    return true;
}

bool find_new_pos(Mesh& mesh, const int v_id, Vector3& x)
{
    auto& tets         = mesh.tets;
    auto& tet_vertices = mesh.tet_vertices;

    std::vector<int> js;
    js.reserve(tet_vertices[v_id].conn_tets.size());
    std::vector<std::array<Scalar, 12>> Ts;
    for (int t_id : tet_vertices[v_id].conn_tets) {
        int j = tets[t_id].find(v_id);
        js.push_back(j);

        std::array<int, 4> loop_ids = {{0, 1, 2, 3}};
        if (is_inverted(tet_vertices[tets[t_id][j]].pos,
                        tet_vertices[tets[t_id][(j + 1) % 4]].pos,
                        tet_vertices[tets[t_id][(j + 2) % 4]].pos,
                        tet_vertices[tets[t_id][(j + 3) % 4]].pos))
            std::swap(loop_ids[2], loop_ids[3]);

        std::array<Scalar, 12> T;
        for (int k = 0; k < loop_ids.size(); k++) {
            const Vector3& p = tet_pos(mesh, t_id, (j + loop_ids[k]) % 4);
            for (int c = 0; c < 3; c++)
                T[k * 3 + c] = p[c];
        }
        Ts.push_back(T);
    }

    const int    max_newton_it = 15;
    const int    max_search_it = 10;
    const Scalar J_delta       = 1e-8;

    x = tet_vertices[v_id].pos;
    Vector3 J;
    Matrix3 H;

    // The vertex being moved is the first corner of every T, so a candidate position is tried by
    // writing it there.
    const auto set_candidate = [&Ts](const Vector3& p) {
        for (auto& T : Ts) {
            T[0] = p(0);
            T[1] = p(1);
            T[2] = p(2);
        }
    };
    const auto total_energy = [&Ts]() {
        Scalar sum = 0;
        for (auto& T : Ts)
            sum += AMIPS_energy(T);
        return sum;
    };

    for (int newton_it = 0; newton_it < max_newton_it; newton_it++) {
        if (newton_it > 0)
            set_candidate(x);

        const Scalar f = total_energy();

        J.setZero();
        for (auto& T : Ts) {
            Vector3 tmp_J;
            AMIPS_jacobian(T, tmp_J);
            J += tmp_J;
        }
        if (!J.allFinite() || (std::abs(J(0)) < J_delta && std::abs(J(1)) < J_delta &&
                               std::abs(J(2)) < J_delta))
            break;

        H.setZero();
        for (auto& T : Ts) {
            Matrix3 tmp_H;
            AMIPS_hessian(T, tmp_H);
            H += tmp_H;
        }
        if (!H.allFinite())
            break;

        bool    found_step = false;
        Scalar  a          = 1;
        Vector3 x_next;
        for (int i = 0; i < max_search_it; i++) {
            x_next = solve_col_piv_householder_qr(H, H * x - a * J);
            if (!x_next.allFinite())
                break;
            set_candidate(x_next);

            bool is_valid = true;
            int  ii       = 0;
            for (int t_id : tet_vertices[v_id].conn_tets) {
                int j = js[ii++];
                if (is_inverted(mesh, t_id, j, x_next)) {
                    is_valid = false;
                    break;
                }
            }
            if (!is_valid) {
                a /= 2;
                continue;
            }

            if (total_energy() >= f) {
                a /= 2;
                continue;
            }

            found_step = true;
            break;
        }

        if (!found_step || !x_next.allFinite())
            break;
        x = x_next;
    }

    return x != tet_vertices[v_id].pos;
}

}  // namespace
}  // namespace floatTetWild

void floatTetWild::vertex_smoothing(Mesh& mesh, const AABBWrapper& tree)
{
    auto& tets         = mesh.tets;
    auto& tet_vertices = mesh.tet_vertices;

    const auto smooth_one = [&](const int v_id) {
        if (tet_vertices[v_id].is_removed)
            return;
        if (tet_vertices[v_id].is_freezed)
            return;
        if (tet_vertices[v_id].is_on_bbox)
            return;

        Vector3 p;
        if (!find_new_pos(mesh, v_id, p))
            return;

        std::vector<Scalar> new_qs;
        if (tet_vertices[v_id].is_on_boundary) {
            if (!project_and_check(mesh, v_id, p, tree, false, new_qs))
                return;
            if (is_out_boundary_envelope(mesh, v_id, p, tree))
                return;
            else if (is_out_envelope(mesh, v_id, p, tree))
                return;
        }
        else if (tet_vertices[v_id].is_on_surface) {
            if (!project_and_check(mesh, v_id, p, tree, true, new_qs))
                return;
            if (is_out_envelope(mesh, v_id, p, tree))
                return;
        }

        tet_vertices[v_id].pos = p;

        int cnt = 0;
        for (int t_id : tet_vertices[v_id].conn_tets) {
            if (!new_qs.empty())
                tets[t_id].quality = new_qs[cnt++];
            else
                tets[t_id].quality = get_quality(mesh, t_id);
        }
    };

    std::vector<std::vector<int>> concurrent_sets;
    std::vector<int>              serial_set;
    // The partition is the one a serial run always used, whatever the thread count: serial_set
    // mixes colours and so changes the order neighbours are smoothed in, which used to make the
    // output depend on how many threads ran.
    one_ring_vertex_sets(mesh, concurrent_sets, serial_set);

    for (const auto& s : concurrent_sets) {
        parallel_for(size_t(0), size_t(s.size()), [&](size_t i) { smooth_one(s[i]); });
    }

    for (size_t v_id : serial_set)
        smooth_one(v_id);
}
