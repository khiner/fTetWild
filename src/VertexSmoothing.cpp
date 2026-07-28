// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/LocalOperations.h>
#include <floattetwild/VertexSmoothing.h>

#include <floattetwild/MeshImprovement.h>
#include <floattetwild/ParallelFor.hpp>

#include <atomic>

void floatTetWild::vertex_smoothing(Mesh& mesh, const AABBWrapper& tree)
{
    auto& tets         = mesh.tets;
    auto& tet_vertices = mesh.tet_vertices;

    // smooth_one runs under a parallel_for, so these have to be atomic even though they only feed
    // the log line below.
    std::atomic<int> counter{0};
    std::atomic<int> suc_counter{0};
    std::atomic<int> suc_counter_sf{0};

    const auto smooth_one = [&](const int v_id) {
        if (tet_vertices[v_id].is_removed)
            return;
        if (tet_vertices[v_id].is_freezed)
            return;
        if (tet_vertices[v_id].is_on_bbox)
            return;
        counter++;

        ////newton
        Vector3 p;
        if (!find_new_pos(mesh, v_id, p))
            return;


        ////check
        // envelope
        std::vector<Scalar> new_qs;
        if (tet_vertices[v_id].is_on_boundary) {
            if (!project_and_check(mesh, v_id, p, tree, false, new_qs))
                return;
            if (is_out_boundary_envelope(mesh, v_id, p, tree))
                return;
            else if (is_out_envelope(mesh, v_id, p, tree))
                return;
            suc_counter_sf++;
        }
        else if (tet_vertices[v_id].is_on_surface) {
            if (!project_and_check(mesh, v_id, p, tree, true, new_qs))
                return;
            if (is_out_envelope(mesh, v_id, p, tree))
                return;
            suc_counter_sf++;
        }
        suc_counter++;

        ////real update
        tet_vertices[v_id].pos = p;

        // quality
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
    // 2 is what params.num_threads * 2 gave at one thread, so this is the partition a serial run
    // always used. serial_set mixes colours and so changes the order neighbours are smoothed in,
    // which made the output depend on the thread count.
    mesh.one_ring_vertex_sets(2, concurrent_sets, serial_set);

    for (const auto& s : concurrent_sets) {
        parallel_for(size_t(0), size_t(s.size()), [&](size_t i) { smooth_one(s[i]); });
    }

    for (size_t v_id : serial_set)
        smooth_one(v_id);

    cout << "success = " << suc_counter << "(" << counter << ")" << endl;
}

bool floatTetWild::project_and_check(Mesh&                mesh,
                                     int                  v_id,
                                     Vector3&             p,
                                     const AABBWrapper&   tree,
                                     bool                 is_sf,
                                     std::vector<Scalar>& new_qs)
{
    // project to surface
    if (is_sf)
        tree.project_to_sf(p);
    else {
        if (mesh.is_input_all_inserted)
            tree.project_to_b(p);
        else
            tree.project_to_tmp_b(p);
    }


    // check inversion & quality
    double max_q = 0;
    for (int t_id : mesh.tet_vertices[v_id].conn_tets) {
        if (mesh.tets[t_id].quality > max_q)
            max_q = mesh.tets[t_id].quality;
    }
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

bool floatTetWild::find_new_pos(Mesh& mesh, const int v_id, Vector3& x)
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
        if (is_inverted(tet_vertices[tets[t_id][j]],
                        tet_vertices[tets[t_id][(j + 1) % 4]],
                        tet_vertices[tets[t_id][(j + 2) % 4]],
                        tet_vertices[tets[t_id][(j + 3) % 4]]))
            std::swap(loop_ids[2], loop_ids[3]);

        std::array<Scalar, 12> T;
        for (int k = 0; k < loop_ids.size(); k++) {
            T[k * 3]     = tet_vertices[tets[t_id][(j + loop_ids[k]) % 4]].pos[0];
            T[k * 3 + 1] = tet_vertices[tets[t_id][(j + loop_ids[k]) % 4]].pos[1];
            T[k * 3 + 2] = tet_vertices[tets[t_id][(j + loop_ids[k]) % 4]].pos[2];
        }
        Ts.push_back(T);
    }

    ////newton
    const int    max_newton_it = 15;
    const int    max_search_it = 10;
    const Scalar f_delta       = 1e-8;
    const Scalar J_delta       = 1e-8;

    x = tet_vertices[v_id].pos;
    Vector3 J;
    Matrix3 H;

    for (int newton_it = 0; newton_it < max_newton_it; newton_it++) {
        if (newton_it > 0) {
            for (auto& T : Ts) {
                T[0] = x(0);
                T[1] = x(1);
                T[2] = x(2);
            }
        }

        // f
        Scalar f = 0;
        for (auto& T : Ts) {
            f += AMIPS_energy(T);
        }


        // J
        J << 0, 0, 0;
        for (auto& T : Ts) {
            Vector3 tmp_J;
            AMIPS_jacobian(T, tmp_J);
            J += tmp_J;
        }
        if (!J.allFinite() || (std::abs(J(0)) < J_delta && std::abs(J(1)) < J_delta &&
                               std::abs(J(2)) < J_delta))  // gradient is also zero
            break;

        // H
        H << 0, 0, 0, 0, 0, 0, 0, 0, 0;
        for (auto& T : Ts) {
            Matrix3 tmp_H;
            AMIPS_hessian(T, tmp_H);
            H += tmp_H;
        }
        if (!H.allFinite())
            break;

        // x
        bool    found_step = false;
        Scalar  a          = 1;
        Vector3 x_next;
        for (int i = 0; i < max_search_it; i++) {
            x_next = solve_col_piv_householder_qr(H, H * x - a * J);
            if (!x_next.allFinite())
                break;
            for (auto& T : Ts) {
                T[0] = x_next(0);
                T[1] = x_next(1);
                T[2] = x_next(2);
            }

            // check inversion
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

            // check energy
            Scalar f_next = 0;
            for (auto& T : Ts) {
                f_next += AMIPS_energy(T);
            }
            if (f_next >= f) {
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

    if (x != tet_vertices[v_id].pos)
        return true;
    return false;
}
