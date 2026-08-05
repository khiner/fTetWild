// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LocalOperations.h"

#include "ParallelFor.hpp"

#include <limits>
#include <utility>

namespace floatTetWild {
namespace {

// Column-pivoting Householder QR for the Newton step's 3x3 system, and the pieces it needs. This
// follows Eigen 3.4's ColPivHouseholderQR::computeInPlace and _solve_impl step for step, including
// the pivot threshold and the LAPACK norm downdate, because the mesh has to come out the same.

inline Scalar inner_product(const Scalar* a, const Scalar* b, int n)
{
    Scalar s = a[0] * b[0];
    for (int i = 1; i < n; i++) s = add_product(s, a[i], b[i]);
    return s;
}

// Householder reflector for the column segment v[0..n-1], stored back in place: v[1..n-1] becomes
// the essential part, beta is the new diagonal coefficient and tau the reflector coefficient.
inline void make_householder_in_place(Scalar* v, int stride, int n, Scalar& tau, Scalar& beta)
{
    Scalar tail_sq_norm = 0;
    if (n > 1) {
        const Scalar x = v[stride];
        tail_sq_norm = x * x;
        for (int i = 2; i < n; i++) tail_sq_norm = add_product(tail_sq_norm, v[i * stride], v[i * stride]);
    }
    const Scalar c0  = v[0];
    const Scalar tol = (std::numeric_limits<Scalar>::min)();

    if (n == 1 || tail_sq_norm <= tol) {
        tau  = 0;
        beta = c0;
        for (int i = 1; i < n; i++) v[i * stride] = 0;
    }
    else {
        beta = std::sqrt(add_product(tail_sq_norm, c0, c0));
        if (c0 >= 0) beta = -beta;
        const Scalar d = c0 - beta;
        for (int i = 1; i < n; i++) v[i * stride] /= d;
        tau = (beta - c0) / beta;
    }
}

// Applies the reflector whose coefficient is tau and whose essential part is `essential` to one
// column, from row k down. Both the decomposition and the solve need this, and the two have to
// stay the same arithmetic.
inline void apply_householder(Scalar* col, int k, const Scalar* essential, int essential_size, Scalar tau)
{
    Scalar tmp = inner_product(essential, col + k + 1, essential_size);
    tmp    = tmp + col[k];
    col[k] = sub_product(col[k], tau, tmp);
    for (int i = 0; i < essential_size; i++)
        col[k + 1 + i] = sub_product(col[k + 1 + i], tau * essential[i], tmp);
}

Vector3 solve_col_piv_householder_qr(const Matrix3& matrix, const Vector3& rhs)
{
    // Column-major, so that a column is a contiguous run and reads like Eigen's m_qr.col(k).
    Scalar qr[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) qr[j * 3 + i] = matrix(i, j);

    Scalar h_coeffs[3];
    Scalar norms_updated[3];
    Scalar norms_direct[3];
    int    transpositions[3];

    for (int k = 0; k < 3; k++) {
        const Scalar* c = qr + k * 3;
        Scalar        s = c[0] * c[0];
        s               = add_product(s, c[1], c[1]);
        s               = add_product(s, c[2], c[2]);
        norms_direct[k]  = std::sqrt(s);
        norms_updated[k] = norms_direct[k];
    }

    const Scalar eps      = std::numeric_limits<Scalar>::epsilon();
    Scalar       max_norm = norms_updated[0];
    for (int k = 1; k < 3; k++)
        if (norms_updated[k] > max_norm) max_norm = norms_updated[k];
    const Scalar threshold_helper        = (max_norm * eps) * (max_norm * eps) / 3;
    const Scalar norm_downdate_threshold = std::sqrt(eps);

    int nonzero_pivots = 3;

    for (int k = 0; k < 3; k++) {
        int    biggest_col_index = k;
        Scalar biggest           = norms_updated[k];
        for (int i = k + 1; i < 3; i++) {
            if (norms_updated[i] > biggest) {
                biggest           = norms_updated[i];
                biggest_col_index = i;
            }
        }
        const Scalar biggest_col_sq_norm = biggest * biggest;

        // Bug 941 in Eigen: the decomposition runs to the end even once the pivots stop being
        // meaningful, so that the original matrix is still reproduced.
        if (nonzero_pivots == 3 && biggest_col_sq_norm < threshold_helper * Scalar(3 - k))
            nonzero_pivots = k;

        transpositions[k] = biggest_col_index;
        if (k != biggest_col_index) {
            for (int i = 0; i < 3; i++) std::swap(qr[k * 3 + i], qr[biggest_col_index * 3 + i]);
            std::swap(norms_updated[k], norms_updated[biggest_col_index]);
            std::swap(norms_direct[k], norms_direct[biggest_col_index]);
        }

        Scalar beta;
        make_householder_in_place(qr + k * 3 + k, 1, 3 - k, h_coeffs[k], beta);
        qr[k * 3 + k] = beta;

        // Apply the reflector to the trailing columns.
        const Scalar tau            = h_coeffs[k];
        const int    essential_size = 3 - k - 1;
        if (essential_size > 0 && tau != 0) {
            const Scalar* essential = qr + k * 3 + k + 1;
            for (int j = k + 1; j < 3; j++)
                apply_householder(qr + j * 3, k, essential, essential_size, tau);
        }

        // LAPACK's stable norm downdate, xGEQP3 lines 278-297.
        for (int j = k + 1; j < 3; j++) {
            if (norms_updated[j] != 0) {
                Scalar temp = std::abs(qr[j * 3 + k]) / norms_updated[j];
                temp        = (1 + temp) * (1 - temp);
                temp        = temp < 0 ? 0 : temp;
                const Scalar ratio = norms_updated[j] / norms_direct[j];
                const Scalar temp2 = temp * (ratio * ratio);
                if (temp2 <= norm_downdate_threshold) {
                    Scalar s = qr[j * 3 + k + 1] * qr[j * 3 + k + 1];
                    for (int i = k + 2; i < 3; i++) s = add_product(s, qr[j * 3 + i], qr[j * 3 + i]);
                    norms_direct[j]  = std::sqrt(s);
                    norms_updated[j] = norms_direct[j];
                }
                else {
                    norms_updated[j] *= std::sqrt(temp);
                }
            }
        }
    }

    int permutation[3] = {0, 1, 2};
    for (int k = 0; k < 3; k++) std::swap(permutation[k], permutation[transpositions[k]]);

    Vector3 dst;
    if (nonzero_pivots == 0) {
        dst.setZero();
        return dst;
    }

    Scalar c[3] = {rhs[0], rhs[1], rhs[2]};

    // c <- Q^T c, reflectors in increasing order.
    for (int k = 0; k < nonzero_pivots; k++) {
        const Scalar tau            = h_coeffs[k];
        const int    essential_size = 3 - k - 1;
        if (essential_size == 0) {
            c[k] *= 1 - tau;
        }
        else if (tau != 0)
            apply_householder(c, k, qr + k * 3 + k + 1, essential_size, tau);
    }

    // Back substitution over the leading nonzero_pivots block.
    for (int kk = 0; kk < nonzero_pivots; kk++) {
        const int i = nonzero_pivots - kk - 1;
        c[i] /= qr[i * 3 + i];
        for (int j = 0; j < i; j++) c[j] = sub_product(c[j], c[i], qr[i * 3 + j]);
    }

    for (int i = 0; i < nonzero_pivots; i++) dst[permutation[i]] = c[i];
    for (int i = nonzero_pivots; i < 3; i++) dst[permutation[i]] = 0;
    return dst;
}

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
        Scalar new_q = get_quality(mesh, t_id, j, p);
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

    std::vector<std::array<Scalar, 12>> Ts;
    for (int t_id : tet_vertices[v_id].conn_tets) {
        int j = tets[t_id].find(v_id);

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
            for (int t_id : tet_vertices[v_id].conn_tets) {
                if (is_inverted(mesh, t_id, tets[t_id].find(v_id), x_next)) {
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
