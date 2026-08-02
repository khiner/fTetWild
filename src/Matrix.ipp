// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Column-pivoting Householder QR for a 3x3 system, and the pieces it needs. This follows Eigen
// 3.4's ColPivHouseholderQR::computeInPlace and _solve_impl step for step, including the pivot
// threshold and the LAPACK norm downdate, because the vertex smoother's Newton step runs through it
// and the mesh has to come out the same.
//
// Included by Matrix.hpp. Not a standalone header.

#pragma once

#include <limits>

namespace floatTetWild {
namespace matrix_detail {

template <typename T>
inline T inner_product(const T* a, const T* b, int n)
{
    T s = a[0] * b[0];
    for (int i = 1; i < n; i++) s = add_product(s, a[i], b[i]);
    return s;
}

// Householder reflector for the column segment v[0..n-1], stored back in place: v[1..n-1] becomes
// the essential part, beta is the new diagonal coefficient and tau the reflector coefficient.
template <typename T>
inline void make_householder_in_place(T* v, int stride, int n, T& tau, T& beta)
{
    T tail_sq_norm = T(0);
    if (n > 1) {
        const T x = v[stride];
        tail_sq_norm = x * x;
        for (int i = 2; i < n; i++) tail_sq_norm = add_product(tail_sq_norm, v[i * stride], v[i * stride]);
    }
    const T c0  = v[0];
    const T tol = (std::numeric_limits<T>::min)();

    if (n == 1 || tail_sq_norm <= tol) {
        tau  = T(0);
        beta = c0;
        for (int i = 1; i < n; i++) v[i * stride] = T(0);
    }
    else {
        beta = std::sqrt(add_product(tail_sq_norm, c0, c0));
        if (c0 >= T(0)) beta = -beta;
        const T d = c0 - beta;
        for (int i = 1; i < n; i++) v[i * stride] /= d;
        tau = (beta - c0) / beta;
    }
}

// Applies the reflector whose coefficient is tau and whose essential part is `essential` to one
// column, from row k down. Both the decomposition and the solve below need this, and the two have
// to stay the same arithmetic.
template <typename T>
inline void apply_householder(T* col, int k, const T* essential, int essential_size, T tau)
{
    T tmp  = inner_product(essential, col + k + 1, essential_size);
    tmp    = tmp + col[k];
    col[k] = sub_product(col[k], tau, tmp);
    for (int i = 0; i < essential_size; i++)
        col[k + 1 + i] = sub_product(col[k + 1 + i], tau * essential[i], tmp);
}

}  // namespace matrix_detail

template <typename T>
Vector<T, 3> solve_col_piv_householder_qr(const Matrix33<T>& matrix, const Vector<T, 3>& rhs)
{
    using matrix_detail::apply_householder;
    using matrix_detail::make_householder_in_place;

    // Column-major, so that a column is a contiguous run and reads like Eigen's m_qr.col(k).
    T qr[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) qr[j * 3 + i] = matrix(i, j);

    T   h_coeffs[3];
    T   norms_updated[3];
    T   norms_direct[3];
    int transpositions[3];

    for (int k = 0; k < 3; k++) {
        const T* c = qr + k * 3;
        T        s = c[0] * c[0];
        s          = add_product(s, c[1], c[1]);
        s          = add_product(s, c[2], c[2]);
        norms_direct[k]  = std::sqrt(s);
        norms_updated[k] = norms_direct[k];
    }

    const T eps      = std::numeric_limits<T>::epsilon();
    T       max_norm = norms_updated[0];
    for (int k = 1; k < 3; k++)
        if (norms_updated[k] > max_norm) max_norm = norms_updated[k];
    const T threshold_helper       = (max_norm * eps) * (max_norm * eps) / T(3);
    const T norm_downdate_threshold = std::sqrt(eps);

    int nonzero_pivots = 3;

    for (int k = 0; k < 3; k++) {
        int biggest_col_index = k;
        T   biggest           = norms_updated[k];
        for (int i = k + 1; i < 3; i++) {
            if (norms_updated[i] > biggest) {
                biggest           = norms_updated[i];
                biggest_col_index = i;
            }
        }
        const T biggest_col_sq_norm = biggest * biggest;

        // Bug 941 in Eigen: the decomposition runs to the end even once the pivots stop being
        // meaningful, so that the original matrix is still reproduced.
        if (nonzero_pivots == 3 && biggest_col_sq_norm < threshold_helper * T(3 - k))
            nonzero_pivots = k;

        transpositions[k] = biggest_col_index;
        if (k != biggest_col_index) {
            for (int i = 0; i < 3; i++) std::swap(qr[k * 3 + i], qr[biggest_col_index * 3 + i]);
            std::swap(norms_updated[k], norms_updated[biggest_col_index]);
            std::swap(norms_direct[k], norms_direct[biggest_col_index]);
        }

        T beta;
        make_householder_in_place(qr + k * 3 + k, 1, 3 - k, h_coeffs[k], beta);
        qr[k * 3 + k] = beta;

        // Apply the reflector to the trailing columns.
        const T tau = h_coeffs[k];
        const int essential_size = 3 - k - 1;
        if (essential_size > 0 && tau != T(0)) {
            const T* essential = qr + k * 3 + k + 1;
            for (int j = k + 1; j < 3; j++)
                apply_householder(qr + j * 3, k, essential, essential_size, tau);
        }

        // LAPACK's stable norm downdate, xGEQP3 lines 278-297.
        for (int j = k + 1; j < 3; j++) {
            if (norms_updated[j] != T(0)) {
                T temp = std::abs(qr[j * 3 + k]) / norms_updated[j];
                temp   = (T(1) + temp) * (T(1) - temp);
                temp   = temp < T(0) ? T(0) : temp;
                const T ratio = norms_updated[j] / norms_direct[j];
                const T temp2 = temp * (ratio * ratio);
                if (temp2 <= norm_downdate_threshold) {
                    T s = qr[j * 3 + k + 1] * qr[j * 3 + k + 1];
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

    Vector<T, 3> dst;
    if (nonzero_pivots == 0) {
        dst.setZero();
        return dst;
    }

    T c[3] = {rhs[0], rhs[1], rhs[2]};

    // c <- Q^T c, reflectors in increasing order.
    for (int k = 0; k < nonzero_pivots; k++) {
        const T   tau            = h_coeffs[k];
        const int essential_size = 3 - k - 1;
        if (essential_size == 0) {
            c[k] *= T(1) - tau;
        }
        else if (tau != T(0))
            apply_householder(c, k, qr + k * 3 + k + 1, essential_size, tau);
    }

    // Back substitution over the leading nonzero_pivots block.
    for (int kk = 0; kk < nonzero_pivots; kk++) {
        const int i = nonzero_pivots - kk - 1;
        c[i] /= qr[i * 3 + i];
        for (int j = 0; j < i; j++) c[j] = sub_product(c[j], c[i], qr[i * 3 + j]);
    }

    for (int i = 0; i < nonzero_pivots; i++) dst[permutation[i]] = c[i];
    for (int i = nonzero_pivots; i < 3; i++) dst[permutation[i]] = T(0);
    return dst;
}

}  // namespace floatTetWild
