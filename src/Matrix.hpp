// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// The slice of dense linear algebra the mesher uses: fixed 2/3/4-vectors, a 3x3 matrix with a
// column-pivoting Householder QR solve, and a dynamic row-major matrix used as a list of rows.
//
// This replaces Eigen. Output has to stay byte identical, so the arithmetic is written to match
// what Eigen 3.4 emitted for these shapes, and the two places that differ from the obvious spelling
// are called out where they appear:
//
//   - dot() sums the products left to right with no fused multiply-add. Eigen reduced a 3-vector as
//     a 2-wide packet plus a scalar tail, and the tail add never contracted, because the multiply
//     and the add reached the compiler from separate inlined functions.
//   - cross() does contract. Eigen wrote each component as one expression, so the compiler fused
//     the first product into the subtraction. std::fma says that outright instead of relying on the
//     compiler making the same choice here.

#pragma once

#include <cassert>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace floatTetWild {

namespace matrix_detail {
// Keeps a scalar argument out of template deduction, so that `v * 2` picks up the vector's scalar
// type the way it did through Eigen's Scalar typedef.
template <typename T>
struct identity
{
    typedef T type;
};
}  // namespace matrix_detail

// `s + a * b`, kept unfused. The multiply is its own statement so that no compiler folds it into
// the add, which is what Eigen's reductions did by accident: the two halves reached the compiler
// from separate inlined calls. Anywhere the fusion is wanted, std::fma says so.
template <typename T>
inline T add_product(T s, T a, T b)
{
    const T p = a * b;
    return s + p;
}

// `s - a * b`, kept unfused for the same reason.
template <typename T>
inline T sub_product(T s, T a, T b)
{
    const T p = a * b;
    return s - p;
}

template <typename T, int N>
struct Vector
{
    T v[N];

    // Left uninitialized, as Eigen's fixed-size default constructor was.
    Vector() {}
    Vector(T x, T y)
    {
        static_assert(N == 2, "2 coefficients for a 2-vector");
        v[0] = x;
        v[1] = y;
    }
    Vector(T x, T y, T z)
    {
        static_assert(N == 3, "3 coefficients for a 3-vector");
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }
    Vector(T x, T y, T z, T w)
    {
        static_assert(N == 4, "4 coefficients for a 4-vector");
        v[0] = x;
        v[1] = y;
        v[2] = z;
        v[3] = w;
    }

    T&       operator()(int i) { return v[i]; }
    const T& operator()(int i) const { return v[i]; }
    T&       operator[](int i) { return v[i]; }
    const T& operator[](int i) const { return v[i]; }

    T*       data() { return v; }
    const T* data() const { return v; }

    void setZero()
    {
        for (int i = 0; i < N; i++) v[i] = T(0);
    }

    Vector& operator+=(const Vector& o)
    {
        for (int i = 0; i < N; i++) v[i] += o.v[i];
        return *this;
    }
    Vector& operator/=(T s)
    {
        for (int i = 0; i < N; i++) v[i] /= s;
        return *this;
    }

    template <typename U>
    Vector<U, N> cast() const
    {
        Vector<U, N> r;
        for (int i = 0; i < N; i++) r.v[i] = U(v[i]);
        return r;
    }

    bool operator==(const Vector& o) const
    {
        for (int i = 0; i < N; i++)
            if (!(v[i] == o.v[i])) return false;
        return true;
    }
    bool operator!=(const Vector& o) const { return !(*this == o); }

    T dot(const Vector& o) const
    {
        // Left to right and unfused. See the note at the top of the file.
        T s = v[0] * o.v[0];
        for (int i = 1; i < N; i++) s = add_product(s, v[i], o.v[i]);
        return s;
    }

    Vector cross(const Vector& o) const
    {
        static_assert(N == 3, "cross product is 3D");
        // Deliberately fused: Eigen wrote these as single expressions and the compiler contracted
        // them. See the note at the top of the file.
        return Vector(std::fma(v[1], o.v[2], -(v[2] * o.v[1])),
                      std::fma(v[2], o.v[0], -(v[0] * o.v[2])),
                      std::fma(v[0], o.v[1], -(v[1] * o.v[0])));
    }

    T squaredNorm() const { return dot(*this); }
    T norm() const { return std::sqrt(squaredNorm()); }

    Vector normalized() const
    {
        const T z = squaredNorm();
        if (z > T(0)) return *this / std::sqrt(z);
        return *this;
    }
    void normalize()
    {
        const T z = squaredNorm();
        if (z > T(0)) *this /= std::sqrt(z);
    }

    T maxCoeff() const
    {
        T m = v[0];
        for (int i = 1; i < N; i++) m = m < v[i] ? v[i] : m;
        return m;
    }

    bool allFinite() const
    {
        for (int i = 0; i < N; i++)
            if (!std::isfinite(v[i])) return false;
        return true;
    }
};

template <typename T, int N>
inline Vector<T, N> operator+(const Vector<T, N>& a, const Vector<T, N>& b)
{
    Vector<T, N> r;
    for (int i = 0; i < N; i++) r.v[i] = a.v[i] + b.v[i];
    return r;
}

template <typename T, int N>
inline Vector<T, N> operator-(const Vector<T, N>& a, const Vector<T, N>& b)
{
    Vector<T, N> r;
    for (int i = 0; i < N; i++) r.v[i] = a.v[i] - b.v[i];
    return r;
}

template <typename T, int N>
inline Vector<T, N> operator-(const Vector<T, N>& a)
{
    Vector<T, N> r;
    for (int i = 0; i < N; i++) r.v[i] = -a.v[i];
    return r;
}

template <typename T, int N>
inline Vector<T, N> operator*(const Vector<T, N>& a, typename matrix_detail::identity<T>::type s)
{
    Vector<T, N> r;
    for (int i = 0; i < N; i++) r.v[i] = a.v[i] * s;
    return r;
}

template <typename T, int N>
inline Vector<T, N> operator*(typename matrix_detail::identity<T>::type s, const Vector<T, N>& a)
{
    Vector<T, N> r;
    for (int i = 0; i < N; i++) r.v[i] = s * a.v[i];
    return r;
}

template <typename T, int N>
inline Vector<T, N> operator/(const Vector<T, N>& a, typename matrix_detail::identity<T>::type s)
{
    Vector<T, N> r;
    for (int i = 0; i < N; i++) r.v[i] = a.v[i] / s;
    return r;
}

// 3x3 matrix. Row-major storage: nothing reads the buffer directly, and the solve below walks
// columns explicitly either way.

template <typename T>
struct Matrix33
{
    T m[9];

    Matrix33() {}

    void setZero()
    {
        for (int i = 0; i < 9; i++) m[i] = T(0);
    }

    T&       operator()(int i, int j) { return m[i * 3 + j]; }
    const T& operator()(int i, int j) const { return m[i * 3 + j]; }

    Matrix33& operator+=(const Matrix33& o)
    {
        for (int i = 0; i < 9; i++) m[i] += o.m[i];
        return *this;
    }

    bool allFinite() const
    {
        for (int i = 0; i < 9; i++)
            if (!std::isfinite(m[i])) return false;
        return true;
    }
};

// Matrix times vector. The two shapes below are Eigen's, and they are genuinely different from
// each other: writing a 3-vector result took one 2-wide packet plus a scalar tail, so rows 0 and 1
// came out of the packet path, which multiply-adds one column of the matrix at a time, and row 2
// came out of coeff(), which is a plain reduction over the row with no fused multiply-add.
template <typename T>
inline Vector<T, 3> operator*(const Matrix33<T>& a, const Vector<T, 3>& b)
{
    Vector<T, 3> r;
    for (int i = 0; i < 2; i++) {
        T s  = a(i, 0) * b[0];
        s    = std::fma(a(i, 1), b[1], s);
        s    = std::fma(a(i, 2), b[2], s);
        r[i] = s;
    }
    const T tail = add_product(a(2, 1) * b[1], a(2, 2), b[2]);
    const T head = a(2, 0) * b[0];
    r[2]         = head + tail;
    return r;
}

// Dynamic matrix, stored row-major and used throughout as a list of rows.

template <typename T>
struct MatrixX;

// A reference to one row, so that assigning to it writes through to the matrix. That is why the
// copy assignment below is not the implicit one, and why the copy constructor has to be spelled
// out alongside it.
template <typename T>
struct RowRef
{
    T*  p;
    int n;

    RowRef(T* q, int m) : p(q), n(m) {}
    RowRef(const RowRef&) = default;

    template <int N>
    operator Vector<T, N>() const
    {
        assert(n == N);
        Vector<T, N> r;
        for (int i = 0; i < N; i++) r.v[i] = p[i];
        return r;
    }

    template <int N>
    const RowRef& operator=(const Vector<T, N>& v) const
    {
        assert(n == N);
        for (int i = 0; i < N; i++) p[i] = v.v[i];
        return *this;
    }

    const RowRef& operator=(const RowRef& o) const
    {
        assert(n == o.n);
        for (int i = 0; i < n; i++) p[i] = o.p[i];
        return *this;
    }

    template <int N>
    const RowRef& operator+=(const Vector<T, N>& v) const
    {
        assert(n == N);
        for (int i = 0; i < N; i++) p[i] += v.v[i];
        return *this;
    }

    const RowRef& operator/=(T s) const
    {
        for (int i = 0; i < n; i++) p[i] /= s;
        return *this;
    }
};

template <typename T>
struct MatrixX
{
    std::vector<T> a;
    int            nrows = 0;
    int            ncols = 0;

    MatrixX() {}
    MatrixX(int r, int c) : a(size_t(r) * size_t(c)), nrows(r), ncols(c) {}
    // A dynamic vector is an n by 1 matrix, matching how Eigen's VectorX behaved here.
    explicit MatrixX(int n) : a(size_t(n)), nrows(n), ncols(1) {}

    int rows() const { return nrows; }
    int cols() const { return ncols; }
    int size() const { return nrows * ncols; }

    void resize(int r, int c)
    {
        a.assign(size_t(r) * size_t(c), T());
        nrows = r;
        ncols = c;
    }
    void resize(int n) { resize(n, 1); }

    // Row-major storage makes this a truncation or a zero-filled extension of the row list.
    void conservativeResize(int r, int c)
    {
        assert(c == ncols || nrows == 0);
        a.resize(size_t(r) * size_t(c));
        nrows = r;
        ncols = c;
    }

    T&       operator()(int i, int j) { return a[size_t(i) * ncols + j]; }
    const T& operator()(int i, int j) const { return a[size_t(i) * ncols + j]; }
    // Linear access, for the n by 1 case.
    T&       operator()(int i) { return a[i]; }
    const T& operator()(int i) const { return a[i]; }

    // Vector-style subscript, as Eigen offered on an n by 1 matrix.
    T&       operator[](int i) { return a[i]; }
    const T& operator[](int i) const { return a[i]; }

    T maxCoeff() const
    {
        assert(!a.empty());
        T m = a[0];
        for (size_t i = 1; i < a.size(); i++) m = m < a[i] ? a[i] : m;
        return m;
    }

    T*       data() { return a.data(); }
    const T* data() const { return a.data(); }

    RowRef<T> row(int i) { return RowRef<T>(a.data() + size_t(i) * ncols, ncols); }
};

// Column-pivoting Householder QR for a 3x3 system, and the pieces it needs. This follows Eigen
// 3.4's ColPivHouseholderQR::computeInPlace and _solve_impl step for step, including the pivot
// threshold and the LAPACK norm downdate, because the vertex smoother's Newton step runs through
// it and the mesh has to come out the same.

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
