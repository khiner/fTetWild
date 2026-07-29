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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
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

    static constexpr int size() { return N; }
    static constexpr int rows() { return N; }
    static constexpr int cols() { return 1; }

    T*       begin() { return v; }
    T*       end() { return v + N; }
    const T* begin() const { return v; }
    const T* end() const { return v + N; }

    void setZero()
    {
        for (int i = 0; i < N; i++) v[i] = T(0);
    }

    Vector& operator+=(const Vector& o)
    {
        for (int i = 0; i < N; i++) v[i] += o.v[i];
        return *this;
    }
    Vector& operator-=(const Vector& o)
    {
        for (int i = 0; i < N; i++) v[i] -= o.v[i];
        return *this;
    }
    Vector& operator*=(T s)
    {
        for (int i = 0; i < N; i++) v[i] *= s;
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

template <typename T, int N>
struct VectorCommaInit
{
    Vector<T, N>& target;
    int           index;

    VectorCommaInit& operator,(typename matrix_detail::identity<T>::type x)
    {
        assert(index < N);
        target.v[index++] = x;
        return *this;
    }
};

template <typename T, int N>
inline VectorCommaInit<T, N> operator<<(Vector<T, N>& target, typename matrix_detail::identity<T>::type x)
{
    target.v[0] = x;
    return VectorCommaInit<T, N>{target, 1};
}

// ------------------------------------------------------------------------------------------------
// 3x3 matrix. Row-major storage: nothing reads the buffer directly, and the solve below walks
// columns explicitly either way.
// ------------------------------------------------------------------------------------------------

template <typename T>
struct Matrix33
{
    T m[9];

    Matrix33() {}

    T&       operator()(int i, int j) { return m[i * 3 + j]; }
    const T& operator()(int i, int j) const { return m[i * 3 + j]; }

    void setZero()
    {
        for (int i = 0; i < 9; i++) m[i] = T(0);
    }

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

template <typename T>
struct Matrix33CommaInit
{
    Matrix33<T>& target;
    int          index;

    Matrix33CommaInit& operator,(typename matrix_detail::identity<T>::type x)
    {
        assert(index < 9);
        target.m[index++] = x;
        return *this;
    }
};

template <typename T>
inline Matrix33CommaInit<T> operator<<(Matrix33<T>& target, typename matrix_detail::identity<T>::type x)
{
    target.m[0] = x;
    return Matrix33CommaInit<T>{target, 1};
}

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

// ------------------------------------------------------------------------------------------------
// Dynamic matrix, stored row-major and used throughout as a list of rows.
// ------------------------------------------------------------------------------------------------

template <typename T>
struct MatrixX;

template <typename T, bool Const>
struct RowRefBase
{
    typedef typename std::conditional<Const, const T, T>::type Elem;

    Elem* p;
    int   n;

    Elem&       operator()(int j) const { return p[j]; }
    Elem&       operator[](int j) const { return p[j]; }
    int         size() const { return n; }

    template <int N>
    operator Vector<T, N>() const
    {
        assert(n == N);
        Vector<T, N> r;
        for (int i = 0; i < N; i++) r.v[i] = p[i];
        return r;
    }

    bool operator==(const RowRefBase& o) const
    {
        assert(n == o.n);
        for (int i = 0; i < n; i++)
            if (!(p[i] == o.p[i])) return false;
        return true;
    }
    bool operator!=(const RowRefBase& o) const { return !(*this == o); }
};

template <typename T>
struct RowRef : RowRefBase<T, false>
{
    RowRef(T* q, int m) : RowRefBase<T, false>{q, m} {}
    // Assigning a row writes through to the matrix, so the copy assignment below is not the
    // implicit one and the copy constructor has to be spelled out alongside it.
    RowRef(const RowRef&) = default;

    template <int N>
    const RowRef& operator=(const Vector<T, N>& v) const
    {
        assert(this->n == N);
        for (int i = 0; i < N; i++) this->p[i] = v.v[i];
        return *this;
    }

    template <bool C>
    const RowRef& operator=(const RowRefBase<T, C>& o) const
    {
        assert(this->n == o.n);
        for (int i = 0; i < this->n; i++) this->p[i] = o.p[i];
        return *this;
    }

    const RowRef& operator=(const RowRef& o) const
    {
        assert(this->n == o.n);
        for (int i = 0; i < this->n; i++) this->p[i] = o.p[i];
        return *this;
    }

    template <int N>
    const RowRef& operator+=(const Vector<T, N>& v) const
    {
        assert(this->n == N);
        for (int i = 0; i < N; i++) this->p[i] += v.v[i];
        return *this;
    }

    const RowRef& operator/=(T s) const
    {
        for (int i = 0; i < this->n; i++) this->p[i] /= s;
        return *this;
    }
};

template <typename T>
struct RowCommaInit
{
    T*  p;
    int n;
    int index;

    RowCommaInit& operator,(typename matrix_detail::identity<T>::type x)
    {
        assert(index < n);
        p[index++] = x;
        return *this;
    }
};

template <typename T>
inline RowCommaInit<T> operator<<(const RowRef<T>& target, typename matrix_detail::identity<T>::type x)
{
    target.p[0] = x;
    return RowCommaInit<T>{target.p, target.n, 1};
}

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

    void setZero()
    {
        std::fill(a.begin(), a.end(), T(0));
    }

    T&       operator()(int i, int j) { return a[size_t(i) * ncols + j]; }
    const T& operator()(int i, int j) const { return a[size_t(i) * ncols + j]; }
    // Linear access, for the n by 1 case.
    T&       operator()(int i) { return a[i]; }
    const T& operator()(int i) const { return a[i]; }

    // Vector-style subscript, as Eigen offered on an n by 1 matrix.
    T&       operator[](int i) { return a[i]; }
    const T& operator[](int i) const { return a[i]; }

    T        coeff(int i, int j) const { return a[size_t(i) * ncols + j]; }

    MatrixX segment(int start, int n) const
    {
        assert(ncols == 1);
        MatrixX r(n);
        for (int i = 0; i < n; i++) r.a[i] = a[start + i];
        return r;
    }

    T*       data() { return a.data(); }
    const T* data() const { return a.data(); }

    RowRef<T>            row(int i) { return RowRef<T>(a.data() + size_t(i) * ncols, ncols); }
    RowRefBase<T, true>  row(int i) const
    {
        return RowRefBase<T, true>{a.data() + size_t(i) * ncols, ncols};
    }

    template <typename U>
    MatrixX<U> cast() const
    {
        MatrixX<U> r(nrows, ncols);
        for (size_t i = 0; i < a.size(); i++) r.a[i] = U(a[i]);
        return r;
    }

    T maxCoeff() const
    {
        assert(!a.empty());
        T m = a[0];
        for (size_t i = 1; i < a.size(); i++) m = m < a[i] ? a[i] : m;
        return m;
    }

    bool allFinite() const
    {
        for (size_t i = 0; i < a.size(); i++)
            if (!std::isfinite(a[i])) return false;
        return true;
    }
};

// ------------------------------------------------------------------------------------------------
// Column-pivoting Householder QR solve for a 3x3 system, ported from Eigen 3.4's
// ColPivHouseholderQR so that the smoother takes the same steps it always has.
// ------------------------------------------------------------------------------------------------

template <typename T>
Vector<T, 3> solve_col_piv_householder_qr(const Matrix33<T>& matrix, const Vector<T, 3>& rhs);

}  // namespace floatTetWild

#include <floattetwild/Matrix.ipp>
