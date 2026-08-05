// The fast winding number for triangle soups.
//
// Copyright (c) 2018 Side Effects Software Inc. MIT licence, reproduced in FastWindingNumber.cpp
// and in THIRD_PARTY.md.

#ifndef FLOATTETWILD_FASTWINDINGNUMBER_H
#define FLOATTETWILD_FASTWINDINGNUMBER_H

#include <cmath>
#include <vector>

namespace floatTetWild {
namespace FastWindingNumber {

// A 3-vector of whatever the solid angle is being accumulated in: plain floats for the exact
// per-triangle part, four-at-a-time lanes for the part that approximates a BVH node's four
// children together. The arithmetic is written the way the paper's implementation wrote it,
// down to reciprocal-multiply for scalar division, because the mesher's output depends on it.
template<typename T>
struct Vec3
{
    T v[3];

    Vec3() {}
    explicit Vec3(const T a[3])
    {
        for (int i = 0; i < 3; ++i)
            v[i] = a[i];
    }

    const T& operator[](int i) const { return v[i]; }
    T&       operator[](int i) { return v[i]; }

    void operator=(T a)
    {
        for (int i = 0; i < 3; ++i)
            v[i] = a;
    }
    void operator+=(const Vec3& r)
    {
        for (int i = 0; i < 3; ++i)
            v[i] += r.v[i];
    }
    void operator-=(const Vec3& r)
    {
        for (int i = 0; i < 3; ++i)
            v[i] -= r.v[i];
    }
    void operator*=(const Vec3& r)
    {
        for (int i = 0; i < 3; ++i)
            v[i] *= r.v[i];
    }
    void operator*=(T r)
    {
        for (int i = 0; i < 3; ++i)
            v[i] *= r;
    }
    void operator/=(T r)
    {
        r = 1 / r;
        for (int i = 0; i < 3; ++i)
            v[i] *= r;
    }

    Vec3 operator+(const Vec3& r) const
    {
        Vec3 o;
        for (int i = 0; i < 3; ++i)
            o.v[i] = v[i] + r.v[i];
        return o;
    }
    Vec3 operator-(const Vec3& r) const
    {
        Vec3 o;
        for (int i = 0; i < 3; ++i)
            o.v[i] = v[i] - r.v[i];
        return o;
    }
    Vec3 operator*(const Vec3& r) const
    {
        Vec3 o;
        for (int i = 0; i < 3; ++i)
            o.v[i] = v[i] * r.v[i];
        return o;
    }
    Vec3 operator*(T r) const
    {
        Vec3 o;
        for (int i = 0; i < 3; ++i)
            o.v[i] = v[i] * r;
        return o;
    }
    Vec3 operator/(T r) const
    {
        Vec3 o;
        r = 1 / r;
        for (int i = 0; i < 3; ++i)
            o.v[i] = v[i] * r;
        return o;
    }

    T length2() const
    {
        T a0(v[0]);
        T result(a0 * a0);
        for (int i = 1; i < 3; ++i) {
            T ai(v[i]);
            result += ai * ai;
        }
        return result;
    }
    T length() const { return std::sqrt(length2()); }

    T dot(const Vec3& r) const
    {
        T result(v[0] * r.v[0]);
        for (int i = 1; i < 3; ++i)
            result += v[i] * r.v[i];
        return result;
    }
};

template<typename T>
Vec3<T> operator*(T a, const Vec3<T>& b)
{
    Vec3<T> o;
    for (int i = 0; i < 3; ++i)
        o[i] = a * b[i];
    return o;
}

template<typename T>
T dot(const Vec3<T>& a, const Vec3<T>& b)
{ return a.dot(b); }

using Vec3f = Vec3<float>;

// The signed solid angle each query point subtends against the triangle soup, approximated by
// "Fast Winding Numbers for Soups and Clouds" [Barill et al. 2018]: a four-wide BVH over the
// triangles, each node carrying a second order Taylor expansion of the angle its triangles
// subtend, descended only where the expansion is not accurate enough. Divide by 4 pi to get a
// winding number.
std::vector<float> solid_angles(const std::vector<Vec3f>& points,
                                const std::vector<int>&   triangles,
                                const std::vector<Vec3f>& queries);

}  // namespace FastWindingNumber
}  // namespace floatTetWild

#endif  // FLOATTETWILD_FASTWINDINGNUMBER_H
