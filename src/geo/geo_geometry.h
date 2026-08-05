// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/basic/determinant.h, geogram/basic/vecg.h, geogram/basic/geometry.h and
// geogram/basic/geometry_nd.h, which were four headers each including the next.
// Copied rather than reimplemented: point_triangle_squared_distance and friends. Its arithmetic is
// unchanged, but the edge clamp its regions 3, 4 and 5 each wrote out, and the interior quadratic
// four of them end at, are now one lambda each.
//
// geogram writes all of this generically: determinants over any T, vectors over a dimension and a
// coordinate type, and the Geom distances over any class with data() and dimension(). Every one of
// those is fixed here -- 3d points of doubles -- so the parameters are gone. The Geom helpers that
// lived here as well, tetra_signed_volume in both its shapes, had no caller and are gone too.

#pragma once

#include "geo_basic.h"

#include <cmath>

namespace floatTetWild {
namespace geo {

    inline double det2x2(
        double a11, double a12,
        double a21, double a22
    ) {
        return a11*a22-a12*a21 ;
    }

    inline double det3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33
    ) {
        return
            a11*det2x2(a22,a23,a32,a33)
            -a21*det2x2(a12,a13,a32,a33)
            +a31*det2x2(a12,a13,a22,a23);
    }

    inline double det4x4(
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44
    ) {
        double m12 = a21*a12 - a11*a22;
        double m13 = a31*a12 - a11*a32;
        double m14 = a41*a12 - a11*a42;
        double m23 = a31*a22 - a21*a32;
        double m24 = a41*a22 - a21*a42;
        double m34 = a41*a32 - a31*a42;

        double m123 = m23*a13 - m13*a23 + m12*a33;
        double m124 = m24*a13 - m14*a23 + m12*a43;
        double m134 = m34*a13 - m14*a33 + m13*a43;
        double m234 = m34*a23 - m24*a33 + m23*a43;

        return (m234*a14 - m134*a24 + m124*a34 - m123*a44);
    }

    /************************************************************************/

    // A point or a vector in 3d. Syntax is (mostly) compatible with GLSL.
    class vec3 {
    public:
        vec3() :
            x(0.0),
            y(0.0),
            z(0.0) {
        }

        vec3(double x_in, double y_in, double z_in) :
            x(x_in),
            y(y_in),
            z(z_in) {
        }

        inline double length2() const {
            return x * x + y * y + z * z;
        }

        inline double length() const {
            return sqrt(x * x + y * y + z * z);
        }

        // As one expression, which is not the same double as summing the three terms in a loop:
        // see Geom::distance2() below.
        inline double distance2(const vec3& v) const {
            double dx = v.x - x;
            double dy = v.y - y;
            double dz = v.z - z;
            return dx * dx + dy * dy + dz * dz;
        }

        inline double distance(const vec3& v) const {
            return sqrt(distance2(v));
        }

        inline vec3 operator+ (const vec3& v) const {
            return vec3(x + v.x, y + v.y, z + v.z);
        }

        inline vec3 operator- (const vec3& v) const {
            return vec3(x - v.x, y - v.y, z - v.z);
        }

        inline vec3 operator* (double s) const {
            return vec3(x * s, y * s, z * s);
        }

        inline vec3 operator/ (double s) const {
            return vec3(x / s, y / s, z / s);
        }

        double* data() {
            return &x;
        }

        const double* data() const {
            return &x;
        }

        inline double& operator[] (index_t i) {
            geo_debug_assert(i < 3);
            return data()[i];
        }

        inline const double& operator[] (index_t i) const {
            geo_debug_assert(i < 3);
            return data()[i];
        }

        double x;
        double y;
        double z;
    };

    inline double length(const vec3& v) {
        return v.length();
    }

    inline double length2(const vec3& v) {
        return v.length2();
    }

    inline double distance2(const vec3& v1, const vec3& v2) {
        return v2.distance2(v1);
    }

    inline double distance(const vec3& v1, const vec3& v2) {
        return v2.distance(v1);
    }

    inline vec3 operator* (double s, const vec3& v) {
        return vec3(s * v.x, s * v.y, s * v.z);
    }

    // Undefined when the norm is 0.
    inline vec3 normalize(const vec3& v) {
        double s = length(v);
        if(s > 1e-30) {
            s = 1.0 / s;
        }
        return s * v;
    }

    inline double dot(const vec3& v1, const vec3& v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    inline vec3 cross(const vec3& v1, const vec3& v2) {
        return vec3(
            det2x2(v1.y, v2.y, v1.z, v2.z),
            det2x2(v1.z, v2.z, v1.x, v2.x),
            det2x2(v1.x, v2.x, v1.y, v2.y)
        );
    }

    /*******************************************************************/

    // Axis-aligned bounding box.
    class Box {
    public:
        double xyz_min[3];
        double xyz_max[3];
    };

    // The smallest Box that encloses \p B1 and \p B2.
    inline void bbox_union(Box& target, const Box& B1, const Box& B2) {
        for(index_t c = 0; c < 3; ++c) {
            target.xyz_min[c] = std::min(B1.xyz_min[c], B2.xyz_min[c]);
            target.xyz_max[c] = std::max(B1.xyz_max[c], B2.xyz_max[c]);
        }
    }

    /******************************************************************/

    namespace Geom {

        // The squared distance between two 3d points, accumulated one coordinate at a time. That
        // is a different double from geo::distance2() above, which sums the three terms as one
        // expression, because the compiler contracts the two differently. Which one a call site
        // wants is not interchangeable, so both are spelled out qualified wherever they are used.
        inline double distance2(const double* p1, const double* p2) {
            double result = 0.0;
            for(index_t i = 0; i < 3; i++) {
                result += geo_sqr(p2[i] - p1[i]);
            }
            return result;
        }

        inline double distance2(const vec3& p1, const vec3& p2) {
            return distance2(p1.data(), p2.data());
        }

        // The point of segment [\p V0, \p V1] closest to \p point, its barycentric coordinates
        // relative to the two extremities, and the squared distance to it.
        inline double point_segment_squared_distance(
            const vec3& point,
            const vec3& V0,
            const vec3& V1,
            vec3& closest_point,
            double& lambda0,
            double& lambda1
        ) {
            double l2 = geo::distance2(V0,V1);
            double t = dot(point - V0, V1 - V0);
            if(t <= 0.0 || l2 == 0.0) {
                closest_point = V0;
                lambda0 = 1.0;
                lambda1 = 0.0;
                return geo::distance2(point, V0);
            } else if(t > l2) {
                closest_point = V1;
                lambda0 = 0.0;
                lambda1 = 1.0;
                return geo::distance2(point, V1);
            }
            lambda1 = t / l2;
            lambda0 = 1.0-lambda1;
            closest_point = lambda0 * V0 + lambda1 * V1;
            return geo::distance2(point, closest_point);
        }

        // The same for the triangle (\p V0, \p V1, \p V2). See
        // http://www.geometrictools.com/LibMathematics/Distance/Distance.html
        inline double point_triangle_squared_distance(
            const vec3& point,
            const vec3& V0,
            const vec3& V1,
            const vec3& V2,
            vec3& closest_point,
            double& lambda0, double& lambda1, double& lambda2
        ) {
            vec3 diff = V0 - point;
            vec3 edge0 = V1 - V0;
            vec3 edge1 = V2 - V0;
            double a00 = length2(edge0);
            double a01 = dot(edge0, edge1);
            double a11 = length2(edge1);
            double b0 = dot(diff, edge0);
            double b1 = dot(diff, edge1);
            double c = length2(diff);
            double det = ::fabs(a00 * a11 - a01 * a01);
            double s = a01 * b1 - a11 * b0;
            double t = a01 * b0 - a00 * b1;
            double sqrDistance;

            // The two shapes the regions below reduce to, spelled once. The arithmetic in each is
            // exactly what those regions wrote out longhand.

            // The nearest point on one of the two edges leaving V0: \p u runs along the edge with
            // squared length \p a and projection \p b, clamped to it, and \p other goes to 0.
            const auto on_edge = [&c](double b, double a, double& u, double& other) {
                other = 0.0;
                if(b >= 0.0) {
                    u = 0.0;
                    return c;
                }
                if(-b >= a) {
                    u = 1.0;
                    return a + 2.0 * b + c;
                }
                u = -b / a;
                return b * u + c;
            };

            // The squared distance at barycentric coordinates (\p su, \p tv).
            const auto sq_dist_at = [&](double su, double tv) {
                return su * (a00 * su + a01 * tv + 2.0 * b0) +
                    tv * (a01 * su + a11 * tv + 2.0 * b1) + c;
            };

            // The nearest point on the edge from V1 to V2, at \p numer along it, clamped to the
            // end where \p u is 1 and \p other is 0. \p a_end and \p b_end are the coefficients
            // at that end. Regions 1, 2 and 6 all finish here.
            const auto on_diagonal = [&](
                double numer, double& u, double& other, double a_end, double b_end
            ) {
                const double denom = a00 - 2.0 * a01 + a11;
                if(numer >= denom) {
                    u = 1.0;
                    other = 0.0;
                    return a_end + 2.0 * b_end + c;
                }
                u = numer / denom;
                other = 1.0 - u;
                return sq_dist_at(s, t);
            };

            // If the triangle is degenerate, the answer is the nearest of its three edges. The
            // first one always takes, so that the winner is one of the three even under a NaN.
            if(det < 1e-30) {
                double result = 0.0;
                bool first = true;
                const auto try_edge = [&](
                    const vec3& A, const vec3& B, double& la, double& lb, double& lz
                ) {
                    double cur_l1, cur_l2;
                    vec3 cur_closest;
                    double cur_dist = point_segment_squared_distance(
                        point, A, B, cur_closest, cur_l1, cur_l2
                    );
                    if(first || cur_dist < result) {
                        result = cur_dist;
                        closest_point = cur_closest;
                        la = cur_l1;
                        lb = cur_l2;
                        lz = 0.0;
                    }
                    first = false;
                };
                try_edge(V0, V1, lambda0, lambda1, lambda2);
                try_edge(V0, V2, lambda0, lambda2, lambda1);
                try_edge(V1, V2, lambda1, lambda2, lambda0);
                return result;
            }

            if(s + t <= det) {
                if(s < 0.0) {
                    if(t < 0.0) {   // region 4
                        // b0 < 0.0 already rules out on_edge()'s first case.
                        sqrDistance = b0 < 0.0 ? on_edge(b0, a00, s, t)
                                               : on_edge(b1, a11, t, s);
                    } else {  // region 3
                        sqrDistance = on_edge(b1, a11, t, s);
                    }
                } else if(t < 0.0) {  // region 5
                    sqrDistance = on_edge(b0, a00, s, t);
                } else {  // region 0
                    // minimum at interior point
                    double invDet = double(1.0) / det;
                    s *= invDet;
                    t *= invDet;
                    sqrDistance = sq_dist_at(s, t);
                }
            } else {
                double tmp0, tmp1, numer;

                if(s < 0.0) {   // region 2
                    tmp0 = a01 + b0;
                    tmp1 = a11 + b1;
                    if(tmp1 > tmp0) {
                        sqrDistance = on_diagonal(tmp1 - tmp0, s, t, a00, b0);
                    } else {
                        s = 0.0;
                        if(tmp1 <= 0.0) {
                            t = 1.0;
                            sqrDistance = a11 + 2.0 * b1 + c;
                        }
                        else if(b1 >= 0.0) {
                            t = 0.0;
                            sqrDistance = c;
                        } else {
                            t = -b1 / a11;
                            sqrDistance = b1 * t + c;
                        }
                    }
                } else if(t < 0.0) {  // region 6
                    tmp0 = a01 + b1;
                    tmp1 = a00 + b0;
                    if(tmp1 > tmp0) {
                        sqrDistance = on_diagonal(tmp1 - tmp0, t, s, a11, b1);
                    } else {
                        t = 0.0;
                        if(tmp1 <= 0.0) {
                            s = 1.0;
                            sqrDistance = a00 + 2.0 * b0 + c;
                        } else if(b0 >= 0.0) {
                            s = 0.0;
                            sqrDistance = c;
                        } else {
                            s = -b0 / a00;
                            sqrDistance = b0 * s + c;
                        }
                    }
                } else { // region 1
                    numer = a11 + b1 - a01 - b0;
                    if(numer <= 0.0) {
                        s = 0.0;
                        t = 1.0;
                        sqrDistance = a11 + 2.0 * b1 + c;
                    } else {
                        sqrDistance = on_diagonal(numer, s, t, a00, b0);
                    }
                }
            }

            // Account for numerical round-off error.
            if(sqrDistance < 0.0) {
                sqrDistance = 0.0;
            }

            closest_point = V0 + s * edge0 + t * edge1;
            lambda0 = 1.0 - s - t;
            lambda1 = s;
            lambda2 = t;
            return sqrDistance;
        }

        // \overload Without the closest point or the barycentric coordinates.
        inline double point_triangle_squared_distance(
            const vec3& p, const vec3& q1, const vec3& q2, const vec3& q3
        ) {
            vec3 closest_point;
            double lambda1, lambda2, lambda3;
            return point_triangle_squared_distance(
                p, q1, q2, q3, closest_point, lambda1, lambda2, lambda3
            );
        }

    }

} }
