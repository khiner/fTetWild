// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/basic/determinant.h, geogram/basic/vecg.h, geogram/basic/geometry.h and
// geogram/basic/geometry_nd.h, which were four headers each including the next.
// Copied rather than reimplemented: point_triangle_squared_distance and friends, unchanged.
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

            // If the triangle is degenerate
            if(det < 1e-30) {
                double cur_l1, cur_l2;
                vec3 cur_closest;
                double result;
                double cur_dist = point_segment_squared_distance(
		    point, V0, V1, cur_closest, cur_l1, cur_l2
		);
                result = cur_dist;
                closest_point = cur_closest;
                lambda0 = cur_l1;
                lambda1 = cur_l2;
                lambda2 = 0.0;
                cur_dist = point_segment_squared_distance(
		    point, V0, V2, cur_closest, cur_l1, cur_l2
		);
                if(cur_dist < result) {
                    result = cur_dist;
                    closest_point = cur_closest;
                    lambda0 = cur_l1;
                    lambda2 = cur_l2;
                    lambda1 = 0.0;
                }
                cur_dist = point_segment_squared_distance(
		    point, V1, V2, cur_closest, cur_l1, cur_l2
		);
                if(cur_dist < result) {
                    result = cur_dist;
                    closest_point = cur_closest;
                    lambda1 = cur_l1;
                    lambda2 = cur_l2;
                    lambda0 = 0.0;
                }
                return result;
            }

            if(s + t <= det) {
                if(s < 0.0) {
                    if(t < 0.0) {   // region 4
                        if(b0 < 0.0) {
                            t = 0.0;
                            if(-b0 >= a00) {
                                s = 1.0;
                                sqrDistance = a00 + 2.0 * b0 + c;
                            } else {
                                s = -b0 / a00;
                                sqrDistance = b0 * s + c;
                            }
                        } else {
                            s = 0.0;
                            if(b1 >= 0.0) {
                                t = 0.0;
                                sqrDistance = c;
                            } else if(-b1 >= a11) {
                                t = 1.0;
                                sqrDistance = a11 + 2.0 * b1 + c;
                            } else {
                                t = -b1 / a11;
                                sqrDistance = b1 * t + c;
                            }
                        }
                    } else {  // region 3
                        s = 0.0;
                        if(b1 >= 0.0) {
                            t = 0.0;
                            sqrDistance = c;
                        } else if(-b1 >= a11) {
                            t = 1.0;
                            sqrDistance = a11 + 2.0 * b1 + c;
                        } else {
                            t = -b1 / a11;
                            sqrDistance = b1 * t + c;
                        }
                    }
                } else if(t < 0.0) {  // region 5
                    t = 0.0;
                    if(b0 >= 0.0) {
                        s = 0.0;
                        sqrDistance = c;
                    } else if(-b0 >= a00) {
                        s = 1.0;
                        sqrDistance = a00 + 2.0 * b0 + c;
                    } else {
                        s = -b0 / a00;
                        sqrDistance = b0 * s + c;
                    }
                } else {  // region 0
                    // minimum at interior point
                    double invDet = double(1.0) / det;
                    s *= invDet;
                    t *= invDet;
                    sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                        t * (a01 * s + a11 * t + 2.0 * b1) + c;
                }
            } else {
                double tmp0, tmp1, numer, denom;

                if(s < 0.0) {   // region 2
                    tmp0 = a01 + b0;
                    tmp1 = a11 + b1;
                    if(tmp1 > tmp0) {
                        numer = tmp1 - tmp0;
                        denom = a00 - 2.0 * a01 + a11;
                        if(numer >= denom) {
                            s = 1.0;
                            t = 0.0;
                            sqrDistance = a00 + 2.0 * b0 + c;
                        } else {
                            s = numer / denom;
                            t = 1.0 - s;
                            sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                                t * (a01 * s + a11 * t + 2.0 * b1) + c;
                        }
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
                        numer = tmp1 - tmp0;
                        denom = a00 - 2.0 * a01 + a11;
                        if(numer >= denom) {
                            t = 1.0;
                            s = 0.0;
                            sqrDistance = a11 + 2.0 * b1 + c;
                        } else {
                            t = numer / denom;
                            s = 1.0 - t;
                            sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                                t * (a01 * s + a11 * t + 2.0 * b1) + c;
                        }
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
                        denom = a00 - 2.0 * a01 + a11;
                        if(numer >= denom) {
                            s = 1.0;
                            t = 0.0;
                            sqrDistance = a00 + 2.0 * b0 + c;
                        } else {
                            s = numer / denom;
                            t = 1.0 - s;
                            sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                                t * (a01 * s + a11 * t + 2.0 * b1) + c;
                        }
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
