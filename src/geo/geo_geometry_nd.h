// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/basic/geometry_nd.h
// Copied rather than reimplemented: point_triangle_squared_distance and friends, unchanged

#pragma once

#include "geo_geometry.h"

/**
 * \file geogram/basic/geometry_nd.h
 * \brief Geometric functions in arbitrary dimension
 */

namespace floatTetWild {
namespace geo {

    namespace Geom {

        /**
         * \brief Computes the squared distance between two nd points.
         * \param[in] p1 a pointer to the coordinates of the first point
         * \param[in] p2 a pointer to the coordinates of the second point
         * \param[in] dim dimension (number of coordinates of the points)
         * \return the squared distance between \p p1 and \p p2
         * \tparam COORD_T the numeric type of the point coordinates
         */
        template <class COORD_T>
        inline double distance2(
            const COORD_T* p1, const COORD_T* p2, coord_index_t dim
        ) {
            double result = 0.0;
            for(coord_index_t i = 0; i < dim; i++) {
                result += geo::geo_sqr(double(p2[i]) - double(p1[i]));
            }
            return result;
        }

        /**
         * \brief Computes the distance between two nd points.
         * \param[in] p1 a pointer to the coordinates of the first point
         * \param[in] p2 a pointer to the coordinates of the second point
         * \param[in] dim dimension (number of coordinates of the points)
         * \return the distance between \p p1 and \p p2
         * \tparam COORD_T the numeric type of the point coordinates
         */
        template <class COORD_T>
        inline double distance(
            const COORD_T* p1, const COORD_T* p2, coord_index_t dim
        ) {
            return ::sqrt(distance2(p1, p2, dim));
        }

        /**
         * \brief Computes the squared distance between two nd points.
         * \param[in] p1 first point
         * \param[in] p2 second point
         * \tparam VEC the class that represents the points. VEC needs to
         *   implement data(), that returns a pointer to the coordinates
         *   of the point.
         * \return the squared distance between \p p1 and \p p2
         */
        template <class VEC>
        inline double distance2(
            const VEC& p1, const VEC& p2
        ) {
            geo_debug_assert(p1.dimension() == p2.dimension());
            return distance2(
                p1.data(), p2.data(), coord_index_t(p1.dimension())
            );
        }

        /**

         * \brief Computes the distance between two nd points.
         * \param[in] p1 first point
         * \param[in] p2 second point
         * \tparam VEC the class that represents the points. VEC needs to
         *   implement data(), that returns a pointer to the coordinates
         *   of the point.
         * \return the distance between \p p1 and \p p2
         */
        template <class VEC>
        inline double distance(
            const VEC& p1, const VEC& p2
        ) {
            geo_debug_assert(p1.dimension() == p2.dimension());
            return distance(p1.data(), p2.data(), coord_index_t(p1.dimension()));
        }








        /**
         * \brief Computes the point closest to a given point in a nd segment
         * \param[in] point the query point
         * \param[in] V0 first extremity of the segment
         * \param[in] V1 second extremity of the segment
         * \param[out] closest_point the point closest to \p point in the
         *  segment [\p V0, \p V1]
         * \param[out] lambda0 barycentric coordinate of the closest point
         *  relative to \p V0
         * \param[out] lambda1 barycentric coordinate of the closest point
         *  relative to \p V1
         * \tparam VEC the class that represents the points.
         * \return the squared distance between the point and
         *  the segment [\p V0, \p V1]
         */
        template <class VEC>
        inline double point_segment_squared_distance(
            const VEC& point,
            const VEC& V0,
            const VEC& V1,
            VEC& closest_point,
            double& lambda0,
            double& lambda1
        ) {
            double l2 = distance2(V0,V1);
            double t = dot(point - V0, V1 - V0);
            if(t <= 0.0 || l2 == 0.0) {
                closest_point = V0;
                lambda0 = 1.0;
                lambda1 = 0.0;
                return distance2(point, V0);
            } else if(t > l2) {
                closest_point = V1;
                lambda0 = 0.0;
                lambda1 = 1.0;
                return distance2(point, V1);
            }
            lambda1 = t / l2;
            lambda0 = 1.0-lambda1;
            closest_point = lambda0 * V0 + lambda1 * V1;
            return distance2(point, closest_point);
        }


        /**
         * \brief Computes the point closest to a given point in a nd segment
         * \param[in] point the query point
         * \param[in] V0 first extremity of the segment
         * \param[in] V1 second extremity of the segment
         * \tparam VEC the class that represents the points.
         * \return the squared distance between the point and
         *  the segment [\p V0, \p V1]
         */
        template <class VEC>
        inline double point_segment_squared_distance(
            const VEC& point,
            const VEC& V0,
            const VEC& V1
        ) {
            VEC closest_point;
            double lambda0;
            double lambda1;
            return point_segment_squared_distance(
                point, V0, V1, closest_point, lambda0, lambda1
            );
        }

        /**
         * \brief Computes the point closest to a given point in a nd triangle
         * \details See
         *  http://www.geometrictools.com/LibMathematics/Distance/Distance.html
         * \param[in] point the query point
         * \param[in] V0 first vertex of the triangle
         * \param[in] V1 second vertex of the triangle
         * \param[in] V2 third vertex of the triangle
         * \param[out] closest_point the point closest to \p point in the
         *  triangle (\p V0, \p V1, \p V2)
         * \param[out] lambda0 barycentric coordinate of the closest point
         *  relative to \p V0
         * \param[out] lambda1 barycentric coordinate of the closest point
         *  relative to \p V1
         * \param[out] lambda2 barycentric coordinate of the closest point
         *  relative to \p V2
         * \tparam VEC the class that represents the points.
         * \return the squared distance between the point and
         *  the triangle (\p V0, \p V1, \p V2)
         */
        template <class VEC>
        inline double point_triangle_squared_distance(
            const VEC& point,
            const VEC& V0,
            const VEC& V1,
            const VEC& V2,
            VEC& closest_point,
            double& lambda0, double& lambda1, double& lambda2
        ) {
            VEC diff = V0 - point;
            VEC edge0 = V1 - V0;
            VEC edge1 = V2 - V0;
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
                VEC cur_closest;
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

        /**
         * \brief Computes the squared distance between a point and a nd
         *  triangle.
         * \details See
         *  http://www.geometrictools.com/LibMathematics/Distance/Distance.html
         * \param[in] p the query point
         * \param[in] q1 first vertex of the triangle
         * \param[in] q2 second vertex of the triangle
         * \param[in] q3 third vertex of the triangle
         * \tparam VEC the class that represents the points.
         * \return the squared distance between the point and
         *  the triangle (\p V0, \p V1, \p V2)
         */

        template <class VEC>
        inline double point_triangle_squared_distance(
            const VEC& p, const VEC& q1, const VEC& q2, const VEC& q3
        ) {
            VEC closest_point;
            double lambda1, lambda2, lambda3;
            return point_triangle_squared_distance(
                p, q1, q2, q3, closest_point, lambda1, lambda2, lambda3
            );
        }




    }
} }

