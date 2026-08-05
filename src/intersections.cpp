// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "intersections.h"

#include "Predicates.hpp"
#include "geo/geo_predicates.h"

#include <algorithm>

bool floatTetWild::line_line_intersection_2d(const std::array<Vector2, 2> &seg1, const std::array<Vector2, 2> &seg2,
                                             Scalar& t1, Scalar& t2) {
    const Scalar& x1 = seg1[0][0];
    const Scalar& y1 = seg1[0][1];
    const Scalar& x2 = seg1[1][0];
    const Scalar& y2 = seg1[1][1];

    const Scalar& x3 = seg2[0][0];
    const Scalar& y3 = seg2[0][1];
    const Scalar& x4 = seg2[1][0];
    const Scalar& y4 = seg2[1][1];

    Scalar d = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);
    if(d == 0)
        return false;
    t1 = ((y3 - y4) * (x1 - x3) + (x4 - x3) * (y1 - y3)) / d;
    t2 = ((y1 - y2) * (x1 - x3) + (x2 - x1) * (y1 - y3)) / d;
    return true;
}

floatTetWild::Scalar floatTetWild::p_seg_squared_dist_3d(const Vector3 &v, const Vector3 &a, const Vector3 &b){
    Vector3 av = v-a;
    Vector3 ab = b-a;
    if(av.dot(ab)<0)
        return av.squaredNorm();
    Vector3 bv = v-b;
    if(bv.dot(-ab)<0)
        return bv.squaredNorm();

    return (ab.cross(-av)).squaredNorm()/ab.squaredNorm();
}

bool floatTetWild::seg_plane_intersection(const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& n,
                                          Vector3& p) {
    Vector3 u = p2 - p1;
    Vector3 w = p1 - a;

    Scalar D = n.dot(u);
    Scalar d1 = -n.dot(w);

    if (D == 0)  // parallel to the plane, whether or not it lies in it
        return false;

    Scalar t = d1 / D;
    if (t <= 0 || t >= 1) {
        return false;
    }

    p = p1 + t * u;
    return true;
}

namespace floatTetWild {
namespace {
// oris[i * 3 + k] is the side of edge i of `edges` that ps[k] lies on. Returns true when some
// point is strictly inside, or when no point is strictly outside, both of which already answer
// the overlap question.
bool sides_of(const std::array<Vector2, 3>& edges, const std::array<Vector2, 3>& ps,
              std::array<int, 9>& oris) {
    std::array<int, 3> cnt_pos = {{0, 0, 0}};
    std::array<int, 3> cnt_neg = {{0, 0, 0}};
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            const int ori = Predicates::orient_2d(edges[i], edges[(i + 1) % 3], ps[k]);
            oris[i * 3 + k] = ori;
            if (ori == Predicates::ORI_POSITIVE)
                cnt_pos[k]++;
            else if (ori == Predicates::ORI_NEGATIVE)
                cnt_neg[k]++;
        }
    }
    for (int k = 0; k < 3; k++) {
        if (cnt_pos[k] == 3 || cnt_neg[k] == 3)
            return true;  // a vertex is strictly inside the other triangle
    }
    // the whole triangle is contained by the other one
    return std::find(oris.begin(), oris.end(), Predicates::ORI_NEGATIVE) == oris.end()
        || std::find(oris.begin(), oris.end(), Predicates::ORI_POSITIVE) == oris.end();
}
}  // namespace
}  // namespace floatTetWild

bool floatTetWild::is_tri_tri_cutted_2d(const std::array<Vector2, 3>& vs_tet, const std::array<Vector2, 3>& vs_tri) {
    std::array<int, 9> tri_tet;
    if (sides_of(vs_tri, vs_tet, tri_tet))
        return true;

    std::array<int, 9> tet_tri;
    if (sides_of(vs_tet, vs_tri, tet_tri))
        return true;

    for (int tri_e_id = 0; tri_e_id < 3; tri_e_id++) {
        for (int tet_e_id = 0; tet_e_id < 3; tet_e_id++) {
            if (is_crossing(tri_tet[tri_e_id * 3 + tet_e_id], tri_tet[tri_e_id * 3 + (tet_e_id + 1) % 3])
                && is_crossing(tet_tri[tet_e_id * 3 + tri_e_id], tet_tri[tet_e_id * 3 + (tri_e_id + 1) % 3]))
                return true;
        }
    }

    return false;
}

floatTetWild::OriCounts floatTetWild::count_orientations(const int* oris, int n) {
    return {int(std::count(oris, oris + n, int(Predicates::ORI_POSITIVE))),
            int(std::count(oris, oris + n, int(Predicates::ORI_NEGATIVE)))};
}

bool floatTetWild::is_p_inside_tri_2d(const Vector2& p, const std::array<Vector2, 3> &tri) {
    int oris[3];
    for (int i = 0; i < 3; i++)
        oris[i] = Predicates::orient_2d(tri[i], tri[(i + 1) % 3], p);

    const OriCounts cnt = count_orientations(oris, 3);
    return cnt.neg == 3 || cnt.pos == 3;  // strictly inside
}

// The axis the triangle's normal leans on most, which is the one to drop when flattening to 2d.
int floatTetWild::get_t(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2) {
    Vector3 n = tri_normal(p0, p1, p2);
    Scalar max = 0;
    int t = 0;
    for (int i = 0; i < 3; i++) {
        Scalar cos_a = abs(n[i]);
        if (cos_a > max) {
            max = cos_a;
            t = i;
        }
    }
    return t;
}

floatTetWild::Vector2 floatTetWild::to_2d(const Vector3 &p, int t) {
    return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
}

floatTetWild::Vector2 floatTetWild::to_2d(const Vector3 &p, const Vector3& n, const Vector3& pp, int t) {
    Scalar dist = n.dot(p - pp);
    return to_2d(p - dist * n, t);
}

bool floatTetWild::is_crossing(int s1, int s2) {
    return (s1 == Predicates::ORI_POSITIVE && s2 == Predicates::ORI_NEGATIVE)
        || (s2 == Predicates::ORI_POSITIVE && s1 == Predicates::ORI_NEGATIVE);
}

bool floatTetWild::is_tri_tri_cut(const Vector3& p1, const Vector3& p2, const Vector3& p3,
                                  const Vector3& q1, const Vector3& q2, const Vector3& q3,
                                  bool is_coplanar) {
    if(is_coplanar){
        int axis = get_t(p1, p2, p3);

        return is_tri_tri_cutted_2d({{to_2d(p1, axis), to_2d(p2, axis), to_2d(p3, axis)}}, {{to_2d(q1, axis), to_2d(q2, axis), to_2d(q3, axis)}});
    }

    std::array<Scalar, 3> s = {{0,0,0}}, t = {{0,0,0}};
    int result = tri_tri_intersection_test_3d(p1.data(), p2.data(), p3.data(),
                                              q1.data(), q2.data(), q3.data(),
                                              &s[0], &t[0]);
    if (result != 1) {
        return false;
    }

    // An intersection segment of zero length is a touch, not a cut.
    if (std::abs(s[0] - t[0]) <= SCALAR_ZERO && std::abs(s[1] - t[1]) <= SCALAR_ZERO && std::abs(s[2] - t[2]) <= SCALAR_ZERO)
        return false;

    return true;
}

void floatTetWild::get_bbox(std::initializer_list<Vector3> ps, Vector3& min, Vector3& max) {
    min = *ps.begin();
    max = *ps.begin();
    for (const Vector3& p : ps) {
        for (int j = 0; j < 3; j++) {
            if (p[j] < min[j])
                min[j] = p[j];
            if (p[j] > max[j])
                max[j] = p[j];
        }
    }
}

bool floatTetWild::is_bbox_intersected(const Vector3& min1, const Vector3& max1, const Vector3& min2, const Vector3& max2) {
    for (int j = 0; j < 3; j++) {
        if (min1[j] > max2[j] || max1[j] < min2[j])
            return false;
    }
    return true;
}

// Inside the tet or on its boundary.
bool floatTetWild::is_point_inside_tet(const Vector3& p, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t) {
    // p against each face, in the order that replaces one corner of the tet at a time.
    const int oris[4] = {Predicates::orient_3d(p, p1t, p2t, p3t),
                         Predicates::orient_3d(p0t, p, p2t, p3t),
                         Predicates::orient_3d(p0t, p1t, p, p3t),
                         Predicates::orient_3d(p0t, p1t, p2t, p)};

    const OriCounts cnt = count_orientations(oris, 4);
    return cnt.pos == 0 || cnt.neg == 0;
}

/*
*  Triangle-Triangle Overlap Test Routines
*  July, 2002
*  Updated December 2003
*
*  This file contains C implementation of algorithms for
*  performing two and three-dimensional triangle-triangle intersection test
*  The algorithms and underlying theory are described in
*
* "Fast and Robust Triangle-Triangle Overlap Test
*  Using Orientation Predicates"  P. Guigue - O. Devillers
*
*  Journal of Graphics Tools, 8(1), 2003
*
*  Other information are available from the Web page
*  http:<i>//www.acm.org/jgt/papers/GuigueDevillers03/
*
*  Trimmed to tri_tri_intersection_test_3d(), which computes the segment of
*  intersection when the triangles overlap and are not coplanar. Upstream also
*  had the 2d test and overlap-only versions of both.
*/

// modified by Aaron to better detect coplanarity

/* The sign of the volume of the tetrahedron pa, pb, pc, pd, inverted: which side of the plane
*  through the first three points the fourth is on. Exact, so a coplanar quadruple gives 0.
*/
static int sub_sub_cross_sub_dot(const double pa[3], const double pb[3], const double pc[3], const double pd[3]) {
    return -floatTetWild::geo::orient_3d(pa, pb, pc, pd);
}

/* some 3D macros */

#define CROSS(dest,v1,v2)                       \
               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; \
                        dest[1]=v1[1]-v2[1]; \
                        dest[2]=v1[2]-v2[2];

#define SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
                             dest[1] = alpha * v[1]; \
                             dest[2] = alpha * v[2];

/* Where the segment from A towards C crosses the plane with normal N through B, written into
*  dest. CONSTRUCT_INTERSECTION below is eight of these: one for each end of the intersection
*  segment, under each of its four cases.
*/
#define MEET_PLANE(dest,A,B,C,N) \
         SUB(v1,A,B) \
         SUB(v2,A,C) \
         alpha = DOT(v1,N) / DOT(v2,N); \
         SCALAR(v1,alpha,v2) \
         SUB(dest,A,v1)

/*
*
*  Three-dimensional Triangle-Triangle Intersection
*
*/

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

#define CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
  if (sub_sub_cross_sub_dot(q1, r2, p1, p2) > 0) {\
    if (sub_sub_cross_sub_dot(r1, r2, p1, p2) <=  0) {\
      if (sub_sub_cross_sub_dot(r1, q2, p1, p2) >  0) {\
         MEET_PLANE(source,p1,p2,r1,N2) \
         MEET_PLANE(target,p2,p1,r2,N1) \
         return 1; \
      }\
      else { \
         MEET_PLANE(source,p2,p1,q2,N1) \
         MEET_PLANE(target,p2,p1,r2,N1) \
         return 1; \
      } \
    } \
    else { \
      return 0; \
    } \
  } \
  else { \
    if (sub_sub_cross_sub_dot(q1, q2, p1, p2) <  0) {\
      return 0; \
    } \
    else { \
      if (sub_sub_cross_sub_dot(r1, q2, p1, p2) >= 0) {\
        MEET_PLANE(source,p1,p2,r1,N2) \
        MEET_PLANE(target,p1,p2,q1,N2) \
        return 1; \
      } \
      else { \
        MEET_PLANE(source,p2,p1,q2,N1) \
        MEET_PLANE(target,p1,p2,q1,N2) \
        return 1; \
      } \
    } \
  } \
}

#define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
  if (dp2 > 0) { \
     if (dq2 > 0) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
     else if (dr2 > 0) CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
     else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
  else if (dp2 < 0) { \
    if (dq2 < 0) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0) CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
  } else { \
    if (dq2 < 0) { \
      if (dr2 >= 0)  CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
      else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0) { \
      if (dr2 > 0) CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
      else  CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
      if (dr2 > 0) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
      else if (dr2 < 0) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
      else { \
       return -1;\
     } \
  }} }

/*
   Computes the segment of intersection of the two triangles if it exists.
   source and target return the endpoints of the line segment of intersection.
   Returns -1 when the triangles are coplanar.
*/

int tri_tri_intersection_test_3d(const double p1[3], const double q1[3], const double r1[3],
                                 const double p2[3], const double q2[3], const double r2[3],
                                 double source[3], double target[3] )

{
    int dp1, dq1, dr1, dp2, dq2, dr2;
    double v1[3], v2[3];
    double N1[3], N2[3];
    double alpha;

    SUB(v1,q1,p1)
    SUB(v2,r1,p1)
    CROSS(N1,v1,v2)

    SUB(v1,p2,r2)
    SUB(v2,q2,r2)
    CROSS(N2,v1,v2)

    // Compute distance signs  of p1, q1 and r1
    // to the plane of triangle(p2,q2,r2)

    dp1 = sub_sub_cross_sub_dot(p2, q2, r2, p1);
    dq1 = sub_sub_cross_sub_dot(p2, q2, r2, q1);
    dr1 = sub_sub_cross_sub_dot(p2, q2, r2, r1);

    if (((dp1 * dq1) > 0) && ((dp1 * dr1) > 0))  return 666;

    // Compute distance signs  of p2, q2 and r2
    // to the plane of triangle(p1,q1,r1)

    dp2 = sub_sub_cross_sub_dot(p1, q1, r1, p2);
    dq2 = sub_sub_cross_sub_dot(p1, q1, r1, q2);
    dr2 = sub_sub_cross_sub_dot(p1, q1, r1, r2);

    if (((dp2 * dq2) > 0) && ((dp2 * dr2) > 0)) return 666;

    // Permutation in a canonical form of T1's vertices

    if (dp1 > 0) {
        if (dq1 > 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
        else if (dr1 > 0) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)

        else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    } else if (dp1 < 0) {
        if (dq1 < 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
        else if (dr1 < 0) TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
        else TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
    } else {
        if (dq1 < 0) {
            if (dr1 >= 0) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
            else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
        }
        else if (dq1 > 0) {
            if (dr1 > 0) TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
            else TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
        }
        else  {
            if (dr1 > 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
            else if (dr1 < 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
            else {
                // triangles are co-planar
                return -1;
            }
        }
    }
}
