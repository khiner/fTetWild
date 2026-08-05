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

    typedef double real;

#include "geo/geo_predicates.h"

/* The sign of the volume of the tetrahedron pa, pb, pc, pd, inverted: which side of the plane
*  through the first three points the fourth is on. Exact, so a coplanar quadruple gives 0.
*/
static int sub_sub_cross_sub_dot(const real pa[3], const real pb[3], const real pc[3], const real pd[3]) {
    return -floatTetWild::geo::PCK::orient_3d(pa, pb, pc, pd);
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

int tri_tri_intersection_test_3d(const real p1[3], const real q1[3], const real r1[3],
                                 const real p2[3], const real q2[3], const real r2[3],
                                 real source[3], real target[3] )

{
    int dp1, dq1, dr1, dp2, dq2, dr2;
    real v1[3], v2[3];
    real N1[3], N2[3];
    real alpha;

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
