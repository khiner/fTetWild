// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/numerics/predicates.cpp
// Copied rather than reimplemented: the symbolic-perturbation tie-breaking is what decides
// Delaunay's output on the cospherical configurations a regular background grid is full of.
//
// Trimmed to the two predicates fTetWild reaches -- orient_3d and in_sphere_3d_SOS -- plus the
// exact cascade they bottom out in. The bodies are unchanged.

#include "geo_predicates.h"
#include "geo_multi_precision.h"

#include <algorithm>

// The filters below rely on the compiler not contracting y = a*x+b into a fused multiply-add,
// which would break them. Under gcc that needs -ffp-contract=off, and under MSVC
// #pragma fp_contract(off).

namespace {

    using namespace floatTetWild::geo;

    constexpr int FPG_UNCERTAIN_VALUE = 0;

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

    // Arithmetic filter for orient_3d(). Generated from source file orient3d.pck, unchanged.
    inline int orient_3d_filter(const double* p0, const double* p1, const double* p2, const double* p3) {
        double a11;
        a11 = (p1[0] - p0[0]);
        double a12;
        a12 = (p1[1] - p0[1]);
        double a13;
        a13 = (p1[2] - p0[2]);
        double a21;
        a21 = (p2[0] - p0[0]);
        double a22;
        a22 = (p2[1] - p0[1]);
        double a23;
        a23 = (p2[2] - p0[2]);
        double a31;
        a31 = (p3[0] - p0[0]);
        double a32;
        a32 = (p3[1] - p0[1]);
        double a33;
        a33 = (p3[2] - p0[2]);
        double Delta;
        Delta = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
        int int_tmp_result;
        double eps;
        double max1 = fabs(a11);
        if((max1 < fabs(a21)))
        {
            max1 = fabs(a21);
        }
        if((max1 < fabs(a31)))
        {
            max1 = fabs(a31);
        }
        double max2 = fabs(a12);
        if((max2 < fabs(a13)))
        {
            max2 = fabs(a13);
        }
        if((max2 < fabs(a22)))
        {
            max2 = fabs(a22);
        }
        if((max2 < fabs(a23)))
        {
            max2 = fabs(a23);
        }
        double max3 = fabs(a22);
        if((max3 < fabs(a23)))
        {
            max3 = fabs(a23);
        }
        if((max3 < fabs(a32)))
        {
            max3 = fabs(a32);
        }
        if((max3 < fabs(a33)))
        {
            max3 = fabs(a33);
        }
        double lower_bound_1;
        double upper_bound_1;
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if((max2 < lower_bound_1))
        {
            lower_bound_1 = max2;
        }
        else
        {
            if((max2 > upper_bound_1))
            {
                upper_bound_1 = max2;
            }
        }
        if((max3 < lower_bound_1))
        {
            lower_bound_1 = max3;
        }
        else
        {
            if((max3 > upper_bound_1))
            {
                upper_bound_1 = max3;
            }
        }
        if((lower_bound_1 < 1.63288018496748314939e-98))
        {
            return FPG_UNCERTAIN_VALUE;
        }
        else
        {
            if((upper_bound_1 > 5.59936185544450928309e+101))
            {
                return FPG_UNCERTAIN_VALUE;
            }
            eps = (5.11071278299732992696e-15 * ((max2 * max3) * max1));
            if((Delta > eps))
            {
                int_tmp_result = 1;
            }
            else
            {
                if((Delta < -eps))
                {
                    int_tmp_result = -1;
                }
                else
                {
                    return FPG_UNCERTAIN_VALUE;
                }
            }
        }
        return int_tmp_result;
    }

    // geogram writes these two with SSE2 intrinsics under __SSE2__, and the plain expressions
    // otherwise. It never included the intrinsics header, so only the plain path was ever
    // compiled here, and min and max of doubles are exact either way.
    inline double max4(double x1, double x2, double x3, double x4) {
        return std::max(std::max(x1,x2),std::max(x3,x4));
    }

    inline void get_minmax3(
        double& m, double& M, double x1, double x2, double x3
    ) {
        m = std::min(std::min(x1,x2), x3);
        M = std::max(std::max(x1,x2), x3);
    }

    // Arithmetic filter for in_sphere_3d_SOS(): +1 when \p t was determined to be outside the
    // circumsphere of \p p, \p q, \p r, \p s, -1 when inside, and FPG_UNCERTAIN_VALUE when the
    // filter could not tell. Optimized by hand by Sylvain Pion, which may beat the FPG/PCK
    // generated filter, and Delaunay_3d uses it massively enough for that to be worth it.
    inline int in_sphere_3d_filter_optim(
        const double* p, const double* q,
        const double* r, const double* s, const double* t
    ) {
        double ptx = p[0] - t[0];
        double pty = p[1] - t[1];
        double ptz = p[2] - t[2];
        double pt2 = geo_sqr(ptx) + geo_sqr(pty) + geo_sqr(ptz);

        double qtx = q[0] - t[0];
        double qty = q[1] - t[1];
        double qtz = q[2] - t[2];
        double qt2 = geo_sqr(qtx) + geo_sqr(qty) + geo_sqr(qtz);

        double rtx = r[0] - t[0];
        double rty = r[1] - t[1];
        double rtz = r[2] - t[2];
        double rt2 = geo_sqr(rtx) + geo_sqr(rty) + geo_sqr(rtz);

        double stx = s[0] - t[0];
        double sty = s[1] - t[1];
        double stz = s[2] - t[2];
        double st2 = geo_sqr(stx) + geo_sqr(sty) + geo_sqr(stz);

        // Compute the semi-static bound.
        double maxx = ::fabs(ptx);
        double maxy = ::fabs(pty);
        double maxz = ::fabs(ptz);

        double aqtx = ::fabs(qtx);
        double artx = ::fabs(rtx);
        double astx = ::fabs(stx);

        double aqty = ::fabs(qty);
        double arty = ::fabs(rty);
        double asty = ::fabs(sty);

        double aqtz = ::fabs(qtz);
        double artz = ::fabs(rtz);
        double astz = ::fabs(stz);

        maxx = max4(maxx, aqtx, artx, astx);
        maxy = max4(maxy, aqty, arty, asty);
        maxz = max4(maxz, aqtz, artz, astz);

        double eps = 1.2466136531027298e-13 * maxx * maxy * maxz;

        double min_max;
        double max_max;
        get_minmax3(min_max, max_max, maxx, maxy, maxz);

        double det = det4x4(
            ptx,pty,ptz,pt2,
            rtx,rty,rtz,rt2,
            qtx,qty,qtz,qt2,
            stx,sty,stz,st2
        );

        if (min_max < 1e-58)  { /* sqrt^5(min_double/eps) */
            // Protect against underflow in the computation of eps.
            return FPG_UNCERTAIN_VALUE;
        } else if (max_max < 1e61)  { /* sqrt^5(max_double/4 [hadamard]) */
            // Protect against overflow in the computation of det.
            eps *= (max_max * max_max);
            // Note: inverted as compared to CGAL
            //   CGAL: in_sphere_3d (called side_of_oriented_sphere())
            //      positive side is outside the sphere.
            //   PCK: in_sphere_3d : positive side is inside the sphere
            if (det > eps)  return -1;
            if (det < -eps) return  1;
        }

        return FPG_UNCERTAIN_VALUE;
    }

    // ================= side4 =========================================

    // The exact side4_3d_SOS() predicate, over the expansion class. Symbolic perturbation is
    // always applied. geogram could be asked to return zero instead, through a parameter no
    // caller here sets.
    Sign side4_3d_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* p4
    ) {

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);
        const expansion& a14 = expansion_sq_dist(p1, p0).negate();

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);
        const expansion& a24 = expansion_sq_dist(p2, p0).negate();

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);
        const expansion& a34 = expansion_sq_dist(p3, p0).negate();

        const expansion& a41 = expansion_diff(p4[0], p0[0]);
        const expansion& a42 = expansion_diff(p4[1], p0[1]);
        const expansion& a43 = expansion_diff(p4[2], p0[2]);
        const expansion& a44 = expansion_sq_dist(p4, p0).negate();

        // The four cofactors, sharing the 2x2 minors.

        const expansion& m12 = expansion_det2x2(a12,a13,a22,a23);
        const expansion& m13 = expansion_det2x2(a12,a13,a32,a33);
        const expansion& m14 = expansion_det2x2(a12,a13,a42,a43);
        const expansion& m23 = expansion_det2x2(a22,a23,a32,a33);
        const expansion& m24 = expansion_det2x2(a22,a23,a42,a43);
        const expansion& m34 = expansion_det2x2(a32,a33,a42,a43);

        const expansion& z11 = expansion_product(a21,m34);
        const expansion& z12 = expansion_product(a31,m24).negate();
        const expansion& z13 = expansion_product(a41,m23);
        const expansion& Delta1 = expansion_sum3(z11,z12,z13);

        const expansion& z21 = expansion_product(a11,m34);
        const expansion& z22 = expansion_product(a31,m14).negate();
        const expansion& z23 = expansion_product(a41,m13);
        const expansion& Delta2 = expansion_sum3(z21,z22,z23);

        const expansion& z31 = expansion_product(a11,m24);
        const expansion& z32 = expansion_product(a21,m14).negate();
        const expansion& z33 = expansion_product(a41,m12);
        const expansion& Delta3 = expansion_sum3(z31,z32,z33);

        const expansion& z41 = expansion_product(a11,m23);
        const expansion& z42 = expansion_product(a21,m13).negate();
        const expansion& z43 = expansion_product(a31,m12);
        const expansion& Delta4 = expansion_sum3(z41,z42,z43);

        Sign Delta4_sign = Delta4.sign();
        geo_assert(Delta4_sign != ZERO);

        const expansion& r_1 = expansion_product(Delta1, a14);
        const expansion& r_2 = expansion_product(Delta2, a24).negate();
        const expansion& r_3 = expansion_product(Delta3, a34);
        const expansion& r_4 = expansion_product(Delta4, a44).negate();
        const expansion& r = expansion_sum4(r_1, r_2, r_3, r_4);
        Sign r_sign = r.sign();

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {

            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;
            // Symbolic perturbation ordered by point address, which is what geogram's
            // SOS_ADDRESS mode does. Its other mode, sorting by the coordinates themselves,
            // was reachable only through set_SOS_mode(), which is not vendored.
            std::sort(p_sort, p_sort + 5);
            for(index_t i = 0; i < 5; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta2, Delta1);
                    const expansion& z2 = expansion_diff(Delta4, Delta3);
                    const expansion& z = expansion_sum(z1, z2);
                    Sign z_sign = z.sign();
                    if(z_sign != ZERO) {
                        return Sign(Delta4_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    Sign Delta1_sign = Delta1.sign();
                    if(Delta1_sign != ZERO) {
                        return Sign(Delta4_sign * Delta1_sign);
                    }
                } else if(p_sort[i] == p2) {
                    Sign Delta2_sign = Delta2.sign();
                    if(Delta2_sign != ZERO) {
                        return Sign(-Delta4_sign * Delta2_sign);
                    }
                } else if(p_sort[i] == p3) {
                    Sign Delta3_sign = Delta3.sign();
                    if(Delta3_sign != ZERO) {
                        return Sign(Delta4_sign * Delta3_sign);
                    }
                } else if(p_sort[i] == p4) {
                    return NEGATIVE;
                }
            }
        }
        return Sign(Delta4_sign * r_sign);
    }

    // ================= orient3d ======================================

    Sign orient_3d_exact(
        const double* p0, const double* p1,
        const double* p2, const double* p3
    ) {

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);

        const expansion& Delta = expansion_det3x3(
            a11, a12, a13, a21, a22, a23, a31, a32, a33
        );

        return Delta.sign();
    }


}

namespace floatTetWild {
namespace geo {

    Sign in_sphere_3d_SOS(
        const double* p0, const double* p1,
        const double* p2, const double* p3,
        const double* p4
    ) {
        // in_sphere_3d is simply implemented using side4_3d.
        // Both predicates are equivalent through duality as can
        // be easily seen:
        // side4_3d(p0,p1,p2,p3,p4) returns POSITIVE if
        //    d(q,p0) < d(q,p4)
        //    where q denotes the circumcenter of (p0,p1,p2,p3)
        // Note that d(q,p0) = R  (radius of circumscribed sphere)
        // In other words, side4_3d(p0,p1,p2,p3,p4) returns POSITIVE if
        //   d(q,p4) > R which means whenever p4 is not in the
        //   circumscribed sphere of (p0,p1,p2,p3).
        // Therefore:
        // in_sphere_3d(p0,p1,p2,p3,p4) = -side4_3d(p0,p1,p2,p3,p4)

        // This specialized filter supposes that orient_3d(p0,p1,p2,p3) > 0

        Sign result = Sign(in_sphere_3d_filter_optim(p0, p1, p2, p3, p4));

        if(result == 0) {
            result = side4_3d_exact_SOS(p0, p1, p2, p3, p4);
        }
        return Sign(-result);
    }

    Sign orient_3d(
        const double* p0, const double* p1,
        const double* p2, const double* p3
    ) {
        Sign result = Sign(orient_3d_filter(p0, p1, p2, p3));
        if(result == 0) {
            result = orient_3d_exact(p0, p1, p2, p3);
        }
        return result;
    }

} }
