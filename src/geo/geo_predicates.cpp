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

// This makes sure the compiler will not optimize y = a*x+b
// with fused multiply-add, this would break the exact
// predicates.
#ifdef GEO_COMPILER_MSVC
#pragma fp_contract(off)
#endif

#define FPG_UNCERTAIN_VALUE 0

#include "predicates_orient3d.h"

namespace {

    using namespace floatTetWild::geo;

    PCK::SOSMode SOS_mode_ = PCK::SOS_ADDRESS;

    /**
     * \brief Comparator class for nD points using lexicographic order.
     * \details Used by symbolic perturbations.
     */
    class LexicoCompare {
    public:

        /**
         * \brief LexicoCompare constructor.
         * \param[in] dim dimension of the points to compare.
         */
        LexicoCompare(index_t dim) : dim_(dim) {
        }

        /**
         * \brief Compares two points with respect to the lexicographic
         *  order.
         * \param[in] x , y pointers to the coordinates of the two points.
         * \retval true if x is strictly before y in the lexicographic order.
         * \retval false otherwise.
         */
        bool operator()(const double* x, const double* y) const {
            for(index_t i=0; i<dim_-1; ++i) {
                if(x[i] < y[i]) {
                    return true;
                }
                if(x[i] > y[i]) {
                    return false;
                }
            }
            return (x[dim_-1] < y[dim_-1]);
        }
    private:
        index_t dim_;
    };

    /**
     * \brief Compares two 3D points with respect to the lexicographic
     *  order.
     * \param[in] x , y pointers to the coordinates of the two 3D points.
     * \retval true if x is strictly before y in the lexicographic order.
     * \retval false otherwise.
     */
    bool lexico_compare_3d(const double* x, const double* y) {
        if(x[0] < y[0]) {
            return true;
        }
        if(x[0] > y[0]) {
            return false;
        }
        if(x[1] < y[1]) {
            return true;
        }
        if(x[1] > y[1]) {
            return false;
        }
        return x[2] < y[2];
    }

    /**
     * \brief Sorts an array of pointers to points.
     * \details set_SOS_mode() alters the behavior of this function.
     *  If set to PCK::SOS_ADDRESS, then just the addresses of the points
     *  are sorted. If set to PCK::SOS_LEXICO, then the points are sorted
     *  in function of the lexicographic order of their coordinates.
     * \param[in] begin a pointer to the first point.
     * \param[in] end one position past the pointer to the last point.
     * \param[in] dim the dimension of the points.
     */
    void SOS_sort(
        const double** begin, const double** end, index_t dim
    ) {
        if(SOS_mode_ == PCK::SOS_ADDRESS) {
            std::sort(begin, end);
        } else {
            if(dim == 3) {
                std::sort(begin, end, lexico_compare_3d);
            } else {
                std::sort(begin, end, LexicoCompare(dim));
            }
        }
    }

    /**
     * \brief Gets the maximum of 4 double precision numbers.
     * \param[in] x1 , x2 , x3 , x4 the four numbers.
     * \return the maximum.
     */
    inline double max4(double x1, double x2, double x3, double x4) {
#ifdef __SSE2__
        double result;
        __m128d X1 =_mm_load_sd(&x1);
        __m128d X2 =_mm_load_sd(&x2);
        __m128d X3 =_mm_load_sd(&x3);
        __m128d X4 =_mm_load_sd(&x4);
        X1 = _mm_max_sd(X1,X2);
        X3 = _mm_max_sd(X3,X4);
        X1 = _mm_max_sd(X1,X3);
        _mm_store_sd(&result, X1);
        return result;
#else
        return std::max(std::max(x1,x2),std::max(x3,x4));
#endif
    }

    /**
     * \brief Gets the minimum and maximum of 3 double precision numbers.
     * \param[in] x1 , x2 , x3 the three numbers.
     * \param[out] m the minimum
     * \param[out] M the maximum
     */
    inline void get_minmax3(
        double& m, double& M, double x1, double x2, double x3
    ) {
#ifdef __SSE2__
        __m128d X1 =_mm_load_sd(&x1);
        __m128d X2 =_mm_load_sd(&x2);
        __m128d X3 =_mm_load_sd(&x3);
        __m128d MIN12 = _mm_min_sd(X1,X2);
        __m128d MAX12 = _mm_max_sd(X1,X2);
        X1 = _mm_min_sd(MIN12, X3);
        X3 = _mm_max_sd(MAX12, X3);
        _mm_store_sd(&m, X1);
        _mm_store_sd(&M, X3);
#else
        m = std::min(std::min(x1,x2), x3);
        M = std::max(std::max(x1,x2), x3);
#endif
    }

    /**
     * \brief Arithmetic filter for the in_sphere_3d_SOS() predicate.
     * \details This filter was optimized by hand by Sylvain Pion
     *  (may be faster than FPG/PCK-generated filter).
     * Since it is used massively by Delaunay_3d, using the
     *   optimized version may be worth it.
     * \param[in] p first vertex of the tetrahedron
     * \param[in] q second vertex of the tetrahedron
     * \param[in] r third vertex of the tetrahedron
     * \param[in] s fourth vertex of the tetrahedron
     * \param[in] t point to be tested
     * \retval +1 if \p t was determined to be outside
     *   the circumsphere of \p p,\p q,\p r,\p s
     * \retval -1 if \p t was determined to be inside
     *   the circumsphere of  \p p,\p q,\p r,\p s
     * \retval 0 if the position of \p t could be be determined
     */
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

    /**
     * \brief Exact implementation of the side4_3d_SOS() predicate
     *  using low-level exact arithmetics API (expansion class).
     * \param[in] sos if true, applies symbolic perturbation when
     *  result is zero, else returns zero
     */
    Sign side4_3d_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* p4, bool sos = true
    ) {

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);
        const expansion& a14 = expansion_sq_dist(p1, p0, 3).negate();

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);
        const expansion& a24 = expansion_sq_dist(p2, p0, 3).negate();

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);
        const expansion& a34 = expansion_sq_dist(p3, p0, 3).negate();

        const expansion& a41 = expansion_diff(p4[0], p0[0]);
        const expansion& a42 = expansion_diff(p4[1], p0[1]);
        const expansion& a43 = expansion_diff(p4[2], p0[2]);
        const expansion& a44 = expansion_sq_dist(p4, p0, 3).negate();

        // This commented-out version does not reuse
        // the 2x2 minors.
/*
  const expansion& Delta1 = expansion_det3x3(
  a21, a22, a23,
  a31, a32, a33,
  a41, a42, a43
  );
  const expansion& Delta2 = expansion_det3x3(
  a11, a12, a13,
  a31, a32, a33,
  a41, a42, a43
  );
  const expansion& Delta3 = expansion_det3x3(
  a11, a12, a13,
  a21, a22, a23,
  a41, a42, a43
  );
  const expansion& Delta4 = expansion_det3x3(
  a11, a12, a13,
  a21, a22, a23,
  a31, a32, a33
  );
*/

        // Optimized version that reuses the 2x2 minors

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
        if(sos && r_sign == ZERO) {

            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;
            SOS_sort(p_sort, p_sort + 5, 3);
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

    namespace PCK {

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

        bool points_are_identical_3d(
            const double* p1,
            const double* p2
        ) {
            return
                (p1[0] == p2[0]) &&
                (p1[1] == p2[1]) &&
                (p1[2] == p2[2])
                ;
        }

        bool points_are_colinear_3d(
            const double* p1,
            const double* p2,
            const double* p3
        ) {
            // Colinearity is tested by using four coplanarity
            // tests with four points that are not coplanar.
            // TODO: use PCK::aligned_3d() instead (to be tested)
            static const double q000[3] = {0.0, 0.0, 0.0};
            static const double q001[3] = {0.0, 0.0, 1.0};
            static const double q010[3] = {0.0, 1.0, 0.0};
            static const double q100[3] = {1.0, 0.0, 0.0};
            return
                PCK::orient_3d(p1, p2, p3, q000) == ZERO &&
                PCK::orient_3d(p1, p2, p3, q001) == ZERO &&
                PCK::orient_3d(p1, p2, p3, q010) == ZERO &&
                PCK::orient_3d(p1, p2, p3, q100) == ZERO
                ;
        }

        namespace {
            // geogram made callers do this through GEO::initialize(). Here the predicates are
            // reached from library entry points that no caller has to set up, so run it during
            // static initialization instead. expansion::initialize() only writes two static
            // doubles, so it does not depend on any other translation unit's statics.
            struct Initializer {
                Initializer() {
                    expansion::initialize();
                }
            };
            Initializer initializer;
        }
    }
} }
