// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/numerics/multi_precision.cpp
// Copied rather than reimplemented: Shewchuk expansion arithmetic, unchanged

#include "geo_multi_precision.h"

namespace {

    using namespace floatTetWild::geo;

    // The helpers below are Shewchuk's, in his naming: two_one means a length 2 expansion and a
    // double, and the trailing numbers on the outputs run from high magnitude to low.

    // Each of the next three turns an exact operation on two doubles into a length 2 expansion:
    // x is the rounded result and y the part that did not fit. two_diff, their sibling, stays in
    // the header for assign_diff().
    inline void two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        double bvirt = x - a;
        double avirt = x - bvirt;
        double bround = b - bvirt;
        double around = a - avirt;
        y = around + bround;
    }

    // fma() is exact, so it recovers the rounding error of the product outright. Shewchuk's
    // splitting version does the same in seven operations and is what geogram falls back to
    // without hardware fused multiply-add; nothing here needs that path.
    //
    // Note: under gcc, automatic generation of fma() for a*b+c has to be deactivated with
    // -ffp-contract=off, else it breaks functions such as fast_expansion_sum_zeroelim().
    inline void two_product(double a, double b, double& x, double& y) {
        x = a*b;
        y = fma(a,b,-x);
    }

    inline void square(double a, double& x, double& y) {
        x = a*a;
        y = fma(a,a,-x);
    }

    // Like two_sum, but only correct when |a| > |b|.
    inline void fast_two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        double bvirt = x - a;
        y = b - bvirt;
    }

    inline void two_one_sum(
        double a1, double a0, double b, double& x2, double& x1, double& x0
    ) {
        double _i;
        two_sum(a0, b, _i, x0);
        two_sum(a1, _i, x2, x1);
    }

    inline void two_two_sum(
        double a1, double a0, double b1, double b0,
        double& x3, double& x2, double& x1, double& x0
    ) {
        double _j, _0;
        two_one_sum(a1, a0, b0, _j, _0, x0);
        two_one_sum(_j, _0, b1, x3, x2, x1);
    }

    // The square of a length 2 expansion, into the six doubles at \p x. Squaring is quicker than
    // multiplying two different length 2 expansions, and needs six components rather than eight.
    inline void two_square(
        double a1, double a0,
        double* x
    ) {
        double _0, _1, _2;
        double _j, _k, _l;
        square(a0, _j, x[0]);
        _0 = a0 + a0;
        two_product(a1, _0, _k, _1);
        two_one_sum(_k, _1, _j, _l, _2, x[1]);
        square(a1, _j, _1);
        two_two_sum(_j, _1, _l, _2, x[5], x[4], x[3], x[2]);
    }

    // (a - b)^2, exactly, in the six components of \p e.
    inline void assign_sq_diff(expansion& e, double a, double b) {
        double d0, d1;
        two_diff(a, b, d1, d0);
        two_square(d1, d0, &e[0]);
        e.set_length(6);
    }

    // The product of two length 2 expansions, into the eight doubles at \p x.
    void two_two_product(
        const double* a,
        const double* b,
        double* x
    ) {
        double _0, _1, _2;
        double _i, _j, _k, _l, _m, _n;

        two_product(a[0],b[0],_i,x[0]);
        two_product(a[1],b[0],_j,_0);
        two_sum(_i, _0, _k, _1);
        fast_two_sum(_j, _k, _l, _2);
        two_product(a[0], b[1], _i, _0);
        two_sum(_1, _0, _k, x[1]);
        two_sum(_2, _k, _j, _1);
        two_sum(_l, _j, _m, _2);
        two_product(a[1], b[1], _j, _0);
        two_sum(_i, _0, _n, _0);
        two_sum(_1, _0, _i, x[2]);
        two_sum(_2, _i, _k, _1);
        two_sum(_m, _k, _l, _2);
        two_sum(_j, _n, _k, _0);
        two_sum(_1, _0, _j, x[3]);
        two_sum(_2, _j, _i, _1);
        two_sum(_l, _i, _m, _2);
        two_sum(_1, _k, _i, x[4]);
        two_sum(_2, _i, _k, x[5]);
        two_sum(_m, _k, x[7], x[6]);
    }

}

namespace floatTetWild {
namespace geo {

    // This and fast_expansion_combine_zeroelim below are adapted from Jonathan Shewchuk's code,
    // see either version of his paper for details. In each, \p h cannot be \p e or \p f, and zero
    // components are dropped from the result.

    // h = b * e. Maintains the nonoverlapping property, and with round-to-even (as with IEEE 754)
    // the strongly nonoverlapping and nonadjacent properties as well.
    void scale_expansion_zeroelim(
        const expansion& e, double b, expansion& h
    ) {
        double Q, sum;
        double hh;
        double product1;
        double product0;
        index_t eindex, hindex;

        index_t elen = e.length();

        // Sanity check: e and h cannot be the same.

        two_product(e[0], b, Q, hh);

        hindex = 0;
        if(hh != 0) {
            h[hindex++] = hh;
        }
        for(eindex = 1; eindex < elen; eindex++) {
            double enow = e[eindex];
            two_product(enow, b,  product1, product0);
            two_sum(Q, product0, sum, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
            fast_two_sum(product1, sum, Q, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    // h = e +/- f. The two are the same walk over the two expansions, merging by magnitude, and
    // differ only in that the difference negates every component of f as it reads it. Negation is
    // exact, so both are Shewchuk's fast_expansion_sum_zeroelim() with his arithmetic unchanged.
    // With round-to-even, the sum maintains the strongly nonoverlapping property, but NOT the
    // nonoverlapping or nonadjacent ones.
    template <bool Negate>
    void fast_expansion_combine_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        double Q;
        double Qnew;
        double hh;
        index_t eindex, findex, hindex;
        double enow, fnow;
        index_t elen = e.length();
        index_t flen = f.length();

        // sanity check: h cannot be e or f

        // Reading f through here is the whole of the difference between the two operations.
        const auto f_at = [&f](index_t i) { return Negate ? -f[i] : f[i]; };

        enow = e[0];
        fnow = f_at(0);
        eindex = findex = 0;
        if((fnow > enow) == (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = f_at(++findex);
        }
        hindex = 0;
        if((eindex < elen) && (findex < flen)) {
            if((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = e[++eindex];
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = f_at(++findex);
            }
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
            while((eindex < elen) && (findex < flen)) {
                if((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = e[++eindex];
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = f_at(++findex);
                }
                Q = Qnew;
                if(hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = e[++eindex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while(findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = f_at(++findex);
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    void fast_expansion_sum_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        fast_expansion_combine_zeroelim<false>(e, f, h);
    }

    void fast_expansion_diff_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        fast_expansion_combine_zeroelim<true>(e, f, h);
    }

    /************************************************************************/

    // Reached only by the intermediates of assign_product() that are too big for the stack.
    // geogram pooled these behind a global lock, for a high-level expansion_nt that allocated
    // every value on the heap and that nothing here uses.
    expansion* expansion::new_expansion_on_heap(index_t capa) {
        return new(new unsigned char[expansion::bytes(capa)])expansion(capa);
    }

    void expansion::delete_expansion_on_heap(expansion* e) {
        delete[] reinterpret_cast<unsigned char*>(e);
    }

    // ====== Initialization from expansion and double ===============

    expansion& expansion::assign_product(const expansion& a, double b) {
        scale_expansion_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion sum and difference =========================

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b
    ) {
        fast_expansion_sum_zeroelim(a, b, *this);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        expansion& ab = expansion_sum(a, b);
        this->assign_sum(ab, c);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    ) {
        expansion& ab = expansion_sum(a, b);
        expansion& cd = expansion_sum(c, d);
        this->assign_sum(ab, cd);
        return *this;
    }

    expansion& expansion::assign_diff(const expansion& a, const expansion& b) {
        fast_expansion_diff_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion product ==================================

    // Recursive helper function for product implementation
    expansion& expansion::assign_sub_product(
        const double* a, index_t a_length, const expansion& b
    ) {
        if(a_length == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.

            const double* a1 = a;
            index_t a1_length = a_length / 2;
            const double* a2 = a1 + a1_length;
            index_t a2_length = a_length - a1_length;

            // Allocate both halves on the stack or on the heap if too large
            // (some platformes, e.g. MacOSX, have a small stack)

            index_t a1b_capa = sub_product_capacity(a1_length, b.length());
            index_t a2b_capa = sub_product_capacity(a2_length, b.length());

            bool a1b_on_heap = (a1b_capa > MAX_CAPACITY_ON_STACK);
            bool a2b_on_heap = (a2b_capa > MAX_CAPACITY_ON_STACK);

            expansion* a1b = a1b_on_heap ?
                new_expansion_on_heap(a1b_capa) :
                new_expansion_on_stack(a1b_capa);

            a1b->assign_sub_product(a1, a1_length, b);

            expansion* a2b = a2b_on_heap ?
                new_expansion_on_heap(a2b_capa) :
                new_expansion_on_stack(a2b_capa);

            a2b->assign_sub_product(a2, a2_length, b);

            this->assign_sum(*a1b, *a2b);

            if(a1b_on_heap) {
                delete_expansion_on_heap(a1b);
            }

            if(a2b_on_heap) {
                delete_expansion_on_heap(a2b);
            }
        }
        return *this;
    }

    expansion& expansion::assign_product(
        const expansion& a, const expansion& b
    ) {
        if(a.length() == 0 || b.length() == 0) {
            x_[0] = 0.0;
            set_length(0);
        } else if(a.length() == 1 && b.length() == 1) {
            two_product(a[0], b[0], x_[1], x_[0]);
            set_length(2);
        } else if(a.length() == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else if(b.length() == 1) {
            scale_expansion_zeroelim(a, b[0], *this);
        } else if(a.length() == 2 && b.length() == 2) {
            two_two_product(a.data(), b.data(), x_);
            set_length(8);
        } else {

            const expansion* pa = &a;
            const expansion* pb = &b;

            if(pa->length() > pb->length()) {
                std::swap(pa, pb);
            }

            // [Shewchuk 97]
            // (https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf)
            // Section 2.8: other operations
            // Distillation: sum of k values.
            //    Worst case: 1/2*k*(k-1)
            //    But O(k log(k)) if the "summing tree" is well balanced
            //      and using fast_expansion_sum().
            // Recommended way of computing a product:
            //    compute a1*b, a2*b ... ak*b using scale_expansion_zeroelim()
            //    sum them using a well-balanced tree
            // However, there is an extra cost for the recursion (and more
            // importantly, for allocating the intermediary sums, especially
            // when they do not fit on the stack). So when there are less than
            // 16 values to add, we simply accumulate them.

            bool use_balanced_distillation = (pa->length() >= 16);

            if(use_balanced_distillation) {
                // assign_sub_product() is a recursive function that
                // creates a balanced distillation tree on the stack.
                assign_sub_product(pa->data(), pa->length(),*pb);
            } else {
                // trivial implementation: compute all the products
                // P = ak*b and accumulate them into S

                index_t P_capa = product_capacity(*pb, 3.0); // 3.0, or any
                                                             // number that is
                                                             // not a power of 2

                index_t S_capa = capacity(); // same capacity as this,
                                             // enough to store sum.

                bool P_on_heap = (P_capa > MAX_CAPACITY_ON_STACK);
                bool S_on_heap = (S_capa > MAX_CAPACITY_ON_STACK);

                expansion* P = P_on_heap ?
                    new_expansion_on_heap(P_capa) :
                    new_expansion_on_stack(P_capa);

                expansion* S = S_on_heap ?
                    new_expansion_on_heap(S_capa) :
                    new_expansion_on_stack(S_capa);

                expansion* S1 = S;
                expansion* S2 = this;

                if((pa->length()%2) == 0) {
                    std::swap(S1,S2);
                }

                for(index_t i=0; i<pa->length(); ++i) {
                    if(i == 0) {
                        S2->assign_product(*pb, (*pa)[i]);
                    } else {
                        P->assign_product(*pb, (*pa)[i]);
                        S2->assign_sum(*S1,*P);
                    }
                    std::swap(S1,S2);
                }

                geo_assert(S1 == this);

                if(S_on_heap) {
                    delete_expansion_on_heap(S);
                }

                if(P_on_heap) {
                    delete_expansion_on_heap(P);
                }
            }
        }
        return *this;
    }

    // =============  determinants ==========================================

    expansion& expansion::assign_det2x2(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    ) {
        const expansion& a11a22 = expansion_product(a11, a22);
        const expansion& a12a21 = expansion_product(a12, a21);
        return this->assign_diff(a11a22, a12a21);
    }

    expansion& expansion::assign_det3x3(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        // Development w.r.t. first row
        const expansion& c11 = expansion_det2x2(a22, a23, a32, a33);
        const expansion& c12 = expansion_det2x2(a23, a21, a33, a31);
        const expansion& c13 = expansion_det2x2(a21, a22, a31, a32);
        const expansion& a11c11 = expansion_product(a11, c11);
        const expansion& a12c12 = expansion_product(a12, c12);
        const expansion& a13c13 = expansion_product(a13, c13);
        return this->assign_sum(a11c11, a12c12, a13c13);
    }

    // =============  geometric operations ==================================

    expansion& expansion::assign_sq_dist(const double* p1, const double* p2) {
        // The distillation tree geogram builds by halving the dimension, written out: the y and z
        // terms are summed first, then the x term is added to that.
        expansion& dx = *new_expansion_on_stack(6);
        expansion& dy = *new_expansion_on_stack(6);
        expansion& dz = *new_expansion_on_stack(6);
        assign_sq_diff(dx, p1[0], p2[0]);
        assign_sq_diff(dy, p1[1], p2[1]);
        assign_sq_diff(dz, p1[2], p2[2]);
        expansion& dyz = expansion_sum(dy, dz);
        return this->assign_sum(dx, dyz);
    }

    /************************************************************************/

} }
