// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/numerics/multi_precision.h
// Copied rather than reimplemented: Shewchuk expansion arithmetic, arithmetic unchanged.
//
// Multi-precision arithmetic on expansions, as described by Jonathan Shewchuk in "Adaptive
// Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates", Discrete &
// Computational Geometry 18(3):305-363, October 1997.
//
// Trimmed to the entry points the mesher reaches: gone are compare(), is_same_as(), optimize(),
// compress_expansion(), sign_of_expansion_determinant(), grow_expansion_zeroelim() and the sum and
// difference of an expansion with a double. Squared distance is 3d rather than nd, since that is
// the only shape asked for. The expansion constructor is explicit here, which upstream's is not --
// see the note on it.

#pragma once

#include "geo_basic.h"
#include <new>
#include <math.h>

namespace floatTetWild {
namespace geo {

    // Shewchuk's primitive: turns an exact difference of two doubles into a length 2 expansion,
    // x the rounded result and y the part that did not fit. Its siblings two_sum, two_product and
    // square live in geo_multi_precision.cpp, the only place that uses them.
    inline void two_diff(double a, double b, double& x, double& y) {
        x = a - b;
        double bvirt = a - x;
        double avirt = x + bvirt;
        double bround = bvirt - b;
        double around = a - avirt;
        y = around + bround;
    }

    /************************************************************************/

    // A number in arbitrary precision: a length and that many doubles that sum to it exactly, in
    // increasing magnitude. Allocated with the macros below rather than constructed, since the
    // storage runs past the end of the object.
    class expansion {
    public:
    index_t length() const {
        return length_;
    }

    index_t capacity() const {
        return capacity_;
    }

    void set_length(index_t new_length) {
        length_ = new_length;
    }

    // Component i. One past capacity is valid to read: the storage carries a sentry because
    // fast_expansion_sum_zeroelim() may look one past the end without using what it finds.
    const double& operator[] (index_t i) const {
        return x_[i];
    }

    double& operator[] (index_t i) {
        return x_[i];
    }

    const double* data() const {
        return x_;
    }

    // Bytes needed to store an expansion of the given capacity: the object itself, less the two
    // doubles x_ is declared with, plus the components and the sentry.
    static size_t bytes(index_t capa) {
        return
            sizeof(expansion) - 2 * sizeof(double) +
            (capa + 1) * sizeof(double);
    }

    // explicit, unlike upstream. Without it a double converts to an expansion of that many
    // components, so dropping any assign_xxx(const expansion&, double) overload silently rebinds
    // its callers to the (expansion, expansion) one and the exact arithmetic quietly goes wrong
    // instead of failing to build.
    explicit expansion(index_t capa) :
    length_(0),
    capacity_(capa) {
    }

    // Allocates an expansion in the calling function's frame. A macro rather than an inline
    // function because alloca() cannot be called from one.
#define new_expansion_on_stack(capa)                                    \
    (new (alloca(expansion::bytes(capa)))expansion(capa))

    static expansion* new_expansion_on_heap(index_t capa);

    static void delete_expansion_on_heap(expansion* e);

    // Each operation comes as a capacity, which the macros below allocate, and an assignment,
    // which fills what they allocated. The capacities never depend on the values, only on the
    // lengths, so they can be computed before the operands are.

    static index_t diff_capacity(double, double) {
        return 2;
    }

    expansion& assign_diff(double a, double b) {
        set_length(2);
        two_diff(a, b, x_[1], x_[0]);
        return *this;
    }

    static index_t product_capacity(const expansion& a, double) {
        return a.length() * 2;
    }

    expansion& assign_product(const expansion& a, double b);

    static index_t sum_capacity(const expansion& a, const expansion& b) {
        return a.length() + b.length();
    }

    expansion& assign_sum(const expansion& a, const expansion& b);

    static index_t sum_capacity(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        return a.length() + b.length() + c.length();
    }

    expansion& assign_sum(
        const expansion& a, const expansion& b, const expansion& c
    );

    static index_t sum_capacity(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    ) {
        return a.length() + b.length() + c.length() + d.length();
    }

    expansion& assign_sum(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    );

    static index_t diff_capacity(const expansion& a, const expansion& b) {
        return a.length() + b.length();
    }

    expansion& assign_diff(const expansion& a, const expansion& b);

    static index_t product_capacity(
        const expansion& a, const expansion& b
    ) {
        return a.length() * b.length() * 2;
    }

    expansion& assign_product(const expansion& a, const expansion& b);

    // a11 * a22 - a21 * a12.
    static index_t det2x2_capacity(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    ) {
        return
            product_capacity(a11, a22) +
            product_capacity(a21, a12);
    }

    expansion& assign_det2x2(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    );

    static index_t det3x3_capacity(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        // Development w.r.t. first row
        index_t c11_capa = det2x2_capacity(a22, a23, a32, a33);
        index_t c12_capa = det2x2_capacity(a21, a23, a31, a33);
        index_t c13_capa = det2x2_capacity(a21, a22, a31, a32);
        return 2 * (
            a11.length() * c11_capa +
            a12.length() * c12_capa +
            a13.length() * c13_capa
        );
    }

    expansion& assign_det3x3(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    );

    // The squared distance between two 3d points. Six components per coordinate.
    static index_t sq_dist_capacity() {
        return 18;
    }

    expansion& assign_sq_dist(const double* p1, const double* p2);

    expansion& negate() {
        for(index_t i = 0; i < length_; ++i) {
            x_[i] = -x_[i];
        }
        return *this;
    }

    // An approximation of the stored value.
    double estimate() const {
        double result = 0.0;
        for(index_t i = 0; i < length(); ++i) {
            result += x_[i];
        }
        return result;
    }

    // The sign, exactly: the components increase in magnitude, so the last one decides.
    Sign sign() const {
        if(length() == 0) {
            return ZERO;
        }
        return geo_sgn(x_[length() - 1]);
    }

    private:
    static index_t sub_product_capacity(
        index_t a_length, index_t b_length
    ) {
        return a_length * b_length * 2;
    }

    // [a1...a_length] * b. Used by assign_product() in balanced distillation mode, where it
    // recursively assembles sub-sums.
    expansion& assign_sub_product(
        const double* a, index_t a_length, const expansion& b
    );

    expansion(const expansion& rhs) = delete;

    expansion& operator= (const expansion& rhs) = delete;

    private:

    // Above this capacity an intermediate goes on the heap rather than the stack. geogram drops
    // to 256 on Apple, where a default pthread stack is 512 KB. The thread pool here gives every
    // worker 64 MB, so fTetWild keeps 1024 on every platform.
    static constexpr index_t MAX_CAPACITY_ON_STACK = 1024;
    index_t length_;
    index_t capacity_;
    double x_[2];  // x_ is in fact of size [capacity_]
    };

    // Each of these allocates a result in the calling function's frame and fills it, so the
    // reference must not outlive the call. \p a and \p b are doubles or expansions where the
    // capacity accepts either.

#define expansion_sum(a, b)                     \
    new_expansion_on_stack(                     \
        expansion::sum_capacity(a, b)           \
    )->assign_sum(a, b)

#define expansion_sum3(a, b, c)                 \
    new_expansion_on_stack(                     \
        expansion::sum_capacity(a, b, c)        \
    )->assign_sum(a, b, c)

#define expansion_sum4(a, b, c, d)              \
    new_expansion_on_stack(                     \
        expansion::sum_capacity(a, b, c, d)     \
    )->assign_sum(a, b, c, d)

#define expansion_diff(a, b)                    \
    new_expansion_on_stack(                     \
        expansion::diff_capacity(a, b)          \
    )->assign_diff(a, b)

#define expansion_product(a, b)                 \
    new_expansion_on_stack(                     \
        expansion::product_capacity(a, b)       \
    )->assign_product(a, b)

#define expansion_det2x2(a11, a12, a21, a22)            \
    new_expansion_on_stack(                             \
        expansion::det2x2_capacity(a11, a12, a21, a22)  \
    )->assign_det2x2(a11, a12, a21, a22)

#define expansion_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33)   \
    new_expansion_on_stack(                                             \
        expansion::det3x3_capacity(a11,a12,a13,a21,a22,a23,a31,a32,a33) \
    )->assign_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33)

    // \p a and \p b are const double*, three coordinates each.
#define expansion_sq_dist(a, b)                 \
    new_expansion_on_stack(                     \
        expansion::sq_dist_capacity()           \
    )->assign_sq_dist(a, b)

    /************************************************************************/
} }
