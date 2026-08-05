// The fast winding number for triangle soups, from "Fast Winding Numbers for Soups and Clouds"
// [Barill et al. 2018], https://github.com/alecjacobson/WindingNumber, by way of libigl.
//
// Upstream is an amalgamation of the Houdini Development Kit: a general array container, a
// general fixed-size vector, two interchangeable SIMD backends behind sixty macros, a BVH
// templated on its width, its index type, its box type and one of seven split heuristics, and
// the 2D subtended angle. The mesher builds one four-wide BVH of float boxes, splits it on box
// area, and asks it for one thing, so that is all that is left here.
//
// Copyright (c) 2018 Side Effects Software Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
// BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "FastWindingNumber.h"

#include "ParallelFor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <thread>
#include <utility>

namespace floatTetWild {
namespace FastWindingNumber {
namespace {

using std::int32_t;
using std::uint32_t;
using std::uint64_t;

constexpr float Inf = std::numeric_limits<float>::infinity();

// Four floats at a time, and a four-lane mask over them. The BVH is four wide, so a node's four
// children are approximated together. Upstream reached these through an SSE backend or a scalar
// one; every operation the approximation uses is exact IEEE, so the two agreed lane for lane.
//
// Plain aggregates behind free functions that take and return whole vectors, because that is the
// shape clang turns into four-lane NEON. The same arithmetic as member operators over references
// comes out two lanes at a time, and so does fusing any one step into a per-lane expression: a
// `? :` in place of the bitwise and below, or std::isfinite in place of the exponent test, each
// cost a third of the query on their own.
struct Float4
{
    float v[4];
};

struct Mask4
{
    // All bits set or none, per lane.
    int32_t v[4];
};

Mask4 as_mask(Float4 a)
{
    Mask4 r;
    std::memcpy(&r, &a, sizeof(r));
    return r;
}

Float4 as_float(Mask4 a)
{
    Float4 r;
    std::memcpy(&r, &a, sizeof(r));
    return r;
}

int32_t lane_mask(bool c) { return c ? int32_t(-1) : 0; }

Float4 splat(float a) { return {{a, a, a, a}}; }

Float4 operator+(Float4 a, Float4 b)
{ return {{a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2], a.v[3] + b.v[3]}}; }
Float4 operator-(Float4 a, Float4 b)
{ return {{a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2], a.v[3] - b.v[3]}}; }
Float4 operator*(Float4 a, Float4 b)
{ return {{a.v[0] * b.v[0], a.v[1] * b.v[1], a.v[2] * b.v[2], a.v[3] * b.v[3]}}; }
Float4 operator/(Float4 a, Float4 b)
{ return {{a.v[0] / b.v[0], a.v[1] / b.v[1], a.v[2] / b.v[2], a.v[3] / b.v[3]}}; }
Float4 operator-(Float4 a) { return {{-a.v[0], -a.v[1], -a.v[2], -a.v[3]}}; }
Float4 operator*(Float4 a, float b) { return a * splat(b); }

Float4& operator+=(Float4& a, Float4 b) { return a = a + b; }
Float4& operator-=(Float4& a, Float4 b) { return a = a - b; }
Float4& operator*=(Float4& a, Float4 b) { return a = a * b; }

Float4 sqrt(Float4 a)
{ return {{std::sqrt(a.v[0]), std::sqrt(a.v[1]), std::sqrt(a.v[2]), std::sqrt(a.v[3])}}; }

Mask4 operator&(Mask4 a, Mask4 b)
{ return {{a.v[0] & b.v[0], a.v[1] & b.v[1], a.v[2] & b.v[2], a.v[3] & b.v[3]}}; }
Mask4 operator~(Mask4 a) { return {{~a.v[0], ~a.v[1], ~a.v[2], ~a.v[3]}}; }

Mask4 operator<=(Float4 a, Float4 b)
{
    return {{lane_mask(a.v[0] <= b.v[0]), lane_mask(a.v[1] <= b.v[1]),
             lane_mask(a.v[2] <= b.v[2]), lane_mask(a.v[3] <= b.v[3])}};
}

// Masked-off lanes go to zero, an all-zero float being +0.
Float4 operator&(Float4 a, Mask4 m) { return as_float(as_mask(a) & m); }

// A maximum exponent means either infinite or NaN.
Mask4 is_finite(Float4 a)
{
    const Mask4 all = {{0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000}};
    const Mask4 exp = as_mask(a) & all;
    return ~Mask4{{lane_mask(exp.v[0] == all.v[0]), lane_mask(exp.v[1] == all.v[1]),
                   lane_mask(exp.v[2] == all.v[2]), lane_mask(exp.v[3] == all.v[3])}};
}

// One bit per lane, low lane first, the way the SSE movemask this replaces was read.
uint32_t mask_bits(Mask4 a)
{
    return uint32_t(a.v[0] < 0) | (uint32_t(a.v[1] < 0) << 1) | (uint32_t(a.v[2] < 0) << 2) |
           (uint32_t(a.v[3] < 0) << 3);
}

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

using Vec3f  = Vec3<float>;
using Vec3x4 = Vec3<Float4>;

Vec3f cross(const Vec3f& v1, const Vec3f& v2)
{
    Vec3f result;
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}

Vec3f max_components(const Vec3f& a, const Vec3f& b)
{
    Vec3f result;
    for (int i = 0; i < 3; ++i)
        result[i] = (a[i] < b[i]) ? b[i] : a[i];
    return result;
}

// ---------------------------------------------------------------------------------------------
// Axis-aligned boxes

struct Box
{
    // Min and max per axis. Left uninitialized by the vector constructors that fill it in
    // afterwards, so it stays an aggregate.
    float vals[3][2];

    void reset()
    {
        for (int axis = 0; axis < 3; ++axis) {
            vals[axis][0] = std::numeric_limits<float>::max();
            vals[axis][1] = -std::numeric_limits<float>::max();
        }
    }
    void set(const Vec3f& p)
    {
        for (int axis = 0; axis < 3; ++axis) {
            vals[axis][0] = p[axis];
            vals[axis][1] = p[axis];
        }
    }
    void enlarge(const Box& src)
    {
        for (int axis = 0; axis < 3; ++axis) {
            float&      minv    = vals[axis][0];
            float&      maxv    = vals[axis][1];
            const float curminv = src.vals[axis][0];
            const float curmaxv = src.vals[axis][1];
            minv                = (minv < curminv) ? minv : curminv;
            maxv                = (maxv > curmaxv) ? maxv : curmaxv;
        }
    }
    void enlarge(const Vec3f& p)
    {
        for (int axis = 0; axis < 3; ++axis) {
            float&      minv = vals[axis][0];
            float&      maxv = vals[axis][1];
            const float cur  = p[axis];
            minv             = (minv < cur) ? minv : cur;
            maxv             = (maxv > cur) ? maxv : cur;
        }
    }

    Vec3f min() const
    {
        Vec3f v;
        for (int axis = 0; axis < 3; ++axis)
            v[axis] = vals[axis][0];
        return v;
    }
    Vec3f max() const
    {
        Vec3f v;
        for (int axis = 0; axis < 3; ++axis)
            v[axis] = vals[axis][1];
        return v;
    }

    // Half the surface area, the only split heuristic the mesher asks for.
    float half_surface_area() const
    {
        const float d0 = (vals[0][1] - vals[0][0]);
        const float d1 = (vals[1][1] - vals[1][0]);
        const float d2 = (vals[2][1] - vals[2][0]);
        return d0 * d1 + d1 * d2 + d2 * d0;
    }

    bool has_nan_or_inf() const
    {
        bool bad = false;
        for (int axis = 0; axis < 3; ++axis) {
            bad |= !std::isfinite(vals[axis][0]);
            bad |= !std::isfinite(vals[axis][1]);
        }
        return bad;
    }
};

// Twice the box centre along one axis. The splitting code works in these doubled coordinates
// throughout, so the factor of two never has to be undone.
float box_centre_x2(const Box& box, uint32_t axis)
{
    const float* v = box.vals[axis];
    return v[0] + v[1];
}

// ---------------------------------------------------------------------------------------------
// A four-wide BVH over the triangle boxes

constexpr uint32_t BvhN = 4;

struct Node
{
    uint32_t child[BvhN];

    // A child is either a triangle index or, with the top bit set, an index into the node array.
    // Empty is every bit set, so it reads as internal and is never a valid node index.
    static constexpr uint32_t Empty       = uint32_t(-1);
    static constexpr uint32_t InternalBit = uint32_t(1) << (sizeof(uint32_t) * 8 - 1);

    static uint32_t mark_internal(uint32_t node_num) { return node_num | InternalBit; }
    static bool     is_internal(uint32_t c) { return (c & InternalBit) != 0; }
    static uint32_t internal_num(uint32_t c) { return c & ~InternalBit; }
};

// 16 equal-length spans (15 evenly-spaced splits) is enough for a decent heuristic, and at least
// 1/16 of the boxes must land on each side, else the tree could get very deep.
constexpr uint32_t NSpans      = 16;
constexpr uint32_t NSplits     = NSpans - 1;
constexpr uint32_t MinFraction = 16;

int num_processors()
{ return std::thread::hardware_concurrency(); }

// An overestimate of the number of nodes needed. At worst we could have only 2 children in every
// leaf, and above that a geometric series with r=1/N and a=(nboxes/2)/N.
uint32_t node_estimate(uint32_t nboxes)
{ return nboxes / 2 + nboxes / (2 * (BvhN - 1)); }

void compute_full_bounding_box(Box&            axes_minmax,
                               const Box*      boxes,
                               const uint32_t  nboxes,
                               const uint32_t* indices)
{
    if (!nboxes) {
        axes_minmax.reset();
        return;
    }
    uint32_t ntasks = 1;
    if (nboxes >= 2 * 4096) {
        const uint32_t nprocessors = num_processors();
        ntasks = (nprocessors > 1) ? std::min(4 * nprocessors, nboxes / 4096) : 1;
    }
    if (ntasks == 1) {
        Box box = boxes[indices[0]];
        for (uint32_t i = 1; i < nboxes; ++i)
            box.enlarge(boxes[indices[i]]);
        axes_minmax = box;
    }
    else {
        // NOTE: the per-chunk boxes start zeroed rather than empty, so the union always contains
        // the origin. That is what upstream did and the mesher's output depends on it.
        const size_t     nslots = detail::chunk_count(0, nboxes);
        std::vector<Box> parallel_boxes(nslots);
        detail::parallel_ranges(0, nboxes, [&](size_t lo, size_t hi, size_t t) {
            for (size_t i = lo; i < hi; ++i)
                parallel_boxes[t].enlarge(boxes[indices[i]]);
        });
        Box box = parallel_boxes[0];
        for (size_t t = 1; t < nslots; ++t)
            box.enlarge(parallel_boxes[t]);
        axes_minmax = box;
    }
}

// Drops the indices of every box with a NaN or infinite bound by shifting the rest down over
// them. Returns whether any were dropped and updates nboxes.
bool exclude_nan_inf_boxes(const Box* boxes, uint32_t* indices, uint32_t& nboxes)
{
    const uint32_t* indices_end = indices + nboxes;

    uint32_t* psrc_index = indices;
    for (; psrc_index != indices_end; ++psrc_index) {
        if (boxes[*psrc_index].has_nan_or_inf())
            break;
    }
    if (psrc_index == indices_end)
        return false;

    // First NaN or infinite box
    uint32_t* nan_start = psrc_index;
    for (++psrc_index; psrc_index != indices_end; ++psrc_index) {
        if (!boxes[*psrc_index].has_nan_or_inf()) {
            *nan_start = *psrc_index;
            ++nan_start;
        }
    }
    nboxes = nan_start - indices;
    return true;
}

// Reorders [indices, indices_end) so that everything with a centre below pivotx2 comes first,
// everything equal to it next, and everything above it last, and reports the two boundaries.
void partition_by_centre(const Box*            boxes,
                         uint32_t* const       indices,
                         const uint32_t* const indices_end,
                         const uint32_t        axis,
                         const float           pivotx2,
                         uint32_t*&            ppivot_start,
                         uint32_t*&            ppivot_end)
{
    // First element >= pivot
    uint32_t* pivot_start = indices;
    // First element > pivot
    uint32_t* pivot_end = indices;

    for (uint32_t* psrc_index = indices; psrc_index != indices_end; ++psrc_index) {
        const float srcsum = box_centre_x2(boxes[*psrc_index], axis);
        if (srcsum < pivotx2) {
            if (psrc_index != pivot_start) {
                if (pivot_start == pivot_end) {
                    // Common case: nothing equal to the pivot
                    std::swap(*psrc_index, *pivot_start);
                }
                else {
                    // Less common case: at least one thing equal to the pivot
                    const uint32_t temp = *psrc_index;
                    *psrc_index         = *pivot_end;
                    *pivot_end          = *pivot_start;
                    *pivot_start        = temp;
                }
            }
            ++pivot_start;
            ++pivot_end;
        }
        else if (srcsum == pivotx2) {
            // Add to the pivot area
            if (psrc_index != pivot_end)
                std::swap(*psrc_index, *pivot_end);
            ++pivot_end;
        }
    }
    ppivot_start = pivot_start;
    ppivot_end   = pivot_end;
}

// Partially sorts by box centre until the entry at nth is the one that belongs there.
void select_nth_by_centre(const Box*      boxes,
                          uint32_t*       indices,
                          const uint32_t* indices_end,
                          const uint32_t  axis,
                          uint32_t* const nth)
{
    while (true) {
        // Choose median of first, middle, and last as the pivot
        float pivots[3] = {box_centre_x2(boxes[indices[0]], axis),
                           box_centre_x2(boxes[indices[(indices_end - indices) / 2]], axis),
                           box_centre_x2(boxes[*(indices_end - 1)], axis)};
        if (pivots[0] < pivots[1])
            std::swap(pivots[0], pivots[1]);
        if (pivots[0] < pivots[2])
            std::swap(pivots[0], pivots[2]);
        if (pivots[1] < pivots[2])
            std::swap(pivots[1], pivots[2]);

        uint32_t* pivot_start;
        uint32_t* pivot_end;
        partition_by_centre(boxes, indices, indices_end, axis, pivots[1], pivot_start, pivot_end);
        if (nth < pivot_start) {
            indices_end = pivot_start;
        }
        else if (nth < pivot_end) {
            // nth is in the middle of the pivot range, which is in the right place.
            return;
        }
        else {
            indices = pivot_end;
        }
        if (indices_end <= indices + 1)
            return;
    }
}

// Splits [indices, indices+nboxes) in two, minimizing the area heuristic. Reports the boundary
// and the two boxes.
void split(const Box& axes_minmax,
           const Box* boxes,
           uint32_t*  indices,
           uint32_t   nboxes,
           uint32_t*& split_indices,
           Box*       split_boxes)
{
    if (nboxes == 2) {
        split_boxes[0] = boxes[indices[0]];
        split_boxes[1] = boxes[indices[1]];
        split_indices  = indices + 1;
        return;
    }

    constexpr uint32_t SmallLimit = 6;
    if (nboxes <= SmallLimit) {
        // Special case for a small number of boxes: check all (2^(n-1))-1 partitions. Without
        // loss of generality, box 0 is in partition 0, and not all boxes are in partition 0.
        Box local_boxes[SmallLimit];
        for (uint32_t box = 0; box < nboxes; ++box)
            local_boxes[box] = boxes[indices[box]];

        const uint32_t partition_limit = (uint32_t(1) << (nboxes - 1));
        uint32_t       best_partition  = uint32_t(-1);
        float          best_heuristic;
        for (uint32_t partition_bits = 1; partition_bits < partition_limit; ++partition_bits) {
            Box sub_boxes[2];
            sub_boxes[0] = local_boxes[0];
            sub_boxes[1].reset();
            uint32_t sub_counts[2] = {1, 0};
            for (uint32_t bit = 0; bit < nboxes - 1; ++bit) {
                uint32_t dest = (partition_bits >> bit) & 1;
                sub_boxes[dest].enlarge(local_boxes[bit + 1]);
                ++sub_counts[dest];
            }
            const float heuristic = sub_boxes[0].half_surface_area() * sub_counts[0] +
                                    sub_boxes[1].half_surface_area() * sub_counts[1];
            if (best_partition == uint32_t(-1) || heuristic < best_heuristic) {
                best_partition = partition_bits;
                best_heuristic = heuristic;
                split_boxes[0] = sub_boxes[0];
                split_boxes[1] = sub_boxes[1];
            }
        }

        // Reorder the indices. Index 0 is always in partition 0, so it can stay put.
        uint32_t local_indices[SmallLimit - 1];
        for (uint32_t box = 0; box < nboxes - 1; ++box)
            local_indices[box] = indices[box + 1];
        uint32_t* dest_indices = indices + 1;
        uint32_t* src_indices  = local_indices;
        // Copy partition 0
        for (uint32_t bit = 0; bit < nboxes - 1; ++bit, ++src_indices) {
            if (!((best_partition >> bit) & 1)) {
                *dest_indices = *src_indices;
                ++dest_indices;
            }
        }
        split_indices = dest_indices;
        // Copy partition 1
        src_indices = local_indices;
        for (uint32_t bit = 0; bit < nboxes - 1; ++bit, ++src_indices) {
            if ((best_partition >> bit) & 1) {
                *dest_indices = *src_indices;
                ++dest_indices;
            }
        }
        return;
    }

    uint32_t max_axis        = 0;
    float    max_axis_length = axes_minmax.vals[0][1] - axes_minmax.vals[0][0];
    for (uint32_t axis = 1; axis < 3; ++axis) {
        const float axis_length = axes_minmax.vals[axis][1] - axes_minmax.vals[axis][0];
        if (axis_length > max_axis_length) {
            max_axis        = axis;
            max_axis_length = axis_length;
        }
    }

    if (!(max_axis_length > 0.0f)) {
        // All boxes are a single point or NaN, so pick an arbitrary split point.
        split_indices  = indices + nboxes / 2;
        split_boxes[0] = axes_minmax;
        split_boxes[1] = axes_minmax;
        return;
    }

    const uint32_t axis = max_axis;

    constexpr uint32_t MidLimit = 2 * NSpans;
    if (nboxes <= MidLimit) {
        // Sort along the axis, and try all possible splits.
        float midpointsx2[MidLimit];
        for (uint32_t i = 0; i < nboxes; ++i)
            midpointsx2[i] = box_centre_x2(boxes[indices[i]], axis);
        uint32_t local_indices[MidLimit];
        for (uint32_t i = 0; i < nboxes; ++i)
            local_indices[i] = i;

        const uint32_t chunk_starts[5] = {
          0, nboxes / 4, nboxes / 2, uint32_t((3 * uint64_t(nboxes)) / 4), nboxes};

        // For sorting, insertion sort 4 chunks and merge them
        for (uint32_t chunk = 0; chunk < 4; ++chunk) {
            const uint32_t start = chunk_starts[chunk];
            const uint32_t end   = chunk_starts[chunk + 1];
            for (uint32_t i = start + 1; i < end; ++i) {
                uint32_t indexi = local_indices[i];
                float    vi     = midpointsx2[indexi];
                for (uint32_t j = start; j < i; ++j) {
                    uint32_t indexj = local_indices[j];
                    float    vj     = midpointsx2[indexj];
                    if (vi < vj) {
                        do {
                            local_indices[j] = indexi;
                            indexi           = indexj;
                            ++j;
                            if (j == i) {
                                local_indices[j] = indexi;
                                break;
                            }
                            indexj = local_indices[j];
                        } while (true);
                        break;
                    }
                }
            }
        }
        // Merge chunks into another buffer
        const auto by_midpoint = [&midpointsx2](const uint32_t a, const uint32_t b) -> bool {
            return midpointsx2[a] < midpointsx2[b];
        };
        uint32_t local_indices_temp[MidLimit];
        std::merge(local_indices,
                   local_indices + chunk_starts[1],
                   local_indices + chunk_starts[1],
                   local_indices + chunk_starts[2],
                   local_indices_temp,
                   by_midpoint);
        std::merge(local_indices + chunk_starts[2],
                   local_indices + chunk_starts[3],
                   local_indices + chunk_starts[3],
                   local_indices + chunk_starts[4],
                   local_indices_temp + chunk_starts[2],
                   by_midpoint);
        std::merge(local_indices_temp,
                   local_indices_temp + chunk_starts[2],
                   local_indices_temp + chunk_starts[2],
                   local_indices_temp + chunk_starts[4],
                   local_indices,
                   by_midpoint);

        // Translate local_indices into indices, then copy back
        for (uint32_t i = 0; i < nboxes; ++i)
            local_indices[i] = indices[local_indices[i]];
        for (uint32_t i = 0; i < nboxes; ++i)
            indices[i] = local_indices[i];

        // Accumulate boxes
        Box            left_boxes[MidLimit - 1];
        Box            right_boxes[MidLimit - 1];
        const uint32_t nsplits         = nboxes - 1;
        Box            box_accumulator = boxes[local_indices[0]];
        left_boxes[0]                  = box_accumulator;
        for (uint32_t i = 1; i < nsplits; ++i) {
            box_accumulator.enlarge(boxes[local_indices[i]]);
            left_boxes[i] = box_accumulator;
        }
        box_accumulator          = boxes[local_indices[nsplits - 1]];
        right_boxes[nsplits - 1] = box_accumulator;
        for (uint32_t i = nsplits - 1; i > 0; --i) {
            box_accumulator.enlarge(boxes[local_indices[i]]);
            right_boxes[i - 1] = box_accumulator;
        }

        uint32_t best_split = 0;
        float    best_local_heuristic =
          left_boxes[0].half_surface_area() + right_boxes[0].half_surface_area() * (nboxes - 1);
        for (uint32_t s = 1; s < nsplits; ++s) {
            const float heuristic = left_boxes[s].half_surface_area() * (s + 1) +
                                    right_boxes[s].half_surface_area() * (nboxes - (s + 1));
            if (heuristic < best_local_heuristic) {
                best_split           = s;
                best_local_heuristic = heuristic;
            }
        }
        split_indices  = indices + best_split + 1;
        split_boxes[0] = left_boxes[best_split];
        split_boxes[1] = right_boxes[best_split];
        return;
    }

    const float axis_min    = axes_minmax.vals[max_axis][0];
    const float axis_length = max_axis_length;
    Box         span_boxes[NSpans];
    for (uint32_t i = 0; i < NSpans; ++i)
        span_boxes[i].reset();
    uint32_t span_counts[NSpans];
    for (uint32_t i = 0; i < NSpans; ++i)
        span_counts[i] = 0;

    const float axis_min_x2 = 2 * axis_min;
    // NOTE: the factor of 0.5 undoes the doubling in the box centres.
    const float        axis_index_scale          = (0.5f * NSpans) / axis_length;
    constexpr uint32_t BoxSpansParallelThreshold = 2048;
    uint32_t           ntasks                    = 1;
    if (nboxes >= BoxSpansParallelThreshold) {
        const uint32_t nprocessors = num_processors();
        ntasks = (nprocessors > 1)
                   ? std::min(4 * nprocessors, nboxes / (BoxSpansParallelThreshold / 2))
                   : 1;
    }
    const auto span_of = [axis, axis_min_x2, axis_index_scale](const Box& box) {
        const int span = int((box_centre_x2(box, axis) - axis_min_x2) * axis_index_scale);
        return uint32_t(span <= 0 ? 0 : (span >= int(NSpans - 1) ? int(NSpans - 1) : span));
    };
    if (ntasks == 1) {
        for (uint32_t indexi = 0; indexi < nboxes; ++indexi) {
            const Box&     box        = boxes[indices[indexi]];
            const uint32_t span_index = span_of(box);
            ++span_counts[span_index];
            span_boxes[span_index].enlarge(box);
        }
    }
    else {
        const size_t          nslots = detail::chunk_count(0, nboxes);
        std::vector<Box>      parallel_boxes(NSpans * nslots);
        std::vector<uint32_t> parallel_counts(NSpans * nslots, 0);
        for (size_t t = 0; t < nslots; t++) {
            for (uint32_t i = 0; i < NSpans; ++i)
                parallel_boxes[t * NSpans + i].reset();
        }
        detail::parallel_ranges(0, nboxes, [&](size_t lo, size_t hi, size_t t) {
            for (size_t j = lo; j < hi; ++j) {
                const Box&     box        = boxes[indices[j]];
                const uint32_t span_index = span_of(box);
                ++parallel_counts[t * NSpans + span_index];
                parallel_boxes[t * NSpans + span_index].enlarge(box);
            }
        });
        for (size_t t = 0; t < nslots; t++) {
            for (uint32_t i = 0; i < NSpans; i++) {
                span_counts[i] += parallel_counts[t * NSpans + i];
                span_boxes[i].enlarge(parallel_boxes[t * NSpans + i]);
            }
        }
    }

    // Accumulate boxes: spans 0 to NSpans-2 on the left, 1 to NSpans-1 on the right.
    Box left_boxes[NSplits];
    Box right_boxes[NSplits];
    Box box_accumulator = span_boxes[0];
    left_boxes[0]       = box_accumulator;
    for (uint32_t i = 1; i < NSplits; ++i) {
        box_accumulator.enlarge(span_boxes[i]);
        left_boxes[i] = box_accumulator;
    }
    box_accumulator          = span_boxes[NSpans - 1];
    right_boxes[NSplits - 1] = box_accumulator;
    for (uint32_t i = NSplits - 1; i > 0; --i) {
        box_accumulator.enlarge(span_boxes[i]);
        right_boxes[i - 1] = box_accumulator;
    }

    // Accumulate counts
    uint32_t left_counts[NSplits];
    uint32_t count_accumulator = span_counts[0];
    left_counts[0]             = count_accumulator;
    for (uint32_t spliti = 1; spliti < NSplits; ++spliti) {
        count_accumulator += span_counts[spliti];
        left_counts[spliti] = count_accumulator;
    }

    // Check which split is optimal, making sure at least 1/MinFraction of the boxes are on each
    // side.
    const uint32_t min_count          = nboxes / MinFraction;
    const uint32_t max_count          = ((MinFraction - 1) * uint64_t(nboxes)) / MinFraction;
    float          smallest_heuristic = Inf;
    uint32_t       split_index        = -1;
    for (uint32_t spliti = 0; spliti < NSplits; ++spliti) {
        const uint32_t left_count = left_counts[spliti];
        if (left_count < min_count || left_count > max_count)
            continue;
        const uint32_t right_count = nboxes - left_count;
        const float    heuristic   = left_count * left_boxes[spliti].half_surface_area() +
                                     right_count * right_boxes[spliti].half_surface_area();
        if (heuristic < smallest_heuristic) {
            smallest_heuristic = heuristic;
            split_index        = spliti;
        }
    }

    uint32_t* const indices_end = indices + nboxes;

    if (split_index == uint32_t(-1)) {
        // No split was anywhere close to balanced, so fall back to searching for one. Find the
        // span containing the balance point, namely where left_counts goes from being less than
        // min_count to more than max_count: if that is span 0 select max_count, if it is the last
        // span select min_count, else select nboxes/2.
        uint32_t* nth_index;
        if (left_counts[0] > max_count)
            nth_index = indices + max_count;
        else if (left_counts[NSplits - 1] < min_count)
            nth_index = indices + min_count;
        else
            nth_index = indices + nboxes / 2;
        select_nth_by_centre(boxes, indices, indices + nboxes, max_axis, nth_index);

        split_indices = nth_index;
        Box left_box  = boxes[indices[0]];
        for (uint32_t* left_indices = indices + 1; left_indices < nth_index; ++left_indices)
            left_box.enlarge(boxes[*left_indices]);
        Box right_box = boxes[nth_index[0]];
        for (uint32_t* right_indices = nth_index + 1; right_indices < indices_end; ++right_indices)
            right_box.enlarge(boxes[*right_indices]);
        split_boxes[0] = left_box;
        split_boxes[1] = right_box;
    }
    else {
        const float pivotx2 = axis_min_x2 + (split_index + 1) * axis_length / (NSpans / 2);
        uint32_t*   ppivot_start;
        uint32_t*   ppivot_end;
        partition_by_centre(
          boxes, indices, indices + nboxes, max_axis, pivotx2, ppivot_start, ppivot_end);

        split_indices = indices + left_counts[split_index];

        // Ignoring roundoff error we would have split_indices in [ppivot_start, ppivot_end], but
        // in practice it may not be.
        if (split_indices >= ppivot_start && split_indices <= ppivot_end) {
            split_boxes[0] = left_boxes[split_index];
            split_boxes[1] = right_boxes[split_index];
            return;
        }

        // Roundoff error changed the split, so recompute the boxes.
        split_indices = (split_indices < ppivot_start) ? ppivot_start : ppivot_end;

        // Emergency checks, just in case
        if (split_indices == indices)
            ++split_indices;
        else if (split_indices == indices_end)
            --split_indices;

        Box left_box = boxes[indices[0]];
        for (uint32_t* left_indices = indices + 1; left_indices < split_indices; ++left_indices)
            left_box.enlarge(boxes[*left_indices]);
        Box right_box = boxes[split_indices[0]];
        for (uint32_t* right_indices = split_indices + 1; right_indices < indices_end;
             ++right_indices)
            right_box.enlarge(boxes[*right_indices]);
        split_boxes[0] = left_box;
        split_boxes[1] = right_box;
    }
}

// Splits into BvhN parts by repeatedly splitting whichever part has the largest area times count.
void multi_split(const Box& axes_minmax,
                 const Box* boxes,
                 uint32_t*  indices,
                 uint32_t   nboxes,
                 uint32_t*  sub_indices[BvhN + 1],
                 Box        sub_boxes[BvhN])
{
    sub_indices[0] = indices;
    sub_indices[2] = indices + nboxes;
    split(axes_minmax, boxes, indices, nboxes, sub_indices[1], &sub_boxes[0]);

    float sub_box_areas[BvhN];
    sub_box_areas[0] = sub_boxes[0].half_surface_area();
    sub_box_areas[1] = sub_boxes[1].half_surface_area();
    for (uint32_t nsub = 2; nsub < BvhN; ++nsub) {
        // Choose which one to split
        uint32_t split_choice = uint32_t(-1);
        float    max_heuristic;
        for (uint32_t i = 0; i < nsub; ++i) {
            const uint32_t index_count = (sub_indices[i + 1] - sub_indices[i]);
            if (index_count > 1) {
                const float heuristic = sub_box_areas[i] * index_count;
                if (split_choice == uint32_t(-1) || heuristic > max_heuristic) {
                    split_choice  = i;
                    max_heuristic = heuristic;
                }
            }
        }

        uint32_t* selected_start = sub_indices[split_choice];
        uint32_t* selected_end   = sub_indices[split_choice + 1];

        // Shift results over; we can skip the one we selected.
        for (uint32_t i = nsub; i > split_choice; --i)
            sub_indices[i + 1] = sub_indices[i];
        for (uint32_t i = nsub - 1; i > split_choice; --i)
            sub_boxes[i + 1] = sub_boxes[i];
        for (uint32_t i = nsub - 1; i > split_choice; --i)
            sub_box_areas[i + 1] = sub_box_areas[i];

        split(sub_boxes[split_choice],
              boxes,
              selected_start,
              selected_end - selected_start,
              sub_indices[split_choice + 1],
              &sub_boxes[split_choice]);
        sub_box_areas[split_choice]     = sub_boxes[split_choice].half_surface_area();
        sub_box_areas[split_choice + 1] = sub_boxes[split_choice + 1].half_surface_area();
    }
}

constexpr uint32_t BuildParallelThreshold = 1024;

// Copies the subtrees that were built in parallel into the node array, renumbering their internal
// child indices to where they landed.
void adjust_parallel_child_nodes(uint32_t                 nparallel,
                                 std::vector<Node>&       nodes,
                                 const Node&              node,
                                 const std::vector<Node>* parallel_nodes,
                                 uint32_t* const*         sub_indices)
{
    uint32_t counted_parallel = 0;
    uint32_t childi           = 0;
    for (uint32_t taski = 0; taski < nparallel; taski++) {
        // First, find which child this is
        for (; childi < BvhN; ++childi) {
            const uint32_t sub_nboxes = sub_indices[childi + 1] - sub_indices[childi];
            if (sub_nboxes >= BuildParallelThreshold) {
                if (counted_parallel == taski)
                    break;
                ++counted_parallel;
            }
        }

        const std::vector<Node>& local_nodes       = parallel_nodes[counted_parallel];
        const uint32_t           n                 = local_nodes.size();
        const uint32_t           local_nodes_start = Node::internal_num(node.child[childi]) + 1;
        ++counted_parallel;
        ++childi;

        for (uint32_t j = 0; j < n; ++j) {
            Node local_node = local_nodes[j];
            for (uint32_t childj = 0; childj < BvhN; ++childj) {
                uint32_t local_child = local_node.child[childj];
                if (Node::is_internal(local_child) && local_child != Node::Empty)
                    local_node.child[childj] = local_child + local_nodes_start;
            }
            nodes[local_nodes_start + j] = local_node;
        }
    }
}

// Fills in node and appends its subtree to nodes. The caller has already made room for node
// itself, and reserved an overestimate of the whole tree, so nodes never reallocates under the
// reference.
void init_node(std::vector<Node>& nodes,
               Node&              node,
               const Box&         axes_minmax,
               const Box*         boxes,
               uint32_t*          indices,
               const uint32_t     nboxes)
{
    if (nboxes <= BvhN) {
        // Fits in one node
        for (uint32_t i = 0; i < nboxes; ++i)
            node.child[i] = indices[i];
        for (uint32_t i = nboxes; i < BvhN; ++i)
            node.child[i] = Node::Empty;
        return;
    }

    uint32_t* sub_indices[BvhN + 1];
    Box       sub_boxes[BvhN];
    multi_split(axes_minmax, boxes, indices, nboxes, sub_indices, sub_boxes);

    // Count the number of nodes to run in parallel and fill in single items in this node
    uint32_t nparallel = 0;
    for (uint32_t i = 0; i < BvhN; ++i) {
        const uint32_t sub_nboxes = sub_indices[i + 1] - sub_indices[i];
        if (sub_nboxes == 1)
            node.child[i] = sub_indices[i][0];
        else if (sub_nboxes >= BuildParallelThreshold)
            ++nparallel;
    }

    // NOTE: child nodes of this node need to be placed just before the nodes in their
    // corresponding subtree, in between the subtrees, because the parallel traversal uses the
    // difference between the child node IDs to determine the number of nodes in the subtree.

    if (nparallel >= 2) {
        std::vector<std::vector<Node>> parallel_nodes(nparallel);
        std::vector<Node>              parallel_parent_nodes(nparallel);
        parallel_for(0, nparallel, [&](size_t taski) {
            // First, find which child this is
            uint32_t counted_parallel = 0;
            uint32_t sub_nboxes       = 0;
            uint32_t childi;
            for (childi = 0; childi < BvhN; ++childi) {
                sub_nboxes = sub_indices[childi + 1] - sub_indices[childi];
                if (sub_nboxes >= BuildParallelThreshold) {
                    if (counted_parallel == taski)
                        break;
                    ++counted_parallel;
                }
            }

            std::vector<Node>& local_nodes = parallel_nodes[taski];
            // Preallocate an overestimate of the number of nodes needed.
            local_nodes.reserve(node_estimate(sub_nboxes));
            // We'll have to fix the internal node numbers in these later.
            init_node(local_nodes,
                      parallel_parent_nodes[taski],
                      sub_boxes[childi],
                      boxes,
                      sub_indices[childi],
                      sub_nboxes);
        });

        uint32_t counted_parallel = 0;
        for (uint32_t i = 0; i < BvhN; ++i) {
            const uint32_t sub_nboxes = sub_indices[i + 1] - sub_indices[i];
            if (sub_nboxes == 1)
                continue;

            uint32_t local_nodes_start = nodes.size();
            node.child[i]              = Node::mark_internal(local_nodes_start);
            if (sub_nboxes >= BuildParallelThreshold) {
                // First, adjust the root child node
                Node child_node = parallel_parent_nodes[counted_parallel];
                ++local_nodes_start;
                for (uint32_t childi = 0; childi < BvhN; ++childi) {
                    const uint32_t child_child = child_node.child[childi];
                    if (Node::is_internal(child_child) && child_child != Node::Empty)
                        child_node.child[childi] = child_child + local_nodes_start;
                }

                // Make space in the array for the sub-child nodes
                const uint32_t n = parallel_nodes[counted_parallel].size();
                ++counted_parallel;
                nodes.resize(local_nodes_start + n);
                nodes[local_nodes_start - 1] = child_node;
            }
            else {
                nodes.resize(local_nodes_start + 1);
                init_node(
                  nodes, nodes[local_nodes_start], sub_boxes[i], boxes, sub_indices[i], sub_nboxes);
            }
        }

        // Now, adjust and copy all sub-child nodes that were made in parallel
        adjust_parallel_child_nodes(nparallel, nodes, node, parallel_nodes.data(), sub_indices);
    }
    else {
        for (uint32_t i = 0; i < BvhN; ++i) {
            const uint32_t sub_nboxes = sub_indices[i + 1] - sub_indices[i];
            if (sub_nboxes == 1)
                continue;
            const uint32_t local_nodes_start = nodes.size();
            node.child[i]                    = Node::mark_internal(local_nodes_start);
            nodes.resize(local_nodes_start + 1);
            init_node(
              nodes, nodes[local_nodes_start], sub_boxes[i], boxes, sub_indices[i], sub_nboxes);
        }
    }
}

std::vector<Node> build_bvh(const Box* boxes, uint32_t nboxes)
{
    if (nboxes == 0)
        return {};

    std::vector<uint32_t> indices(nboxes);
    std::iota(indices.begin(), indices.end(), 0u);

    Box axes_minmax;
    compute_full_bounding_box(axes_minmax, boxes, nboxes, indices.data());

    // Exclude any boxes with NaNs or infinities, which shifts the good indices down over them.
    if (exclude_nan_inf_boxes(boxes, indices.data(), nboxes)) {
        if (nboxes == 0)
            return {};
        compute_full_bounding_box(axes_minmax, boxes, nboxes, indices.data());
    }

    std::vector<Node> nodes;
    // Preallocate an overestimate, which is also what keeps init_node's node reference valid.
    nodes.reserve(node_estimate(nboxes));
    nodes.resize(1);
    init_node(nodes, nodes[0], axes_minmax, boxes, indices.data(), nboxes);
    return nodes;
}

// ---------------------------------------------------------------------------------------------
// The solid angle of one triangle, exactly

// The signed solid angle triangle abc subtends from query.
//
// WARNING: this uses the right-handed normal convention, so negate the output, or swap b and c,
// if you want it positive inside and negative outside.
float signed_solid_angle_tri(const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& query)
{
    // Make a, b, and c relative to query
    Vec3f qa = a - query;
    Vec3f qb = b - query;
    Vec3f qc = c - query;

    const float alength = qa.length();
    const float blength = qb.length();
    const float clength = qc.length();

    // If any triangle vertex is coincident with query, query is on the surface, which we treat as
    // no solid angle.
    if (alength == 0 || blength == 0 || clength == 0)
        return 0.0f;

    qa /= alength;
    qb /= blength;
    qc /= clength;

    // The formula on Wikipedia has roughly dot(qa,cross(qb,qc)), but that is unstable when qa, qb
    // and qc are very close, e.g. if the input triangle was very far away. This is equivalent but
    // more stable.
    const float numerator = dot(qa, cross(qb - qa, qc - qa));

    // If the numerator is 0, regardless of denominator, query is on the surface.
    if (numerator == 0)
        return 0.0f;

    const float denominator = 1.0f + dot(qa, qb) + dot(qa, qc) + dot(qb, qc);

    return 2.0f * std::atan2(numerator, denominator);
}

// ---------------------------------------------------------------------------------------------
// Per-node Taylor expansion of the solid angle, to second order

// One node's four children, four lanes wide.
struct BoxData
{
    BoxData() {}

    /// An upper bound on the squared distance from average_p to the farthest point in the box.
    Float4 max_p_dist2;

    /// Centre of mass of the mesh surface in this box
    Vec3x4 average_p;

    /// Unnormalized, area-weighted normal of the mesh in this box
    Vec3x4 n;

    /// Values for Omega_1
    Vec3x4 nij_diag;  // Nxx, Nyy, Nzz
    Float4 nxy_nyx;   // Nxy+Nyx
    Float4 nyz_nzy;   // Nyz+Nzy
    Float4 nzx_nxz;   // Nzx+Nxz

    /// Values for Omega_2
    Vec3x4 nijk_diag;         // Nxxx, Nyyy, Nzzz
    Float4 sum_permute_nxyz;  // (Nxyz+Nxzy+Nyzx+Nyxz+Nzxy+Nzyx) = 2*(Nxyz+Nyzx+Nzxy)
    Float4 n2xxy_nyxx;        // Nxxy+Nxyx+Nyxx = 2Nxxy+Nyxx
    Float4 n2xxz_nzxx;        // Nxxz+Nxzx+Nzxx = 2Nxxz+Nzxx
    Float4 n2yyz_nzyy;        // Nyyz+Nyzy+Nzyy = 2Nyyz+Nzyy
    Float4 n2yyx_nxyy;        // Nyyx+Nyxy+Nxyy = 2Nyyx+Nxyy
    Float4 n2zzx_nxzz;        // Nzzx+Nzxz+Nxzz = 2Nzzx+Nxzz
    Float4 n2zzy_nyzz;        // Nzzy+Nzyz+Nyzz = 2Nzzy+Nyzz
};

// What a node hands up to its parent, one lane wide, and only needed while precomputing.
struct LocalData
{
    Box box;

    // P and N are needed from each child for computing Nij.
    Vec3f average_p;
    Vec3f area_p;
    Vec3f n;

    // Unsigned area is needed for computing the average position.
    float area;

    // These are needed for computing Nijk.
    Vec3f nij_diag;
    float nxy, nyx;
    float nyz, nzy;
    float nzx, nxz;

    Vec3f nijk_diag;  // Nxxx, Nyyy, Nzzz
    float sum_permute_nxyz;
    float n2xxy_nyxx;
    float n2xxz_nzxx;
    float n2yyz_nzyy;
    float n2yyx_nxyy;
    float n2zzx_nxzz;
    float n2zzy_nyzz;
};

// Subtrees of at least this many nodes are precomputed in parallel with their siblings.
constexpr uint32_t PrecomputeParallelThreshold = 4096;

struct Precompute
{
    BoxData*     box_data;
    const Node*  nodes;
    const Box*   triangle_boxes;
    const int*   triangle_points;
    const Vec3f* positions;

    void triangle(uint32_t itemi, LocalData& out) const;
    bool child_data(const Node& n, uint32_t s, LocalData* child) const;
    void combine(uint32_t nodei, LocalData* out, uint32_t nchildren, const LocalData* child) const;
    void node(uint32_t nodei, LocalData* out) const;
    void node_parallel(uint32_t nodei, uint32_t next_node_id, LocalData* out) const;
};

void Precompute::triangle(uint32_t itemi, LocalData& out) const
{
    const int* const cur_triangle_points = triangle_points + 3 * itemi;
    const Vec3f      a                   = positions[cur_triangle_points[0]];
    const Vec3f      b                   = positions[cur_triangle_points[1]];
    const Vec3f      c                   = positions[cur_triangle_points[2]];
    const Vec3f      ab                  = b - a;
    const Vec3f      ac                  = c - a;

    out.box = triangle_boxes[itemi];

    // Area-weighted normal (unnormalized)
    const Vec3f N     = 0.5f * cross(ab, ac);
    const float area2 = N.length2();
    const float area  = std::sqrt(area2);
    const Vec3f P     = (a + b + c) / 3.0f;
    out.average_p     = P;
    out.area_p        = P * area;
    out.n             = N;
    out.area          = area;

    // NOTE: because P is at the centroid, a triangle contributes Nij = 0.
    out.nij_diag = 0.0f;
    out.nxy      = 0;
    out.nyx      = 0;
    out.nyz      = 0;
    out.nzy      = 0;
    out.nzx      = 0;
    out.nxz      = 0;

    // If it's zero-length, the results are zero, so we can skip.
    if (area == 0) {
        out.nijk_diag        = 0.0f;
        out.sum_permute_nxyz = 0;
        out.n2xxy_nyxx       = 0;
        out.n2xxz_nzxx       = 0;
        out.n2yyz_nzyy       = 0;
        out.n2yyx_nxyy       = 0;
        out.n2zzx_nxzz       = 0;
        out.n2zzy_nyzz       = 0;
        return;
    }

    // We need to use the NORMALIZED normal to multiply the integrals by.
    Vec3f n = N / area;

    // Figure out the order of a, b, and c in x, y, and z for computing the integrals for Nijk.
    Vec3f values[3] = {a, b, c};

    int order_x[3] = {0, 1, 2};
    if (a[0] > b[0])
        std::swap(order_x[0], order_x[1]);
    if (values[order_x[0]][0] > c[0])
        std::swap(order_x[0], order_x[2]);
    if (values[order_x[1]][0] > values[order_x[2]][0])
        std::swap(order_x[1], order_x[2]);
    float dx = values[order_x[2]][0] - values[order_x[0]][0];

    int order_y[3] = {0, 1, 2};
    if (a[1] > b[1])
        std::swap(order_y[0], order_y[1]);
    if (values[order_y[0]][1] > c[1])
        std::swap(order_y[0], order_y[2]);
    if (values[order_y[1]][1] > values[order_y[2]][1])
        std::swap(order_y[1], order_y[2]);
    float dy = values[order_y[2]][1] - values[order_y[0]][1];

    int order_z[3] = {0, 1, 2};
    if (a[2] > b[2])
        std::swap(order_z[0], order_z[1]);
    if (values[order_z[0]][2] > c[2])
        std::swap(order_z[0], order_z[2]);
    if (values[order_z[1]][2] > values[order_z[2]][2])
        std::swap(order_z[1], order_z[2]);
    float dz = values[order_z[2]][2] - values[order_z[0]][2];

    auto&& compute_integrals = [](const Vec3f& a,
                                  const Vec3f& b,
                                  const Vec3f& c,
                                  const Vec3f& P,
                                  float*       integral_ii,
                                  float*       integral_ij,
                                  float*       integral_ik,
                                  const int    i) {
        // NOTE: a, b, and c must be in order of the i axis. We're splitting the triangle at the
        // middle i coordinate.
        const Vec3f oab = b - a;
        const Vec3f oac = c - a;
        const Vec3f ocb = b - c;
        const float t   = oab[i] / oac[i];

        const int   j     = (i == 2) ? 0 : (i + 1);
        const int   k     = (j == 2) ? 0 : (j + 1);
        const float jdiff = t * oac[j] - oab[j];
        const float kdiff = t * oac[k] - oab[k];
        Vec3f       cross_a;
        cross_a[0] = (jdiff * oab[k] - kdiff * oab[j]);
        cross_a[1] = kdiff * oab[i];
        cross_a[2] = jdiff * oab[i];
        Vec3f cross_c;
        cross_c[0]               = (jdiff * ocb[k] - kdiff * ocb[j]);
        cross_c[1]               = kdiff * ocb[i];
        cross_c[2]               = jdiff * ocb[i];
        const float area_scale_a = cross_a.length();
        const float area_scale_c = cross_c.length();
        const float Pai          = a[i] - P[i];
        const float Pci          = c[i] - P[i];

        // Integral over the area of the triangle of (pi^2)dA, by splitting the triangle into two
        // at b, the a side and the c side.
        const float int_ii_a = area_scale_a * (0.5f * Pai * Pai + (2.0f / 3.0f) * Pai * oab[i] +
                                               0.25f * oab[i] * oab[i]);
        const float int_ii_c = area_scale_c * (0.5f * Pci * Pci + (2.0f / 3.0f) * Pci * ocb[i] +
                                               0.25f * ocb[i] * ocb[i]);
        *integral_ii         = int_ii_a + int_ii_c;

        int    jk       = j;
        float* integral = integral_ij;
        float  diff     = jdiff;
        while (true)  // This only does 2 iterations, one for j and one for k
        {
            if (integral) {
                float obmidj  = b[jk] + 0.5f * diff;
                float oabmidj = obmidj - a[jk];
                float ocbmidj = obmidj - c[jk];
                float Paj     = a[jk] - P[jk];
                float Pcj     = c[jk] - P[jk];
                // Integral over the area of the triangle of (pi*pj)dA
                const float int_ij_a =
                  area_scale_a * (0.5f * Pai * Paj + (1.0f / 3.0f) * Pai * oabmidj +
                                  (1.0f / 3.0f) * Paj * oab[i] + 0.25f * oab[i] * oabmidj);
                const float int_ij_c =
                  area_scale_c * (0.5f * Pci * Pcj + (1.0f / 3.0f) * Pci * ocbmidj +
                                  (1.0f / 3.0f) * Pcj * ocb[i] + 0.25f * ocb[i] * ocbmidj);
                *integral = int_ij_a + int_ij_c;
            }
            if (jk == k)
                break;
            jk       = k;
            integral = integral_ik;
            diff     = kdiff;
        }
    };

    float integral_xx = 0;
    float integral_xy = 0;
    float integral_yy = 0;
    float integral_yz = 0;
    float integral_zz = 0;
    float integral_zx = 0;
    // Note that if the span of any axis is zero, the integral must be zero, since there's a factor
    // of (p_i-P_i), i.e. value minus average, and every value must be equal to the average.
    if (dx > 0) {
        compute_integrals(values[order_x[0]],
                          values[order_x[1]],
                          values[order_x[2]],
                          P,
                          &integral_xx,
                          ((dx >= dy && dy > 0) ? &integral_xy : nullptr),
                          ((dx >= dz && dz > 0) ? &integral_zx : nullptr),
                          0);
    }
    if (dy > 0) {
        compute_integrals(values[order_y[0]],
                          values[order_y[1]],
                          values[order_y[2]],
                          P,
                          &integral_yy,
                          ((dy >= dz && dz > 0) ? &integral_yz : nullptr),
                          ((dx < dy && dx > 0) ? &integral_xy : nullptr),
                          1);
    }
    if (dz > 0) {
        compute_integrals(values[order_z[0]],
                          values[order_z[1]],
                          values[order_z[2]],
                          P,
                          &integral_zz,
                          ((dx < dz && dx > 0) ? &integral_zx : nullptr),
                          ((dy < dz && dy > 0) ? &integral_yz : nullptr),
                          2);
    }

    Vec3f niii;
    niii[0] = integral_xx;
    niii[1] = integral_yy;
    niii[2] = integral_zz;
    niii *= n;
    out.nijk_diag        = niii;
    out.sum_permute_nxyz = 2 * (n[0] * integral_yz + n[1] * integral_zx + n[2] * integral_xy);
    float nxxy           = n[0] * integral_xy;
    float nxxz           = n[0] * integral_zx;
    float nyyz           = n[1] * integral_yz;
    float nyyx           = n[1] * integral_xy;
    float nzzx           = n[2] * integral_zx;
    float nzzy           = n[2] * integral_yz;
    out.n2xxy_nyxx       = 2 * nxxy + n[1] * integral_xx;
    out.n2xxz_nzxx       = 2 * nxxz + n[2] * integral_xx;
    out.n2yyz_nzyy       = 2 * nyyz + n[2] * integral_yy;
    out.n2yyx_nxyy       = 2 * nyyx + n[0] * integral_yy;
    out.n2zzx_nxzz       = 2 * nzzx + n[0] * integral_zz;
    out.n2zzy_nyzz       = 2 * nzzy + n[1] * integral_zz;
}

void Precompute::combine(uint32_t         nodei,
                         LocalData*       out,
                         uint32_t         nchildren,
                         const LocalData* child) const
{
    BoxData& data = box_data[nodei];

    Vec3f N                = child[0].n;
    data.n[0].v[0]         = N[0];
    data.n[1].v[0]         = N[1];
    data.n[2].v[0]         = N[2];
    Vec3f areaP            = child[0].area_p;
    float area             = child[0].area;
    Vec3f local_P          = child[0].average_p;
    data.average_p[0].v[0] = local_P[0];
    data.average_p[1].v[0] = local_P[1];
    data.average_p[2].v[0] = local_P[2];
    for (uint32_t i = 1; i < nchildren; ++i) {
        const Vec3f local_N = child[i].n;
        N += local_N;
        data.n[0].v[i] = local_N[0];
        data.n[1].v[i] = local_N[1];
        data.n[2].v[i] = local_N[2];
        areaP += child[i].area_p;
        area += child[i].area;
        const Vec3f local_P2   = child[i].average_p;
        data.average_p[0].v[i] = local_P2[0];
        data.average_p[1].v[i] = local_P2[1];
        data.average_p[2].v[i] = local_P2[2];
    }
    for (uint32_t i = nchildren; i < BvhN; ++i) {
        // Set to zero, just to avoid false positives for uses of uninitialized memory.
        data.n[0].v[i]         = 0;
        data.n[1].v[i]         = 0;
        data.n[2].v[i]         = 0;
        data.average_p[0].v[i] = 0;
        data.average_p[1].v[i] = 0;
        data.average_p[2].v[i] = 0;
    }
    out->n      = N;
    out->area_p = areaP;
    out->area   = area;

    Box box = child[0].box;
    for (uint32_t i = 1; i < nchildren; ++i)
        box.enlarge(child[i].box);

    // Normalize P
    Vec3f averageP;
    if (area > 0)
        averageP = areaP / area;
    else
        averageP = 0.5f * (box.min() + box.max());
    out->average_p = averageP;
    out->box       = box;

    for (uint32_t i = 0; i < nchildren; ++i) {
        const Box&   local_box = child[i].box;
        const Vec3f& local_P2  = child[i].average_p;
        const Vec3f  maxPDiff =
          max_components(local_P2 - local_box.min(), local_box.max() - local_P2);
        data.max_p_dist2.v[i] = maxPDiff.length2();
    }
    for (uint32_t i = nchildren; i < BvhN; ++i) {
        // This child is non-existent. Infinity means the approximation is never used for it, and
        // the query traversal can check for Empty.
        data.max_p_dist2.v[i] = Inf;
    }

    // We now have the current box's P, so we can adjust Nij and Nijk
    out->nij_diag         = child[0].nij_diag;
    out->nxy              = 0;
    out->nyx              = 0;
    out->nyz              = 0;
    out->nzy              = 0;
    out->nzx              = 0;
    out->nxz              = 0;
    out->nijk_diag        = child[0].nijk_diag;
    out->sum_permute_nxyz = child[0].sum_permute_nxyz;
    out->n2xxy_nyxx       = child[0].n2xxy_nyxx;
    out->n2xxz_nzxx       = child[0].n2xxz_nzxx;
    out->n2yyz_nzyy       = child[0].n2yyz_nzyy;
    out->n2yyx_nxyy       = child[0].n2yyx_nxyy;
    out->n2zzx_nxzz       = child[0].n2zzx_nxzz;
    out->n2zzy_nyzz       = child[0].n2zzy_nyzz;

    for (uint32_t i = 1; i < nchildren; ++i) {
        out->nij_diag += child[i].nij_diag;
        out->nijk_diag += child[i].nijk_diag;
        out->sum_permute_nxyz += child[i].sum_permute_nxyz;
        out->n2xxy_nyxx += child[i].n2xxy_nyxx;
        out->n2xxz_nzxx += child[i].n2xxz_nzxx;
        out->n2yyz_nzyy += child[i].n2yyz_nzyy;
        out->n2yyx_nxyy += child[i].n2yyx_nxyy;
        out->n2zzx_nxzz += child[i].n2zzx_nxzz;
        out->n2zzy_nyzz += child[i].n2zzy_nyzz;
    }
    for (uint32_t i = 0; i < nchildren; ++i) {
        for (int j = 0; j < 3; ++j)
            data.nij_diag[j].v[i] = child[i].nij_diag[j];
        data.nxy_nyx.v[i] = child[i].nxy + child[i].nyx;
        data.nyz_nzy.v[i] = child[i].nyz + child[i].nzy;
        data.nzx_nxz.v[i] = child[i].nzx + child[i].nxz;
        for (int j = 0; j < 3; ++j)
            data.nijk_diag[j].v[i] = child[i].nijk_diag[j];
        data.sum_permute_nxyz.v[i] = child[i].sum_permute_nxyz;
        data.n2xxy_nyxx.v[i]       = child[i].n2xxy_nyxx;
        data.n2xxz_nzxx.v[i]       = child[i].n2xxz_nzxx;
        data.n2yyz_nzyy.v[i]       = child[i].n2yyz_nzyy;
        data.n2yyx_nxyy.v[i]       = child[i].n2yyx_nxyy;
        data.n2zzx_nxzz.v[i]       = child[i].n2zzx_nxzz;
        data.n2zzy_nyzz.v[i]       = child[i].n2zzy_nyzz;
    }
    for (uint32_t i = nchildren; i < BvhN; ++i) {
        // Set to zero, just to avoid false positives for uses of uninitialized memory.
        for (int j = 0; j < 3; ++j)
            data.nij_diag[j].v[i] = 0;
        data.nxy_nyx.v[i] = 0;
        data.nyz_nzy.v[i] = 0;
        data.nzx_nxz.v[i] = 0;
        for (int j = 0; j < 3; ++j)
            data.nijk_diag[j].v[i] = 0;
        data.sum_permute_nxyz.v[i] = 0;
        data.n2xxy_nyxx.v[i]       = 0;
        data.n2xxz_nzxx.v[i]       = 0;
        data.n2yyz_nzyy.v[i]       = 0;
        data.n2yyx_nxyy.v[i]       = 0;
        data.n2zzx_nxzz.v[i]       = 0;
        data.n2zzy_nyzz.v[i]       = 0;
    }

    for (uint32_t i = 0; i < nchildren; ++i) {
        const LocalData& child_data   = child[i];
        Vec3f            displacement = child_data.average_p - out->average_p;
        Vec3f            N2           = child_data.n;

        // Adjust Nij for the change in centre P
        out->nij_diag += N2 * displacement;
        float nxy = child_data.nxy + N2[0] * displacement[1];
        float nyx = child_data.nyx + N2[1] * displacement[0];
        float nyz = child_data.nyz + N2[1] * displacement[2];
        float nzy = child_data.nzy + N2[2] * displacement[1];
        float nzx = child_data.nzx + N2[2] * displacement[0];
        float nxz = child_data.nxz + N2[0] * displacement[2];

        out->nxy += nxy;
        out->nyx += nyx;
        out->nyz += nyz;
        out->nzy += nzy;
        out->nzx += nzx;
        out->nxz += nxz;

        // Adjust Nijk for the change in centre P
        out->nijk_diag +=
          2.0f * displacement * child_data.nij_diag + displacement * displacement * child_data.n;
        out->sum_permute_nxyz += (displacement[0] * (nyz + nzy) + displacement[1] * (nzx + nxz) +
                                  displacement[2] * (nxy + nyx));
        out->n2xxy_nyxx +=
          2 * (displacement[1] * child_data.nij_diag[0] + displacement[0] * child_data.nxy +
               N2[0] * displacement[0] * displacement[1]) +
          2 * child_data.nyx * displacement[0] + N2[1] * displacement[0] * displacement[0];
        out->n2xxz_nzxx +=
          2 * (displacement[2] * child_data.nij_diag[0] + displacement[0] * child_data.nxz +
               N2[0] * displacement[0] * displacement[2]) +
          2 * child_data.nzx * displacement[0] + N2[2] * displacement[0] * displacement[0];
        out->n2yyz_nzyy +=
          2 * (displacement[2] * child_data.nij_diag[1] + displacement[1] * child_data.nyz +
               N2[1] * displacement[1] * displacement[2]) +
          2 * child_data.nzy * displacement[1] + N2[2] * displacement[1] * displacement[1];
        out->n2yyx_nxyy +=
          2 * (displacement[0] * child_data.nij_diag[1] + displacement[1] * child_data.nyx +
               N2[1] * displacement[1] * displacement[0]) +
          2 * child_data.nxy * displacement[1] + N2[0] * displacement[1] * displacement[1];
        out->n2zzx_nxzz +=
          2 * (displacement[0] * child_data.nij_diag[2] + displacement[2] * child_data.nzx +
               N2[2] * displacement[2] * displacement[0]) +
          2 * child_data.nxz * displacement[2] + N2[0] * displacement[2] * displacement[2];
        out->n2zzy_nyzz +=
          2 * (displacement[1] * child_data.nij_diag[2] + displacement[2] * child_data.nzy +
               N2[2] * displacement[2] * displacement[1]) +
          2 * child_data.nyz * displacement[2] + N2[1] * displacement[2] * displacement[2];
    }
}

// Fills child[s] from the s-th child of n: recursing for an internal node, or the one triangle
// for a leaf entry. False at an empty slot, and anything after one is empty too, so the callers
// below all stop there.
bool Precompute::child_data(const Node& n, uint32_t s, LocalData* child) const
{
    const uint32_t c = n.child[s];
    if (!Node::is_internal(c)) {
        triangle(c, child[s]);
        return true;
    }
    if (c == Node::Empty)
        return false;
    node(Node::internal_num(c), &child[s]);
    return true;
}

void Precompute::node(uint32_t nodei, LocalData* out) const
{
    const Node& n = nodes[nodei];
    LocalData   child[BvhN];
    uint32_t    s;
    for (s = 0; s < BvhN; ++s) {
        if (!child_data(n, s, child))
            break;
    }
    // NOTE: s is now the number of non-empty entries in this node.
    combine(nodei, out, s, child);
}

void Precompute::node_parallel(uint32_t nodei, uint32_t next_node_id, LocalData* out) const
{
    const Node& n = nodes[nodei];

    // To determine the number of nodes in a child's subtree, take the next node ID minus the
    // current child's node ID.
    uint32_t next_nodes[BvhN];
    uint32_t nnodes[BvhN] = {0, 0, 0, 0};
    uint32_t nchildren    = BvhN;
    uint32_t nparallel    = 0;
    for (uint32_t s = BvhN - 1; s < BvhN; --s) {
        const uint32_t c = n.child[s];
        if (c == Node::Empty) {
            --nchildren;
            continue;
        }
        next_nodes[s] = next_node_id;
        if (Node::is_internal(c)) {
            // NOTE: this depends on init_node appending the child nodes in between their content,
            // instead of all at once.
            const uint32_t child_node_id = Node::internal_num(c);
            nnodes[s]                    = next_node_id - child_node_id;
            next_node_id                 = child_node_id;
        }
        nparallel += (nnodes[s] >= PrecomputeParallelThreshold);
    }

    LocalData child[BvhN];
    if (nparallel >= 2) {
        // Do any non-parallel ones first
        for (uint32_t s = 0; s < BvhN; ++s) {
            if (nnodes[s] >= PrecomputeParallelThreshold)
                continue;
            if (!child_data(n, s, child))
                break;
        }
        parallel_for(0, nparallel, [&](size_t taski) {
            // First, find which child this is
            uint32_t parallel_count = 0;
            uint32_t s;
            for (s = 0; s < BvhN; ++s) {
                if (nnodes[s] < PrecomputeParallelThreshold)
                    continue;
                if (parallel_count == taski)
                    break;
                ++parallel_count;
            }
            const uint32_t c = n.child[s];
            if (Node::is_internal(c))
                node_parallel(Node::internal_num(c), next_nodes[s], &child[s]);
            else
                triangle(c, child[s]);
        });
        combine(nodei, out, nchildren, child);
    }
    else {
        // Empty children sit at the tail, so node() combines the same nchildren entries.
        node(nodei, out);
    }
}

// ---------------------------------------------------------------------------------------------
// Querying

// The value of (maxP/q) beyond which a box's approximation is accurate enough to use, squared.
constexpr float AccuracyScale2 = 2.0f * 2.0f;

constexpr uint32_t AllChildBits = (uint32_t(1) << BvhN) - 1;

struct Query
{
    const BoxData* box_data;
    const Node*    nodes;
    const int*     triangle_points;
    const Vec3f*   positions;
    Vec3f          query_point;

    float approximate(uint32_t nodei, uint32_t& descend_bits) const;
    float node(uint32_t nodei) const;
};

// The approximation of the node's four children, summed over those it is accurate enough for.
// The rest come back in descend_bits, to be descended into instead.
float Query::approximate(uint32_t nodei, uint32_t& descend_bits) const
{
    const BoxData& data  = box_data[nodei];
    const Float4   maxP2 = data.max_p_dist2;
    Vec3x4         q;
    q[0] = splat(query_point[0]);
    q[1] = splat(query_point[1]);
    q[2] = splat(query_point[2]);
    q -= data.average_p;
    const Float4 qlength2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];

    // If the query point is within a factor of accuracy_scale of the box radius, the
    // approximation is assumed not good enough, so it needs to descend.
    const Mask4 descend_mask = (qlength2 <= maxP2 * AccuracyScale2);
    descend_bits             = mask_bits(descend_mask);
    if (descend_bits == AllChildBits)
        return 0;

    // qlength2 must be non-zero, since it's strictly greater than something. We still need to be
    // careful of NaNs, though, because the 4th power might cause problems.
    const Float4 qlength_m2 = splat(1.0f) / qlength2;
    const Float4 qlength_m1 = sqrt(qlength_m2);

    // Normalize q to reduce issues with overflow/underflow, since we'd need the 7th power if we
    // didn't normalize, and (1e-6)^-7 = 1e42, which overflows single precision.
    q *= qlength_m1;

    Float4 Omega_approx = -qlength_m2 * dot(q, data.n);

    const Vec3x4 q2         = q * q;
    const Float4 qlength_m3 = qlength_m2 * qlength_m1;
    const Float4 Omega_1 =
      qlength_m3 * (data.nij_diag[0] + data.nij_diag[1] + data.nij_diag[2] -
                    splat(3.0f) * (dot(q2, data.nij_diag) + q[0] * q[1] * data.nxy_nyx +
                                    q[0] * q[2] * data.nzx_nxz + q[1] * q[2] * data.nyz_nzy));
    Omega_approx += Omega_1;

    const Vec3x4 q3         = q2 * q;
    const Float4 qlength_m4 = qlength_m2 * qlength_m2;
    Float4       temp0[3]   = {data.n2yyx_nxyy + data.n2zzx_nxzz,
                               data.n2zzy_nyzz + data.n2xxy_nyxx,
                               data.n2xxz_nzxx + data.n2yyz_nzyy};
    Float4       temp1[3]   = {q[1] * data.n2xxy_nyxx + q[2] * data.n2xxz_nzxx,
                               q[2] * data.n2yyz_nzyy + q[0] * data.n2yyx_nxyy,
                               q[0] * data.n2zzx_nxzz + q[1] * data.n2zzy_nyzz};
    const Float4 Omega_2 =
      qlength_m4 *
      (splat(1.5f) * dot(q, splat(3) * data.nijk_diag + Vec3x4(temp0)) -
       splat(7.5f) * (dot(q3, data.nijk_diag) + q[0] * q[1] * q[2] * data.sum_permute_nxyz +
                       dot(q2, Vec3x4(temp1))));
    Omega_approx += Omega_2;

    // If q is so small that we got NaNs and we just have a small bounding box, it needs to
    // descend.
    const Mask4 mask = is_finite(Omega_approx) & ~descend_mask;
    Omega_approx     = Omega_approx & mask;
    descend_bits     = (~mask_bits(mask)) & AllChildBits;

    float sum = Omega_approx.v[0];
    for (uint32_t i = 1; i < BvhN; ++i)
        sum += Omega_approx.v[i];
    return sum;
}

float Query::node(uint32_t nodei) const
{
    uint32_t    descend;
    const float approximated = approximate(nodei, descend);
    if (!descend)
        return approximated;

    const Node& n = nodes[nodei];
    float       child[BvhN];
    uint32_t    s;
    for (s = 0; s < BvhN; ++s) {
        if (!((descend >> s) & 1))
            continue;
        const uint32_t c = n.child[s];
        if (Node::is_internal(c)) {
            // NOTE: anything after an empty child is empty too, so we can break.
            if (c == Node::Empty) {
                descend &= (uint32_t(1) << s) - 1;
                break;
            }
            child[s] = node(Node::internal_num(c));
        }
        else {
            const int* const cur_triangle_points = triangle_points + 3 * c;
            child[s] = signed_solid_angle_tri(positions[cur_triangle_points[0]],
                                              positions[cur_triangle_points[1]],
                                              positions[cur_triangle_points[2]],
                                              query_point);
        }
    }
    // NOTE: s is BvhN unless the loop broke on an Empty child, in which case everything at s and
    // beyond is empty and its descend bits were just cleared.
    float sum = (descend & 1) ? child[0] : 0;
    for (uint32_t i = 1; i < s; ++i)
        sum += ((descend >> i) & 1) ? child[i] : 0;
    return approximated + sum;
}

}  // namespace

std::vector<float> solid_angles(const std::vector<float>& points,
                                const std::vector<int>&   triangles,
                                const std::vector<float>& queries)
{
    const int* triangle_points = triangles.data();
    // Vec3f is three floats and nothing else, so the flat arrays are read in place.
    const Vec3f*   positions    = reinterpret_cast<const Vec3f*>(points.data());
    const Vec3f*   query_points = reinterpret_cast<const Vec3f*>(queries.data());
    const size_t   nqueries     = queries.size() / 3;
    const uint32_t ntriangles   = triangles.size() / 3;

    std::vector<Box> triangle_boxes(ntriangles);
    const auto       fill_box = [&](size_t i) {
        const int* cur_triangle_points = triangle_points + i * 3;
        Box&       box                 = triangle_boxes[i];
        box.set(positions[cur_triangle_points[0]]);
        box.enlarge(positions[cur_triangle_points[1]]);
        box.enlarge(positions[cur_triangle_points[2]]);
    };
    if (ntriangles < 16 * 1024) {
        for (uint32_t i = 0; i < ntriangles; ++i)
            fill_box(i);
    }
    else {
        parallel_for(0, ntriangles, fill_box);
    }

    const std::vector<Node> nodes = build_bvh(triangle_boxes.data(), ntriangles);

    std::vector<float> out(nqueries, 0.0f);
    if (nodes.empty())
        return out;

    std::vector<BoxData> box_data(nodes.size());
    const Precompute     precompute {
      box_data.data(), nodes.data(), triangle_boxes.data(), triangle_points, positions};
    LocalData root;
    precompute.node_parallel(0, nodes.size(), &root);

    parallel_for(0, nqueries, [&](size_t p) {
        const Query query {
          box_data.data(), nodes.data(), triangle_points, positions, query_points[p]};
        out[p] = query.node(0);
    });
    return out;
}

}  // namespace FastWindingNumber
}  // namespace floatTetWild
