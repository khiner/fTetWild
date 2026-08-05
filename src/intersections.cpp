// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "intersections.h"

#include "Predicates.hpp"

namespace floatTetWild {
namespace {
// Where the lines through seg1 and seg2 meet, as t1 along seg1 and t2 along seg2. False when the
// two are parallel.
//
// Assumptions: the segments are not degenerate and not coplanar.
bool line_line_intersection_2d(const std::array<Vector2, 2> &seg1, const std::array<Vector2, 2> &seg2,
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
}  // namespace
}  // namespace floatTetWild

bool floatTetWild::seg_line_intersection_2d(const std::array<Vector2, 2> &seg, const std::array<Vector2, 2> &line, Scalar& t_seg){
    Scalar _;
    if (!line_line_intersection_2d(seg, line, t_seg, _))
        return false;
    return t_seg >= 0 && t_seg <= 1;
}

bool floatTetWild::seg_seg_intersection_2d(const std::array<Vector2, 2> &seg1, const std::array<Vector2, 2> &seg2, Scalar& t2){
    Scalar t1;
    if (!line_line_intersection_2d(seg1, seg2, t1, t2))
        return false;
    return t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1;
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

    if (fabs(D) == 0)  // parallel to the plane, whether or not it lies in it
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

void floatTetWild::count_orientations(const int* oris, int n, int& cnt_pos, int& cnt_neg) {
    cnt_pos = 0;
    cnt_neg = 0;
    for (int i = 0; i < n; i++) {
        if (oris[i] == Predicates::ORI_POSITIVE)
            cnt_pos++;
        else if (oris[i] == Predicates::ORI_NEGATIVE)
            cnt_neg++;
    }
}

bool floatTetWild::is_p_inside_tri_2d(const Vector2& p, const std::array<Vector2, 3> &tri) {
    int oris[3];
    for (int i = 0; i < 3; i++)
        oris[i] = Predicates::orient_2d(tri[i], tri[(i + 1) % 3], p);

    int cnt_pos, cnt_neg;
    count_orientations(oris, 3, cnt_pos, cnt_neg);
    return cnt_neg == 3 || cnt_pos == 3;  // strictly inside
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
    Vector3 proj_p = p - dist * n;
    return Vector2(proj_p[(t + 1) % 3], proj_p[(t + 2) % 3]);
}

bool floatTetWild::is_crossing(int s1, int s2) {
    return (s1 == Predicates::ORI_POSITIVE && s2 == Predicates::ORI_NEGATIVE)
        || (s2 == Predicates::ORI_POSITIVE && s1 == Predicates::ORI_NEGATIVE);
}

bool floatTetWild::is_tri_tri_cut(const Vector3& p1, const Vector3& p2, const Vector3& p3,
                                  const Vector3& q1, const Vector3& q2, const Vector3& q3, int hint) {
    if(hint == CUT_COPLANAR){
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

    int cnt_pos, cnt_neg;
    count_orientations(oris, 4, cnt_pos, cnt_neg);
    return cnt_pos == 0 || cnt_neg == 0;
}