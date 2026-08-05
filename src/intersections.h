// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_INTERSECTIONS_H
#define FLOATTETWILD_INTERSECTIONS_H

#include "Types.hpp"

#include <initializer_list>

namespace floatTetWild {
    // The kind of cut is_tri_tri_cut is asked about. The three edge cases are 0..2 so that they
    // double as the local edge index.
    constexpr int CUT_EDGE_0 = 0;
    constexpr int CUT_EDGE_2 = 2;
    constexpr int CUT_FACE = 3;
    constexpr int CUT_COPLANAR = 4;

    Scalar p_seg_squared_dist_3d(const Vector3 &p, const Vector3 &a, const Vector3 &b);

    bool is_p_inside_tri_2d(const Vector2& p, const std::array<Vector2, 3> &tri);
    bool is_tri_tri_cutted_2d(const std::array<Vector2, 3> &p_tet, const std::array<Vector2, 3> &p_tri);

    bool seg_seg_intersection_2d(const std::array<Vector2, 2> &seg1, const std::array<Vector2, 2> &seg2, Scalar& t2);
    bool seg_line_intersection_2d(const std::array<Vector2, 2> &seg, const std::array<Vector2, 2> &line, Scalar& t_seg);
    bool seg_plane_intersection(const Vector3 &p1, const Vector3 &p2, const Vector3 &a, const Vector3 &n,
                                Vector3 &p);

    int get_t(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2);
    Vector2 to_2d(const Vector3 &p, int t);
    Vector2 to_2d(const Vector3 &p, const Vector3& n, const Vector3& pp, int t);

    bool is_crossing(int s1, int s2);

    bool is_tri_tri_cut(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,//cutting tri
                        const Vector3 &q1, const Vector3 &q2, const Vector3 &q3, int hint);//face of tet

    // The corner-wise minimum and maximum of the given points.
    void get_bbox(std::initializer_list<Vector3> ps, Vector3& min, Vector3& max);

    bool is_bbox_intersected(const Vector3& min1, const Vector3& max1, const Vector3& min2, const Vector3& max2);

    bool is_point_inside_tet(const Vector3& p, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t);

    // How many of the orientations are strictly positive and how many strictly negative. A zero is
    // neither, so n - cnt_pos - cnt_neg is how many lie on the plane they were taken against.
    void count_orientations(const int* oris, int n, int& cnt_pos, int& cnt_neg);
}

// Guigue-Devillers, in triangle_triangle_intersection.cpp, which has no header of its own. Returns
// 1 when the triangles overlap, and then source and target are the ends of the shared segment.
int tri_tri_intersection_test_3d(const floatTetWild::Scalar p1[3], const floatTetWild::Scalar q1[3],
                                 const floatTetWild::Scalar r1[3], const floatTetWild::Scalar p2[3],
                                 const floatTetWild::Scalar q2[3], const floatTetWild::Scalar r2[3],
                                 floatTetWild::Scalar source[3], floatTetWild::Scalar target[3]);

#endif //FLOATTETWILD_INTERSECTIONS_H
