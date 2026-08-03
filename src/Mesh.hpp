// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Parameters.h>
#include <floattetwild/Types.hpp>

#include <vector>
#include <array>
#include <cassert>

namespace floatTetWild {

#define NOT_SURFACE SCHAR_MAX
#define KNOWN_NOT_SURFACE -SCHAR_MAX/2
#define KNOWN_SURFACE SCHAR_MAX/2

#define NO_SURFACE_TAG 0
#define NOT_BBOX -1
#define OPP_T_ID_BOUNDARY -1

    class MeshVertex {
    public:
        MeshVertex() {}

        Vector3 pos;

        inline Scalar &operator[](const int index) {
            assert(index >= 0 && index < 3);
            return pos[index];
        }

        inline Scalar operator[](const int index) const {
            assert(index >= 0 && index < 3);
            return pos[index];
        }

        std::vector<int> conn_tets;

        bool is_on_surface = false;
        bool is_on_boundary = false;
        bool is_on_cut = false;
        bool is_on_bbox = false;

        bool is_removed = false;
        bool is_freezed = false;

        Scalar sizing_scalar = 1;
    };

    class MeshTet {
    public:
        Vector4i indices;

        MeshTet() {}
        MeshTet(const Vector4i &idx) : indices(idx) {}
        MeshTet(int v0, int v1, int v2, int v3) : indices(v0, v1, v2, v3) {}

        inline int &operator[](const int index) {
            assert(index >= 0 && index < 4);
            return indices[index];
        }

        inline int operator[](const int index) const {
            assert(index >= 0 && index < 4);
            return indices[index];
        }

        inline int find(int ele) const {
            for (int j = 0; j < 4; j++)
                if (indices[j] == ele)
                    return j;
            return -1;
        }

        std::array<char, 4> is_surface_fs = {{NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE}};
        std::array<char, 4> is_bbox_fs = {{NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX}};
        std::array<char, 4> surface_tags = {{0, 0, 0, 0}};

        Scalar quality = 0;
        Scalar scalar = 0;
        bool is_removed = false;
        bool is_outside = false;
    };

    // Vertices and tets are both kept in place and marked, never erased, so the two lists below
    // are walked the same way.
    template<typename T>
    inline int count_live(const std::vector<T> &v) {
        int cnt = 0;
        for (const auto &e:v) {
            if (!e.is_removed)
                cnt++;
        }
        return cnt;
    }

    // The first slot a removal left free, or the end when there is none.
    template<typename T>
    inline int first_removed(const std::vector<T> &v) {
        for (int i = 0; i < int(v.size()); i++) {
            if (v[i].is_removed)
                return i;
        }
        return int(v.size());
    }

    class Mesh {
    public:
        std::vector<MeshVertex> tet_vertices;
        std::vector<MeshTet> tets;

        Parameters params;

        int t_empty_start = 0;
        int v_empty_start = 0;
        bool is_closed = true;

        bool is_input_all_inserted = false;
        bool is_coarsening = false;

        inline int get_v_num() const { return count_live(tet_vertices); }
        inline int get_t_num() const { return count_live(tets); }

        inline void reset_t_empty_start() { t_empty_start = first_removed(tets); }
        inline void reset_v_empty_start() { v_empty_start = first_removed(tet_vertices); }

        // The largest and the mean tet quality. Both in one pass, because the optimization loop
        // asks for the pair several times per iteration and the tet array is the big one.
        inline void get_max_avg_energy(Scalar &max_energy, Scalar &avg_energy) const {
            max_energy = 0;
            avg_energy = 0;
            int cnt = 0;
            for (auto &t: tets) {
                if (t.is_removed)
                    continue;
                if (t.quality > max_energy)
                    max_energy = t.quality;
                avg_energy += t.quality;
                cnt++;
            }
            avg_energy /= cnt;
        }

        inline Scalar get_max_energy() const {
            Scalar max_energy, _;
            get_max_avg_energy(max_energy, _);
            return max_energy;
        }

        inline Scalar get_avg_energy() const {
            Scalar _, avg_energy;
            get_max_avg_energy(_, avg_energy);
            return avg_energy;
        }
    };
}

