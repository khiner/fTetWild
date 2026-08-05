// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/delaunay/cavity.h
// Copied rather than reimplemented: Delaunay3d cavity bookkeeping, unchanged

#pragma once

#include "geo_basic.h"
#include <string.h>

namespace floatTetWild {
namespace geo {

    // The facets on the boundary of the cavity left by removing the tetrahedra in conflict with a
    // point, collected while that conflict zone is traversed. Fixed capacity, in open-addressed
    // arrays: overflow is not an error, it sets OK() false and Delaunay3d falls back to walking
    // the border instead.
    class Cavity {

    public:

        typedef uint8_t local_index_t;

        Cavity() {
            clear();
        }

        void clear() {
            nb_f_ = 0;
            OK_ = true;
            ::memset(h2t_, END_OF_LIST, sizeof(h2t_));
        }

        // False once capacity has been exceeded, which makes the whole collection unusable.
        bool OK() const {
            return OK_;
        }

        // Records that facet \p boundary_f of tetrahedron \p tglobal, with vertices \p v0, \p v1
        // and \p v2, is on the boundary.
        void new_facet(
            index_t tglobal, index_t boundary_f,
            index_t v0, index_t v1, index_t v2
        ) {
            if(!OK_) {
                return;
            }

            geo_debug_assert(v0 != v1);
            geo_debug_assert(v1 != v2);
            geo_debug_assert(v2 != v0);

            local_index_t new_t = local_index_t(nb_f_);

            if(nb_f_ == MAX_F) {
                OK_ = false;
                return;
            }

            set_vv2t(v0, v1, new_t);
            set_vv2t(v1, v2, new_t);
            set_vv2t(v2, v0, new_t);

            if(!OK_) {
                return;
            }

            ++nb_f_;
            tglobal_[new_t] = tglobal;
            boundary_f_[new_t] = boundary_f;
            f2v_[new_t][0] = v0;
            f2v_[new_t][1] = v1;
            f2v_[new_t][2] = v2;
        }

        index_t nb_facets() const {
            return nb_f_;
        }

        index_t facet_tet(index_t f) const {
            geo_debug_assert(f < nb_facets());
            return tglobal_[f];
        }

        void set_facet_tet(index_t f, index_t t) {
            geo_debug_assert(f < nb_facets());
            tglobal_[f] = t;
        }

        // The local index, in 0..3, of the tetrahedron facet that \p f corresponds to.
        index_t facet_facet(index_t f) const {
            geo_debug_assert(f < nb_facets());
            return boundary_f_[f];
        }

        // The global index of vertex \p lv, in 0..2, of facet \p f.
        index_t facet_vertex(index_t f, index_t lv) const {
            geo_debug_assert(f < nb_facets());
            geo_debug_assert(lv < 3);
            return f2v_[f][lv];
        }

        // The tetrahedra of the three facets that share an edge with \p f, found by looking up
        // each of its edges reversed.
        void get_facet_neighbor_tets(
            index_t f, index_t& t0, index_t& t1, index_t& t2
        ) const {
            index_t v0 = f2v_[f][0];
            index_t v1 = f2v_[f][1];
            index_t v2 = f2v_[f][2];
            t0 = tglobal_[get_vv2t(v2,v1)];
            t1 = tglobal_[get_vv2t(v0,v2)];
            t2 = tglobal_[get_vv2t(v1,v0)];
        }

    private:
        static constexpr index_t        MAX_H = 1033;
        static constexpr local_index_t  END_OF_LIST = 255;
        static constexpr index_t        MAX_F = 128;

        index_t hash(index_t v1, index_t v2) const {
            return (
		((index_t(v1+1) * 73856093) ^
		 (index_t(v2+1) * 83492791)) % MAX_H
	    );
        }

        // Associates local facet \p f with the oriented edge \p v1, \p v2, by linear probing.
        void set_vv2t(
            index_t v1, index_t v2, local_index_t f
        ) {
            index_t h = hash(v1,v2);
            index_t cur = h;
            do {
                if(h2t_[cur] == END_OF_LIST) {
                    h2t_[cur] = f;
                    h2v_[cur] = (uint64_t(v1+1) << 32) |
                        uint64_t(v2+1);
                    return;
                }
                cur = (cur+1)%MAX_H;
            } while(cur != h);
            OK_ = false;
        }

        local_index_t get_vv2t(index_t v1, index_t v2) const {
            uint64_t K = (uint64_t(v1+1) << 32) |
                uint64_t(v2+1);
            index_t h = hash(v1,v2);
            index_t cur = h;
            do {
                if(h2v_[cur] == K) {
                    return h2t_[cur];
                }
                cur = (cur+1)%MAX_H;
            } while(cur != h);
            geo_assert(false);
        }

        // The hash table: an oriented edge to the local facet it belongs to, and the edge itself
        // so a probe can tell a match from a collision.
        local_index_t  h2t_[MAX_H];
        uint64_t h2v_[MAX_H];

        index_t nb_f_;

        // Per local facet: its tetrahedron, which of that tetrahedron's four facets it is, and
        // its three global vertex indices.
        index_t tglobal_[MAX_F];
        index_t boundary_f_[MAX_F];
        index_t f2v_[MAX_F][3];

        bool OK_;
    };

    }

}
