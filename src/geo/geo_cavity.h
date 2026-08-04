// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/delaunay/cavity.h
// Copied rather than reimplemented: Delaunay3d cavity bookkeeping, unchanged

#pragma once

#include "geo_basic.h"
#include <string.h>

namespace floatTetWild {
namespace geo {

    /**
     * \brief Represents the set of tetrahedra on the boundary
     *  of the cavity in a 3D Delaunay triangulation.
     */
    class Cavity {

    public:

        /**
         * \brief Type used for local indices.
         */
        typedef Numeric::uint8 local_index_t;

        /**
         * \brief Cavity constructor.
         */
        Cavity() {
            clear();
        }

        /**
         * \brief Clears this cavity.
         */
        void clear() {
            nb_f_ = 0;
            OK_ = true;
            ::memset(h2t_, END_OF_LIST, sizeof(h2t_));
        }

        /**
         * \brief Tests whether this Cavity is valid.
         * \retval true if this Cavity is valid.
         * \retval false otherwise. A Cavity is not valid
         *  when there was overflow.
         */
        bool OK() const {
            return OK_;
        }

        /**
         * \brief Inserts a new boundary facet in the structure.
         * \param[in] tglobal global tetrahedron index
         * \param[in] boundary_f index of the facet that is on the boundary
         * \param[in] v0 , v1 , v2 the three vertices of the facet that
         *  is on the boundary
         */
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

        /**
         * \brief Gets the number of facets.
         * \return the number of facets.
         */
        index_t nb_facets() const {
            return nb_f_;
        }

        /**
         * \brief Gets the tetrahedron associated with a facet.
         * \param[in] f the facet
         * \return the tetrahedron associated with \p f.
         */
        index_t facet_tet(index_t f) const {
            geo_debug_assert(f < nb_facets());
            return tglobal_[f];
        }

        /**
         * \brief Sets the tetrahedron associated with a facet.
         * \param[in] f the facet.
         * \param[in] t the tetrahedron to be associated with \p f.
         */
        void set_facet_tet(index_t f, index_t t) {
            geo_debug_assert(f < nb_facets());
            tglobal_[f] = t;
        }

        /**
         * \brief Gets the local tetrahedron facet that corresponds
         *  to a facet.
         * \param[in] f the facet.
         * \return the local index of the tetrahedron facet associated
         *  with \p f, in 0..3
         */
        index_t facet_facet(index_t f) const {
            geo_debug_assert(f < nb_facets());
            return boundary_f_[f];
        }

        /**
         * \brief Gets the vertex of a facet.
         * \param[in] f a facet.
         * \param[in] lv local index of the vertex, in 0..2.
         * \return the global vertex index.
         */
        index_t facet_vertex(index_t f, index_t lv) const {
            geo_debug_assert(f < nb_facets());
            geo_debug_assert(lv < 3);
            return f2v_[f][lv];
        }

        /**
         * \brief Gets the neighbors of a facet.
         * \param[in] f a facet
         * \param[out] t0 , t1 , t2 the global tetrahedron
         *  indices that correspond to the neighbors of \p f.
         */
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

        /**
         * \brief Computes the hash code associated with an oriented
         *  edge.
         * \param[in] v1 , v2 the global indices of the two vertices
         * \return the hash code, in 0 .. MAX_H -1
         */
        index_t hash(index_t v1, index_t v2) const {
            return (
		((index_t(v1+1) * 73856093) ^
		 (index_t(v2+1) * 83492791)) % MAX_H
	    );
        }

        /**
         * \brief Sets the local facet associated with an oriented
         *  edge.
         * \param[in] v1 , v2 the global indices of the two vertices
         * \param[in] f the local face index.
         */
        void set_vv2t(
            index_t v1, index_t v2, local_index_t f
        ) {
            index_t h = hash(v1,v2);
            index_t cur = h;
            do {
                if(h2t_[cur] == END_OF_LIST) {
                    h2t_[cur] = f;
                    h2v_[cur] = (Numeric::uint64(v1+1) << 32) |
                        Numeric::uint64(v2+1);
                    return;
                }
                cur = (cur+1)%MAX_H;
            } while(cur != h);
            OK_ = false;
        }

        /**
         * \brief gets the local facet associated with an oriented
         *  edge.
         * \param[in] v1 , v2 the global indices of the two vertices
         * \return the local facet index.
         */
        local_index_t get_vv2t(index_t v1, index_t v2) const {
            Numeric::uint64 K = (Numeric::uint64(v1+1) << 32) |
                Numeric::uint64(v2+1);
            index_t h = hash(v1,v2);
            index_t cur = h;
            do {
                if(h2v_[cur] == K) {
                    return h2t_[cur];
                }
                cur = (cur+1)%MAX_H;
            } while(cur != h);
            geo_assert_not_reached;
        }

        /** \brief Hash index to local facet id. */
        local_index_t  h2t_[MAX_H];

        /** \brief Hash index to global vertex id. */
        Numeric::uint64 h2v_[MAX_H];

        /** \brief Number of facets. */
        index_t nb_f_;

        /** \brief Local facet index to tetrahedra index. */
        index_t tglobal_[MAX_F];

        /** \brief Local facet index to facet on border index. */
        index_t boundary_f_[MAX_F];

        /** \brief Local facet index to three global vertex indices. */
        index_t f2v_[MAX_F][3];

        /**
         * \brief True if the structure is correct, false
         *  otherwise, if capacity was exceeded.
         */
        bool OK_;
    };

    }

}
