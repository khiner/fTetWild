// A stand-in for the slice of geogram/delaunay/delaunay.h that Delaunay3d needs.
//
// Reimplemented rather than copied: this base class only holds state and hands out accessors, so
// it decides nothing. geogram's version is a Counted, is built through a Factory, and carries a
// PackedArrays neighbour service, constraint and region support and a nearest_vertex() search --
// none of which fTetWild reaches, and all of which would pull in another four headers. The
// declarations Delaunay3d does use are unchanged.

#pragma once

#include "geo_basic.h"

namespace floatTetWild {
namespace geo {

    /**
     * \brief Abstract interface for a Delaunay triangulation in arbitrary dimension.
     */
    class Delaunay {
    public:
        virtual ~Delaunay() {
        }

        coord_index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Number of vertices in a cell: dimension + 1.
         */
        index_t cell_size() const {
            return cell_size_;
        }

        /**
         * \brief Sets the vertices and recomputes the cells.
         * \param[in] nb_vertices number of vertices
         * \param[in] vertices a pointer to the coordinates, dimension() doubles per vertex
         */
        virtual void set_vertices(index_t nb_vertices, const double* vertices) {
            nb_vertices_ = nb_vertices;
            vertices_ = vertices;
        }

        index_t nb_vertices() const {
            return nb_vertices_;
        }

        const double* vertex_ptr(index_t i) const {
            geo_debug_assert(i < nb_vertices());
            return vertices_ + vertex_stride_ * i;
        }

        index_t nb_cells() const {
            return nb_cells_;
        }

        /**
         * \brief Number of finite cells.
         * \details Only meaningful when keeps_infinite() is set; then the finite cells are
         *  0 .. nb_finite_cells() - 1.
         */
        index_t nb_finite_cells() const {
            geo_debug_assert(keep_infinite_);
            return nb_finite_cells_;
        }

        /**
         * \brief The cell-to-vertex array, cell_size() indices per cell.
         */
        const index_t* cell_to_v() const {
            return cell_to_v_;
        }

        /**
         * \brief The cell-to-cell array, cell_size() indices per cell.
         */
        const index_t* cell_to_cell() const {
            return cell_to_cell_;
        }

        index_t cell_vertex(index_t c, index_t lv) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lv < cell_size());
            return cell_to_v_[c * cell_v_stride_ + lv];
        }

        index_t cell_adjacent(index_t c, index_t lf) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lf < cell_size());
            return cell_to_cell_[c * cell_neigh_stride_ + lf];
        }

        /**
         * \brief Specifies whether the vertices should be reordered with BRIO before insertion.
         * \details Reordering is on by default and makes point location much faster. Turning it
         *  off changes the insertion order, and with it the triangulation of cospherical input.
         */
        void set_reorder(bool x) {
            do_reorder_ = x;
        }



    protected:
        explicit Delaunay(coord_index_t dimension) {
            set_dimension(dimension);
        }

        /**
         * \brief Publishes the computed triangulation.
         */
        void set_arrays(
            index_t nb_cells, const index_t* cell_to_v, const index_t* cell_to_cell
        ) {
            nb_cells_ = nb_cells;
            cell_to_v_ = cell_to_v;
            cell_to_cell_ = cell_to_cell;
        }

        /**
         * \brief Sets the dimension and everything derived from it: the number of doubles
         *  between two consecutive vertices, the number of vertices in a cell, and the strides
         *  of the cell vertex and cell adjacency arrays.
         */
        void set_dimension(coord_index_t dim) {
            dimension_ = dim;
            vertex_stride_ = dim;
            cell_size_ = index_t(dim) + 1;
            cell_v_stride_ = cell_size_;
            cell_neigh_stride_ = cell_size_;
        }

        coord_index_t dimension_ = 3;
        index_t vertex_stride_ = 3;
        index_t cell_size_ = 4;
        index_t cell_v_stride_ = 4;
        index_t cell_neigh_stride_ = 4;

        const double* vertices_ = nullptr;
        index_t nb_vertices_ = 0;
        index_t nb_cells_ = 0;
        const index_t* cell_to_v_ = nullptr;
        const index_t* cell_to_cell_ = nullptr;

        /** \brief If true, uses BRIO reordering. */
        bool do_reorder_ = true;

        /** \brief If true, the infinite vertex and infinite simplices are kept. */
        bool keep_infinite_ = false;

        /**
         * \brief If keep_infinite_ is true, then finite cells are 0 .. nb_finite_cells_ - 1 and
         *  infinite cells are nb_finite_cells_ .. nb_cells_.
         */
        index_t nb_finite_cells_ = 0;
    };
}
}
