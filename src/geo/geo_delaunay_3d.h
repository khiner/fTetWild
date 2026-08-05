// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/delaunay/delaunay_3d.h
// Copied rather than reimplemented: create_first_tetrahedron scans raw indices and SOS breaks ties by point address, so insertion order decides the mesh

#pragma once

#include "geo_cavity.h"
#include "geo_predicates.h"

#include <stack>

namespace floatTetWild {
namespace geo {

    // Incremental Delaunay triangulation in 3d, by Bowyer-Watson: each new point's conflict zone
    // is traversed from the inside, as CGAL does, since that crosses fewer tetrahedra than
    // traversing from the outside. Vertex deletion is not supported, and neither is degenerate
    // input with all points coplanar or all colinear.
    //
    // Points are inserted in BRIO order, which speeds up point location dramatically; see
    // compute_BRIO_order(). locate() itself randomizes which neighbour it walks to next, and
    // starts from the result of an inexact walk ("structural filtering"), both of which are what
    // CGAL and tetgen do.
    //
    // References: Boissonnat, Devillers, Teillaud and Yvinec, "Triangulations in CGAL", Proc. 16th
    // Annu. ACM Sympos. Comput. Geom. 11-18, 2000. Si, "Constrained Delaunay tetrahedral mesh
    // generation and refinement", Finite elements in Analysis and Design 46(1-2):33-46, 2010.
    // Bowyer, "Computing Dirichlet tessellations", Comput. J. 24(2):162-166, 1981. Watson,
    // "Computing the n-dimensional Delaunay tessellation with application to Voronoi polytopes",
    // Comput. J. 24(2):167-172, 1981. Amenta, Choi and Rote, "Incremental constructions con brio",
    // ACM Sympos. Comput. Geom. 2003. Delage and Devillers, "Spatial Sorting", CGAL User and
    // Reference Manual, 3.9 edition, 2011. Devillers, Pion and Teillaud, "Walking in a
    // triangulation", 17th Annu. Sympos. Comput. Geom. 106-114. Funke, Mehlhorn and Naher,
    // "Structural filtering, a paradigm for efficient and exact geometric programs", Comput.
    // Geom., 1999.
    class Delaunay3d {
    public:
        Delaunay3d();

        // Sets the vertices, three doubles each, and computes the cells.
        void set_vertices(index_t nb_vertices, const double* vertices);

        index_t nb_vertices() const {
            return nb_vertices_;
        }

        const double* vertex_ptr(index_t i) const {
            geo_debug_assert(i < nb_vertices());
            return vertices_ + 3 * i;
        }

        // Also the number of slots the two stores below are divided into while the triangulation
        // is being built, where some of them are free or virtual.
        index_t nb_cells() const {
            return cell_to_v_store_.size() / 4;
        }

        // The cell-to-vertex array, four indices per cell.
        const index_t* cell_to_v() const {
            return cell_to_v_store_.data();
        }

    private:
        // No hint given for point location.
        static constexpr index_t NO_TETRAHEDRON = NO_INDEX;

        // Four non-coplanar points to start the incremental construction from, or false if the
        // points are all coplanar.
        bool create_first_tetrahedron(
            index_t& iv0, index_t& iv1, index_t& iv2, index_t& iv3
        );

        // The tetrahedron that contains \p p, or one incident to it when \p p is on a face, edge
        // or vertex. If \p orient is non-null it receives the orientation of \p p with respect to
        // the four facets. Outside the convex hull of the points inserted so far, the answer is a
        // virtual tetrahedron, or NO_TETRAHEDRON once those have been removed.
        index_t locate(
            const double* p, index_t hint = NO_TETRAHEDRON,
            Sign* orient = nullptr
        ) const;

        // The same walk with inexact predicates, stopping after \p max_iter tetrahedra. Its
        // result is a hint for locate(), which is what makes locate() faster than walking exactly
        // the whole way ("structural filtering").
        index_t locate_inexact(
            const double* p, index_t hint, index_t max_iter
        ) const;

        // Inserts point \p v, and returns one of the tetrahedra incident to it.
        index_t insert(index_t v, index_t hint = NO_TETRAHEDRON);

        // The tetrahedra in conflict with point \p v, given \p t containing it and the \p orient
        // locate() returned alongside. \p t_bndry and \p f_bndry come back as a tetrahedron
        // adjacent to the boundary of the conflict zone and the facet it is adjacent along. The
        // conflict tetrahedra are chained from \p first to \p last through tet_next(). The chain
        // is empty when \p v already exists in the triangulation.
        void find_conflict_zone(
            index_t v,
            index_t t, const Sign* orient,
            index_t& t_bndry, index_t& f_bndry,
            index_t& first, index_t& last
        );

        // Propagates find_conflict_zone() outwards from \p t, which must already be marked as in
        // conflict.
        void find_conflict_zone_iterative(
            const double* p, index_t t,
            index_t& t_bndry, index_t& f_bndry,
            index_t& first, index_t& last
        );

        // Fills the conflict zone with a star of tetrahedra around \p v, one per facet on its
        // border, and returns one of them. Used when the Cavity collected while traversing the
        // conflict zone did not overflow its arrays.
        index_t stellate_cavity(index_t v);

        // The same star, built by walking the border instead: starting from facet \p f_bndry of
        // \p t_bndry and recursing until the zone is filled. \p prev_f is the facet of \p t_bndry
        // that it was reached from, or NO_INDEX for the first one.
        index_t stellate_conflict_zone_iterative(
            index_t v,
            index_t t_bndry, index_t f_bndry,
            index_t prev_f=NO_INDEX
        );

        // Used by stellate_conflict_zone_iterative(): given \p t1 on the border of the conflict
        // zone along facet \p t1fborder, finds the tetrahedron \p t2 across facet \p t1ft2 that
        // is also on the border and shares an edge with both facets. \p t2fborder and \p t2ft1
        // come back as t2's facets on the border and towards t1. Returns true if \p t2 is newly
        // created, false if it is an old tetrahedron still in conflict.
        //
        // Long for an inline function, but geogram measured a modest gain from keeping it one.
        bool get_neighbor_along_conflict_zone_border(
            index_t t1,
            index_t t1fborder,
            index_t t1ft2,
            index_t& t2,
            index_t& t2fborder,
            index_t& t2ft1
        ) const {

            //   Find two vertices that are both on facets new_f and f1
            //  (the edge around which we are turning)
            //  This uses duality as follows:
            //  Primal form (not used here):
            //    halfedge_facet_[v1][v2] returns a facet that is incident
            //    to both v1 and v2.
            //  Dual form (used here):
            //    halfedge_facet_[f1][f2] returns a vertex that both
            //    f1 and f2 are incident to.
            index_t ev1 = tet_vertex(t1, halfedge_facet_[t1ft2][t1fborder]);
            index_t ev2 = tet_vertex(t1, halfedge_facet_[t1fborder][t1ft2]);

            //   Turn around edge [ev1,ev2] inside the conflict zone
            // until we reach again the boundary of the conflict zone.
            // Traversing inside the conflict zone is faster (as compared
            // to outside) since it traverses a smaller number of tets.
            index_t cur_t = t1;
            index_t cur_f = t1ft2;
            index_t next_t = tet_adjacent(cur_t,cur_f);
            while(tet_is_in_list(next_t)) {
                geo_debug_assert(next_t != t1);
                cur_t = next_t;
                cur_f = get_facet_by_halfedge(cur_t,ev1,ev2);
                next_t = tet_adjacent(cur_t, cur_f);
            }

            //  At this point, cur_t is in conflict zone and
            // next_t is outside the conflict zone.
            index_t f12,f21;
            get_facets_by_halfedge(next_t, ev1, ev2, f12, f21);
            t2 = tet_adjacent(next_t,f21);
            index_t v_neigh_opposite = tet_vertex(next_t,f12);
            t2ft1 = find_tet_vertex(t2, v_neigh_opposite);
            t2fborder = cur_f;

            //  Test whether the found neighboring tet was created
            //  (then return true) or is an old tet in conflict
            //  (then return false).
            return(t2 != cur_t);
        }

        /****** Combinatorics - new and delete ***************************/

        //   Tetrahedra are chained into lists, which is how the free list that recycles deleted
        // tetrahedra, the conflict region and the list of newly created tetrahedra are all
        // managed. A tetrahedron that is not in a list can instead be marked, using the same
        // space: the index of the point being inserted serves as the stamp. So a tetrahedron is
        // in exactly one of three states, read off cell_next_[t]:
        //   - in a list                      (NOT_IN_LIST_BIT clear)
        //   - not in a list and marked       (== cur_stamp_, which has NOT_IN_LIST_BIT set)
        //   - not in a list and not marked   (anything else with NOT_IN_LIST_BIT set)

        static constexpr index_t NOT_IN_LIST = ~index_t(0);

        static constexpr index_t NOT_IN_LIST_BIT =
	    index_t(1) << (sizeof(index_t)*8-1) ;

        static constexpr index_t END_OF_LIST = ~NOT_IN_LIST_BIT;

        bool tet_is_in_list(index_t t) const {
            geo_debug_assert(t < nb_cells());
            return (cell_next_[t] & NOT_IN_LIST_BIT) == 0;
        }

        // \pre tet_is_in_list(t)
        index_t tet_next(index_t t) const {
            geo_debug_assert(t < nb_cells());
            geo_debug_assert(tet_is_in_list(t));
            return cell_next_[t];
        }

        // Prepends \p t to the list running from \p first to \p last, both END_OF_LIST when it is
        // empty.
        void add_tet_to_list(index_t t, index_t& first, index_t& last) {
            geo_debug_assert(t < nb_cells());
            geo_debug_assert(!tet_is_in_list(t));
            if(last == END_OF_LIST) {
                geo_debug_assert(first == END_OF_LIST);
                first = last = t;
                cell_next_[t] = END_OF_LIST;
            } else {
                cell_next_[t] = first;
                first = t;
            }
        }

        bool tet_is_marked(index_t t) const {
            return cell_next_[t] == cur_stamp_;
        }

        void mark_tet(index_t t) {
            cell_next_[t] = cur_stamp_;
        }

        // The first vertex of a virtual tetrahedron, whose three other vertices are then a
        // triangle on the convex hull of the points.
        static constexpr index_t VERTEX_AT_INFINITY = NO_INDEX;

        bool tet_is_finite(index_t t) const {
            return
                cell_to_v_store_[4 * t]     != NO_INDEX &&
                cell_to_v_store_[4 * t + 1] != NO_INDEX &&
                cell_to_v_store_[4 * t + 2] != NO_INDEX &&
                cell_to_v_store_[4 * t + 3] != NO_INDEX;
        }

        // A real tetrahedron is incident to four user-specified vertices, so neither recycled nor
        // incident to the vertex at infinity.
        bool tet_is_real(index_t t) const {
            return !tet_is_free(t) && tet_is_finite(t);
        }

        bool tet_is_virtual(index_t t) const {
            return
                !tet_is_free(t) && (
                    cell_to_v_store_[4 * t] == VERTEX_AT_INFINITY ||
                    cell_to_v_store_[4 * t + 1] == VERTEX_AT_INFINITY ||
                    cell_to_v_store_[4 * t + 2] == VERTEX_AT_INFINITY ||
                    cell_to_v_store_[4 * t + 3] == VERTEX_AT_INFINITY) ;
        }

        // Deleted tetrahedra are chained in the free list, so being in a list is being free --
        // outside find_conflict_zone(), which chains the conflict region through the same field.
        bool tet_is_free(index_t t) const {
            return tet_is_in_list(t);
        }

        // Recycles a tetrahedron from the free list, or grows the stores by one. Adjacencies are
        // cleared, vertices are not.
        index_t new_tetrahedron(
            index_t v1, index_t v2,
            index_t v3, index_t v4
        ) {
            index_t result;
            if(first_free_ == END_OF_LIST) {
                cell_to_v_store_.resize(
		    cell_to_v_store_.size() + 4, NO_INDEX
		);
                cell_to_cell_store_.resize(
		    cell_to_cell_store_.size() + 4, NO_INDEX
		);
                // index_t(NOT_IN_LIST) is necessary, else with
                // NOT_IN_LIST alone the compiler tries to generate a
                // reference to NOT_IN_LIST resulting in a link error.
                cell_next_.push_back(index_t(NOT_IN_LIST));
                result = nb_cells() - 1;
            } else {
                result = first_free_;
                first_free_ = tet_next(first_free_);
                cell_next_[result] = NOT_IN_LIST;
            }

            cell_to_cell_store_[4 * result] = NO_INDEX;
            cell_to_cell_store_[4 * result + 1] = NO_INDEX;
            cell_to_cell_store_[4 * result + 2] = NO_INDEX;
            cell_to_cell_store_[4 * result + 3] = NO_INDEX;

            cell_to_v_store_[4 * result] = v1;
            cell_to_v_store_[4 * result + 1] = v2;
            cell_to_v_store_[4 * result + 2] = v3;
            cell_to_v_store_[4 * result + 3] = v4;
            return result;
        }

        /********* Combinatorics ******************************************/

        // The global index of the \p lv%th vertex of tetrahedron \p t, or NO_INDEX if that vertex
        // is at infinity.
        index_t tet_vertex(index_t t, index_t lv) const {
            geo_debug_assert(t < nb_cells());
            geo_debug_assert(lv < 4);
            return cell_to_v_store_[4 * t + lv];
        }

        // lv such that tet_vertex(t,lv)==v.
        // \pre \p t is incident to \p v
        index_t find_tet_vertex(index_t t, index_t v) const {
            geo_debug_assert(t < nb_cells());
            return find_4(&(cell_to_v_store_[4 * t]),v);
        }

        void set_tet_vertex(index_t t, index_t lv, index_t v) {
            geo_debug_assert(t < nb_cells());
            geo_debug_assert(lv < 4);
            cell_to_v_store_[4 * t + lv] = v;
        }

        // The tetrahedron adjacent to \p t across its facet \p lf.
        index_t tet_adjacent(index_t t, index_t lf) const {
            geo_debug_assert(t < nb_cells());
            geo_debug_assert(lf < 4);
            return cell_to_cell_store_[4 * t + lf];
        }

        void set_tet_adjacent(index_t t1, index_t lf1, index_t t2) {
            geo_debug_assert(t1 < nb_cells());
            geo_debug_assert(t2 < nb_cells());
            geo_debug_assert(lf1 < 4);
            cell_to_cell_store_[4 * t1 + lf1] = t2;
        }

        // f such that tet_adjacent(t1,f)==t2.
        // \pre \p t1 and \p t2 are adjacent
        index_t find_tet_adjacent(index_t t1, index_t t2) const {
            geo_debug_assert(t1 < nb_cells());
            geo_debug_assert(t2 < nb_cells());
            geo_debug_assert(t1 != t2);

            index_t result = find_4(&(cell_to_cell_store_[4 * t1]),t2);

            // Sanity check: make sure that t1 is adjacent to t2
            // only once!
            geo_debug_assert(tet_adjacent(t1,(result+1)%4) != t2);
            geo_debug_assert(tet_adjacent(t1,(result+2)%4) != t2);
            geo_debug_assert(tet_adjacent(t1,(result+3)%4) != t2);
            return result;
        }

        /****** Combinatorics - traversals ************************/

        // The local index of the facet of \p t incident to the oriented edge \p v1, \p v2.
        index_t get_facet_by_halfedge(index_t t, index_t v1, index_t v2) const {
            geo_debug_assert(t < nb_cells());
            geo_debug_assert(v1 != v2);
            const index_t* T = &(cell_to_v_store_[4 * t]);
            index_t lv1 = find_4(T,v1);
            index_t lv2 = find_4(T,v2);
            geo_debug_assert(lv1 != lv2);
            return halfedge_facet_[lv1][lv2];
        }

        // The same for both orientations at once: \p f12 for [v1,v2] and \p f21 for [v2,v1].
        void get_facets_by_halfedge(
            index_t t, index_t v1, index_t v2,
            index_t& f12, index_t& f21
        ) const {
            geo_debug_assert(t < nb_cells());
            geo_debug_assert(v1 != v2);

            //   Find local index of v1 and v2 in tetrahedron t
            // The following expression is 10% faster than using
            // if() statements (multiply by boolean result of test).
            // Thank to Laurent Alonso for this idea.
            const index_t* T = &(cell_to_v_store_[4 * t]);

            index_t lv1 = index_t(
                (T[1] == v1) | ((T[2] == v1) * 2) | ((T[3] == v1) * 3)
	    );

            index_t lv2 = index_t(
                (T[1] == v2) | ((T[2] == v2) * 2) | ((T[3] == v2) * 3)
	    );

            geo_debug_assert(lv1 != 0 || T[0] == v1);
            geo_debug_assert(lv2 != 0 || T[0] == v2);
            geo_debug_assert(lv1 != lv2);

            f12 = halfedge_facet_[lv1][lv2];
            f21 = halfedge_facet_[lv2][lv1];
        }


        /****** Predicates **********************************************/

        // Whether \p t is in conflict with the 3d point \p p: for a real tetrahedron, whether its
        // circumscribed sphere contains \p p; for a virtual one, whether the tetrahedron formed
        // by its real face and \p p is oriented positively.
        bool tet_is_conflict(index_t t, const double* p) const {

            // Lookup tetrahedron vertices
            const double* pv[4];
            for(index_t i=0; i<4; ++i) {
                index_t v = tet_vertex(t,i);
                pv[i] = (v == NO_INDEX) ? nullptr : vertex_ptr(v);
            }

            // Check for virtual tetrahedra (then in_sphere()
            // is replaced with orient3d())
            for(index_t lf = 0; lf < 4; ++lf) {

                if(pv[lf] == nullptr) {

                    // Facet of a virtual tetrahedron opposite to
                    // infinite vertex corresponds to
                    // the triangle on the convex hull of the points.
                    // Orientation is obtained by replacing vertex lf
                    // with p.
                    pv[lf] = p;
                    Sign sign = PCK::orient_3d(pv[0],pv[1],pv[2],pv[3]);

                    if(sign > 0) {
                        return true;
                    }

                    if(sign < 0) {
                        return false;
                    }

                    // If sign is zero, we check the real tetrahedron
                    // adjacent to the facet on the convex hull.
                    geo_debug_assert(tet_adjacent(t, lf) != NO_INDEX);
                    index_t t2 = tet_adjacent(t, lf);
                    geo_debug_assert(!tet_is_virtual(t2));

                    //  If t2 is already chained in the conflict list,
                    // then it is conflict
                    if(tet_is_in_list(t2)) {
                        return true;
                    }

                    //  If t2 is marked, then it is not in conflict.
                    if(tet_is_marked(t2)) {
                        return false;
                    }

                    return tet_is_conflict(t2, p);
                }
            }

            //   If the tetrahedron is a finite one, it is in conflict
            // if its circumscribed sphere contains the point (this is
            // the standard case).

            return (PCK::in_sphere_3d_SOS(pv[0], pv[1], pv[2], pv[3], p) > 0);
        }

        // The index (0,1,2 or 3) of \p v in the four different integers at \p T, one of which is
        // \p v.
        static index_t find_4(const index_t* T, index_t v) {
            // The following expression is 10% faster than using
            // if() statements. This uses the C++ norm, that
            // ensures that the 'true' boolean value converted to
            // an int is always 1. With most compilers, this avoids
            // generating branching instructions.
            // Thank to Laurent Alonso for this idea.
            // Note: Laurent also has this version:
            //    (T[0] != v)+(T[2]==v)+2*(T[3]==v)
            // that avoids a *3 multiply, but it is not faster in
            // practice.
            index_t result = index_t(
                (T[1] == v) | ((T[2] == v) * 2) | ((T[3] == v) * 3)
            );
            // Sanity check, important if it was T[0], not explicitly
            // tested (detects input that does not meet the precondition).
            geo_debug_assert(T[result] == v);
            return result;
        }

        const double* vertices_ = nullptr;
        index_t nb_vertices_ = 0;

        // The triangulation, four vertex indices per tetrahedron once set_vertices() has
        // compacted the free and virtual tetrahedra out of it.
        vector<index_t> cell_to_v_store_;
        vector<index_t> cell_to_cell_store_;
        vector<index_t> cell_next_;
        vector<index_t> reorder_;
        index_t cur_stamp_; // used for marking
        index_t first_free_;

        // tet_facet_vertex_[lf][lv] is the local vertex index, in 0..3, of the lv%th vertex of
        // local facet lf. The tetrahedron formed by a vertex lv and tet_facet_vertex_[lv][0..2]
        // has the same orientation as the original one, for any lv.
        static const index_t tet_facet_vertex_[4][3];

        // A local facet index by the local indices of the two extremities of a halfedge.
        static const index_t halfedge_facet_[4][4];

        // Used by the de-recursified find_conflict_zone_iterative().
        std::stack<index_t> S_;

        // One frame of the de-recursified stellate_conflict_zone_iterative(): the parameters it
        // was called with, then the locals it restores when a nested call returns. The narrow
        // types keep the frame small, because degenerate input can make the stack very deep.
        struct StellateFrame {
            // Parameters
            index_t t1;
            uint8_t t1fbord;
            uint8_t t1fprev;

            // Local variables
            index_t new_t;
            uint8_t t1ft2;
            uint8_t t2ft1;
        };

        std::vector<StellateFrame> S2_;

        Cavity cavity_;
    };

} }
