// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/delaunay/delaunay_3d.cpp and delaunay_3d.h, and cavity.h
// Copied rather than reimplemented: create_first_tetrahedron scans raw indices and SOS breaks ties by point address, so insertion order decides the mesh

#include "geo_delaunay_3d.h"
#include "geo_geometry.h"
#include "geo_mesh_reorder.h"
#include "geo_predicates.h"

#include <cstring>
#include <stack>

namespace floatTetWild {
namespace geo {

    namespace {

        // orient_3d from the plain floating point determinant, for the inexact walk.
        inline Sign orient_3d_inexact(
            const double* p0, const double* p1,
            const double* p2, const double* p3
        ) {
            double a11 = p1[0] - p0[0] ;
            double a12 = p1[1] - p0[1] ;
            double a13 = p1[2] - p0[2] ;

            double a21 = p2[0] - p0[0] ;
            double a22 = p2[1] - p0[1] ;
            double a23 = p2[2] - p0[2] ;

            double a31 = p3[0] - p0[0] ;
            double a32 = p3[1] - p0[1] ;
            double a33 = p3[2] - p0[2] ;

            double Delta =
                a11*det2x2(a22,a23,a32,a33)
                -a21*det2x2(a12,a13,a32,a33)
                +a31*det2x2(a12,a13,a22,a23);

            return geo_sgn(Delta);
        }

        bool points_are_identical_3d(const double* p1, const double* p2) {
            return
                (p1[0] == p2[0]) &&
                (p1[1] == p2[1]) &&
                (p1[2] == p2[2])
                ;
        }

        // Colinearity tested by four coplanarity tests against four points that are themselves
        // not coplanar.
        bool points_are_colinear_3d(
            const double* p1, const double* p2, const double* p3
        ) {
            static const double q000[3] = {0.0, 0.0, 0.0};
            static const double q001[3] = {0.0, 0.0, 1.0};
            static const double q010[3] = {0.0, 1.0, 0.0};
            static const double q100[3] = {1.0, 0.0, 0.0};
            return
                PCK::orient_3d(p1, p2, p3, q000) == ZERO &&
                PCK::orient_3d(p1, p2, p3, q001) == ZERO &&
                PCK::orient_3d(p1, p2, p3, q010) == ZERO &&
                PCK::orient_3d(p1, p2, p3, q100) == ZERO
                ;
        }

        // A local facet index by the local indices of the two extremities of a halfedge.
        const index_t halfedge_facet_[4][4] = {
            {4, 2, 3, 1},
            {3, 4, 0, 2},
            {1, 3, 4, 0},
            {2, 0, 1, 4}
        };

        // tet_facet_vertex_[lf][lv] is the local vertex index, in 0..3, of the lv%th vertex of
        // local facet lf. The tetrahedron formed by a vertex lv and tet_facet_vertex_[lv][0..2]
        // has the same orientation as the original one, for any lv.
        const index_t tet_facet_vertex_[4][3] = {
            {1, 2, 3},
            {0, 3, 2},
            {3, 0, 1},
            {1, 0, 2}
        };

        // The index (0,1,2 or 3) of \p v in the four different integers at \p T, one of which is
        // \p v. Branchless: the boolean tests convert to 0 or 1, and if none of T[1..3] matches
        // the result is 0.
        index_t find_4(const index_t* T, index_t v) {
            return index_t(
                (T[1] == v) | ((T[2] == v) * 2) | ((T[3] == v) * 3)
            );
        }

        // The facets on the boundary of the cavity left by removing the tetrahedra in conflict
        // with a point, collected while that conflict zone is traversed. Fixed capacity, in
        // open-addressed arrays: overflow is not an error, it sets OK false and Delaunay3d falls
        // back to walking the border instead. From geogram/delaunay/cavity.h.
        struct Cavity {
            typedef uint8_t local_index_t;

            static constexpr index_t        MAX_H = 1033;
            static constexpr local_index_t  END_OF_LIST = 255;
            static constexpr index_t        MAX_F = 128;

            Cavity() {
                clear();
            }

            void clear() {
                nb_f_ = 0;
                OK_ = true;
                ::memset(h2t_, END_OF_LIST, sizeof(h2t_));
            }

            // Records that facet \p boundary_f of tetrahedron \p tglobal, with vertices \p v0,
            // \p v1 and \p v2, is on the boundary.
            void new_facet(
                index_t tglobal, index_t boundary_f,
                index_t v0, index_t v1, index_t v2
            ) {
                if(!OK_) {
                    return;
                }

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

            // The tetrahedra of the three facets that share an edge with \p f, found by looking
            // up each of its edges reversed.
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

            // The hash table: an oriented edge to the local facet it belongs to, and the edge
            // itself so a probe can tell a match from a collision.
            local_index_t  h2t_[MAX_H];
            uint64_t h2v_[MAX_H];

            index_t nb_f_;

            // Per local facet: its tetrahedron, which of that tetrahedron's four facets it is,
            // and its three global vertex indices.
            index_t tglobal_[MAX_F];
            index_t boundary_f_[MAX_F];
            index_t f2v_[MAX_F][3];

            // False once capacity has been exceeded, which makes the whole collection unusable.
            bool OK_;
        };

        // Incremental Delaunay triangulation in 3d, by Bowyer-Watson: each new point's conflict
        // zone is traversed from the inside, as CGAL does, since that crosses fewer tetrahedra
        // than traversing from the outside. Vertex deletion is not supported, and neither is
        // degenerate input with all points coplanar or all colinear.
        //
        // Points are inserted in BRIO order, which speeds up point location dramatically; see
        // compute_BRIO_order(). locate() itself randomizes which neighbour it walks to next, and
        // starts from the result of an inexact walk ("structural filtering"), both of which are
        // what CGAL and tetgen do.
        //
        // geogram also had a weighted mode, entered by asking for dimension 4, where the vertices
        // are 4d and the combinatorics stay 3d, and a keep_infinite mode that published the
        // virtual tetrahedra incident to the vertex at infinity. Nothing here asks for either, so
        // this is plain 3d and publishes only the real tetrahedra.
        //
        // References: Boissonnat, Devillers, Teillaud and Yvinec, "Triangulations in CGAL", Proc.
        // 16th Annu. ACM Sympos. Comput. Geom. 11-18, 2000. Si, "Constrained Delaunay tetrahedral
        // mesh generation and refinement", Finite elements in Analysis and Design 46(1-2):33-46,
        // 2010. Bowyer, "Computing Dirichlet tessellations", Comput. J. 24(2):162-166, 1981.
        // Watson, "Computing the n-dimensional Delaunay tessellation with application to Voronoi
        // polytopes", Comput. J. 24(2):167-172, 1981. Amenta, Choi and Rote, "Incremental
        // constructions con brio", ACM Sympos. Comput. Geom. 2003. Delage and Devillers, "Spatial
        // Sorting", CGAL User and Reference Manual, 3.9 edition, 2011. Devillers, Pion and
        // Teillaud, "Walking in a triangulation", 17th Annu. Sympos. Comput. Geom. 106-114.
        // Funke, Mehlhorn and Naher, "Structural filtering, a paradigm for efficient and exact
        // geometric programs", Comput. Geom., 1999.
        struct Delaunay3d {

            /****** Combinatorics - new and delete ***************************/

            //   Tetrahedra are chained into lists, which is how the free list that recycles
            // deleted tetrahedra, the conflict region and the list of newly created tetrahedra
            // are all managed. A tetrahedron that is not in a list can instead be marked, using
            // the same space: the index of the point being inserted serves as the stamp. So a
            // tetrahedron is in exactly one of three states, read off cell_next_[t]:
            //   - in a list                      (NOT_IN_LIST_BIT clear)
            //   - not in a list and marked       (== cur_stamp_, which has NOT_IN_LIST_BIT set)
            //   - not in a list and not marked   (anything else with NOT_IN_LIST_BIT set)

            static constexpr index_t NOT_IN_LIST = ~index_t(0);

            static constexpr index_t NOT_IN_LIST_BIT =
                index_t(1) << (sizeof(index_t)*8-1) ;

            static constexpr index_t END_OF_LIST = ~NOT_IN_LIST_BIT;

            // The first vertex of a virtual tetrahedron, whose three other vertices are then a
            // triangle on the convex hull of the points.
            static constexpr index_t VERTEX_AT_INFINITY = NO_INDEX;

            bool tet_is_in_list(index_t t) const {
                return (cell_next_[t] & NOT_IN_LIST_BIT) == 0;
            }

            // \pre tet_is_in_list(t)
            index_t tet_next(index_t t) const {
                return cell_next_[t];
            }

            // Prepends \p t to the list running from \p first to \p last, both END_OF_LIST when
            // it is empty.
            void add_tet_to_list(index_t t, index_t& first, index_t& last) {
                if(last == END_OF_LIST) {
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

            //   Deleted tetrahedra are chained in the free list, so outside find_conflict_zone(),
            // which chains the conflict region through the same field, being in a list is being
            // free. A real tetrahedron is incident to four user-specified vertices, so neither
            // free nor incident to the vertex at infinity.
            bool tet_is_real(index_t t) const {
                return
                    !tet_is_in_list(t) &&
                    cell_to_v_store_[4 * t]     != NO_INDEX &&
                    cell_to_v_store_[4 * t + 1] != NO_INDEX &&
                    cell_to_v_store_[4 * t + 2] != NO_INDEX &&
                    cell_to_v_store_[4 * t + 3] != NO_INDEX;
            }

            bool tet_is_virtual(index_t t) const {
                return
                    !tet_is_in_list(t) && (
                        cell_to_v_store_[4 * t] == VERTEX_AT_INFINITY ||
                        cell_to_v_store_[4 * t + 1] == VERTEX_AT_INFINITY ||
                        cell_to_v_store_[4 * t + 2] == VERTEX_AT_INFINITY ||
                        cell_to_v_store_[4 * t + 3] == VERTEX_AT_INFINITY) ;
            }

            // Recycles a tetrahedron from the free list, or grows the stores by one. Adjacencies
            // are cleared, vertices are not.
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
                    // The cast keeps NOT_IN_LIST from being odr-used, which would need an
                    // out-of-class definition.
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

            // Also the number of slots the two stores are divided into while the triangulation is
            // being built, where some of them are free or virtual.
            index_t nb_cells() const {
                return cell_to_v_store_.size() / 4;
            }

            const double* vertex_ptr(index_t i) const {
                return vertices_ + 3 * i;
            }

            // The global index of the \p lv%th vertex of tetrahedron \p t, or NO_INDEX if that
            // vertex is at infinity.
            index_t tet_vertex(index_t t, index_t lv) const {
                return cell_to_v_store_[4 * t + lv];
            }

            // The tetrahedron adjacent to \p t across its facet \p lf.
            index_t tet_adjacent(index_t t, index_t lf) const {
                return cell_to_cell_store_[4 * t + lf];
            }

            void set_tet_adjacent(index_t t1, index_t lf1, index_t t2) {
                cell_to_cell_store_[4 * t1 + lf1] = t2;
            }

            // f such that tet_adjacent(t1,f)==t2.
            // \pre \p t1 and \p t2 are adjacent
            index_t find_tet_adjacent(index_t t1, index_t t2) const {
                return find_4(&(cell_to_cell_store_[4 * t1]), t2);
            }

            /****** Predicates **********************************************/

            // Whether \p t is in conflict with the 3d point \p p: for a real tetrahedron, whether
            // its circumscribed sphere contains \p p; for a virtual one, whether the tetrahedron
            // formed by its real face and \p p is oriented positively.
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
                        index_t t2 = tet_adjacent(t, lf);

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

            /****** Insertion ***********************************************/

            // A tetrahedron the walks below can start from: \p hint if it is usable, else a
            // random one, and in either case stepping off a virtual tetrahedron onto its real
            // neighbour.
            index_t walkable_hint(index_t hint) const {
                // If no hint specified, find a tetrahedron randomly. Both walks below draw from
                // the same process-wide stream, so which tetrahedra they start from depends on
                // the order they run in, and that is what decides the triangulation on
                // cospherical input.
                while(hint == NO_INDEX) {
                    hint = index_t(random_int32()) % nb_cells();
                    if(tet_is_in_list(hint)) {
                        hint = NO_INDEX;
                    }
                }

                //  Always start from a real tet. If the tet is virtual,
                // find its real neighbor (always opposite to the
                // infinite vertex)
                if(tet_is_virtual(hint)) {
                    for(index_t lf = 0; lf < 4; ++lf) {
                        if(tet_vertex(hint, lf) == VERTEX_AT_INFINITY) {
                            hint = tet_adjacent(hint, lf);
                            break;
                        }
                    }
                }
                return hint;
            }

            // The same walk as locate() with inexact predicates, stopping after \p max_iter
            // tetrahedra. Its result is a hint for locate(), which is what makes locate() faster
            // than walking exactly the whole way ("structural filtering").
            index_t locate_inexact(
                const double* p, index_t hint, index_t max_iter
            ) const {

                hint = walkable_hint(hint);

                index_t t = hint;
                index_t t_pred = NO_INDEX;

            still_walking:
                {
                    const double* pv[4];
                    pv[0] = vertex_ptr(tet_vertex(t,0));
                    pv[1] = vertex_ptr(tet_vertex(t,1));
                    pv[2] = vertex_ptr(tet_vertex(t,2));
                    pv[3] = vertex_ptr(tet_vertex(t,3));

                    for(index_t f = 0; f < 4; ++f) {

                        index_t t_next = tet_adjacent(t,f);

                        //  If the opposite tet is -1, then it means that
                        // we are trying to locate() within a tetrahedralization
                        // from which the infinite tets were removed.
                        if(t_next == NO_INDEX) {
                            return NO_INDEX;
                        }

                        //   If the candidate next tetrahedron is the
                        // one we came from, then we know already that
                        // the orientation is positive, thus we examine
                        // the next candidate (or exit the loop if they
                        // are exhausted).
                        if(t_next == t_pred) {
                            continue ;
                        }

                        //   To test the orientation of p w.r.t. the facet f of
                        // t, we replace vertex number f with p in t (same
                        // convention as in CGAL).
                        const double* pv_bkp = pv[f];
                        pv[f] = p;
                        Sign ori = orient_3d_inexact(pv[0], pv[1], pv[2], pv[3]);

                        //   If the orientation is not negative, then we cannot
                        // walk towards t_next, and examine the next candidate
                        // (or exit the loop if they are exhausted).
                        if(ori != NEGATIVE) {
                            pv[f] = pv_bkp;
                            continue;
                        }

                        //  If the opposite tet is a virtual tet, then
                        // the point has a positive orientation relative
                        // to the facet on the border of the convex hull,
                        // thus t_next is a tet in conflict and we are
                        // done.
                        if(tet_is_virtual(t_next)) {
                            return t_next;
                        }

                        //   If we reach this point, then t_next is a valid
                        // successor, thus we are still walking.
                        t_pred = t;
                        t = t_next;
                        if(--max_iter != 0) {
                            goto still_walking;
                        }
                    }
                }

                //   If we reach this point, we did not find a valid successor
                // for walking (a face for which p has negative orientation),
                // thus we reached the tet for which p has all positive
                // face orientations (i.e. the tet that contains p).

                return t;
            }

            // The tetrahedron that contains \p p, or one incident to it when \p p is on a face,
            // edge or vertex. \p orient receives the orientation of \p p with respect to the four
            // facets. Outside the convex hull of the points inserted so far, the answer is a
            // virtual tetrahedron.
            index_t locate(
                const double* p, index_t hint, Sign* orient
            ) const {

                //   Try improving the hint by using the
                // inexact locate function. This gains
                // (a little bit) performance (a few
                // percent in total Delaunay computation
                // time), but it is better than nothing...
                //   Note: there is a maximum number of tets
                // traversed by locate_inexact()  (2500)
                // since there exists configurations in which
                // locate_inexact() loops forever !
                hint = locate_inexact(p, hint, 2500);

                hint = walkable_hint(hint);

                index_t t = hint;
                index_t t_pred = NO_INDEX;

            still_walking:
                {
                    const double* pv[4];
                    pv[0] = vertex_ptr(tet_vertex(t,0));
                    pv[1] = vertex_ptr(tet_vertex(t,1));
                    pv[2] = vertex_ptr(tet_vertex(t,2));
                    pv[3] = vertex_ptr(tet_vertex(t,3));

                    // Start from a random facet
                    index_t f0 = index_t(random_int32()) % 4;
                    for(index_t df = 0; df < 4; ++df) {
                        index_t f = (f0 + df) % 4;

                        index_t t_next = tet_adjacent(t,f);

                        //  If the opposite tet is -1, then it means that
                        // we are trying to locate() within a tetrahedralization
                        // from which the infinite tets were removed.
                        if(t_next == NO_INDEX) {
                            return NO_INDEX;
                        }

                        //   If the candidate next tetrahedron is the
                        // one we came from, then we know already that
                        // the orientation is positive, thus we examine
                        // the next candidate (or exit the loop if they
                        // are exhausted).
                        if(t_next == t_pred) {
                            orient[f] = POSITIVE ;
                            continue ;
                        }

                        //   To test the orientation of p w.r.t. the facet f of
                        // t, we replace vertex number f with p in t (same
                        // convention as in CGAL).
                        // This is equivalent to tet_facet_point_orient3d(t,f,p)
                        // (but less costly, saves a couple of lookups)
                        const double* pv_bkp = pv[f];
                        pv[f] = p;
                        orient[f] = PCK::orient_3d(pv[0], pv[1], pv[2], pv[3]);

                        //   If the orientation is not negative, then we cannot
                        // walk towards t_next, and examine the next candidate
                        // (or exit the loop if they are exhausted).
                        if(orient[f] != NEGATIVE) {
                            pv[f] = pv_bkp;
                            continue;
                        }

                        //  If the opposite tet is a virtual tet, then
                        // the point has a positive orientation relative
                        // to the facet on the border of the convex hull,
                        // thus t_next is a tet in conflict and we are
                        // done.
                        if(tet_is_virtual(t_next)) {
                            for(index_t lf = 0; lf < 4; ++lf) {
                                orient[lf] = POSITIVE;
                            }
                            return t_next;
                        }

                        //   If we reach this point, then t_next is a valid
                        // successor, thus we are still walking.
                        t_pred = t;
                        t = t_next;
                        goto still_walking;
                    }
                }

                //   If we reach this point, we did not find a valid successor
                // for walking (a face for which p has negative orientation),
                // thus we reached the tet for which p has all positive
                // face orientations (i.e. the tet that contains p).

                return t;
            }

            // Hands facet \p lf of \p t, which is on the boundary of the conflict zone, to the
            // cavity along with the three vertices it is spanned by.
            void record_cavity_facet(index_t t, index_t lf) {
                cavity_.new_facet(
                    t, lf,
                    tet_vertex(t, tet_facet_vertex_[lf][0]),
                    tet_vertex(t, tet_facet_vertex_[lf][1]),
                    tet_vertex(t, tet_facet_vertex_[lf][2])
                );
            }

            // Propagates find_conflict_zone() outwards from \p t_in, which must already be
            // marked as in conflict.
            void find_conflict_zone_iterative(
                const double* p, index_t t_in,
                index_t& t_bndry, index_t& f_bndry,
                index_t& first, index_t& last
            ) {

                S_.push(t_in);

                while(!S_.empty()) {

                    index_t t = S_.top();
                    S_.pop();

                    for(index_t lf = 0; lf < 4; ++lf) {
                        index_t t2 = tet_adjacent(t, lf);

                        if(
                            tet_is_in_list(t2)  // known as conflict
                        ) {
                            continue;
                        }

                        if(
                            tet_is_marked(t2)     // known as non-conflict
                        ) {
                            record_cavity_facet(t, lf);
                            continue;
                        }


                        if(tet_is_conflict(t2, p)) {
                            // Chain t2 in conflict list
                            add_tet_to_list(t2, first, last);
                            S_.push(t2);
                            continue;
                        }

                        //   At this point, t is in conflict
                        // and t2 is not in conflict.
                        // We keep a reference to a tet on the boundary
                        t_bndry = t;
                        f_bndry = lf;
                        // Mark t2 as visited (but not conflict)
                        mark_tet(t2);

                        record_cavity_facet(t, lf);

                    }
                }
            }

            // The tetrahedra in conflict with point \p v, given \p t containing it and the
            // \p orient locate() returned alongside. \p t_bndry and \p f_bndry come back as a
            // tetrahedron adjacent to the boundary of the conflict zone and the facet it is
            // adjacent along. The conflict tetrahedra are chained from \p first to \p last
            // through tet_next(). The chain is empty when \p v already exists in the
            // triangulation.
            void find_conflict_zone(
                index_t v,
                index_t t, const Sign* orient,
                index_t& t_bndry, index_t& f_bndry,
                index_t& first, index_t& last
            ) {
                cavity_.clear();

                first = last = END_OF_LIST;

                //  Generate a unique stamp from current vertex index, used for marking
                // tetrahedra. Storage is shared with list chaining, so the stamp carries
                // NOT_IN_LIST_BIT.
                cur_stamp_ = (v | NOT_IN_LIST_BIT);

                // Pointer to the coordinates of the point to be inserted
                const double* p = vertex_ptr(v);


                // Test whether the point already exists in
                // the triangulation. The point already exists
                // if it's located on three faces of the
                // tetrahedron returned by locate().
                int nb_zero =
                    (orient[0] == ZERO) +
                    (orient[1] == ZERO) +
                    (orient[2] == ZERO) +
                    (orient[3] == ZERO) ;

                if(nb_zero >= 3) {
                    return;
                }

                // Note: points on edges and on facets are
                // handled by the way tet_is_in_conflict()
                // is implemented, that naturally inserts
                // the correct tetrahedra in the conflict list.


                // Mark t as conflict
                add_tet_to_list(t, first, last);

                // A small optimization: if the point to be inserted
                // is on some faces of the located tetrahedron, insert
                // the neighbors accros those faces in the conflict list.
                // It saves a couple of calls to the predicates in this
                // specific case (combinatorics are in general less
                // expensive than the predicates).
                if(nb_zero != 0) {
                    for(index_t lf = 0; lf < 4; ++lf) {
                        if(orient[lf] == ZERO) {
                            index_t t2 = tet_adjacent(t, lf);
                            add_tet_to_list(t2, first, last);
                        }
                    }
                    for(index_t lf = 0; lf < 4; ++lf) {
                        if(orient[lf] == ZERO) {
                            index_t t2 = tet_adjacent(t, lf);
                            find_conflict_zone_iterative(
                                p,t2,t_bndry,f_bndry,first,last
                            );
                        }
                    }
                }

                // Determine the conflict list by greedy propagation from t.
                find_conflict_zone_iterative(p,t,t_bndry,f_bndry,first,last);
            }

            // Used by stellate_conflict_zone_iterative(): given \p t1 on the border of the
            // conflict zone along facet \p t1fborder, finds the tetrahedron \p t2 across facet
            // \p t1ft2 that is also on the border and shares an edge with both facets.
            // \p t2fborder and \p t2ft1 come back as t2's facets on the border and towards t1.
            // Returns true if \p t2 is newly created, false if it is an old tetrahedron still in
            // conflict.
            //
            // Long for an inline function, but geogram measured a modest gain from keeping it
            // one.
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
                    cur_t = next_t;
                    // The local facet of cur_t incident to the oriented edge [ev1,ev2].
                    const index_t* T = &(cell_to_v_store_[4 * cur_t]);
                    cur_f = halfedge_facet_[find_4(T,ev1)][find_4(T,ev2)];
                    next_t = tet_adjacent(cur_t, cur_f);
                }

                //  At this point, cur_t is in conflict zone and
                // next_t is outside the conflict zone.
                const index_t* T = &(cell_to_v_store_[4 * next_t]);
                index_t lv1 = find_4(T, ev1);
                index_t lv2 = find_4(T, ev2);
                index_t f12 = halfedge_facet_[lv1][lv2];
                index_t f21 = halfedge_facet_[lv2][lv1];
                t2 = tet_adjacent(next_t,f21);
                index_t v_neigh_opposite = tet_vertex(next_t,f12);
                t2ft1 = find_4(&(cell_to_v_store_[4 * t2]), v_neigh_opposite);
                t2fborder = cur_f;

                //  Test whether the found neighboring tet was created
                //  (then return true) or is an old tet in conflict
                //  (then return false).
                return(t2 != cur_t);
            }

            // The star of tetrahedra around \p v filling the conflict zone, built by walking the
            // border: starting from facet \p t1fbord of \p t1 and recursing until the zone is
            // filled. \p t1fprev is the facet of \p t1 that it was reached from, or NO_INDEX for
            // the first one. Used when the Cavity overflowed its arrays.
            index_t stellate_conflict_zone_iterative(
                index_t v, index_t t1, index_t t1fbord, index_t t1fprev = NO_INDEX
            ) {
                //   This function is de-recursified because some degenerate
                // inputs can cause stack overflow (system stack is limited to
                // a few megs). For instance, it can happen when a large number
                // of points are on the same sphere exactly.

                S2_.push_back({t1, uint8_t(t1fbord), uint8_t(t1fprev), 0, 0, 0});

                index_t new_t;   // the newly created tetrahedron.

                index_t t1ft2;   // traverses the 4 facets of t1.

                index_t t2;      // the tetrahedron on the border of
                                 // the conflict zone that shares an
                                 // edge with t1 along t1ft2.

                index_t t2fbord; // the facet of t2 on the border of
                                 // the conflict zone.

                index_t t2ft1;   // the facet of t2 that is incident to t1.

            entry_point:
                t1 = S2_.back().t1;
                t1fbord = S2_.back().t1fbord;
                t1fprev = S2_.back().t1fprev;


                // Create new tetrahedron with same vertices as t_bndry
                new_t = new_tetrahedron(
                    tet_vertex(t1,0),
                    tet_vertex(t1,1),
                    tet_vertex(t1,2),
                    tet_vertex(t1,3)
                );

                // Replace in new_t the vertex opposite to t1fbord with v
                cell_to_v_store_[4 * new_t + t1fbord] = v;

                // Connect new_t with t1's neighbor accros t1fbord
                {
                    index_t tbord = tet_adjacent(t1,t1fbord);
                    set_tet_adjacent(new_t, t1fbord, tbord);
                    set_tet_adjacent(tbord, find_tet_adjacent(tbord,t1), new_t);
                }

                //  Lookup new_t's neighbors accros its three other
                // facets and connect them
                for(t1ft2=0; t1ft2<4; ++t1ft2) {

                    if(t1ft2 == t1fprev || tet_adjacent(new_t,t1ft2) != NO_INDEX) {
                        continue;
                    }

                    // Get t1's neighbor along the border of the conflict zone
                    if(!get_neighbor_along_conflict_zone_border(
                           t1,t1fbord,t1ft2, t2,t2fbord,t2ft1
                       )) {
                        //   If t1's neighbor is not a new tetrahedron,
                        // create a new tetrahedron through a recursive call.
                        S2_.back().new_t = new_t;
                        S2_.back().t1ft2 = uint8_t(t1ft2);
                        S2_.back().t2ft1 = uint8_t(t2ft1);
                        S2_.push_back({t2, uint8_t(t2fbord), uint8_t(t2ft1), 0, 0, 0});
                        goto entry_point;

                    return_point:
                        // This is the return value of the called function.
                        index_t result = new_t;
                        S2_.pop_back();

                        // Special case: we were in the outermost frame,
                        // then we (truly) return from the function.
                        if(S2_.empty()) {
                            return result;
                        }

                        t1 = S2_.back().t1;
                        t1fbord = S2_.back().t1fbord;
                        t1fprev = S2_.back().t1fprev;
                        new_t = S2_.back().new_t;
                        t1ft2 = S2_.back().t1ft2;
                        t2ft1 = S2_.back().t2ft1;
                        t2 = result;
                    }

                    set_tet_adjacent(t2, t2ft1, new_t);
                    set_tet_adjacent(new_t, t1ft2, t2);
                }

                // Except for the initial call (see "Special case" above),
                // the nested calls all come from the same location,
                // thus there is only one possible return point
                // (no need to push any return address).
                goto return_point;
            }

            // Fills the conflict zone with a star of tetrahedra around \p v, one per facet on
            // its border, and returns one of them. Used when the Cavity collected while
            // traversing the conflict zone did not overflow its arrays.
            index_t stellate_cavity(index_t v) {

                index_t new_tet = NO_INDEX;

                for(index_t f=0; f<cavity_.nb_f_; ++f) {
                    index_t old_tet = cavity_.tglobal_[f];
                    index_t lf = cavity_.boundary_f_[f];
                    index_t t_neigh = tet_adjacent(old_tet, lf);
                    index_t v1 = cavity_.f2v_[f][0];
                    index_t v2 = cavity_.f2v_[f][1];
                    index_t v3 = cavity_.f2v_[f][2];
                    new_tet = new_tetrahedron(v, v1, v2, v3);
                    set_tet_adjacent(new_tet, 0, t_neigh);
                    set_tet_adjacent(
                        t_neigh, find_tet_adjacent(t_neigh,old_tet), new_tet
                    );
                    cavity_.tglobal_[f] = new_tet;
                }

                for(index_t f=0; f<cavity_.nb_f_; ++f) {
                    new_tet = cavity_.tglobal_[f];
                    index_t neigh1, neigh2, neigh3;
                    cavity_.get_facet_neighbor_tets(f, neigh1, neigh2, neigh3);
                    set_tet_adjacent(new_tet, 1, neigh1);
                    set_tet_adjacent(new_tet, 2, neigh2);
                    set_tet_adjacent(new_tet, 3, neigh3);
                }

                return new_tet;
            }

            // Inserts point \p v, and returns one of the tetrahedra incident to it, or NO_INDEX
            // when \p v already exists in the triangulation.
            index_t insert(index_t v, index_t hint) {
                index_t t_bndry = NO_INDEX;
                index_t f_bndry = NO_INDEX;
                index_t first_conflict = NO_INDEX;
                index_t last_conflict = NO_INDEX;

                const double* p = vertex_ptr(v);

                Sign orient[4];
                index_t t = locate(p, hint, orient);
                find_conflict_zone(
                    v,t,orient,t_bndry,f_bndry,first_conflict,last_conflict
                );

                if(first_conflict == END_OF_LIST) {
                    return NO_INDEX;
                }

                index_t new_tet = NO_INDEX;
                if(cavity_.OK_) {
                    new_tet = stellate_cavity(v);
                } else {
                    new_tet = stellate_conflict_zone_iterative(v,t_bndry,f_bndry);
                }

                // Recycle the tetrahedra of the conflict zone.
                cell_next_[last_conflict] = first_free_;
                first_free_ = first_conflict;

                // Return one of the newly created tets
                return new_tet;
            }

            // Four non-coplanar points to start the incremental construction from, or false if
            // the points are all coplanar.
            bool create_first_tetrahedron(
                index_t& iv0, index_t& iv1, index_t& iv2, index_t& iv3
            ) {
                if(nb_vertices_ < 4) {
                    return false;
                }

                iv0 = 0;

                iv1 = 1;
                while(
                    iv1 < nb_vertices_ &&
                    points_are_identical_3d(
                        vertex_ptr(iv0), vertex_ptr(iv1)
                    )
                ) {
                    ++iv1;
                }
                if(iv1 == nb_vertices_) {
                    return false;
                }

                iv2 = iv1 + 1;
                while(
                    iv2 < nb_vertices_ &&
                    points_are_colinear_3d(
                        vertex_ptr(iv0), vertex_ptr(iv1), vertex_ptr(iv2)
                    )
                ) {
                    ++iv2;
                }
                if(iv2 == nb_vertices_) {
                    return false;
                }

                iv3 = iv2 + 1;
                Sign s = ZERO;
                while(
                    iv3 < nb_vertices_ &&
                    (s = PCK::orient_3d(
                        vertex_ptr(iv0), vertex_ptr(iv1),
                        vertex_ptr(iv2), vertex_ptr(iv3)
                    )) == ZERO
                ) {
                    ++iv3;
                }

                if(iv3 == nb_vertices_) {
                    return false;
                }


                if(s == NEGATIVE) {
                    std::swap(iv2, iv3);
                }

                // Create the first tetrahedron
                index_t t0 = new_tetrahedron(iv0, iv1, iv2, iv3);

                // Create the first four virtual tetrahedra surrounding it
                index_t t[4];
                for(index_t f = 0; f < 4; ++f) {
                    // In reverse order since it is an adjacent tetrahedron
                    index_t v1 = tet_vertex(t0, tet_facet_vertex_[f][2]);
                    index_t v2 = tet_vertex(t0, tet_facet_vertex_[f][1]);
                    index_t v3 = tet_vertex(t0, tet_facet_vertex_[f][0]);
                    t[f] = new_tetrahedron(VERTEX_AT_INFINITY, v1, v2, v3);
                }

                // Connect the virtual tetrahedra to the real one
                for(index_t f=0; f<4; ++f) {
                    set_tet_adjacent(t[f], 0, t0);
                    set_tet_adjacent(t0, f, t[f]);
                }

                // Interconnect the four virtual tetrahedra along their common
                // faces
                for(index_t f = 0; f < 4; ++f) {
                    // In reverse order since it is an adjacent tetrahedron
                    index_t lv1 = tet_facet_vertex_[f][2];
                    index_t lv2 = tet_facet_vertex_[f][1];
                    index_t lv3 = tet_facet_vertex_[f][0];
                    set_tet_adjacent(t[f], 1, t[lv1]);
                    set_tet_adjacent(t[f], 2, t[lv2]);
                    set_tet_adjacent(t[f], 3, t[lv3]);
                }

                return true;
            }

            // Sets the vertices, three doubles each, and computes the cells into
            // cell_to_v_store_.
            void set_vertices(index_t nb_vertices, const double* vertices) {
                nb_vertices_ = nb_vertices;
                vertices_ = vertices;

                index_t expected_tetra = nb_vertices * 7;

                cell_to_v_store_.reserve(expected_tetra * 4);
                cell_to_cell_store_.reserve(expected_tetra * 4);
                cell_next_.reserve(expected_tetra);

                //   Sort the vertices spatially. This makes localisation
                // faster. geogram could be asked to skip this, which changes the insertion order
                // and with it the triangulation of cospherical input.
                compute_BRIO_order(nb_vertices, vertex_ptr(0), reorder_);

                // The indices of the vertices of the first tetrahedron.
                index_t v0, v1, v2, v3;
                if(!create_first_tetrahedron(v0, v1, v2, v3)) {
                    // All the points are coplanar.
                    return;
                }

                index_t hint = NO_INDEX;
                // Insert all the vertices incrementally.
                for(index_t i = 0; i < nb_vertices; ++i) {
                    index_t v = reorder_[i];
                    // Do not re-insert the first four vertices.
                    if(v != v0 && v != v1 && v != v2 && v != v3) {
                        index_t new_hint = insert(v, hint);
                        if(new_hint != NO_INDEX) {
                            hint = new_hint;
                        }
                    }
                }


                //   Compress cell_to_v_store_: remove the free and virtual tetrahedra, keeping
                // the real ones in order. The adjacencies are not needed once construction is
                // over, so cell_to_cell_store_ is left behind. geogram compacted and remapped it
                // here too.
                index_t nb_tets = 0;
                for(index_t t = 0; t < nb_cells(); ++t) {
                    if(tet_is_real(t)) {
                        if(t != nb_tets) {
                            std::memcpy(
                                &cell_to_v_store_[nb_tets * 4],
                                &cell_to_v_store_[t * 4],
                                4 * sizeof(index_t)
                            );
                        }
                        ++nb_tets;
                    }
                }
                cell_to_v_store_.resize(4 * nb_tets);
            }

            const double* vertices_ = nullptr;
            index_t nb_vertices_ = 0;

            // The triangulation, four vertex indices per tetrahedron once set_vertices() has
            // compacted the free and virtual tetrahedra out of it.
            std::vector<index_t> cell_to_v_store_;
            std::vector<index_t> cell_to_cell_store_;
            std::vector<index_t> cell_next_;
            std::vector<index_t> reorder_;
            index_t cur_stamp_ = 0; // used for marking
            index_t first_free_ = END_OF_LIST;

            // Used by the de-recursified find_conflict_zone_iterative().
            std::stack<index_t> S_;

            // One frame of the de-recursified stellate_conflict_zone_iterative(): the parameters
            // it was called with, then the locals it restores when a nested call returns. The
            // narrow types keep the frame small, because degenerate input can make the stack
            // very deep.
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

    }

    std::vector<index_t> delaunay_3d(index_t nb_vertices, const double* vertices) {
        Delaunay3d T;
        T.set_vertices(nb_vertices, vertices);
        return std::move(T.cell_to_v_store_);
    }

} }
