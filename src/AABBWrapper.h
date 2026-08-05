#pragma once

#include "geo/geo_geometry.h"
#include "geo/geo_mesh.h"
#include "Mesh.hpp"

#include <memory>

namespace floatTetWild {

    namespace geo {
        // Axis-aligned bounding box.
        struct Box {
            double xyz_min[3];
            double xyz_max[3];
        };
    }

    // Axis-aligned bounding box tree of the facets of a mesh, for locating the facet nearest to a
    // 3d query point.
    //
    // Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
    // Original licence: BSD 3-clause, see geo/LICENSE.geogram. Source: geogram/mesh/mesh_AABB.
    // Copied rather than reimplemented: the traversal order decides which of two equidistant
    // facets a query lands on, and the mesher's envelope tests carry that facet forward as a hint.
    class MeshFacetsAABBWithEps {
    public:
        // \p M must already be triangulated and Morton-ordered by mesh_reorder(): geogram did both
        // here, and every caller now does them beforehand.
        MeshFacetsAABBWithEps(const geo::Mesh& M);

        // The facet nearest to \p p, along with the point on it and the squared distance.
        geo::index_t nearest_facet(
            const geo::vec3& p, geo::vec3& nearest_point, double& sq_dist
        ) const;

        // The same, but stopping as soon as a point within \p sq_epsilon is found. Starting from a
        // facet the caller already has saves the search for a hint, so \p nearest_facet is read as
        // well as written, and geo::NO_FACET means there is none yet. \p nearest_point and
        // \p sq_dist are overwritten before they are read, whatever they arrive holding.
        void facet_in_envelope_with_hint(
            const geo::vec3& p, double sq_epsilon,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist
        ) const;

    private:
        // A facet near \p p to start a search from. Descending to whichever child's bounding box
        // has its centre nearer the query point gains up to 10% on a large mesh (20M facets), as
        // compared to picking the starting facet randomly, by letting the walk prune earlier.
        void get_nearest_facet_hint(
            const geo::vec3& p,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist
        ) const;

        std::vector<geo::Box> bboxes_;
        const geo::Mesh& mesh_;
    };

    inline geo::vec3 to_geo_p(const Vector3 &p) { return geo::vec3(p[0], p[1], p[2]); }

    // Three envelopes over the same kind of query: sf is the input surface, b is its open boundary
    // and tmp_b is the boundary as it stands while triangles are still being inserted.
    class AABBWrapper {
    public:
        // The only piece the callers look at directly, to read back the facet a query landed on.
        const geo::Mesh &sf_mesh;

        AABBWrapper(const geo::Mesh &sf_mesh) : sf_mesh(sf_mesh), sf_tree(sf_mesh) {}

        // The length of the diagonal of the bounding box of the surface mesh.
        Scalar get_sf_diag() const;

        void init_b_mesh_and_tree(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, Mesh &mesh);

        void init_tmp_b_mesh_and_tree(const std::vector<Vector3> &input_vertices,
                                      const std::vector<std::array<int, 2>> &b_edges1,
                                      const Mesh &mesh, const std::vector<std::array<int, 2>> &b_edges2);

        // Moves p onto the nearest point of the envelope and returns how far it moved, squared.
        inline Scalar project_to_sf(Vector3 &p) const { return project(sf_tree, p); }
        inline Scalar project_to_b(Vector3 &p) const { return project(*b_tree, p); }
        inline Scalar project_to_tmp_b(Vector3 &p) const { return project(*tmp_b_tree, p); }

        inline int get_nearest_face_sf(const Vector3 &p) const {
            geo::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max();
            return sf_tree.nearest_facet(to_geo_p(p), nearest_p, sq_dist);
        }

        // A set of sample points is outside when any one of them is.
        inline bool is_out_b_envelope(const std::vector<geo::vec3> &ps, const Scalar eps_2,
                                      geo::index_t prev_facet = geo::NO_FACET) const {
            return is_out(*b_tree, ps, eps_2, prev_facet);
        }

        inline bool is_out_tmp_b_envelope(const std::vector<geo::vec3> &ps, const Scalar eps_2,
                                          geo::index_t prev_facet = geo::NO_FACET) const {
            return is_out(*tmp_b_tree, ps, eps_2, prev_facet);
        }

        inline bool is_out_sf_envelope(const Vector3 &p, const Scalar eps_2) const {
            geo::index_t prev_facet;
            return is_out(sf_tree, p, eps_2, prev_facet);
        }

        // The overloads taking prev_facet hand back the facet the query landed on, which the
        // caller passes to the next query as a hint.
        inline bool is_out_sf_envelope(const Vector3 &p, const Scalar eps_2, geo::index_t &prev_facet) const {
            return is_out(sf_tree, p, eps_2, prev_facet);
        }

        inline bool is_out_tmp_b_envelope(const Vector3 &p, const Scalar eps_2, geo::index_t &prev_facet) const {
            return is_out(*tmp_b_tree, p, eps_2, prev_facet);
        }

        // Same query, but carrying the previous facet and its distance across calls. Walking a
        // triangle's samples in order means the next one usually lands on the same facet.
        inline bool is_out_sf_envelope(const geo::vec3& geo_p, const Scalar eps_2,
                                       geo::index_t& prev_facet, double& sq_dist, geo::vec3& nearest_p) const {
            return is_out(sf_tree, geo_p, eps_2, prev_facet, sq_dist, nearest_p);
        }

    private:
        static Scalar project(const MeshFacetsAABBWithEps &tree, Vector3 &p) {
            geo::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max();
            tree.nearest_facet(to_geo_p(p), nearest_p, sq_dist);
            p[0] = nearest_p[0];
            p[1] = nearest_p[1];
            p[2] = nearest_p[2];
            return sq_dist;
        }

        // nearest_p and sq_dist are scratch the query overwrites; the batch overload below hands
        // the same pair back in so that a run of samples reuses one facet.
        static bool is_out(const MeshFacetsAABBWithEps &tree, const geo::vec3 &p,
                           const Scalar eps_2, geo::index_t &prev_facet,
                           double &sq_dist, geo::vec3 &nearest_p) {
            tree.facet_in_envelope_with_hint(p, eps_2, prev_facet, nearest_p, sq_dist);
            return Scalar(sq_dist) > eps_2;
        }

        static bool is_out(const MeshFacetsAABBWithEps &tree, const Vector3 &p, const Scalar eps_2,
                           geo::index_t &prev_facet) {
            geo::vec3 nearest_p;
            double sq_dist;
            prev_facet = geo::NO_FACET;
            return is_out(tree, to_geo_p(p), eps_2, prev_facet, sq_dist, nearest_p);
        }

        static bool is_out(const MeshFacetsAABBWithEps &tree, const std::vector<geo::vec3> &ps,
                           const Scalar eps_2, geo::index_t prev_facet) {
            geo::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max();
            for (const geo::vec3 &p : ps) {
                if (is_out(tree, p, eps_2, prev_facet, sq_dist, nearest_p))
                    return true;
            }
            return false;
        }

        // The b trees are built after construction, once the boundary edges are known, and tmp_b
        // is rebuilt every time more triangles go in.
        geo::Mesh b_mesh;
        geo::Mesh tmp_b_mesh;
        std::unique_ptr<MeshFacetsAABBWithEps> b_tree;
        std::unique_ptr<MeshFacetsAABBWithEps> tmp_b_tree;
        MeshFacetsAABBWithEps sf_tree;
    };
}
