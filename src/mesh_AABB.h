// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see geo/LICENSE.geogram.
// Source: geogram/mesh/mesh_AABB.h
// Copied rather than reimplemented: the traversal order decides which of two equidistant facets a
// query lands on, and the mesher's envelope tests carry that facet forward as a hint.

#pragma once

#include <floattetwild/geo_mesh.h>
#include <floattetwild/geo_geometry.h>

namespace floatTetWild {

    // Axis-aligned bounding box tree of the facets of a mesh, for locating the facet nearest to a
    // 3d query point.
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
        // well as written, and geo::NO_FACET means there is none yet.
        void facet_in_envelope_with_hint(
            const geo::vec3& p, double sq_epsilon,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist
        ) const;

        // \overload With no facet to start from.
        geo::index_t facet_in_envelope(
            const geo::vec3& p, double sq_epsilon, geo::vec3& nearest_point, double& sq_dist
        ) const {
            geo::index_t nearest_facet = geo::NO_FACET;
            facet_in_envelope_with_hint(p, sq_epsilon, nearest_facet, nearest_point, sq_dist);
            return nearest_facet;
        }

    private:
        // A facet near \p p to start a search from. Descending to whichever child's bounding box
        // has its centre nearer the query point gains up to 10% on a large mesh (20M facets), as
        // compared to picking the starting facet randomly, by letting the walk prune earlier.
        void get_nearest_facet_hint(
            const geo::vec3& p,
            geo::index_t& nearest_facet, geo::vec3& nearest_point, double& sq_dist
        ) const;

        geo::vector<geo::Box> bboxes_;
        const geo::Mesh& mesh_;
    };

inline void get_point_facet_nearest_point(
    const geo::Mesh& M,
    const geo::vec3& p,
    geo::index_t f,
    geo::vec3& nearest_p,
    double& squared_dist
) {
    using namespace floatTetWild::geo;
    const vec3& p1 = M.facet_point(f, 0);
    const vec3& p2 = M.facet_point(f, 1);
    const vec3& p3 = M.facet_point(f, 2);
    double lambda1, lambda2, lambda3;  // barycentric coords, not used.
    squared_dist = Geom::point_triangle_squared_distance(
        p, p1, p2, p3, nearest_p, lambda1, lambda2, lambda3
    );
}

}
