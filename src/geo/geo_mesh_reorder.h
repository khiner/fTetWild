// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/mesh/mesh_reorder.h
// Copied rather than reimplemented: spatial sort; nth_element is not stable so the incoming arrangement decides ties

#pragma once

#include "geo_mesh.h"

namespace floatTetWild {
namespace geo {

    // Reorders all the elements of a mesh in Morton order, which improves data locality. If
    // \p facet_permutation is non-null it receives the permutation the facets were reordered by:
    // facet k is the one that was at facet_permutation[k]. It is left alone when the mesh has no
    // facets.
    //
    // geogram also has a Hilbert ordering here, which has a continuous mapping between indices
    // and space. No caller asked for it. The Hilbert curve itself is still what
    // compute_BRIO_order() sorts along.
    void mesh_reorder(Mesh& M, vector<index_t>* facet_permutation = nullptr);

    // Sorts \p nb_vertices points, three coordinates each, into BRIO order, which is what makes
    // incremental insertion into a Delaunay triangulation fast. See Amenta, Choi and Rote,
    // "Incremental constructions con brio", Sympos. Comput. Geom. 2003.
    void compute_BRIO_order(
        index_t nb_vertices, const double* vertices, vector<index_t>& sorted_indices
    );

}
}
