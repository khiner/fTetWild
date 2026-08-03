// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "MeshIO.hpp"

#include <floattetwild/MshSaver.h>
#include <floattetwild/Logger.hpp>

#include <floattetwild/Timer.h>

#include <floattetwild/SurfaceMeshLoad.hpp>
#include <floattetwild/geo_mesh.h>
#include <floattetwild/geo_mesh_reorder.h>

namespace floatTetWild {

void reorder_and_read_back(geo::Mesh&             mesh,
                           std::vector<Vector3>&  points,
                           std::vector<Vector3i>& faces,
                           std::vector<int>&      tags)
{
    const bool has_tags = (tags.size() == mesh.facets.nb());

    geo::vector<geo::index_t> facet_permutation;
    geo::mesh_reorder(mesh, &facet_permutation);

    if (has_tags) {
        std::vector<int> reordered(mesh.facets.nb(), 0);
        for (geo::index_t f = 0; f < facet_permutation.size(); ++f)
            reordered[f] = tags[facet_permutation[f]];
        tags.swap(reordered);
    }

    points.resize(mesh.vertices.nb());
    for (size_t i = 0; i < points.size(); i++)
        points[i] << mesh.vertices.point(i)[0], mesh.vertices.point(i)[1],
          mesh.vertices.point(i)[2];

    faces.resize(mesh.facets.nb());
    for (size_t i = 0; i < faces.size(); i++)
        faces[i] << int(mesh.facets.vertex(i, 0)), int(mesh.facets.vertex(i, 1)),
          int(mesh.facets.vertex(i, 2));
}

bool load_mesh(const std::string&     path,
                       std::vector<Vector3>&  points,
                       std::vector<Vector3i>& faces,
                       geo::Mesh&             input,
                       std::vector<int>&      flags)
{
    logger().debug("Loading mesh at {}...", path);

    if (!load_surface_mesh(path, input))
        return false;

    reorder_and_read_back(input, points, faces, flags);
    return true;
}

void write_mesh(const std::string&         path,
                        const Mesh&                mesh,
                        const std::vector<Scalar>& color,
                        const bool                 binary,
                        const bool                 separate_components)
{
    assert(color.empty() || color.size() == mesh.tet_vertices.size() ||
           color.size() == mesh.tets.size());

    logger().info("Writing mesh to {}...", path);
    Timer timer;
    timer.start();

    PyMesh::MshSaver mesh_saver(path, binary);

    std::vector<int> old_2_new(mesh.tet_vertices.size(), -1);
    int              cnt_v = 0;
    for (size_t i = 0; i < mesh.tet_vertices.size(); i++) {
        if (!mesh.tet_vertices[i].is_removed)
            old_2_new[i] = cnt_v++;
    }
    const int cnt_t = mesh.get_t_num();

    PyMesh::VectorF V_flat(cnt_v * 3);
    PyMesh::VectorI T_flat(cnt_t * 4);
    PyMesh::VectorI C_flat;
    if (separate_components)
        C_flat.resize(cnt_t);

    size_t index = 0;
    for (size_t i = 0; i < mesh.tet_vertices.size(); ++i) {
        if (mesh.tet_vertices[i].is_removed)
            continue;
        for (int j = 0; j < 3; j++)
            V_flat[index * 3 + j] = mesh.tet_vertices[i].pos[j];
        index++;
    }

    index = 0;
    for (size_t i = 0; i < mesh.tets.size(); ++i) {
        if (mesh.tets[i].is_removed)
            continue;
        // The saver wants the opposite orientation, so the last two corners are swapped.
        T_flat[index * 4 + 0] = old_2_new[mesh.tets[i][0]];
        T_flat[index * 4 + 1] = old_2_new[mesh.tets[i][1]];
        T_flat[index * 4 + 2] = old_2_new[mesh.tets[i][3]];
        T_flat[index * 4 + 3] = old_2_new[mesh.tets[i][2]];
        if (separate_components)
            C_flat[index] = mesh.tets[i].scalar;
        index++;
    }

    mesh_saver.save_mesh(V_flat, T_flat, C_flat);

    // One colour per tet or one per vertex, whichever the caller sized the array to.
    const bool per_tet = color.size() == mesh.tets.size();
    if (per_tet || color.size() == mesh.tet_vertices.size()) {
        PyMesh::VectorF color_flat(per_tet ? cnt_t : cnt_v);
        index = 0;
        for (size_t i = 0; i < color.size(); i++) {
            if (per_tet ? mesh.tets[i].is_removed : mesh.tet_vertices[i].is_removed)
                continue;
            color_flat[index++] = color[i];
        }
        if (per_tet)
            mesh_saver.save_elem_scalar_field("color", color_flat);
        else
            mesh_saver.save_scalar_field("color", color_flat);
    }

    logger().info(" took {}s", timer.getElapsedTimeInSec());
}

}  // namespace floatTetWild
