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
#include <floattetwild/writeOBJ.h>

#include <floattetwild/SurfaceMeshLoad.hpp>
#include <floattetwild/geo_mesh.h>
#include <floattetwild/geo_mesh_reorder.h>

#include <numeric>

namespace floatTetWild {
namespace {

void extract_volume_mesh(const Mesh&                     mesh,
                         const std::function<bool(int)>& skip_tet,
                         const std::function<bool(int)>& skip_vertex,
                         MatrixXs&                       V,
                         MatrixXi&                       T)
{
    const auto& points = mesh.tet_vertices;
    const auto& tets   = mesh.tets;

    V.resize(points.size(), 3);
    T.resize(tets.size(), 4);

    size_t           index = 0;
    std::vector<int> old_2_new(points.size(), -1);
    for (size_t i = 0; i < points.size(); ++i) {
        if (skip_vertex(i)) {
            continue;
        }
        old_2_new[i] = index;
        V.row(index) = points[i].pos;
        ++index;
    }

    V.conservativeResize(index, 3);

    index = 0;
    for (size_t i = 0; i < tets.size(); ++i) {
        if (skip_tet(i))
            continue;
        for (int j = 0; j < 4; j++) {
            T(index, j) = old_2_new[tets[i][j]];
        }
        ++index;
    }
    T.conservativeResize(index, 4);
}

void write_mesh_aux(const std::string&              path,
                    const Mesh&                     mesh,
                    const std::vector<int>&         t_ids,
                    const std::vector<Scalar>&      color,
                    const bool                      binary,
                    const bool                      separate_components,
                    const std::function<bool(int)>& skip_tet,
                    const std::function<bool(int)>& skip_vertex)
{
    assert(color.empty() || color.size() == mesh.tet_vertices.size() ||
           color.size() == mesh.tets.size());

    PyMesh::MshSaver mesh_saver(path, binary);

    std::map<int, int> old_2_new;
    int                cnt_v = 0;
    for (int i = 0; i < mesh.tet_vertices.size(); i++) {
        if (!skip_vertex(i)) {
            old_2_new[i] = cnt_v;
            cnt_v++;
        }
    }
    int cnt_t = 0;
    for (const int i : t_ids) {
        if (!skip_tet(i))
            cnt_t++;
    }
    PyMesh::VectorF V_flat(cnt_v * 3);
    PyMesh::VectorI T_flat(cnt_t * 4);
    PyMesh::VectorI C_flat;

    if(separate_components)
        C_flat.resize(cnt_t);

    size_t index = 0;
    for (size_t i = 0; i < mesh.tet_vertices.size(); ++i) {
        if (skip_vertex(i))
            continue;
        for (int j = 0; j < 3; j++)
            V_flat[index * 3 + j] = mesh.tet_vertices[i][j];
        index++;
    }

    index = 0;
    for (const int i : t_ids) {
        if (skip_tet(i))
            continue;
        T_flat[index * 4 + 0] = old_2_new[mesh.tets[i][0]];
        T_flat[index * 4 + 1] = old_2_new[mesh.tets[i][1]];
        T_flat[index * 4 + 2] = old_2_new[mesh.tets[i][3]];
        T_flat[index * 4 + 3] = old_2_new[mesh.tets[i][2]];

        if (separate_components)
            C_flat[index] = mesh.tets[i].scalar;

        index++;
    }

    mesh_saver.save_mesh(V_flat, T_flat, C_flat, 3, mesh_saver.TET);

    if (color.size() == mesh.tets.size()) {
        PyMesh::VectorF color_flat(cnt_t);
        index = 0;
        for (const int i : t_ids) {
            if (skip_tet(i))
                continue;
            color_flat[index++] = color[i];
        }
        mesh_saver.save_elem_scalar_field("color", color_flat);
    }
    else if (color.size() == mesh.tet_vertices.size()) {
        PyMesh::VectorF color_flat(cnt_v);
        index = 0;
        for (int i = 0; i < mesh.tet_vertices.size(); i++) {
            if (skip_vertex(i))
                continue;
            color_flat[index++] = color[i];
        }
        mesh_saver.save_scalar_field("color", color_flat);
    }
}
}  // namespace

bool MeshIO::load_mesh(const std::string&     path,
                       std::vector<Vector3>&  points,
                       std::vector<Vector3i>& faces,
                       geo::Mesh&             input,
                       std::vector<int>&      flags)
{
    logger().debug("Loading mesh at {}...", path);
    Timer timer;
    timer.start();

    if (!load_surface_mesh(path, input))
        return false;

    bool is_valid = (flags.size() == input.facets.nb());
    if (is_valid) {
        assert(flags.size() == input.facets.nb());
        geo::Attribute<int> bflags(input.facets.attributes(), "bbflags");
        for (int index = 0; index < (int)input.facets.nb(); ++index) {
            bflags[index] = flags[index];
        }
    }

    // #ifdef FLOAT_TETWILD_USE_FLOAT
    // #else
    // #endif

    geo::mesh_reorder(input, geo::MESH_ORDER_MORTON);

    if (is_valid) {
        flags.clear();
        flags.resize(input.facets.nb());
        geo::Attribute<int> bflags(input.facets.attributes(), "bbflags");
        for (int index = 0; index < (int)input.facets.nb(); ++index) {
            flags[index] = bflags[index];
        }
    }

    points.resize(input.vertices.nb());
    for (size_t i = 0; i < points.size(); i++)
        points[i] << (input.vertices.point(i))[0], (input.vertices.point(i))[1],
          (input.vertices.point(i))[2];

    faces.resize(input.facets.nb());
    for (size_t i = 0; i < faces.size(); i++)
        faces[i] << input.facets.vertex(i, 0), input.facets.vertex(i, 1), input.facets.vertex(i, 2);

    return true;
}

void MeshIO::write_mesh(const std::string&      path,
                        const Mesh&             mesh,
                        const std::vector<int>& t_ids,
                        const bool              only_interior,
                        const bool              binary,
                        const bool              separate_components)
{
    logger().info("Writing mesh to {}...", path);
    Timer timer;
    timer.start();

    if (only_interior) {
        const auto skip_tet    = [&mesh](const int i) { return mesh.tets[i].is_outside; };
        const auto skip_vertex = [&mesh](const int i) { return mesh.tet_vertices[i].is_outside; };
        write_mesh_aux(path, mesh, t_ids, std::vector<Scalar>(), binary, separate_components, skip_tet, skip_vertex);
    }
    else {
        timer.start();
        const auto skip_tet    = [&mesh](const int i) { return mesh.tets[i].is_removed; };
        const auto skip_vertex = [&mesh](const int i) { return mesh.tet_vertices[i].is_removed; };
        write_mesh_aux(path, mesh, t_ids, std::vector<Scalar>(), binary, separate_components, skip_tet, skip_vertex);
    }

    timer.stop();
    logger().info(" took {}s", timer.getElapsedTime());
}

void MeshIO::write_mesh(const std::string&         path,
                        const Mesh&                mesh,
                        const bool                 only_interior,
                        const std::vector<Scalar>& color,
                        const bool                 binary,
                        const bool                 separate_components)
{
    logger().info("Writing mesh to {}...", path);
    Timer timer;
    timer.start();

    std::vector<int> t_ids(mesh.tets.size());
    std::iota(std::begin(t_ids), std::end(t_ids), 0);  // Fill with 0, 1, ..., n.

    if (only_interior) {
        const auto skip_tet    = [&mesh](const int i) { return mesh.tets[i].is_outside; };
        const auto skip_vertex = [&mesh](const int i) { return mesh.tet_vertices[i].is_outside; };
        write_mesh_aux(path, mesh, t_ids, color, binary, separate_components, skip_tet, skip_vertex);
    }
    else {
        const auto skip_tet    = [&mesh](const int i) { return mesh.tets[i].is_removed; };
        const auto skip_vertex = [&mesh](const int i) { return mesh.tet_vertices[i].is_removed; };
        write_mesh_aux(path, mesh, t_ids, color, binary, separate_components, skip_tet, skip_vertex);
    }

    timer.stop();
    logger().info(" took {}s", timer.getElapsedTime());
}

void MeshIO::extract_volume_mesh(const Mesh&      mesh,
                                 MatrixXs& V,
                                 MatrixXi& T,
                                 bool      only_interior)
{
    if (only_interior) {
        const auto skip_tet    = [&mesh](const int i) { return mesh.tets[i].is_outside; };
        const auto skip_vertex = [&mesh](const int i) { return mesh.tet_vertices[i].is_outside; };
        floatTetWild::extract_volume_mesh(mesh, skip_tet, skip_vertex, V, T);
    }
    else {
        const auto skip_tet    = [&mesh](const int i) { return mesh.tets[i].is_removed; };
        const auto skip_vertex = [&mesh](const int i) { return mesh.tet_vertices[i].is_removed; };
        floatTetWild::extract_volume_mesh(mesh, skip_tet, skip_vertex, V, T);
    }
}
}  // namespace floatTetWild
