// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "MeshIO.hpp"

#include <floattetwild/Logger.hpp>

#include <floattetwild/Timer.h>

#include <floattetwild/SurfaceMeshLoad.hpp>
#include <floattetwild/geo_mesh.h>
#include <floattetwild/geo_mesh_reorder.h>

#include <array>
#include <fstream>
#include <stdexcept>

namespace floatTetWild {
namespace {

// Writing a Gmsh 2.2 file, from PyMesh, Copyright (c) 2015 by Qingnan Zhou. Nodes are 3D and
// elements are tets, which is all fTetWild writes. PyMesh also handled 2D nodes and triangle,
// quad and hex elements. Coordinates and indices arrive flat, three and four per element.

// Gmsh's element type for a 4-node tetrahedron.
constexpr int    TetElementType = 4;
constexpr size_t NodesPerTet    = 4;
constexpr size_t Dim            = 3;

void save_header(std::ofstream& fout, bool binary)
{
    fout << "$MeshFormat" << std::endl;
    fout << "2.2 " << (binary ? 1 : 0) << " " << sizeof(Scalar) << std::endl;
    if (binary) {
        int one = 1;
        fout.write((char*)&one, sizeof(int));
    }
    fout << "$EndMeshFormat" << std::endl;
    fout.flush();
}

void save_nodes(std::ofstream& fout, bool binary, const MatrixXs& nodes)
{
    fout << "$Nodes" << std::endl;
    fout << nodes.size() / Dim << std::endl;
    for (size_t i = 0; i < size_t(nodes.size()); i += Dim) {
        const Scalar* v        = &nodes[int(i)];
        int           node_idx = int(i / Dim) + 1;

        if (!binary) {
            fout << node_idx << " " << v[0] << " " << v[1] << " " << v[2] << std::endl;
        } else {
            fout.write((char*)&node_idx, sizeof(int));
            fout.write(reinterpret_cast<const char*>(v), sizeof(Scalar) * Dim);
        }
    }
    fout << "$EndNodes" << std::endl;
    fout.flush();
}

void save_elements(std::ofstream& fout, bool binary, const MatrixXi& elements,
                   const MatrixXi& components)
{
    const size_t num_elements = elements.size() / NodesPerTet;

    fout << "$Elements" << std::endl;
    fout << num_elements << std::endl;

    if (num_elements > 0) {
        int elem_type = TetElementType;
        int num_elems = int(num_elements);
        int tags      = components.size() > 0 ? 2 : 0;
        if (binary) {
            fout.write((char*)&elem_type, sizeof(int));
            fout.write((char*)&num_elems, sizeof(int));
            fout.write((char*)&tags, sizeof(int));
        }
        for (size_t i = 0; i < size_t(elements.size()); i += NodesPerTet) {
            int                          elem_num = int(i / NodesPerTet) + 1;
            std::array<int, NodesPerTet> elem;
            for (size_t j = 0; j < NodesPerTet; j++)
                elem[j] = elements[int(i + j)] + 1;

            if (!binary) {
                fout << elem_num << " " << elem_type << " " << tags << " ";
                if (components.size() > 0)
                    fout << components[elem_num - 1] << " " << components[elem_num - 1] << " ";

                for (size_t j = 0; j < NodesPerTet; j++)
                    fout << elem[j] << " ";
                fout << std::endl;
            } else {
                fout.write((char*)&elem_num, sizeof(int));
                if (components.size() > 0) {
                    std::array<int, 2> comps = {
                      {components[elem_num - 1], components[elem_num - 1]}};
                    fout.write((char*)comps.data(), sizeof(int) * 2);
                }
                fout.write((char*)elem.data(), sizeof(int) * NodesPerTet);
            }
        }
    }
    fout << "$EndElements" << std::endl;
    fout.flush();
}

// One value per node or per element, which is the same section apart from its name and length.
void save_field(std::ofstream& fout, bool binary, const char* section,
                const std::string& fieldname, const MatrixXs& field)
{
    fout << "$" << section << std::endl;
    fout << "1" << std::endl;            // num string tags.
    fout << "\"" << fieldname << "\"" << std::endl;
    fout << "1" << std::endl;            // num real tags.
    fout << "0.0" << std::endl;          // time value.
    fout << "3" << std::endl;            // num int tags.
    fout << "0" << std::endl;            // the time step
    fout << "1" << std::endl;            // 1-component scalar field.
    fout << field.size() << std::endl;   // number of nodes or elements

    for (int i = 0; i < field.size(); i++) {
        int idx = i + 1;
        if (binary) {
            fout.write((char*)&idx, sizeof(int));
            fout.write(reinterpret_cast<const char*>(&field[i]), sizeof(Scalar));
        } else {
            fout << idx << " " << field[i] << std::endl;
        }
    }
    fout << "$End" << section << std::endl;
    fout.flush();
}

}  // namespace

void reorder_and_read_back(geo::Mesh&             mesh,
                           std::vector<Vector3>&  points,
                           std::vector<Vector3i>& faces,
                           std::vector<int>&      tags)
{
    const bool has_tags = (tags.size() == mesh.nb_facets());

    geo::vector<geo::index_t> facet_permutation;
    geo::mesh_reorder(mesh, &facet_permutation);

    if (has_tags) {
        std::vector<int> reordered(mesh.nb_facets(), 0);
        for (geo::index_t f = 0; f < facet_permutation.size(); ++f)
            reordered[f] = tags[facet_permutation[f]];
        tags.swap(reordered);
    }

    points.resize(mesh.nb_vertices());
    for (size_t i = 0; i < points.size(); i++)
        points[i] << mesh.point(i)[0], mesh.point(i)[1],
          mesh.point(i)[2];

    faces.resize(mesh.nb_facets());
    for (size_t i = 0; i < faces.size(); i++)
        faces[i] << int(mesh.facet_vertex(i, 0)), int(mesh.facet_vertex(i, 1)),
          int(mesh.facet_vertex(i, 2));
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

    std::ofstream fout(path.c_str(), binary ? std::fstream::binary : std::fstream::out);
    if (!fout)
        throw std::runtime_error("Error opening " + path + " to write msh file.");

    std::vector<int> old_2_new(mesh.tet_vertices.size(), -1);
    int              cnt_v = 0;
    for (size_t i = 0; i < mesh.tet_vertices.size(); i++) {
        if (!mesh.tet_vertices[i].is_removed)
            old_2_new[i] = cnt_v++;
    }
    const int cnt_t = mesh.get_t_num();

    MatrixXs V_flat(cnt_v * 3);
    MatrixXi T_flat(cnt_t * 4);
    MatrixXi C_flat;
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

    save_header(fout, binary);
    save_nodes(fout, binary, V_flat);
    save_elements(fout, binary, T_flat, C_flat);

    // One colour per tet or one per vertex, whichever the caller sized the array to.
    const bool per_tet = color.size() == mesh.tets.size();
    if (per_tet || color.size() == mesh.tet_vertices.size()) {
        MatrixXs color_flat(per_tet ? cnt_t : cnt_v);
        index = 0;
        for (size_t i = 0; i < color.size(); i++) {
            if (per_tet ? mesh.tets[i].is_removed : mesh.tet_vertices[i].is_removed)
                continue;
            color_flat[index++] = color[i];
        }
        save_field(fout, binary, per_tet ? "ElementData" : "NodeData", "color", color_flat);
    }

    logger().info(" took {}s", timer.getElapsedTimeInSec());
}

}  // namespace floatTetWild
