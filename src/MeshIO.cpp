// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "MeshIO.hpp"

#include "FloatTetwild.h"
#include "Logger.hpp"

#include "Timer.h"

#include "geo/geo_mesh.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace floatTetWild {
namespace {

// Reading a surface mesh in one of the formats -i advertises, chosen by extension.
//
// This replaces GEO::mesh_load. geogram's readers for these formats merge nothing and sort
// nothing: they append each vertex as it appears in the file and each facet in file order, and
// mesh_load's only post-pass fills facet adjacency, which nothing downstream of here reads.
// Reproducing the file order is therefore all that is required for the mesh handed to
// mesh_reorder to be identical, and the corpus checks that for .stl. Duplicate vertices are
// expected here and are resolved later by remove_duplicates() in Simplification.cpp.
//
// Polygons with more than three vertices are fanned into triangles as they are read, so the
// result is always simplicial. geogram left them as polygons and fTetWild then called
// GEO::mesh_repair(TRIANGULATE|QUIET), which also dropped degenerate facets, reoriented facets
// anti-Moebius, split non-manifold vertices and oriented normals. Those passes are gone: no
// input the corpus or the tests exercise is polygonal, and none of them ran for triangle input
// in the first place.

void add_vertex(geo::Mesh& mesh, double x, double y, double z)
{
    const geo::index_t v = mesh.create_vertices(1);
    double*            p = mesh.point_ptr(v);
    p[0]                 = x;
    p[1]                 = y;
    p[2]                 = z;
}

// Fan the polygon around its first vertex, the same split GEO::MeshFacets::triangulate() makes.
void add_polygon(geo::Mesh& mesh, const std::vector<geo::index_t>& vertices)
{
    for (size_t i = 1; i + 1 < vertices.size(); ++i) {
        const geo::index_t f = mesh.create_triangles(1);
        mesh.set_facet_vertex(f, 0, vertices[0]);
        mesh.set_facet_vertex(f, 1, vertices[i]);
        mesh.set_facet_vertex(f, 2, vertices[i + 1]);
    }
}

// A binary STL is 80 header bytes, a uint32 triangle count and 50 bytes per triangle. Many
// binary files still start with "solid", so the size is what decides, as in geogram.
bool stl_is_binary(const std::string& path, uint32_t& nb_triangles)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        return false;
    in.seekg(80);
    char count[4];
    if (!in.read(count, 4))
        return false;
    nb_triangles = uint32_t(uint8_t(count[0])) | uint32_t(uint8_t(count[1])) << 8 |
                   uint32_t(uint8_t(count[2])) << 16 | uint32_t(uint8_t(count[3])) << 24;
    in.seekg(0, std::ios::end);
    return in.tellg() == std::streamoff(uint64_t(nb_triangles) * 50 + 84);
}

bool load_stl_binary(const std::string& path, geo::Mesh& mesh, uint32_t nb_triangles)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        return false;
    in.seekg(84);

    mesh.create_vertices(nb_triangles * 3);
    mesh.create_triangles(nb_triangles);

    for (uint32_t t = 0; t < nb_triangles; ++t) {
        float record[12];  // normal, then the three corners
        if (!in.read(reinterpret_cast<char*>(record), sizeof(record)))
            return false;
        in.seekg(2, std::ios::cur);  // attribute byte count

        for (geo::index_t lv = 0; lv < 3; ++lv) {
            double* p = mesh.point_ptr(3 * t + lv);
            for (int c = 0; c < 3; ++c)
                p[c] = double(record[3 + lv * 3 + c]);
            mesh.set_facet_vertex(t, lv, 3 * t + lv);
        }
    }
    return true;
}

bool load_stl_ascii(const std::string& path, geo::Mesh& mesh)
{
    std::ifstream in(path);
    if (!in)
        return false;

    std::vector<geo::index_t> facet_vertices;
    std::string               line, keyword;
    bool                      facet_opened = false;

    while (std::getline(in, line)) {
        std::istringstream fields(line);
        if (!(fields >> keyword))
            continue;
        if (keyword == "outer") {
            facet_vertices.clear();
            facet_opened = true;
        }
        else if (keyword == "endloop") {
            facet_opened = false;
            add_polygon(mesh, facet_vertices);
        }
        else if (keyword == "vertex") {
            double xyz[3];
            if (!(fields >> xyz[0] >> xyz[1] >> xyz[2])) {
                logger().error("{}: vertex line has fewer than 3 coordinates", path);
                return false;
            }
            facet_vertices.push_back(mesh.nb_vertices());
            add_vertex(mesh, xyz[0], xyz[1], xyz[2]);
        }
    }

    if (facet_opened) {
        logger().error("{}: last facet is not closed", path);
        return false;
    }
    if (mesh.nb_facets() == 0) {
        logger().error("{}: no facet", path);
        return false;
    }
    return true;
}

bool load_off(const std::string& path, geo::Mesh& mesh)
{
    std::ifstream in(path);
    if (!in)
        return false;

    std::string line;
    auto        next_line = [&]() {
        while (std::getline(in, line)) {
            const size_t first = line.find_first_not_of(" \t\r");
            if (first != std::string::npos && line[first] != '#')
                return true;
        }
        return false;
    };

    if (!next_line() || line.compare(0, 3, "OFF") != 0) {
        logger().error("{}: unrecognized header", path);
        return false;
    }

    geo::index_t nb_vertices = 0, nb_edges = 0, nb_facets = 0;
    if (!next_line() ||
        !(std::istringstream(line) >> nb_vertices >> nb_edges >> nb_facets)) {
        logger().error("{}: unrecognized header", path);
        return false;
    }

    mesh.create_vertices(nb_vertices);
    for (geo::index_t v = 0; v < nb_vertices; ++v) {
        double xyz[3];
        if (!next_line() || !(std::istringstream(line) >> xyz[0] >> xyz[1] >> xyz[2])) {
            logger().error("{}: expected {} vertices", path, nb_vertices);
            return false;
        }
        double* p = mesh.point_ptr(v);
        std::copy(xyz, xyz + 3, p);
    }

    // Facet lines can carry trailing fields such as a per-facet colour, so read exactly the
    // number of corners the leading count announces and ignore the rest of the line.
    while (next_line()) {
        std::istringstream fields(line);
        geo::index_t       nb_facet_vertices = 0;
        if (!(fields >> nb_facet_vertices) || nb_facet_vertices < 3)
            continue;

        std::vector<geo::index_t> facet_vertices(nb_facet_vertices);
        for (geo::index_t j = 0; j < nb_facet_vertices; ++j) {
            if (!(fields >> facet_vertices[j]) || facet_vertices[j] >= mesh.nb_vertices()) {
                logger().error("{}: facet references an invalid vertex", path);
                return false;
            }
        }
        add_polygon(mesh, facet_vertices);
    }
    return true;
}

bool load_obj(const std::string& path, geo::Mesh& mesh)
{
    std::ifstream in(path);
    if (!in)
        return false;

    std::string line, keyword, field;
    while (std::getline(in, line)) {
        std::istringstream fields(line);
        if (!(fields >> keyword))
            continue;

        if (keyword == "v") {
            double xyz[3] = {0, 0, 0};
            fields >> xyz[0] >> xyz[1] >> xyz[2];
            add_vertex(mesh, xyz[0], xyz[1], xyz[2]);
        }
        else if (keyword == "f") {
            std::vector<geo::index_t> facet_vertices;
            while (fields >> field) {
                // A corner is v, v/vt, v//vn or v/vt/vn. Only the vertex index matters here.
                const int index = std::atoi(field.c_str());
                if (index == 0)
                    continue;
                const geo::index_t v = index < 0
                                         ? geo::index_t(int(mesh.nb_vertices()) + index)
                                         : geo::index_t(index - 1);
                if (v >= mesh.nb_vertices()) {
                    logger().error("{}: facet references an invalid vertex: {}", path, index);
                    return false;
                }
                facet_vertices.push_back(v);
            }
            add_polygon(mesh, facet_vertices);
        }
    }
    return true;
}

// Only what a surface mesh needs: the x, y and z properties of the vertex element and the
// vertex-index list of the face element. Other elements and properties are skipped by size.
struct PlyProperty {
    std::string name;
    std::string type;        // scalar type, or the value type of a list
    std::string count_type;  // empty unless the property is a list
};

struct PlyElement {
    std::string              name;
    size_t                   count = 0;
    std::vector<PlyProperty> properties;
};

int ply_type_size(const std::string& type)
{
    if (type == "char" || type == "uchar" || type == "int8" || type == "uint8")
        return 1;
    if (type == "short" || type == "ushort" || type == "int16" || type == "uint16")
        return 2;
    if (type == "int" || type == "uint" || type == "int32" || type == "uint32" ||
        type == "float" || type == "float32")
        return 4;
    if (type == "double" || type == "float64")
        return 8;
    return 0;
}

double ply_read_binary(std::istream& in, const std::string& type, bool big_endian)
{
    const int size = ply_type_size(type);
    char      bytes[8];
    if (size == 0 || !in.read(bytes, size))
        return 0.0;
    if (big_endian)
        std::reverse(bytes, bytes + size);

    if (type == "float" || type == "float32") {
        float value;
        std::memcpy(&value, bytes, 4);
        return double(value);
    }
    if (type == "double" || type == "float64") {
        double value;
        std::memcpy(&value, bytes, 8);
        return value;
    }
    uint64_t value = 0;
    std::memcpy(&value, bytes, size_t(size));
    if (type == "char" || type == "int8")
        return double(int8_t(value));
    if (type == "short" || type == "int16")
        return double(int16_t(value));
    if (type == "int" || type == "int32")
        return double(int32_t(value));
    return double(value);
}

bool load_ply(const std::string& path, geo::Mesh& mesh)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        return false;

    std::string line, keyword;
    if (!std::getline(in, line) || line.compare(0, 3, "ply") != 0) {
        logger().error("{}: unrecognized header", path);
        return false;
    }

    std::string             format;
    std::vector<PlyElement> elements;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        std::istringstream fields(line);
        fields >> keyword;
        if (keyword == "format") {
            fields >> format;
        }
        else if (keyword == "element") {
            elements.emplace_back();
            fields >> elements.back().name >> elements.back().count;
        }
        else if (keyword == "property" && !elements.empty()) {
            PlyProperty property;
            std::string type;
            fields >> type;
            if (type == "list") {
                fields >> property.count_type >> property.type >> property.name;
            }
            else {
                property.type = type;
                fields >> property.name;
            }
            elements.back().properties.push_back(property);
        }
        else if (keyword == "end_header") {
            break;
        }
    }

    if (format != "ascii" && format != "binary_little_endian" && format != "binary_big_endian") {
        logger().error("{}: unsupported ply format '{}'", path, format);
        return false;
    }
    const bool ascii      = format == "ascii";
    const bool big_endian = format == "binary_big_endian";

    for (const PlyElement& element : elements) {
        const bool is_vertex = element.name == "vertex";
        const bool is_face   = element.name == "face";

        for (size_t i = 0; i < element.count; ++i) {
            std::istringstream ascii_fields;
            if (ascii) {
                if (!std::getline(in, line))
                    return false;
                ascii_fields.str(line);
            }

            double                    xyz[3] = {0, 0, 0};
            std::vector<geo::index_t> facet_vertices;

            // One value of the given type, from the text line or from the binary stream.
            const auto read = [&](const std::string& type) {
                double value = 0;
                if (ascii)
                    ascii_fields >> value;
                else
                    value = ply_read_binary(in, type, big_endian);
                return value;
            };

            for (const PlyProperty& property : element.properties) {
                if (!property.count_type.empty()) {
                    const size_t n = size_t(read(property.count_type));
                    for (size_t k = 0; k < n; ++k) {
                        const double value = read(property.type);
                        if (is_face && property.name == "vertex_indices")
                            facet_vertices.push_back(geo::index_t(value));
                    }
                    continue;
                }

                const double value = read(property.type);
                if (is_vertex) {
                    if (property.name == "x")
                        xyz[0] = value;
                    else if (property.name == "y")
                        xyz[1] = value;
                    else if (property.name == "z")
                        xyz[2] = value;
                }
            }

            if (is_vertex)
                add_vertex(mesh, xyz[0], xyz[1], xyz[2]);
            if (is_face && facet_vertices.size() >= 3) {
                for (geo::index_t v : facet_vertices) {
                    if (v >= mesh.nb_vertices()) {
                        logger().error("{}: facet references an invalid vertex", path);
                        return false;
                    }
                }
                add_polygon(mesh, facet_vertices);
            }
        }
    }
    return true;
}

bool load_surface_mesh(const std::string& path, geo::Mesh& mesh)
{
    // Drop any attributes too: the --csg path loads every operand into the same mesh object.
    mesh.clear();

    const size_t dot       = path.find_last_of('.');
    std::string  extension = dot == std::string::npos ? "" : path.substr(dot + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
        return char(std::tolower(c));
    });
    bool ok = false;
    if (extension == "stl") {
        uint32_t nb_triangles = 0;
        ok = stl_is_binary(path, nb_triangles) ? load_stl_binary(path, mesh, nb_triangles)
                                               : load_stl_ascii(path, mesh);
    }
    else if (extension == "off")
        ok = load_off(path, mesh);
    else if (extension == "obj")
        ok = load_obj(path, mesh);
    else if (extension == "ply")
        ok = load_ply(path, mesh);
    else {
        logger().error("Unsupported file format: {}", extension);
        return false;
    }

    if (!ok) {
        logger().error("Could not load file: {}", path);
        return false;
    }

    // geogram zeroed NaNs on load and warned. Keep that, so a file with NaNs still produces a
    // mesh rather than poisoning every predicate downstream.
    bool has_nan = false;
    for (geo::index_t v = 0; v < mesh.nb_vertices(); ++v) {
        double* p = mesh.point_ptr(v);
        for (int c = 0; c < 3; ++c) {
            if (std::isnan(p[c])) {
                has_nan = true;
                p[c]    = 0.0;
            }
        }
    }
    if (has_nan)
        logger().warn("Found NaNs in {}", path);

    return true;
}

}  // namespace

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

// Writing a Gmsh 2.2 file, following PyMesh, Copyright (c) 2015 by Qingnan Zhou. Nodes are 3D and
// elements are tets, which is all fTetWild writes. PyMesh also handled 2D nodes and triangle,
// quad and hex elements.
void write_mesh(const std::string&         path,
                const Mesh&                mesh,
                const std::vector<Scalar>& color,
                const bool                 binary,
                const bool                 separate_components)
{
    assert(color.empty() || color.size() == mesh.tets.size());

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

    fout << "$MeshFormat" << std::endl;
    fout << "2.2 " << (binary ? 1 : 0) << " " << sizeof(Scalar) << std::endl;
    if (binary) {
        int one = 1;
        fout.write((char*)&one, sizeof(int));
    }
    fout << "$EndMeshFormat" << std::endl;
    fout.flush();

    fout << "$Nodes" << std::endl;
    fout << cnt_v << std::endl;
    int node_idx = 1;
    for (const MeshVertex& v : mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        if (!binary) {
            fout << node_idx << " " << v.pos[0] << " " << v.pos[1] << " " << v.pos[2] << std::endl;
        } else {
            fout.write((char*)&node_idx, sizeof(int));
            fout.write(reinterpret_cast<const char*>(v.pos.data()), sizeof(Scalar) * 3);
        }
        node_idx++;
    }
    fout << "$EndNodes" << std::endl;
    fout.flush();

    fout << "$Elements" << std::endl;
    fout << cnt_t << std::endl;
    if (cnt_t > 0) {
        int elem_type = 4;  // Gmsh's element type for a 4-node tetrahedron
        int tags      = separate_components ? 2 : 0;
        if (binary) {
            fout.write((char*)&elem_type, sizeof(int));
            fout.write((const char*)&cnt_t, sizeof(int));
            fout.write((char*)&tags, sizeof(int));
        }
        int elem_num = 1;
        for (const MeshTet& t : mesh.tets) {
            if (t.is_removed)
                continue;
            // The msh reader wants the opposite orientation, so the last two corners are swapped,
            // and its node indices count from 1.
            const std::array<int, 4> elem = {{old_2_new[t[0]] + 1,
                                              old_2_new[t[1]] + 1,
                                              old_2_new[t[3]] + 1,
                                              old_2_new[t[2]] + 1}};
            const int                comp = int(t.scalar);  // the CSG component tag

            if (!binary) {
                fout << elem_num << " " << elem_type << " " << tags << " ";
                if (separate_components)
                    fout << comp << " " << comp << " ";
                for (int j = 0; j < 4; j++)
                    fout << elem[j] << " ";
                fout << std::endl;
            } else {
                fout.write((char*)&elem_num, sizeof(int));
                if (separate_components) {
                    std::array<int, 2> comps = {{comp, comp}};
                    fout.write((char*)comps.data(), sizeof(int) * 2);
                }
                fout.write(reinterpret_cast<const char*>(elem.data()), sizeof(int) * 4);
            }
            elem_num++;
        }
    }
    fout << "$EndElements" << std::endl;
    fout.flush();

    // One value per tet, named "color", which is what the viewers the msh is written for expect.
    // msh also has a per-node form, which nothing writes any more.
    if (color.size() == mesh.tets.size()) {
        fout << "$ElementData" << std::endl;
        fout << "1" << std::endl;          // num string tags.
        fout << "\"color\"" << std::endl;
        fout << "1" << std::endl;          // num real tags.
        fout << "0.0" << std::endl;        // time value.
        fout << "3" << std::endl;          // num int tags.
        fout << "0" << std::endl;          // the time step
        fout << "1" << std::endl;          // 1-component scalar field.
        fout << cnt_t << std::endl;        // number of nodes or elements

        int idx = 1;
        for (size_t i = 0; i < color.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            if (binary) {
                fout.write((char*)&idx, sizeof(int));
                fout.write(reinterpret_cast<const char*>(&color[i]), sizeof(Scalar));
            } else {
                fout << idx << " " << color[i] << std::endl;
            }
            idx++;
        }
        fout << "$EndElementData" << std::endl;
        fout.flush();
    }

    logger().info(" took {}s", timer.getElapsedTimeInSec());
}

}  // namespace floatTetWild
