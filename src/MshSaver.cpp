/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
#include <floattetwild/MshSaver.h>

#include <array>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

namespace floatTetWild {
namespace PyMesh {

// Gmsh's element type for a 4-node tetrahedron.
static const int TetElementType = 4;
static const size_t NodesPerTet = 4;
static const size_t Dim = 3;

MshSaver::MshSaver(const std::string& filename, bool binary) :
        m_binary(binary), m_num_nodes(0), m_num_elements(0) {
    if (!m_binary) {
        fout.open(filename.c_str(), std::fstream::out);
    } else {
        fout.open(filename.c_str(), std::fstream::binary);
    }
    if (!fout)
        throw std::runtime_error("Error opening " + filename + " to write msh file.");
}

MshSaver::~MshSaver() {
    fout.close();
}

void MshSaver::save_mesh(const VectorF& nodes, const VectorI& elements, const VectorI& components) {
    save_header();
    save_nodes(nodes);
    save_elements(elements, components);
}

void MshSaver::save_header() {
    if (!m_binary) {
        fout << "$MeshFormat" << std::endl;
        fout << "2.2 0 " << sizeof(Float) << std::endl;
        fout << "$EndMeshFormat" << std::endl;
    } else {
        fout << "$MeshFormat" << std::endl;
        fout << "2.2 1 " << sizeof(Float) << std::endl;
        int one = 1;
        fout.write((char*)&one, sizeof(int));
        fout << "$EndMeshFormat" << std::endl;
    }
    fout.flush();
}

void MshSaver::save_nodes(const VectorF& nodes) {
    m_num_nodes = nodes.size() / Dim;
    fout << "$Nodes" << std::endl;
    fout << m_num_nodes << std::endl;
    for (size_t i=0; i<nodes.size(); i+=Dim) {
        const Float* v = &nodes[i];
        int node_idx = i/Dim+1;

        if (!m_binary) {
            fout << node_idx << " " << v[0] << " " << v[1] << " " << v[2] << std::endl;
        } else {
            fout.write((char*)&node_idx, sizeof(int));
            fout.write((char*)v, sizeof(Float)*Dim);
        }
    }
    fout << "$EndNodes" << std::endl;
    fout.flush();
}

void MshSaver::save_elements(const VectorI& elements, const VectorI& components) {
    m_num_elements = elements.size() / NodesPerTet;

    fout << "$Elements" << std::endl;
    fout << m_num_elements << std::endl;

    if (m_num_elements > 0) {
        int elem_type = TetElementType;
        int num_elems = m_num_elements;
        int tags = components.size() > 0 ? 2 : 0;
        if (m_binary) {
            fout.write((char*)&elem_type, sizeof(int));
            fout.write((char*)&num_elems, sizeof(int));
            fout.write((char*)&tags, sizeof(int));
        }
        for (size_t i=0; i<elements.size(); i+=NodesPerTet) {
            int elem_num = i/NodesPerTet + 1;
            std::array<int, NodesPerTet> elem;
            for (size_t j=0; j<NodesPerTet; j++) elem[j] = elements[i+j] + 1;

            if (!m_binary) {
                fout << elem_num << " " << elem_type << " " << tags << " ";
                if(components.size() > 0)
                    fout << components[elem_num-1] << " " << components[elem_num-1] << " ";

                for (size_t j=0; j<NodesPerTet; j++) {
                    fout << elem[j] << " ";
                }
                fout << std::endl;
            } else {
                fout.write((char*)&elem_num, sizeof(int));
                if(components.size() > 0){
                    std::array<int, 2> comps = {{components[elem_num-1], components[elem_num-1]}};
                    fout.write((char*)comps.data(), sizeof(int) * 2);
                }
                fout.write((char*)elem.data(), sizeof(int)*NodesPerTet);
            }
        }
    }
    fout << "$EndElements" << std::endl;
    fout.flush();
}

// One value per node or per element, which is the same section apart from its name and length.
void MshSaver::save_field(const char* section, size_t count,
                          const std::string& fieldname, const VectorF& field) {
    assert(field.size() == count);
    fout << "$" << section << std::endl;
    fout << "1" << std::endl; // num string tags.
    fout << "\"" << fieldname << "\"" << std::endl;
    fout << "1" << std::endl; // num real tags.
    fout << "0.0" << std::endl; // time value.
    fout << "3" << std::endl; // num int tags.
    fout << "0" << std::endl; // the time step
    fout << "1" << std::endl; // 1-component scalar field.
    fout << count << std::endl; // number of nodes or elements

    for (size_t i=0; i<count; i++) {
        int idx = i+1;
        if (m_binary) {
            fout.write((char*)&idx, sizeof(int));
            fout.write((char*)&field[i], sizeof(Float));
        } else {
            fout << idx << " " << field[i] << std::endl;
        }
    }
    fout << "$End" << section << std::endl;
    fout.flush();
}

void MshSaver::save_scalar_field(const std::string& fieldname, const VectorF& field) {
    save_field("NodeData", m_num_nodes, fieldname, field);
}

void MshSaver::save_elem_scalar_field(const std::string& fieldname, const VectorF& field) {
    save_field("ElementData", m_num_elements, fieldname, field);
}

}  // namespace PyMesh
}  // namespace floatTetWild
