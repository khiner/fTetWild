// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/CSGTreeParser.hpp>

#include <floattetwild/Logger.hpp>
#include <floattetwild/Timer.h>

#include <geogram/mesh/mesh_reorder.h>


namespace floatTetWild {

namespace {
    // Build a GEO::Mesh from vertex and face vectors, carry the tags through as a facet
// attribute, Morton reorder it, and read the permuted result back out. Named load_mesh
// on MeshIO, but it touches no files.
void build_reordered_sf_mesh(std::vector<Vector3>&  points,
                           std::vector<Vector3i>& faces,
                           GEO::Mesh&             input,
                           std::vector<int>&      flags)
    {
        logger().debug("Loading mesh from data...");
        Timer timer;
        timer.start();
        input.clear(false, false);

        input.clear();
        input.vertices.create_vertices((int)points.size());
        for (int i = 0; i < (int)input.vertices.nb(); ++i)
        {
            GEO::vec3 &p = input.vertices.point(i);
            p[0] = points[i](0);
            p[1] = points[i](1);
            p[2] = points[i](2);
        }
        // Setup faces
        input.facets.create_triangles((int)faces.size());

        for (int c = 0; c < (int)input.facets.nb(); ++c)
        {
            for (int lv = 0; lv < 3; ++lv)
            {
                input.facets.set_vertex(c, lv, faces[c](lv));
            }
        }

        bool is_valid = (flags.size() == input.facets.nb());
        if (is_valid) {
            assert(flags.size() == input.facets.nb());
            GEO::Attribute<int> bflags(input.facets.attributes(), "bbflags");
            for (int index = 0; index < (int)input.facets.nb(); ++index) {
                bflags[index] = flags[index];
            }
        }

        GEO::mesh_reorder(input, GEO::MESH_ORDER_MORTON);

        if (is_valid) {
            flags.clear();
            flags.resize(input.facets.nb());
            GEO::Attribute<int> bflags(input.facets.attributes(), "bbflags");
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
    }
}  // namespace

    void CSGTreeParser::get_meshes_aux(const json &csg_tree_node, std::vector<std::string> &meshes, std::map<std::string, int> &existings, int &index, json &current_node)
    {
        current_node["operation"] = csg_tree_node["operation"];
        if(csg_tree_node["left"].is_string())
        {
            int id = -1;
            const std::string current = csg_tree_node["left"];
            const auto iter = existings.find(current);
            if(iter == existings.end()){
                meshes.push_back(current);
                existings[current] = index;
                id = index;
                ++index;
            }
            else
                id = iter->second;

            current_node["left"] = id;
        }
        else
        {
            json left_node;
            get_meshes_aux(csg_tree_node["left"], meshes, existings, index, left_node);
            current_node["left"] = left_node;
        }


        if(csg_tree_node["right"].is_string())
        {
            int id = -1;
            const std::string current = csg_tree_node["right"];
            const auto iter = existings.find(current);
            if(iter == existings.end()){
                meshes.push_back(current);
                existings[current] = index;
                id = index;
                ++index;
            }
            else
                id = iter->second;

            current_node["right"] = id;
        }
        else
        {
            json right_node;
            get_meshes_aux(csg_tree_node["right"], meshes, existings, index, right_node);
            current_node["right"] = right_node;
        }
    }

    void CSGTreeParser::merge(const std::vector<std::vector<Vector3>> &Vs, const std::vector<std::vector<Vector3i>> &Fs, std::vector<Vector3> &V, std::vector<Vector3i> &F, GEO::Mesh &sf_mesh, std::vector<int> &tags)
    {
        V.clear();
        F.clear();

        for(const auto &vv : Vs)
            V.insert(V.end(), vv.begin(), vv.end());

        int offset = 0;
        int size = 0;
        for(int id = 0; id < Fs.size(); ++id) {
            const auto &ff = Fs[id];

            for(const auto fid : ff)
                F.push_back(Vector3i(fid(0)+offset, fid(1)+offset, fid(2)+offset));

            tags.insert(tags.begin()+size, Fs[id].size(), id);
            size+=Fs[id].size();
            offset += Vs[id].size();
        }



        build_reordered_sf_mesh(V, F, sf_mesh, tags);

    }

    void CSGTreeParser::get_max_id_aux(const json &csg_tree_node, int &max)
    {
        if(csg_tree_node["left"].is_number()){
            const int id = csg_tree_node["left"];
            max = std::max(max, id);
        }
        else
            get_max_id_aux(csg_tree_node["left"], max);


        if(csg_tree_node["right"].is_number()){
            const int id = csg_tree_node["right"];
            max = std::max(max, id);
        }
        else
            get_max_id_aux(csg_tree_node["right"], max);
    }


    bool CSGTreeParser::keep_tet(const json &csg_tree_with_ids, const int t_id, const std::vector<Eigen::VectorXd> &w)
    {
        const std::string op = csg_tree_with_ids["operation"];

        bool left_inside;
        if(csg_tree_with_ids["left"].is_number()){
            int id = csg_tree_with_ids["left"];
            left_inside = w[id][t_id] > 0.5;
        }
        else
        {
            left_inside = keep_tet(csg_tree_with_ids["left"], t_id, w);
        }

        bool right_inside;
        if(csg_tree_with_ids["right"].is_number()){
            int id = csg_tree_with_ids["right"];
            right_inside = w[id][t_id] > 0.5;
        }
        else
        {
            right_inside = keep_tet(csg_tree_with_ids["right"], t_id, w);
        }


        if(op == "union")
            return left_inside || right_inside;
        if(op == "intersection")
            return left_inside && right_inside;
        if(op == "difference")
            return left_inside && !right_inside;

        assert(false);
        return false;
    }

}
