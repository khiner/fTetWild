// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/FloatTetwild.h>

#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/Timer.h>
#include <floattetwild/TriangleInsertion.h>

#include <algorithm>

namespace floatTetWild {

int tetrahedralization(AABBWrapper&           tree,
                       std::vector<Vector3>&  input_vertices,
                       std::vector<Vector3i>& input_faces,
                       std::vector<int>&      input_tags,
                       Mesh&                  mesh,
                       bool                   skip_simplify)
{
    Parameters& params = mesh.params;

    if (input_vertices.empty() || input_faces.empty())
        return EXIT_FAILURE;

    params.init(tree.get_sf_diag());

    stats().push_back(
      {StateInfo::init_id, 0, int(input_vertices.size()), int(input_faces.size()), -1, -1});

    Timer timer;

    // Closes off a step: logs how long it took and records the state it left behind. v_num and
    // t_num count the input for preprocessing, which runs before there is a mesh, and the mesh
    // itself after that. The energies are -1 until there are tets to measure.
    const auto finish_step = [&](int         id,
                                 const char* label,
                                 int         v_num,
                                 int         t_num,
                                 Scalar      max_energy             = -1,
                                 Scalar      avg_energy             = -1,
                                 int         cnt_fail_inserted_face = -1) {
        const double elapsed = timer.getElapsedTimeInSec();
        logger().info("{} {}s", label, elapsed);
        logger().info("");
        stats().push_back(
          {id, elapsed, v_num, t_num, max_energy, avg_energy, cnt_fail_inserted_face});
    };

    timer.start();
    simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);
    finish_step(StateInfo::preprocessing_id,
                "preprocessing",
                input_vertices.size(),
                input_faces.size());

    timer.start();
    std::vector<bool> is_face_inserted(input_faces.size(), false);
    tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    finish_step(StateInfo::tetrahedralization_id,
                "tetrahedralizing",
                mesh.get_v_num(),
                mesh.get_t_num());

    timer.start();
    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
    finish_step(StateInfo::cutting_id,
                "cutting",
                mesh.get_v_num(),
                mesh.get_t_num(),
                mesh.get_max_energy(),
                mesh.get_avg_energy(),
                std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

    timer.start();
    optimization(
      input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
    finish_step(StateInfo::optimization_id,
                "mesh optimization",
                mesh.get_v_num(),
                mesh.get_t_num(),
                mesh.get_max_energy(),
                mesh.get_avg_energy());

    correct_tracked_surface_orientation(mesh, tree);
    logger().info("correct_tracked_surface_orientation done");

    return 0;
}

}  // namespace floatTetWild
