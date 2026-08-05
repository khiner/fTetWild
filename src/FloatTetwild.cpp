// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "FloatTetwild.h"

#include "Logger.hpp"
#include "LocalOperations.h"
#include "MeshImprovement.h"
#include "Simplification.h"
#include "Statistics.h"
#include "Timer.h"
#include "TriangleInsertion.h"
#include "geo/geo_delaunay_3d.h"
#include "geo/geo_mesh_reorder.h"

#include <algorithm>

namespace floatTetWild {
namespace {

void get_bb_corners(const Parameters &params, const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
    min = vertices.front();
    max = vertices.front();

    for (const Vector3 &v : vertices) {
        for (int i = 0; i < 3; i++) {
            min(i) = std::min(min(i), v(i));
            max(i) = std::max(max(i), v(i));
        }
    }

    const Scalar dis = std::max(params.ideal_edge_length, params.eps_input * 2);
    for (int j = 0; j < 3; j++) {
        min[j] -= dis;
        max[j] += dis;
    }

    logger().debug("min = {} {} {}", min[0], min[1], min[2]);
    logger().debug("max = {} {} {}", max[0], max[1], max[2]);
}

void match_bbox_fs(Mesh &mesh, const Vector3 &min, const Vector3 &max) {
    auto get_bbox_fs = [&](const MeshTet &t, int j) {
        std::array<int, 6> cnts = {{0, 0, 0, 0, 0, 0}};
        for (int k = 0; k < 3; k++) {
            Vector3 &pos = mesh.tet_vertices[t[(j + k + 1) % 4]].pos;
            for (int n = 0; n < 3; n++) {
                if (pos[n] == min[n])
                    cnts[n * 2]++;
                else if (pos[n] == max[n])
                    cnts[n * 2 + 1]++;
            }
        }
        for (int i = 0; i < cnts.size(); i++) {
            if (cnts[i] == 3)
                return i;
        }
        return NOT_BBOX;
    };

    for (auto &t: mesh.tets) {
        for (int j = 0; j < 4; j++) {
            t.is_bbox_fs[j] = get_bbox_fs(t, j);
        }
    }
}

void compute_voxel_points(const Vector3 &min, const Vector3 &max, const Parameters &params, const AABBWrapper &tree,
                          std::vector<Vector3> &voxels) {
    const Vector3 diag = max - min;
    // The voxel grid's target cell size, as a fraction of the bounding box diagonal.
    constexpr Scalar box_scale = 1 / 15.0;
    Vector3i n_voxels = (diag / (params.bbox_diag_length * box_scale)).cast<int>();

    for (int d = 0; d < 3; ++d)
        n_voxels(d) = std::max(n_voxels(d), 1);

    Vector3 delta;
    for (int d = 0; d < 3; ++d)
        delta(d) = diag(d) / Scalar(n_voxels(d));

    voxels.clear();
    voxels.reserve((n_voxels(0) + 1) * (n_voxels(1) + 1) * (n_voxels(2) + 1));

    const double sq_distg = 100 * params.eps_2;

    for (int i = 0; i <= n_voxels(0); ++i) {
        const Scalar px = (i == n_voxels(0)) ? max(0) : (min(0) + delta(0) * i);
        for (int j = 0; j <= n_voxels(1); ++j) {
            const Scalar py = (j == n_voxels(1)) ? max(1) : (min(1) + delta(1) * j);
            for (int k = 0; k <= n_voxels(2); ++k) {
                const Scalar pz = (k == n_voxels(2)) ? max(2) : (min(2) + delta(2) * k);
                Vector3      q(px, py, pz);
                if (tree.project_to_sf(q) > sq_distg)
                    voxels.emplace_back(px, py, pz);
            }
        }
    }
}

// The Delaunay mesh of the model vertices and bounding box points.
void tetrahedralize(const std::vector<Vector3>& input_vertices, const AABBWrapper &tree, Mesh &mesh) {
    const Parameters &params = mesh.params;
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    Vector3 min, max;
    get_bb_corners(params, input_vertices, min, max);
    mesh.params.bbox_min = min;
    mesh.params.bbox_max = max;

    std::vector<Vector3> voxel_points;
    compute_voxel_points(min, max, params, tree, voxel_points);

    const int n_pts = input_vertices.size() + voxel_points.size();
    tet_vertices.resize(n_pts);
    for (int i = 0; i < input_vertices.size(); i++)
        tet_vertices[i].pos = input_vertices[i];
    for (int i = 0; i < voxel_points.size(); i++)
        tet_vertices[input_vertices.size() + i].pos = voxel_points[i];

    std::vector<double> V_d(n_pts * 3);
    for (int i = 0; i < tet_vertices.size(); i++) {
        for (int j = 0; j < 3; j++)
            V_d[i * 3 + j] = tet_vertices[i].pos[j];
    }

    const std::vector<geo::index_t> tet2v = geo::delaunay_3d(n_pts, V_d.data());
    tets.resize(tet2v.size() / 4);
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; ++j) {
            const int v_id = tet2v[i * 4 + j];

            tets[i][j] = v_id;
            tet_vertices[v_id].conn_tets.push_back(i);
        }
        std::swap(tets[i][1], tets[i][3]);
    }

    for (const auto &t : mesh.tets) {
        if (is_inverted(mesh.tet_vertices[t[0]].pos, mesh.tet_vertices[t[1]].pos,
                        mesh.tet_vertices[t[2]].pos, mesh.tet_vertices[t[3]].pos)) {
            logger().error("EXIT_INV");
            exit(0);
        }
    }

    match_bbox_fs(mesh, min, max);
}

// Completes the caller-filled Parameters with everything derived from the bounding box diagonal.
void init_params(Parameters& params, Scalar bbox_diag_l)
{
    if (params.stage > 5)
        params.stage = 5;

    params.bbox_diag_length = bbox_diag_l;

    if (params.ideal_edge_length_abs > 0.0) {
        params.ideal_edge_length = params.ideal_edge_length_abs;
        params.ideal_edge_length_rel = params.ideal_edge_length / params.bbox_diag_length;
    }
    else {
        params.ideal_edge_length = params.bbox_diag_length * params.ideal_edge_length_rel;
    }

    params.eps_input = params.bbox_diag_length * params.eps_rel;
    params.dd = params.eps_input / 1.5;

    // The sampled envelope check is conservative by the sampling error bound, the circumradius
    // dd/sqrt(3) of the triangle the samples sit on, so that budget comes out of the epsilon the
    // user asked for.
    double eps_usable = params.eps_input - params.dd / std::sqrt(3);
    params.eps_delta = eps_usable * 0.1;
    params.eps = eps_usable - params.eps_delta * (params.stage - 1);

    params.eps_2 = params.eps * params.eps;

    params.eps_coplanar = params.eps * 0.2;
    if (params.eps_coplanar > params.bbox_diag_length * 1e-6)
        params.eps_coplanar = params.bbox_diag_length * 1e-6;
    params.eps_2_coplanar = params.eps_coplanar * params.eps_coplanar;

    params.eps_simplification = params.eps * 0.8;
    params.eps_2_simplification = params.eps_simplification * params.eps_simplification;
    params.dd_simplification = params.dd / params.eps * params.eps_simplification;

    if (params.min_edge_len_rel < 0)
        params.min_edge_len_rel = params.eps_rel;

    const Scalar split_threshold    = params.ideal_edge_length * (4 / 3.0);
    const Scalar collapse_threshold = params.ideal_edge_length * (4 / 5.0);
    params.split_threshold_2 = split_threshold * split_threshold;
    params.collapse_threshold_2 = collapse_threshold * collapse_threshold;

    logger().debug("bbox_diag_length = {}", params.bbox_diag_length);
    logger().debug("ideal_edge_length = {}", params.ideal_edge_length);

    logger().debug("stage = {}", params.stage);
    logger().debug("eps_input = {}", params.eps_input);
    logger().debug("eps = {}", params.eps);
    logger().debug("eps_simplification = {}", params.eps_simplification);
    logger().debug("eps_coplanar = {}", params.eps_coplanar);

    logger().debug("dd = {}", params.dd);
    logger().debug("dd_simplification = {}", params.dd_simplification);
}

}  // namespace

void reorder_and_read_back(geo::Mesh&             mesh,
                           std::vector<Vector3>&  points,
                           std::vector<Vector3i>& faces,
                           std::vector<int>&      tags)
{
    const bool has_tags = (tags.size() == mesh.nb_facets());

    std::vector<geo::index_t> facet_permutation;
    geo::mesh_reorder(mesh, &facet_permutation);

    if (has_tags) {
        std::vector<int> reordered(mesh.nb_facets(), 0);
        for (geo::index_t f = 0; f < facet_permutation.size(); ++f)
            reordered[f] = tags[facet_permutation[f]];
        tags.swap(reordered);
    }

    points.resize(mesh.nb_vertices());
    for (size_t i = 0; i < points.size(); i++)
        points[i] = Vector3(mesh.point(i)[0], mesh.point(i)[1], mesh.point(i)[2]);

    faces.resize(mesh.nb_facets());
    for (size_t i = 0; i < faces.size(); i++)
        faces[i] = Vector3i(int(mesh.facet_vertex(i, 0)), int(mesh.facet_vertex(i, 1)),
                            int(mesh.facet_vertex(i, 2)));
}

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

    init_params(params, tree.get_sf_diag());

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
    mesh.is_closed = tree.init_b_mesh_and_tree(input_vertices, input_faces);
    finish_step(StateInfo::preprocessing_id,
                "preprocessing",
                input_vertices.size(),
                input_faces.size());

    timer.start();
    std::vector<bool> is_face_inserted(input_faces.size(), false);
    tetrahedralize(input_vertices, tree, mesh);
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    finish_step(StateInfo::tetrahedralization_id,
                "tetrahedralizing",
                mesh.get_v_num(),
                mesh.get_t_num());

    timer.start();
    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
    Scalar max_energy, avg_energy;
    mesh.get_max_avg_energy(max_energy, avg_energy);
    finish_step(StateInfo::cutting_id,
                "cutting",
                mesh.get_v_num(),
                mesh.get_t_num(),
                max_energy,
                avg_energy,
                std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

    timer.start();
    optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree);
    mesh.get_max_avg_energy(max_energy, avg_energy);
    finish_step(StateInfo::optimization_id,
                "mesh optimization",
                mesh.get_v_num(),
                mesh.get_t_num(),
                max_energy,
                avg_energy);

    correct_tracked_surface_orientation(mesh, tree);
    logger().info("correct_tracked_surface_orientation done");

    return 0;
}

}  // namespace floatTetWild
