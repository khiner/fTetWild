// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/FloatTetDelaunay.h>

#include <floattetwild/Logger.hpp>

#include <algorithm>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/geo_delaunay_3d.h>

namespace floatTetWild {
	namespace {
        void
        get_bb_corners(const Parameters &params, const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
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

        void
        compute_voxel_points(const Vector3 &min, const Vector3 &max, const Parameters &params, const AABBWrapper &tree,
                             std::vector<Vector3> &voxels) {
            const Vector3 diag = max - min;
            Vector3i n_voxels = (diag / (params.bbox_diag_length * params.box_scale)).cast<int>();

            for (int d = 0; d < 3; ++d)
                n_voxels(d) = std::max(n_voxels(d), 1);

            Vector3 delta;
            for (int d = 0; d < 3; ++d)
                delta(d) = diag(d) / Scalar(n_voxels(d));

            voxels.clear();
            voxels.reserve((n_voxels(0) + 1) * (n_voxels(1) + 1) * (n_voxels(2) + 1));

            const double sq_distg = 100 * params.eps_2;
            geo::vec3 nearest_point;

            for (int i = 0; i <= n_voxels(0); ++i) {
                const Scalar px = (i == n_voxels(0)) ? max(0) : (min(0) + delta(0) * i);
                for (int j = 0; j <= n_voxels(1); ++j) {
                    const Scalar py = (j == n_voxels(1)) ? max(1) : (min(1) + delta(1) * j);
                    for (int k = 0; k <= n_voxels(2); ++k) {
                        const Scalar pz = (k == n_voxels(2)) ? max(2) : (min(2) + delta(2) * k);
                        if (tree.get_sq_dist_to_sf(Vector3(px, py, pz)) > sq_distg)
                            voxels.emplace_back(px, py, pz);
                    }
                }
            }
        }
    }

	void tetrahedralize(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const AABBWrapper &tree,
	        Mesh &mesh, std::vector<bool> &is_face_inserted) {
        const Parameters &params = mesh.params;
        auto &tet_vertices = mesh.tet_vertices;
        auto &tets = mesh.tets;

        is_face_inserted.resize(input_faces.size(), false);

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

        geo::Delaunay3d T;
        T.set_vertices(n_pts, V_d.data());
        tets.resize(T.nb_cells());
        const auto &tet2v = T.cell_to_v();
        for (int i = 0; i < T.nb_cells(); i++) {
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
}
