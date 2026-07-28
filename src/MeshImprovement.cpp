// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/MeshImprovement.h>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/EdgeSplitting.h>
#include <floattetwild/EdgeCollapsing.h>
#include <floattetwild/EdgeSwapping.h>
#include <floattetwild/VertexSmoothing.h>
#include <floattetwild/Parameters.h>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/CSGTreeParser.hpp>

//#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/Statistics.h>

#include <floattetwild/Logger.hpp>

#include <floattetwild/Timer.h>
#include <floattetwild/fast_winding_number.h>


//#define USE_FWN true

void floatTetWild::init(Mesh &mesh, AABBWrapper& tree) {
    cout << "initializing..." << endl;

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.is_on_surface = false;
        v.is_on_bbox = false;
    }
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] <= 0) {
                mesh.tet_vertices[t[mod4(j + 1)]].is_on_surface = true;
                mesh.tet_vertices[t[mod4(j + 2)]].is_on_surface = true;
                mesh.tet_vertices[t[mod4(j + 3)]].is_on_surface = true;
            }
            if (t.is_bbox_fs[j] != NOT_BBOX) {
                mesh.tet_vertices[t[mod4(j + 1)]].is_on_bbox = true;
                mesh.tet_vertices[t[mod4(j + 2)]].is_on_bbox = true;
                mesh.tet_vertices[t[mod4(j + 3)]].is_on_bbox = true;
            }
        }
    }

    //todo: after all faces are inserted, mark true is_on_boundary

//    //vertices: mark is_on_boundary
//    //rebuild b_tree


    if (mesh.params.log_level<3) {
        output_info(mesh, tree);
        int v_num, t_num;
        double max_energy, avg_energy;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
    }
}

void floatTetWild::optimization(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
        Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops) {
    init(mesh, tree);

    ////pre-processing
    mesh.is_limit_length = false;
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, std::array<int, 5>({{0, 1, 0, 0, 0}}));
    mesh.is_limit_length = true;
    cleanup_empty_slots(mesh);
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, std::array<int, 5>({{0, 0, 0, 0, 1}}));

    const int M = 5;
    const int N = 5;

    ////optimization
    int it_after_al_inserted = 0;
    bool is_just_after_update = false;
    bool is_hit_min_edge_length = false;
    std::vector<std::array<Scalar, 2>> quality_queue;
    int cnt_increase_epsilon = mesh.params.stage - 1;
    for (int it = 0; it < mesh.params.max_its; it++) {
        if (mesh.is_input_all_inserted)
            it_after_al_inserted++;

        Scalar max_energy, avg_energy;
        get_max_avg_energy(mesh, max_energy, avg_energy);
        if (max_energy <= mesh.params.stop_energy && mesh.is_input_all_inserted)
            break;

        if (mesh.params.stop_p > 0) {
            int p = get_max_p(mesh);
            cout << "p = " << p << endl;
            if (p <= mesh.params.stop_p && mesh.is_input_all_inserted)
                break;
        }

        cout << "//////////////// pass " << it << " ////////////////" << endl;
        std::array<int, 5> it_ops;
        if (it % 3 == 2)
            it_ops = {{ops[0], ops[1], ops[2], ops[3], 1}};
        else
            it_ops = {{ops[0], ops[1], ops[2], ops[3], 0}};
        operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, it_ops);

        if (it > mesh.params.max_its / 4 && max_energy > 1e3) {//Scalar check
            if (cnt_increase_epsilon > 0 && cnt_increase_epsilon == mesh.params.stage - 1) {
                mesh.params.eps += mesh.params.eps_delta;
                mesh.params.eps_2 = mesh.params.eps * mesh.params.eps;
                cnt_increase_epsilon--;
                cout << "enlarge envelope, eps = " << mesh.params.eps << endl;
            }
        }

        Scalar new_max_energy, new_avg_energy;
        get_max_avg_energy(mesh, new_max_energy, new_avg_energy);
        if (!is_just_after_update) {
            if (max_energy - new_max_energy < 5e-1 && (avg_energy - new_avg_energy) / avg_energy < 0.1) {
                is_hit_min_edge_length = update_scaling_field(mesh, new_max_energy) || is_hit_min_edge_length;
                is_just_after_update = true;
                if (cnt_increase_epsilon > 0) {
                    mesh.params.eps += mesh.params.eps_delta;
                    mesh.params.eps_2 = mesh.params.eps * mesh.params.eps;
                    cnt_increase_epsilon--;
                    cout << "enlarge envelope, eps = " << mesh.params.eps << endl;
                }
            }
        } else
            is_just_after_update = false;

        quality_queue.push_back(std::array<Scalar, 2>({{new_max_energy, new_avg_energy}}));
        if (is_hit_min_edge_length && mesh.is_input_all_inserted && it_after_al_inserted > M && it > M + N) {
            if (quality_queue[it][0] - quality_queue[it - N][0] >= SCALAR_ZERO
                && quality_queue[it][1] - quality_queue[it - N][1] >= SCALAR_ZERO)
                break;


        }
    }

    ////postprocessing
    cout << "//////////////// postprocessing ////////////////" << endl;
    for (auto &v:mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.sizing_scalar = 1;
    }
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, std::array<int, 5>({{0, 1, 0, 0, 0}}));


    ///apply sizing field
    if(mesh.params.coarsen){
        apply_coarsening(mesh, tree);
    }

    if(mesh.params.apply_sizing_field){
        apply_sizingfield(mesh, tree);
    }

}

void floatTetWild::cleanup_empty_slots(Mesh &mesh, double percentage) {
    if (mesh.tets.size() < 9e5)
        return;
    cout<<mesh.tets.size()<<" ==> ";
    ///
    const int v_end_id = mesh.tet_vertices.size() * percentage;
    const int t_end_id = mesh.tets.size() * percentage;
    //
    std::vector<int> map_v_ids(mesh.tet_vertices.size(), -1);
    int cnt = 0;
    for (int i = 0; i < v_end_id; i++) {
        if (mesh.tet_vertices[i].is_removed)
            continue;
        map_v_ids[i] = cnt++;
    }
    for (int i = v_end_id; i < mesh.tet_vertices.size(); i++) {
        map_v_ids[i] = cnt++;
    }
    //
    std::vector<int> map_t_ids(mesh.tets.size(), -1);
    cnt = 0;
    for (int i = 0; i < t_end_id; i++) {
        if (mesh.tets[i].is_removed)
            continue;
        map_t_ids[i] = cnt++;
    }
    for (int i = t_end_id; i < mesh.tets.size(); i++) {
        map_t_ids[i] = cnt++;
    }

    ///
    mesh.tet_vertices.erase(std::remove_if(mesh.tet_vertices.begin(), mesh.tet_vertices.begin() + v_end_id,
                                           [](const MeshVertex &v) { return v.is_removed; }),
                            mesh.tet_vertices.begin() + v_end_id);
    mesh.tets.erase(std::remove_if(mesh.tets.begin(), mesh.tets.begin() + t_end_id,
                                   [](const MeshTet &t) { return t.is_removed; }),
                    mesh.tets.begin() + t_end_id);

    ///
    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        for (auto &t_id: v.conn_tets)
            t_id = map_t_ids[t_id];
    }
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            t[j] = map_v_ids[t[j]];
    }
    cout<<mesh.tets.size()<<endl;
}

void floatTetWild::operation(Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops){
    Timer igl_timer;
    int v_num, t_num;
    double max_energy, avg_energy;
    double time;

    for (int i = 0; i < ops[0]; i++) {
        igl_timer.start();
        cout << "edge splitting..." << endl;
        untangle(mesh);
        edge_splitting(mesh, tree);
        time = igl_timer.getElapsedTime();
        cout << "edge splitting done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::splitting_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }

    for (int i = 0; i < ops[1]; i++) {
        igl_timer.start();
        cout << "edge collapsing..." << endl;
        untangle(mesh);
        edge_collapsing(mesh, tree);
        time = igl_timer.getElapsedTime();
        cout << "edge collapsing done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::collapsing_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }

    for (int i = 0; i < ops[2]; i++) {
        igl_timer.start();
        cout << "edge swapping..." << endl;
        untangle(mesh);
        edge_swapping(mesh);
        time = igl_timer.getElapsedTime();
        cout << "edge swapping done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::swapping_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }

    for (int i = 0; i < ops[3]; i++) {
        igl_timer.start();
        cout << "vertex smoothing..." << endl;
        vertex_smoothing(mesh, tree);
        time = igl_timer.getElapsedTime();
        cout << "vertex smoothing done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::smoothing_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }
}

void floatTetWild::operation(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
        Mesh &mesh, AABBWrapper& tree, const std::array<int, 5> &ops) {
    operation(mesh, tree, {{ops[0], ops[1], ops[2], ops[3]}});

    Timer igl_timer;

    if (!mesh.is_input_all_inserted) {
        for (int i = 0; i < ops[4]; i++) {
//            //reset boundary points
            igl_timer.start();
            insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, true);
            init(mesh, tree);
            stats().record(StateInfo::cutting_id, igl_timer.getElapsedTimeInSec(),
                           mesh.get_v_num(), mesh.get_t_num(),
                           mesh.get_max_energy(), mesh.get_avg_energy(),
                           std::count(is_face_inserted.begin(), is_face_inserted.end(),
                                      false));

            if(mesh.is_input_all_inserted && mesh.is_closed){
                for (int v_id = 0; v_id < mesh.tet_vertices.size(); v_id++) {
                    if (mesh.tet_vertices[v_id].is_removed)
                        continue;
                    mesh.tet_vertices[v_id].is_on_boundary = false;
                    mesh.tet_vertices[v_id].is_on_cut = false;
                }
            } else {
                for (int v_id = 0; v_id < mesh.tet_vertices.size(); v_id++) {
                    if (mesh.tet_vertices[v_id].is_removed)
                        continue;
                    if (!mesh.tet_vertices[v_id].is_on_boundary)
                        continue;
                    geo::index_t prev_facet;
                    if (tree.is_out_tmp_b_envelope(mesh.tet_vertices[v_id].pos, mesh.params.eps_2, prev_facet)) {
                        mesh.tet_vertices[v_id].is_on_boundary = false;
                        mesh.tet_vertices[v_id].is_on_cut = false;
                    }
                }
            }
        }

//                //check isolate boundary points
//            //todo: keep track of opp_t_ids
//                                                           mesh.get_v_num(), mesh.get_t_num(),
//                                                           mesh.get_max_energy(), mesh.get_avg_energy(),
    }
}

#include <floattetwild/geo_kd_tree.h>
bool floatTetWild::update_scaling_field(Mesh &mesh, Scalar max_energy) {

    cout << "updating sclaing field ..." << endl;
    bool is_hit_min_edge_length = false;

    Scalar radius0 = mesh.params.ideal_edge_length * 1.8;//increasing the radius would increase the #v in output

    static const Scalar stop_filter_energy = mesh.params.stop_energy * 0.8;

    Scalar filter_energy = max_energy / 100 > stop_filter_energy ? max_energy / 100 : stop_filter_energy;

    if (filter_energy > 100) {
        filter_energy = 100;
    }

    cout << "filter_energy = " << filter_energy << endl;
    Scalar recover = 1.5;
    std::vector<Scalar> scale_multipliers(mesh.tet_vertices.size(), recover);
    Scalar refine_scale = 0.5;
    Scalar min_refine_scale = mesh.params.min_edge_len_rel / mesh.params.ideal_edge_length_rel;

    const int N = -int(std::log2(min_refine_scale) - 1);
    std::vector<std::vector<int>> v_ids(N, std::vector<int>());
    for (int i = 0; i < mesh.tet_vertices.size(); i++) {
        auto &v = mesh.tet_vertices[i];
        if (v.is_removed)
            continue;

        bool is_refine = false;
        for (int t_id: v.conn_tets) {
            if (mesh.tets[t_id].quality > filter_energy)
                is_refine = true;
        }
        if (!is_refine)
            continue;

        int n = -int(std::log2(v.sizing_scalar) - 0.5);
        if (n >= N)
            n = N - 1;
        v_ids[n].push_back(i);
    }

    for (int n = 0; n < N; n++) {
        if (v_ids[n].size() == 0)
            continue;

        Scalar radius = radius0 / std::pow(2, n);

        std::unordered_set<int> is_visited;
        std::queue<int> v_queue;

        std::vector<double> pts;//geogram needs double []
        pts.reserve(v_ids[n].size() * 3);
        for (int i = 0; i < v_ids[n].size(); i++) {
            pts.push_back(mesh.tet_vertices[v_ids[n][i]].pos[0]);
            pts.push_back(mesh.tet_vertices[v_ids[n][i]].pos[1]);
            pts.push_back(mesh.tet_vertices[v_ids[n][i]].pos[2]);

            v_queue.push(v_ids[n][i]);
            is_visited.insert(v_ids[n][i]);
            scale_multipliers[v_ids[n][i]] = refine_scale;
        }
        // construct the kdtree
        geo::BalancedKdTree nnsearch(3);
        nnsearch.set_points(int(v_ids[n].size()), pts.data());

        while (!v_queue.empty()) {
            int v_id = v_queue.front();
            v_queue.pop();

            for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
                for (int j = 0; j < 4; j++) {
                    if (is_visited.find(mesh.tets[t_id][j]) != is_visited.end())
                        continue;
                    geo::index_t _;
                    double sq_dist;
                    const double p[3] = {mesh.tet_vertices[mesh.tets[t_id][j]].pos[0],
                                         mesh.tet_vertices[mesh.tets[t_id][j]].pos[1],
                                         mesh.tet_vertices[mesh.tets[t_id][j]].pos[2]};
                    nnsearch.get_nearest_neighbors(1, p, &_, &sq_dist);
                    Scalar dis = sqrt(sq_dist);

                    if (dis < radius) {
                        v_queue.push(mesh.tets[t_id][j]);
                        Scalar new_ss = (dis / radius) * (1 - refine_scale) + refine_scale;
                        if (new_ss < scale_multipliers[mesh.tets[t_id][j]])
                            scale_multipliers[mesh.tets[t_id][j]] = new_ss;
                    }
                    is_visited.insert(mesh.tets[t_id][j]);
                }
            }
        }
    }

    // update scalars
    for (int i=0;i< mesh.tet_vertices.size();i++) {
        auto& v = mesh.tet_vertices[i];
        if (v.is_removed)
            continue;
        Scalar new_scale = v.sizing_scalar * scale_multipliers[i];
        if (new_scale > 1)
            v.sizing_scalar = 1;
        else if (new_scale < min_refine_scale) {
            is_hit_min_edge_length = true;
            v.sizing_scalar = min_refine_scale;
        } else
            v.sizing_scalar = new_scale;
    }

    cout << "is_hit_min_edge_length = " << is_hit_min_edge_length << endl;
    return is_hit_min_edge_length;
}

void floatTetWild::output_info(Mesh& mesh, const AABBWrapper& tree) {
    if(mesh.params.is_quiet)
        return;

    if(mesh.params.log_level >= 2)
        return;

    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    //count
    int cnt_v = mesh.get_v_num();
    int cnt_t = mesh.get_t_num();

//    //quality


    if(mesh.params.log_level > 1) {
        return;
    }

    //euler
    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);
    std::vector<std::array<int, 3>> faces;
    for (auto &t: tets) {
        if(t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[j], t[mod4(j + 1)], t[mod4(j + 2)]}};
            std::sort(f.begin(), f.end());
            faces.push_back(f);
        }
    }
    vector_unique(faces);
    int euler = cnt_v - edges.size() + faces.size() - cnt_t;
    if (euler != 1) {
        cout << "euler error " << euler << endl;
    }

    //inversion
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++){
            if(tet_vertices[tets[i][j]].is_removed){
                cout<<"tet_vertices[tets[i][j]].is_removed"<<endl;
            }
        }
    }
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        if (is_inverted(mesh, i)) {
            cout << "tet " << i << " inverted: " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
                 << tets[i][3] << endl;
        }
    }

    //conn_tets
//            if (std::find(tet_vertices[tets[i][j]].conn_tets.begin(), tet_vertices[tets[i][j]].conn_tets.end(), i)

    //check conn_tets
    std::vector<std::vector<int>> tmp_conn_tets(tet_vertices.size());
    for (int i = 0; i < tets.size(); i++) {
        if(tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            tmp_conn_tets[tets[i][j]].push_back(i);
    }
    for (int i = 0; i < tet_vertices.size(); i++) {
        if(tet_vertices[i].is_removed)
            continue;
        std::vector<int> conn_tets = tet_vertices[i].conn_tets;
        std::sort(conn_tets.begin(), conn_tets.end());
        if (conn_tets != tmp_conn_tets[i]) {
            cout << "conn_tets error" << endl;
            for (auto &ii:tet_vertices[i].conn_tets)
                cout << ii << " ";
            cout << endl;
            for (auto &ii:tmp_conn_tets[i])
                cout << ii << " ";
            cout << endl;
        }
    }


    //surface tracking
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;

        for (int j = 0; j < 4; j++) {
            auto &t = tets[i];
            int opp_t_id = get_opp_t_id(mesh, i, j);
            if (opp_t_id >= 0) {
                int k = get_local_f_id(opp_t_id, t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4], mesh);
                if (tets[opp_t_id].is_surface_fs[k] == NOT_SURFACE && tets[i].is_surface_fs[j] == NOT_SURFACE);
                else if (tets[opp_t_id].is_surface_fs[k] + tets[i].is_surface_fs[j] == 0);
                else
                    cout << "surface faces are not matched" << endl;
            }
        }

        for (int j = 0; j < 4; j++) {
            if (tets[i].is_surface_fs[j] != NOT_SURFACE && tets[i].is_bbox_fs[j] != NOT_BBOX){
                cout<<"tets[i].is_surface_fs[j] != NOT_SURFACE && tets[i].is_bbox_fs[j] != NOT_BBOX"<<endl;
                cout<<i<<" "<<j<<endl;
            }
            if (tets[i].is_surface_fs[j] != NOT_SURFACE) {

                for (int k = 0; k < 3; k++) {
                    if (!tet_vertices[tets[i][(j + 1 + k) % 4]].is_on_surface) {
                        cout << "is_surface_fs error" << endl;
                        cout << "t " << i << ": " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
                             << tets[i][3] << endl;
                        cout << tets[i].is_surface_fs[0] << " " << tets[i].is_surface_fs[1] << " "
                             << tets[i].is_surface_fs[2] << " " << tets[i].is_surface_fs[3] << endl;
                        cout << tet_vertices[tets[i][0]].is_on_surface << " " << tet_vertices[tets[i][1]].is_on_surface
                             << " "
                             << tet_vertices[tets[i][2]].is_on_surface << " " << tet_vertices[tets[i][3]].is_on_surface
                             << endl;
                    }
                }
            }

            if (tets[i].is_bbox_fs[j] != NOT_BBOX) {
                for (int k = 0; k < 3; k++) {
                    if (!tet_vertices[tets[i][(j + 1 + k) % 4]].is_on_bbox) {
                        cout<<"is_bbox_fs error"<<endl;
                        cout << "t " << i << ": " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
                             << tets[i][3] << endl;
                        cout << (int)tets[i].is_bbox_fs[0] << " " << (int)tets[i].is_bbox_fs[1] << " "
                             << (int)tets[i].is_bbox_fs[2] << " " << (int)tets[i].is_bbox_fs[3] << endl;
                        cout << tet_vertices[tets[i][0]].is_on_bbox << " " << tet_vertices[tets[i][1]].is_on_bbox << " "
                             << tet_vertices[tets[i][2]].is_on_bbox << " " << tet_vertices[tets[i][3]].is_on_bbox
                             << endl;
                    }
                }
                if(get_opp_t_id(mesh, i, j) != OPP_T_ID_BOUNDARY){
                    cout<<"wrong-marked bbox face"<<endl;
                }
            } else {
                if(get_opp_t_id(mesh, i, j) == OPP_T_ID_BOUNDARY){
                    cout<<"unmarked bbox face"<<endl;
                }
            }
        }
    }

    for(int i=0;i<tet_vertices.size();i++) {
        if (tet_vertices[i].is_removed)
            continue;

        if(tet_vertices[i].is_on_bbox && tet_vertices[i].is_on_surface){
            cout<<"error tet_vertices[i].is_on_bbox && tet_vertices[i].is_on_surface"<<endl;
            cout<<"v "<<i<<endl;
        }

        if (tet_vertices[i].is_on_bbox) {
            bool is_found = false;
            for (int t_id:tet_vertices[i].conn_tets) {
                int j = tets[t_id].find(i);
                for (int k = 0; k < 3; k++) {
                    if (tets[t_id].is_bbox_fs[(j + 1 + k) % 4] != NOT_BBOX) {
                        is_found = true;
                        break;
                    }
                }
            }
            if (!is_found) {
                cout << "is_on_bbox error" << endl;
                for (int t_id:tet_vertices[i].conn_tets) {
                    cout << "t " << t_id << ": " << tets[t_id][0] << " " << tets[t_id][1] << " " << tets[t_id][2] << " "
                         << tets[t_id][3] << endl;
                    cout << tets[t_id].is_bbox_fs[0] << " " << tets[t_id].is_bbox_fs[1] << " "
                         << tets[t_id].is_bbox_fs[2] << " " << tets[t_id].is_bbox_fs[3] << endl;
                    cout << tet_vertices[tets[t_id][0]].is_on_bbox << " " << tet_vertices[tets[t_id][1]].is_on_bbox << " "
                         << tet_vertices[tets[t_id][2]].is_on_bbox << " " << tet_vertices[tets[t_id][3]].is_on_bbox
                         << endl;
                }
            }
        }
        if (tet_vertices[i].is_on_surface) {
            bool is_found = false;
            for (int t_id:tet_vertices[i].conn_tets) {
                int j = tets[t_id].find(i);
                for (int k = 0; k < 3; k++) {
                    if (tets[t_id].is_bbox_fs[(j + 1 + k) % 4] != NOT_SURFACE) {
                        is_found = true;
                        break;
                    }
                }
            }
            if (!is_found) {
                cout << "is_on_surface error" << endl;
                for (int t_id:tet_vertices[i].conn_tets) {
                    cout << "t " << t_id << ": " << tets[t_id][0] << " " << tets[t_id][1] << " " << tets[t_id][2] << " "
                         << tets[t_id][3] << endl;
                    cout << tets[t_id].is_surface_fs[0] << " " << tets[t_id].is_surface_fs[1] << " "
                         << tets[t_id].is_surface_fs[2] << " " << tets[t_id].is_surface_fs[3] << endl;
                    cout << tet_vertices[tets[t_id][0]].is_on_surface << " " << tet_vertices[tets[t_id][1]].is_on_surface
                         << " "
                         << tet_vertices[tets[t_id][2]].is_on_surface << " " << tet_vertices[tets[t_id][3]].is_on_surface
                         << endl;
                }
            }
        }
    }
    cout<<endl;


    return;
}

void floatTetWild::check_envelope(Mesh& mesh, const AABBWrapper& tree) {//for debug only

    Scalar check_eps = mesh.params.eps_2;

    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] <= 0) {
                std::vector<geo::vec3> ps;
                sample_triangle({{mesh.tet_vertices[t[(j + 1) % 4]].pos, mesh.tet_vertices[t[(j + 2) % 4]].pos,
                                         mesh.tet_vertices[t[(j + 3) % 4]].pos}}, ps, mesh.params.dd);
                if(tree.is_out_sf_envelope(ps, mesh.params.eps_2)){
                    cout << "out of envelope!" << endl;
                }
            }
        }
    }

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed || !v.is_on_surface)
            continue;

        std::vector<geo::vec3> ps = {geo::vec3(v.pos[0], v.pos[1], v.pos[2])};
        Scalar d = tree.dist_sf_envelope(ps, check_eps);
        if (d > mesh.params.eps_2) {
            cout << "v out of envelope!" << endl;
            cout << d << ", eps_input = " << check_eps << endl;
        }
    }

    cout<<"envelope check done"<<endl;
}

int floatTetWild::get_max_p(const Mesh &mesh)
{
    const Scalar scaling = 1.0 / (mesh.params.bbox_max - mesh.params.bbox_min).maxCoeff();
    const Scalar B = 3;
    const int p_ref = 1;

    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);
    Scalar h_ref = 0;
    for(const auto &e : edges){
        const Vector3 edge = (mesh.tet_vertices[e[0]].pos - mesh.tet_vertices[e[1]].pos)*scaling;
        h_ref += edge.norm();
    }
    h_ref /= edges.size();

    const Scalar rho_ref = sqrt(6.)/12.*h_ref;
    const Scalar sigma_ref = rho_ref / h_ref;

    int max_p = 1;

    for(const auto &t : mesh.tets)
    {
        if(t.is_removed)
            continue;

        const auto &v0 = mesh.tet_vertices[t[0]].pos;
        const auto &v1 = mesh.tet_vertices[t[1]].pos;
        const auto &v2 = mesh.tet_vertices[t[2]].pos;
        const auto &v3 = mesh.tet_vertices[t[3]].pos;

        std::array<Vector3, 6> e;
        e[0] = (v0 - v1) * scaling;
        e[1] = (v1 - v2) * scaling;
        e[2] = (v2 - v0) * scaling;

        e[3] = (v0 - v3) * scaling;
        e[4] = (v1 - v3) * scaling;
        e[5] = (v2 - v3) * scaling;

        const Scalar S = (e[0].cross(e[1]).norm() + e[0].cross(e[4]).norm() + e[4].cross(e[1]).norm() + e[2].cross(e[5]).norm()) / 2;
        const Scalar V = std::abs(e[3].dot(e[2].cross(-e[0])))/6;
        const Scalar rho = 3 * V / S;
        Scalar h = e[0].norm();
        for (int i = 1; i < 6; i++)
            h = std::max(h, e[i].norm());

        const Scalar sigma = rho / h;

        const Scalar ptmp = (std::log(B*std::pow(h_ref, p_ref + 1)*sigma*sigma/sigma_ref/sigma_ref) - std::log(h))/std::log(h);
        const int p = (int)std::round(ptmp);
        max_p = std::max(p, max_p);
    }

    return max_p;
}

#include <floattetwild/Predicates.hpp>

//void floatTetWild::apply_sizingfield(const Eigen::VectorXd& V_in, const Eigen::VectorXi& T_in, const Eigen::VectorXd& values,
//        //compute barycenter
//            && max_energy < mesh.params.stop_energy) //refinement and quality enough

void floatTetWild::apply_sizingfield(Mesh& mesh, AABBWrapper& tree) {
    logger().debug("Applying sizing field...");

    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    for (auto &p: tet_vertices) {
        if (p.is_removed)
            continue;
        p.sizing_scalar = 1; //reset
        double value = mesh.params.get_sizing_field_value(p.pos);
        if (value > 0) {
            p.sizing_scalar = value / mesh.params.ideal_edge_length;
        }
    }

    int num_tets = mesh.get_t_num();
    for (int i = 0; i < 20; i++) {
        operation(mesh, tree);
        double tmp_num_tets = mesh.get_t_num();
        double max_energy = mesh.get_max_energy();
        cout<<"/////////"<<i<<" "<<max_energy<<endl;
        if ((tmp_num_tets - num_tets) / num_tets < 0.02
            && max_energy < mesh.params.stop_energy) //refinement and quality enough
            break;
        num_tets = tmp_num_tets;
    }
}

void floatTetWild::apply_coarsening(Mesh& mesh, AABBWrapper& tree) {
    mesh.is_coarsening = true;

    for (auto &v:mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.sizing_scalar = 1;
    }

    int tets_size = mesh.get_t_num();
    int stop_size = tets_size * 0.001;
    for (int i = 0; i < 20; i++) {
        operation(mesh, tree, {{0, 1, 1, 0}});
        int new_size = mesh.get_t_num();
        if (abs(new_size - tets_size) < stop_size)
            break;
        tets_size = new_size;
    }

    mesh.is_coarsening = false;
}

#include <floattetwild/bfs_orient.h>
#include <floattetwild/unique_rows.h>
#include <floattetwild/remove_duplicate_vertices.h>
#include <floattetwild/TriangleInsertion.h>
void floatTetWild::get_tracked_surface(Mesh& mesh, MatrixXs &V_sf, MatrixXi &F_sf, int c_id) {
#define SF_CONDITION t.is_surface_fs[j]<=0&&t.surface_tags[j]==c_id

    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    int cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (SF_CONDITION)
                cnt++;
        }
    }

    V_sf.resize(cnt * 3, 3);
    F_sf.resize(cnt, 3);
    cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (SF_CONDITION) {
                for (int k = 0; k < 3; k++)
                    V_sf.row(cnt * 3 + k) = tet_vertices[t[mod4(j + k + 1)]].pos;
                if (Predicates::orient_3d(tet_vertices[t[mod4(j + 1)]].pos, tet_vertices[t[mod4(j + 2)]].pos,
                                          tet_vertices[t[mod4(j + 3)]].pos, tet_vertices[t[j]].pos) ==
                    Predicates::ORI_POSITIVE)
                    F_sf.row(cnt) << cnt * 3, cnt * 3 + 2, cnt * 3 + 1;
                else
                    F_sf.row(cnt) << cnt * 3, cnt * 3 + 1, cnt * 3 + 2;
                cnt++;
            }
        }
    }

    MatrixXd V;
    MatrixXi F;
    VectorXi _1, _2;
    remove_duplicate_vertices(V_sf, F_sf, -1, V, _1, _2, F);
    V_sf = V;
    F_sf.resize(0, 3);
    bfs_orient(F, F_sf, _1);
}

void floatTetWild::correct_tracked_surface_orientation(Mesh &mesh, AABBWrapper& tree){
    std::vector<std::array<bool, 4>> is_visited(mesh.tets.size(), {{false, false, false, false}});
    for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
        auto &t = mesh.tets[t_id];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] == NOT_SURFACE || is_visited[t_id][j])
                continue;
            is_visited[t_id][j] = true;
            int opp_t_id = get_opp_t_id(t_id, j, mesh);
            if (opp_t_id < 0) {
                t.is_surface_fs[j] = NOT_SURFACE;
                continue;
            }
            int k = get_local_f_id(opp_t_id, t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4], mesh);
            is_visited[opp_t_id][k] = true;
            //
            Vector3 c = (mesh.tet_vertices[t[(j + 1) % 4]].pos + mesh.tet_vertices[t[(j + 2) % 4]].pos +
                         mesh.tet_vertices[t[(j + 3) % 4]].pos) / 3;
            int f_id = tree.get_nearest_face_sf(c);
            const auto &fv1 = tree.sf_mesh.vertices.point(tree.sf_mesh.facets.vertex(f_id, 0));
            const auto &fv2 = tree.sf_mesh.vertices.point(tree.sf_mesh.facets.vertex(f_id, 1));
            const auto &fv3 = tree.sf_mesh.vertices.point(tree.sf_mesh.facets.vertex(f_id, 2));
            auto nf = geo::cross((fv2 - fv1), (fv3 - fv1));
            Vector3 n, nt;
            n << nf[0], nf[1], nf[2];
            //
            auto &tv1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
            auto &tv2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
            auto &tv3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
            if (Predicates::orient_3d(tv1, tv2, tv3, mesh.tet_vertices[mesh.tets[t_id][j]].pos)
                == Predicates::ORI_POSITIVE)
                nt = (tv2 - tv1).cross(tv3 - tv1);
            else
                nt = (tv3 - tv1).cross(tv2 - tv1);
            //
            if (nt.dot(n) > 0) {
                t.is_surface_fs[j] = 1;
                mesh.tets[opp_t_id].is_surface_fs[k] = -1;
            } else {
                t.is_surface_fs[j] = -1;
                mesh.tets[opp_t_id].is_surface_fs[k] = 1;
            }
        }
    }
}


void floatTetWild::boolean_operation(Mesh&                                     mesh,
                                    const CSGTree&                            csg_tree_with_ids,
                                    const std::vector<std::vector<Vector3>>&  Vs,
                                    const std::vector<std::vector<Vector3i>>& Fs)
{
    MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }

    int max_id = CSGTreeParser::get_max_id(csg_tree_with_ids);
    std::vector<VectorXd> w(max_id + 1);

    if(Vs.empty())
    {
        MatrixXs vs;
        MatrixXi fs;


        for (int i = 0; i <= max_id; ++i) {
            get_tracked_surface(mesh, vs, fs, i);

            fast_winding_number(MatrixXd(vs.cast<double>()), MatrixXi(fs), C, w[i]);
        }
    }
    else {
        for (int i = 0; i <= max_id; ++i) {
            MatrixXs vs(Vs[i].size(), 3);
            MatrixXi fs(Fs[i].size(), 3);
            for (int k = 0; k < vs.rows(); ++k)
                vs.row(k) = Vs[i][k];
            for (int k = 0; k < fs.rows(); ++k)
                fs.row(k) = Fs[i][k];

            fast_winding_number(MatrixXd(vs.cast<double>()), MatrixXi(fs), C, w[i]);
            }
        }

    boolean_operation(mesh, csg_tree_with_ids, w);
}

void floatTetWild::boolean_operation(Mesh& mesh, const CSGTree& csg_tree_with_ids){
    boolean_operation(mesh,
                      csg_tree_with_ids,
                      std::vector<std::vector<Vector3>>(),
                      std::vector<std::vector<Vector3i>>());
}

void floatTetWild::boolean_operation(Mesh& mesh, const CSGTree& csg_tree_with_ids, const std::vector<VectorXd> &w)
{
    int max_id = CSGTreeParser::get_max_id(csg_tree_with_ids);

    int cnt = 0;
    for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
        auto &t = mesh.tets[t_id];
        if(t.is_removed)
            continue;
        
        bool keep = CSGTreeParser::keep_tet(csg_tree_with_ids, cnt, w);
        t.is_removed = !keep;
        int tid = 0;
        for (int id = 0; id <= max_id; ++id) {
            bool inside = w[id][cnt] > 0.5;
            if(inside)
                tid = std::max(id+1, tid);
        }
        t.scalar = tid;
        cnt++;
        }

}

void floatTetWild::boolean_operation(Mesh& mesh, int op){
    const int OP_UNION = 0;
    const int OP_INTERSECTION = 1;
    const int OP_DIFFERENCE = 2;

    MatrixXs v1;
    MatrixXi f1;
    get_tracked_surface(mesh, v1, f1, 1);

    MatrixXs v2;
    MatrixXi f2;
    get_tracked_surface(mesh, v2, f2, 2);

    MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }

    VectorXd w1, w2;
    fast_winding_number(MatrixXd(v1.cast<double>()), MatrixXi(f1), C, w1);
    fast_winding_number(MatrixXd(v2.cast<double>()), MatrixXi(f2), C, w2);


    int cnt = 0;
    for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
        auto &t = mesh.tets[t_id];
        if(t.is_removed)
            continue;
        switch(op){
            case OP_UNION:
                if(w1(cnt)<=0.5 && w2(cnt)<=0.5)
                    t.is_removed = true;
                break;
            case OP_INTERSECTION:
                if(w1(cnt)<=0.5 || w2(cnt)<=0.5)
                    t.is_removed = true;
                break;
            case OP_DIFFERENCE:
                if(w1(cnt)<=0.5 || w2(cnt)>0.5)
                    t.is_removed = true;
                break;
            default:
                break;
        }
        cnt++;
    }

}

void floatTetWild::filter_outside(Mesh& mesh, const MatrixXs& V,
                                  const MatrixXi& F_in, bool invert_faces) {
    MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }

    MatrixXi F = F_in;
    VectorXd W;
    if (invert_faces) {
        for (int i = 0; i < F.rows(); i++)
            std::swap(F(i, 1), F(i, 2));
    }
    fast_winding_number(MatrixXd(V.cast<double>()), MatrixXi(F), C, W);

    index = 0;
    int n_tets = 0;
    std::vector<bool> old_flags(mesh.tets.size());
    for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
        auto &t = mesh.tets[t_id];
        old_flags[t_id] = t.is_removed;

        if (t.is_removed)
            continue;
        if (W(index) <= 0.5) {
            t.is_removed = true;
        } else
            n_tets++;
        index++;
    }

    if (n_tets <= 0) {
        if (invert_faces)
            logger().error("Empty mesh, problem with inverted faces");
        else {
            for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
                auto &t = mesh.tets[t_id];
                t.is_removed = old_flags[t_id];
            }
            logger().debug("Empty mesh trying to reverse the faces");
            filter_outside(mesh, V, F_in, true);
        }
    }

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        bool is_remove = true;
        for (int t_id: v.conn_tets) {
            if (!mesh.tets[t_id].is_removed) {
                is_remove = false;
                break;
            }
        }
        v.is_removed = is_remove;
    }
}

void floatTetWild::filter_outside(Mesh& mesh, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces){
    MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }

    MatrixXs V(input_vertices.size(), 3);
    MatrixXi F(input_faces.size(), 3);
    for(int i=0;i<input_vertices.size();i++)
        V.row(i) = input_vertices[i];
    for(int i=0;i<input_faces.size();i++)
        F.row(i) = input_faces[i];
    ///
    VectorXd W;
    fast_winding_number(MatrixXd(V.cast<double>()), MatrixXi(F), C, W);

    index = 0;
    std::vector<bool> old_flags(mesh.tets.size());
    for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
        auto &t = mesh.tets[t_id];
        old_flags[t_id] = t.is_removed;

        if (t.is_removed)
            continue;
        if (W(index) <= 0.5) {
            t.is_removed = true;
        } else
        index++;
    }


    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        bool is_remove = true;
        for (int t_id: v.conn_tets) {
            if (!mesh.tets[t_id].is_removed) {
                is_remove = false;
                break;
            }
        }
        v.is_removed = is_remove;
    }
}

void floatTetWild::filter_outside_floodfill(Mesh& mesh, bool invert_faces) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    std::queue<int> t_queue;
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (tets[i].is_bbox_fs[j] != NOT_BBOX) {
                t_queue.push(i);
                tets[i].is_removed = true;
                break;
            }
        }
    }

    while (!t_queue.empty()) {
        int t_id = t_queue.front();
        t_queue.pop();

        for (int j = 0; j < 4; j++) {
            if (tets[t_id].is_bbox_fs[j] != NOT_BBOX || tets[t_id].is_surface_fs[j] != NOT_SURFACE)
                continue;
            int n_t_id = get_opp_t_id(mesh, t_id, j);
            if (n_t_id < 0 || tets[n_t_id].is_removed)
                continue;
            tets[n_t_id].is_removed = true;
            t_queue.push(n_t_id);
        }
        for (int j = 0; j < 4; j++) {
            vector_erase(tet_vertices[tets[t_id][j]].conn_tets, t_id);
        }
    }

    for (int i = 0; i < tet_vertices.size(); i++) {
        if (tet_vertices[i].is_removed)
            continue;
        if (tet_vertices[i].conn_tets.empty())
            tet_vertices[i].is_removed = true;
    }
}

void floatTetWild::mark_outside(Mesh& mesh, bool invert_faces){
    MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }

    MatrixXs V;
    MatrixXi F;
    get_tracked_surface(mesh, V, F);
    if(invert_faces){
        for (int i = 0; i < F.rows(); i++)
            std::swap(F(i, 1), F(i, 2));
    }
    VectorXd W;
    fast_winding_number(MatrixXd(V.cast<double>()), MatrixXi(F), C, W);

    index = 0;
    int n_tets = 0;
    for (auto& t: mesh.tets) {
        if (t.is_removed){
            t.is_outside = true;
            continue;
        }
        if (W(index) <= 0.5)
            t.is_outside = true;
        else
            n_tets++;
        index++;
    }

    if(n_tets <= 0)
    {
        if(invert_faces)
            logger().error("Empty mesh, problem with inverted faces");
        else{
            for (auto& t: mesh.tets)
            {
                t.is_outside = false;
            }
            logger().debug("Empty mesh trying to reverse the faces");
            mark_outside(mesh, true);
        }
    }

    for (auto& v: mesh.tet_vertices) {
        if (v.is_removed){
            v.is_outside = true;
            continue;
        }
        bool is_outside = true;
        for(int t_id: v.conn_tets) {
            if (!mesh.tets[t_id].is_outside) {
                is_outside = false;
                break;
            }
        }
        v.is_outside = is_outside;
    }
}

void floatTetWild::untangle(Mesh &mesh) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;
    static const Scalar zero_area = 1e2 * SCALAR_ZERO_2;
    static const std::vector<std::array<int, 4>> face_pairs = {{{0, 1, 2, 3}},
                                                               {{0, 2, 1, 3}},
                                                               {{0, 3, 1, 2}}};


    int cnt = 0;
    for (int t_id = 0; t_id < tets.size(); t_id++) {
        auto &t = tets[t_id];
        if (t.is_removed)
            continue;
        if (t.quality < 1e10)
            continue;
        int cnt_on_surface = 0;
        bool has_degenerate_face = false;
        std::array<double, 4> areas;
        double max_area = 0;
        int max_j = -1;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] != NOT_SURFACE)
                cnt_on_surface++;
            areas[j] = get_area(tet_vertices[t[(j + 1) % 4]].pos,
                                tet_vertices[t[(j + 2) % 4]].pos,
                                tet_vertices[t[(j + 3) % 4]].pos);
            if (areas[j] < zero_area)
                has_degenerate_face = true;
            if (areas[j] > max_area) {
                max_area = areas[j];
                max_j = j;
            }
        }
        if (cnt_on_surface == 0)
            continue;

//
//

        if (has_degenerate_face) {
            if (t.is_surface_fs[max_j] != NOT_SURFACE && max_area > zero_area) {
                for (int j = 0; j < 4; j++) {
                    if (j != max_j) {
                        t.is_surface_fs[j] = NOT_SURFACE;
                        int opp_t_id = get_opp_t_id(mesh, t_id, j);
                        if (opp_t_id >= 0) {
                            int k = get_local_f_id(opp_t_id, t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4], mesh);
                            tets[opp_t_id].is_surface_fs[k] = NOT_SURFACE;
                        }
                    }
                }
                cnt++;
            }
        } else {
            if (cnt_on_surface < 2)
                continue;
            if (std::abs(max_area - areas[(max_j + 1) % 4] - areas[(max_j + 2) % 4] - areas[(max_j + 3) % 4]) <
                zero_area) {
                if (t.is_surface_fs[max_j] == NOT_SURFACE)
                    continue;
                for (int j = 0; j < 4; j++) {
                    if (j != max_j && t.is_surface_fs[j] != NOT_SURFACE) {
                        t.is_surface_fs[j] = NOT_SURFACE;
                        int opp_t_id = get_opp_t_id(mesh, t_id, j);
                        if (opp_t_id >= 0) {
                            int k = get_local_f_id(opp_t_id, t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4], mesh);
                            tets[opp_t_id].is_surface_fs[k] = NOT_SURFACE;
                        }
                    }
                }
                cnt++;
            } else {
                for (const auto &fp: face_pairs) {
                    std::array<Vector3, 2> ns;
                    auto &p1 = tet_vertices[tets[t_id][fp[2]]].pos;
                    auto &p2 = tet_vertices[tets[t_id][fp[3]]].pos;
                    Vector3 v = (p2 - p1).normalized();
                    for (int j = 0; j < 2; j++) {
                        auto &p = tet_vertices[tets[t_id][fp[j]]].pos;
                        Vector3 q = p1 + ((p - p1).dot(v)) * v;
                        ns[j] = p - q;
                    }
                    if (ns[0].dot(ns[1]) > 0)
                        continue;

                    if (std::abs(areas[fp[0]] + areas[fp[1]] - areas[fp[2]] - areas[fp[3]]) > zero_area)
                        continue;

                    std::array<int, 2> js = {{-1, -1}};
                    if (t.is_surface_fs[fp[0]] != NOT_SURFACE && t.is_surface_fs[fp[1]] != NOT_SURFACE)
                        js = {{fp[2], fp[3]}};
                    else if (t.is_surface_fs[fp[2]] != NOT_SURFACE && t.is_surface_fs[fp[3]] != NOT_SURFACE)
                        js = {{fp[0], fp[1]}};

                    for (int j: js) {
                        if (j < 0)
                            continue;
                        if (t.is_surface_fs[j] == NOT_SURFACE)
                            continue;
                        t.is_surface_fs[j] = NOT_SURFACE;
                        int opp_t_id = get_opp_t_id(mesh, t_id, j);
                        if (opp_t_id >= 0) {
                            int k = get_local_f_id(opp_t_id, t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4], mesh);
                            tets[opp_t_id].is_surface_fs[k] = NOT_SURFACE;
                        }
                    }
                    if (js[0] >= 0)
                        cnt++;
                    break;
                }
            }
        }

//
//
    }
    cout << "fixed " + std::to_string(cnt) + " tangled element" << endl;
}

void floatTetWild::smooth_open_boundary(Mesh& mesh, const AABBWrapper& tree) {
    mark_outside(mesh);
    smooth_open_boundary_aux(mesh, tree);

    return;
    for(int i=0;i<10;i++) {
        mark_outside(mesh);
        smooth_open_boundary_aux(mesh, tree);
        for(auto& t: mesh.tets)
            t.is_outside = false;
    }
    mark_outside(mesh);
}

void floatTetWild::smooth_open_boundary_aux(Mesh& mesh, const AABBWrapper& tree) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    std::vector<std::array<int, 3>> faces;
    for (auto &t: tets) {
        if (t.is_outside)
            continue;
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[j], t[(j + 1) % 4], t[(j + 2) % 4]}};
            std::sort(f.begin(), f.end());
            faces.push_back(f);
        }
    }
    std::sort(faces.begin(), faces.end());
    if (faces.empty())
        return;

    bool has_open_boundary = false;
    std::vector<bool> is_b_vs(tet_vertices.size(), false);
    std::vector<std::vector<int>> conn_b_fs(tet_vertices.size());
    bool is_boundary = true;
    for (int i = 0; i < faces.size() - 1; i++) {
        if (faces[i] == faces[i + 1]) {
            is_boundary = false;
        } else {
            if (is_boundary) {
                has_open_boundary = true;
                for (int j = 0; j < 3; j++) {
                    if (!tet_vertices[faces[i][j]].is_on_surface) {
                        conn_b_fs[faces[i][j]].push_back(i);
                        is_b_vs[faces[i][j]] = true;
                    }
                }
            }
            is_boundary = true;
        }
    }
    if (is_boundary) {
        has_open_boundary = true;
        for (int j = 0; j < 3; j++) {
            if (!tet_vertices[faces.back()[j]].is_on_surface) {
                conn_b_fs[faces.back()[j]].push_back(faces.size() - 1);
                is_b_vs[faces.back()[j]] = true;
            }
        }
    }
    if (!has_open_boundary)
        return;

    const int IT = 8;
    for (int it = 0; it < IT; it++) {
        ///laplacian
        int cnt = 0;
        int cnt_s = 0;
        for (int v_id = 0; v_id < tet_vertices.size(); v_id++) {
            if (conn_b_fs[v_id].empty())
                continue;

            tet_vertices[v_id].is_freezed = true;

            cnt++;
            std::vector<int> n_v_ids;
            for (auto &f_id: conn_b_fs[v_id]) {
                for (int j = 0; j < 3; j++) {
                    if (faces[f_id][j] != v_id)
                        n_v_ids.push_back(faces[f_id][j]);
                }
            }
            vector_unique(n_v_ids);

            Vector3 c(0, 0, 0);
            for (int n_v_id: n_v_ids) {
                c += tet_vertices[n_v_id].pos;
            }
            c /= n_v_ids.size();

            double dis = (c - tet_vertices[v_id].pos).norm();
            Vector3 v = (c - tet_vertices[v_id].pos).normalized();
            static const int N = 7;
            Vector3 p;
            for (int n = 0; n < N; n++) {
                p = tet_vertices[v_id].pos + dis / pow(2, n) * v;
                bool is_valid = true;
                for (int t_id: tet_vertices[v_id].conn_tets) {
                    int j = tets[t_id].find(v_id);
                    if (is_inverted(mesh, t_id, j, p)) {
                        is_valid = false;
                        break;
                    }
                    if (get_quality(p, tet_vertices[tets[t_id][(j + 1) % 4]].pos,
                                    tet_vertices[tets[t_id][(j + 2) % 4]].pos,
                                    tet_vertices[tets[t_id][(j + 3) % 4]].pos) > mesh.params.stop_energy) {
                        is_valid = false;
                        break;
                    }
                }
                if (!is_valid)
                    continue;

                cnt_s++;
                tet_vertices[v_id].pos = p;
                for (int t_id: tet_vertices[v_id].conn_tets)
                    tets[t_id].quality = get_quality(mesh, t_id);
                break;
            }
        }
        cout<<cnt<<"/"<<cnt_s<<endl;

        ///regular optimization
        for(auto& v: tet_vertices){
            if(v.is_on_surface)
                v.is_freezed = true;
        }
        edge_collapsing(mesh, tree);
        vertex_smoothing(mesh, tree);

        ///unfreeze
        for (int v_id=0; v_id < tet_vertices.size(); v_id++) {
                tet_vertices[v_id].is_freezed = false;
        }
    }
}

void floatTetWild::manifold_edges(Mesh& mesh) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    auto split = [&](int v1_id, int v2_id, std::vector<int> &old_t_ids, std::vector<int> &new_t_ids) {
        ////create new vertex
        MeshVertex new_v;
        new_v.pos = (tet_vertices[v1_id].pos + tet_vertices[v2_id].pos) / 2;
        bool is_found = false;
        for (int i = mesh.v_empty_start; i < tet_vertices.size(); i++) {
            mesh.v_empty_start = i;
            if (tet_vertices[i].is_removed) {
                is_found = true;
                break;
            }
        }
        if (!is_found)
            mesh.v_empty_start = tet_vertices.size();

        int v_id = mesh.v_empty_start;
        if (v_id < tet_vertices.size())
            tet_vertices[v_id] = new_v;
        else
            tet_vertices.push_back(new_v);


        ////check inversion
        set_intersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, old_t_ids);
        for (int t_id: old_t_ids) {
            for (int j = 0; j < 4; j++) {
                if (tets[t_id][j] == v1_id || tets[t_id][j] == v2_id) {
                    if (is_inverted(mesh, t_id, j, new_v.pos)) {
                        tet_vertices[v_id].is_removed = true;
                        return -1;
                    }
                }
            }
        }

        ////real update
        //update tets
        get_new_tet_slots(mesh, old_t_ids.size(), new_t_ids);
        for (int t_id: new_t_ids)
            tets[t_id].reset();

        //update indices & tags
        for (int i = 0; i < old_t_ids.size(); i++) {
            tets[new_t_ids[i]] = tets[old_t_ids[i]];
            for (int j = 0; j < 4; j++) {
                if (tets[old_t_ids[i]][j] == v1_id)
                    tets[old_t_ids[i]][j] = v_id;

                if (tets[new_t_ids[i]][j] == v2_id)
                    tets[new_t_ids[i]][j] = v_id;
            }
            //update quality
            tets[new_t_ids[i]].quality = get_quality(mesh, new_t_ids[i]);
            tets[old_t_ids[i]].quality = get_quality(mesh, old_t_ids[i]);
        }

        tet_vertices[v_id].conn_tets.insert(tet_vertices[v_id].conn_tets.end(), old_t_ids.begin(), old_t_ids.end());
        tet_vertices[v_id].conn_tets.insert(tet_vertices[v_id].conn_tets.end(), new_t_ids.begin(), new_t_ids.end());
        for (int i = 0; i < old_t_ids.size(); i++) {
            for (int j = 0; j < 4; j++) {
                if (tets[old_t_ids[i]][j] != v_id && tets[old_t_ids[i]][j] != v2_id)
                    tet_vertices[tets[old_t_ids[i]][j]].conn_tets.push_back(new_t_ids[i]);
            }
            tet_vertices[v1_id].conn_tets.erase(
                    std::find(tet_vertices[v1_id].conn_tets.begin(), tet_vertices[v1_id].conn_tets.end(),
                              old_t_ids[i]));
            tet_vertices[v1_id].conn_tets.push_back(new_t_ids[i]);
        }

        return v_id;
    };


    std::vector<std::array<int, 3>> faces;
    for (auto &t: tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[j], t[(j + 1) % 4], t[(j + 2) % 4]}};
            std::sort(f.begin(), f.end());
            faces.push_back(f);
        }
    }
    std::sort(faces.begin(), faces.end());
    if (faces.empty())
        return;

    ///
    std::vector<std::array<int, 3>> b_faces;
    bool is_boundary = true;
    for (int i = 0; i < faces.size() - 1; i++) {
        if (faces[i] == faces[i + 1]) {
            is_boundary = false;
        } else {
            if (is_boundary) {
                b_faces.push_back(faces[i]);
            }
            is_boundary = true;
        }
    }
    if (is_boundary) {
        b_faces.push_back(faces.back());
    }

    ///
    std::vector<std::array<int, 2>> b_edges;
    for (int i = 0; i < b_faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            if (b_faces[i][j] < b_faces[i][(j + 1) % 3])
                b_edges.push_back({{b_faces[i][j], b_faces[i][(j + 1) % 3]}});
            else
                b_edges.push_back({{b_faces[i][(j + 1) % 3], b_faces[i][j]}});
        }
    }
    vector_unique(b_edges);

    ///
    std::queue<std::array<int, 2>> edge_queue;
    for (auto &e: b_edges)
        edge_queue.push(e);

    while (!edge_queue.empty()) {
        auto e = edge_queue.front();
        edge_queue.pop();

        std::vector<int> n_t_ids;
        set_intersection(tet_vertices[e[0]].conn_tets, tet_vertices[e[1]].conn_tets, n_t_ids);
        if (n_t_ids.empty())
            continue;

        std::map<int, bool> is_visited;
        for (int t_id: n_t_ids) {
            is_visited[t_id] = false;
        }

        std::vector<std::vector<int>> tet_groups;
        for (int t_id: n_t_ids) {
            if (is_visited.find(t_id) == is_visited.end())
                continue;
            if (is_visited[t_id])
                continue;
            is_visited[t_id] = true;

            tet_groups.emplace_back();
            std::queue<int> tet_queue;
            tet_queue.push(t_id);
            while (!tet_queue.empty()) {
                int t0_id = tet_queue.front();
                tet_queue.pop();
                tet_groups.back().push_back(t0_id);

                for (int j = 0; j < 4; j++) {
                    if (tets[t0_id][j] == e[0] || tets[t0_id][j] == e[1])
                        continue;
                    int opp_t_id = get_opp_t_id(mesh, t0_id, j);
                    if (is_visited.find(opp_t_id) != is_visited.end() && !is_visited[opp_t_id]) {
                        tet_queue.push(opp_t_id);
                        is_visited[opp_t_id] = true;
                    } else {
                    }
                }
            }
        }
        if (tet_groups.size() < 2)
            continue;

        cout<<"find non-manifold edge "<<e[0]<<" "<<e[1]<<endl;
        cout<<tet_groups.size()<<"/"<<n_t_ids.size()<<endl;
        cout<<e[0]<<" "<<e[1]<<endl;
        for (int t_id: n_t_ids) {
            cout<<t_id<<": ";
            tets[t_id].print();
        }

        //split
        std::vector<int> new_t_ids;
        std::vector<int> old_t_ids;
        int v_id = split(e[0], e[1], old_t_ids, new_t_ids);
        if (v_id < 0)
            continue;
        std::map<int, int> old_t_ids_map;
        for (int i = 0; i < old_t_ids.size(); i++) {
            old_t_ids_map[old_t_ids[i]] = i;
        }

        //duplicate v_id
        for (int i = 0; i < tet_groups.size(); i++) {
            int dup_v_id = v_id;
            if(i > 0) {
                tet_vertices.push_back(tet_vertices[v_id]);
                dup_v_id = tet_vertices.size() - 1;
            }
            tet_vertices[dup_v_id].conn_tets.clear();
            for (int old_t_id: tet_groups[i]) {
                int new_t_id = new_t_ids[old_t_ids_map[old_t_id]];
                int j = tets[old_t_id].find(v_id);
                tets[old_t_id][j] = dup_v_id;
                j = tets[new_t_id].find(v_id);
                tets[new_t_id][j] = dup_v_id;
                tet_vertices[dup_v_id].conn_tets.push_back(old_t_id);
                tet_vertices[dup_v_id].conn_tets.push_back(new_t_id);
            }
        }


        //push new edges into the queue
        old_t_ids.insert(old_t_ids.end(), new_t_ids.begin(), new_t_ids.end());
        std::vector<std::array<int, 2>> new_edges;
        static const std::array<std::array<int, 2>, 6> t_es = {{{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
        for (int t_id: old_t_ids) {
            for (auto &le: t_es) {
                if (tets[t_id][le[0]] < tets[t_id][le[1]])
                    new_edges.push_back({{tets[t_id][le[0]], tets[t_id][le[1]]}});
                else
                    new_edges.push_back({{tets[t_id][le[1]], tets[t_id][le[0]]}});
            }
        }
        vector_unique(new_edges);
        for (auto &ne : new_edges)
            edge_queue.push(ne);
    }
}

void floatTetWild::manifold_vertices(Mesh& mesh){
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    std::vector<std::array<int, 3>> faces;
    for (auto &t: tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[j], t[(j + 1) % 4], t[(j + 2) % 4]}};
            std::sort(f.begin(), f.end());
            faces.push_back(f);
        }
    }
    std::sort(faces.begin(), faces.end());
    if (faces.empty())
        return;

    ///
    std::vector<std::array<int, 3>> b_faces;
    bool is_boundary = true;
    for (int i = 0; i < faces.size() - 1; i++) {
        if (faces[i] == faces[i + 1]) {
            is_boundary = false;
        } else {
            if (is_boundary) {
                b_faces.push_back(faces[i]);
            }
            is_boundary = true;
        }
    }
    if (is_boundary) {
        b_faces.push_back(faces.back());
    }

    ///
    std::vector<int> b_v_ids;
    for (int i = 0; i < b_faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            b_v_ids.push_back(b_faces[i][j]);
        }
    }
    vector_unique(b_v_ids);

    ///
    for (int b_v_id: b_v_ids) {
        std::map<int, bool> is_visited;
        for (int t_id: tet_vertices[b_v_id].conn_tets) {
            if (!tets[t_id].is_removed)
                is_visited[t_id] = false;
        }

        std::vector<std::vector<int>> tet_groups;
        for (int t_id: tet_vertices[b_v_id].conn_tets) {
            if (is_visited.find(t_id) == is_visited.end())
                continue;
            if (is_visited[t_id])
                continue;
            is_visited[t_id] = true;

            tet_groups.emplace_back();
            std::queue<int> tet_queue;
            tet_queue.push(t_id);
            while (!tet_queue.empty()) {
                int t0_id = tet_queue.front();
                tet_queue.pop();
                tet_groups.back().push_back(t0_id);

                int j = tets[t0_id].find(b_v_id);
                for (int k = 0; k < 3; k++) {
                    int opp_t_id = get_opp_t_id(mesh, t0_id, (j + 1 + k) % 4);
                    if (is_visited.find(opp_t_id) != is_visited.end() && !is_visited[opp_t_id]) {
                        tet_queue.push(opp_t_id);
                        is_visited[opp_t_id] = true;
                    }
                }
            }
        }
        //
        if (tet_groups.size() < 2) {
            continue;

        }


        cout << "find non-manifold vertex " << b_v_id << endl;

        for (int i = 1; i < tet_groups.size(); i++) {
            tet_vertices.push_back(tet_vertices[b_v_id]);
            for (int t_id: tet_groups[i]) {
                int j = tets[t_id].find(b_v_id);
                tets[t_id][j] = tet_vertices.size() - 1;
            }
        }
    }
}

void floatTetWild::get_surface(Mesh& mesh, MatrixXd& V, MatrixXi& F) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    std::vector<std::array<int, 5>> faces;
    for (int i=0;i<tets.size();i++) {
        auto &t = tets[i];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4]}};
            std::sort(f.begin(), f.end());
            faces.push_back({{f[0], f[1], f[2], i, j}});
        }
    }
    std::sort(faces.begin(), faces.end(), [](const std::array<int, 5>& a, const std::array<int, 5>& b){
        return std::make_tuple(a[0], a[1], a[2]) < std::make_tuple(b[0], b[1], b[2]);
    });
    if (faces.empty())
        return;
    //
    std::vector<std::array<int, 3>> b_faces;
    bool is_boundary = true;
    for (int i = 0; i < faces.size() - 1; i++) {
        if (std::make_tuple(faces[i][0], faces[i][1], faces[i][2])
            == std::make_tuple(faces[i + 1][0], faces[i + 1][1], faces[i + 1][2])) {
            is_boundary = false;
        } else {
            if (is_boundary) {
                b_faces.push_back({{faces[i][0], faces[i][1], faces[i][2]}});
                bool is_inv = is_inverted(tet_vertices[tets[faces[i][3]][faces[i][4]]],
                                          tet_vertices[faces[i][0]],
                                          tet_vertices[faces[i][1]],
                                          tet_vertices[faces[i][2]]);
                if (!is_inv)
                    std::swap(b_faces.back()[1], b_faces.back()[2]);
            }
            is_boundary = true;
        }
    }
    if (is_boundary) {
        b_faces.push_back({{faces.back()[0], faces.back()[1], faces.back()[2]}});
        bool is_inv = is_inverted(tet_vertices[tets[faces.back()[3]][faces.back()[4]]],
                                  tet_vertices[faces.back()[0]],
                                  tet_vertices[faces.back()[1]],
                                  tet_vertices[faces.back()[2]]);
        if(!is_inv)
            std::swap(b_faces.back()[1], b_faces.back()[2]);
    }
    //
    std::vector<int> b_v_ids;
    for (int i = 0; i < b_faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            b_v_ids.push_back(b_faces[i][j]);
        }
    }
    vector_unique(b_v_ids);

    ///
    V.resize(b_v_ids.size(), 3);
    F.resize(b_faces.size(), 3);
    std::map<int, int> map_v_ids;
    for (int i = 0; i < b_v_ids.size(); i++) {
        map_v_ids[b_v_ids[i]] = i;
        V.row(i) = tet_vertices[b_v_ids[i]].pos;
    }
    for (int i = 0; i < b_faces.size(); i++) {
        F.row(i) << map_v_ids[b_faces[i][0]], map_v_ids[b_faces[i][1]], map_v_ids[b_faces[i][2]];
    }
}

void floatTetWild::manifold_surface(Mesh& mesh, MatrixXd& V, MatrixXi& F) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    for (auto &v: tet_vertices) {
        if (v.is_removed)
            continue;

        for (int i = 0; i < v.conn_tets.size(); i++) {
            if (tets[v.conn_tets[i]].is_removed) {
                v.conn_tets.erase(v.conn_tets.begin() + i);
                i--;
            }
        }
        if (v.conn_tets.empty())
            v.is_removed = true;
    }

    manifold_edges(mesh);
    manifold_vertices(mesh);

    get_surface(mesh, V, F);
    //fix pinched-pie
    std::vector<std::vector<int>> conn_f4v(V.rows());
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++)
            conn_f4v[F(i, j)].push_back(i);
    }
    int V_size = V.rows();
    for (int v_id = 0; v_id < V_size; v_id++) {
        if (conn_f4v[v_id].empty())
            continue;
        //
        std::map<int, bool> is_visited;
        for (int f_id: conn_f4v[v_id])
            is_visited[f_id] = false;
        //
        std::queue<int> f_queue;
        f_queue.push(conn_f4v[v_id][0]);
        is_visited[conn_f4v[v_id][0]] = true;
        std::vector<int> f_group;
        while (!f_queue.empty()) {
            int f_id = f_queue.front();
            f_group.push_back(f_id);
            f_queue.pop();
            //
            for (int j = 0; j < 3; j++) {
                if (F(f_id, j) == v_id)
                    continue;
                std::vector<int> tmp;
                set_intersection(conn_f4v[F(f_id, (j + 1) % 3)], conn_f4v[F(f_id, (j + 2) % 3)], tmp);
                if (tmp.size() != 2)
                    continue;
                int n_f_id = tmp[0] == f_id ? tmp[1] : tmp[0];
                if (is_visited.find(n_f_id) == is_visited.end() || is_visited[n_f_id])
                    continue;
                is_visited[n_f_id] = true;
                f_queue.push(n_f_id);
            }
        }
        if (f_group.size() == conn_f4v[v_id].size())
            continue;
        //
        cout << "HHHHHHHHHHHH" << endl;
        V.conservativeResize(V.rows() + 1, V.cols());
        V.row(V.rows() - 1) = V.row(v_id);
        for (int f_id:f_group) {
            for (int j = 0; j < 3; j++) {
                if (F(f_id, j) == v_id) {
                    F(f_id, j) = V.rows() - 1;
                    break;
                }
            }
        }
        conn_f4v.push_back(f_group);
    }

}
