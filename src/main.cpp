// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <CLI/CLI.hpp>

#include <thread>

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetwild.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>

#include <floattetwild/Logger.hpp>

#include <floattetwild/Timer.h>
#include <floattetwild/writeOBJ.h>
#include <floattetwild/ParallelFor.hpp>


#include <floattetwild/geo_mesh.h>

using namespace floatTetWild;



int main(int argc, char** argv)
{
    bool skip_simplify = false;
    bool nobinary      = false;
    bool nocolor       = false;
    bool export_raw    = false;

    Mesh        mesh;
    Parameters& params = mesh.params;

    CLI::App command_line {"float-tetwild"};
    command_line
      .add_option("-i,--input",
                  params.input_path,
                  "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")
      ->check(CLI::ExistingFile);
    command_line.add_option("-o,--output",
                            params.output_path,
                            "Output tetmesh OUTPUT in .msh format. (string, optional, default: "
                            "input_file+postfix+'.msh')");

    command_line.add_option("--tag", params.tag_path, "Tag input faces for Boolean operation.")
      ->check(CLI::ExistingFile);
    int         boolean_op = -1;
    std::string csg_file   = "";
    command_line.add_option(
      "--op", boolean_op, "Boolean operation: 0: union, 1: intersection, 2: difference.");

    auto absolute_op = command_line.add_option(
      "-a,--la",
      params.ideal_edge_length_abs,
      "Ideal edge length not scaled by diag_of_bbox. (double, optional)");
    auto relative_op = command_line.add_option(
      "-l,--lr",
      params.ideal_edge_length_rel,
      "ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)");
    relative_op->excludes(absolute_op);

    command_line.add_option("-e,--epsr",
                            params.eps_rel,
                            "epsilon = diag_of_bbox * EPS. (double, optional, default: 1e-3)");

    command_line.add_option("--max-its", params.max_its, "(for debugging usage only)");
    command_line.add_option(
      "--stop-energy", params.stop_energy, "Stop optimization when max energy is lower than this.");
    command_line.add_option("--stage", params.stage, "(for debugging usage only)");
    command_line.add_option("--stop-p", params.stop_p, "(for debugging usage only)");

    command_line.add_option("--postfix", params.postfix, "(for debugging usage only)");
    command_line.add_option("--log", params.log_path, "Log info to given file.");
    command_line.add_option("--level", params.log_level, "Log level (0 = most verbose, 6 = off).");

    command_line.add_flag("-q,--is-quiet", params.is_quiet, "Mute console output. (optional)");
    command_line.add_flag("--skip-simplify", skip_simplify, "skip preprocessing.");
    command_line.add_flag("--no-binary", nobinary, "export meshes as ascii");
    command_line.add_flag("--no-color", nocolor, "don't export color");
    command_line.add_flag("--not-sort-input", params.not_sort_input, "(for debugging usage only)");

    command_line.add_flag(
      "--smooth-open-boundary", params.smooth_open_boundary, "Smooth the open boundary.");
    command_line.add_flag("--export-raw", export_raw, "Export raw output.");
    command_line.add_flag(
      "--manifold-surface", params.manifold_surface, "Force the output to be manifold.");
    command_line.add_flag("--coarsen", params.coarsen, "Coarsen the output as much as possible.");
    command_line.add_option("--csg", csg_file, "File containing a csg tree.")
      ->check(CLI::ExistingFile);

    command_line.add_flag(
      "--use-old-energy", floatTetWild::use_old_energy, "(for debugging usage only)");  // tmp

    command_line.add_flag(
      "--disable-filtering", params.disable_filtering, "Disable filtering out outside elements.");
    command_line.add_flag(
      "--use-floodfill", params.use_floodfill, "Use flood-fill to extract interior volume.");
    command_line.add_flag(
      "--use-input-for-wn", params.use_input_for_wn, "Use input surface for winding number.");

    unsigned int max_threads = std::numeric_limits<unsigned int>::max();
    command_line.add_option("--max-threads", max_threads, "Maximum number of threads used");

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return command_line.exit(e);
    }

    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads              = std::min(max_threads, num_threads);
    params.num_threads       = num_threads;
    floatTetWild::set_num_threads(num_threads);

    Logger::init(!params.is_quiet, params.log_path);
    params.log_level = std::max(0, std::min(6, params.log_level));
    logger().set_level(params.log_level);
    logger().info("threads {}", num_threads);

    if (params.output_path.empty())
        params.output_path = params.input_path;
    if (params.log_path.empty())
        params.log_path = params.output_path;

    // Everything but the main mesh is named after the output path plus the postfix.
    const std::string out_prefix = params.output_path + "_" + params.postfix;

    // An output path that already names a mesh file is taken as is.
    const auto ends_with = [&](const std::string& suffix) {
        return params.output_path.size() > suffix.size() &&
               params.output_path.compare(params.output_path.size() - suffix.size(),
                                          suffix.size(),
                                          suffix) == 0;
    };
    const std::string output_mesh_name =
      ends_with("msh") || ends_with("mesh") ? params.output_path : out_prefix + ".msh";

    /// set input tage
    std::vector<Vector3>  input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int>      input_tags;

    if (!params.tag_path.empty()) {
        input_tags.reserve(input_faces.size());
        std::string   line;
        std::ifstream fin(params.tag_path);
        if (fin.is_open()) {
            while (getline(fin, line)) {
                input_tags.push_back(std::stoi(line));
            }
            fin.close();
        }
    }

    /// set envelope
    Timer               timer;
    geo::Mesh                sf_mesh;
    CSGTree                  tree_with_ids;
    std::vector<std::string>                 meshes;
    std::vector<std::vector<Vector3>>        csg_Vs;
    std::vector<std::vector<Vector3i>>       csg_Fs;
    if (!csg_file.empty()) {
        CSGTree       csg_tree;
        std::ifstream file(csg_file);

        if (!file.is_open()) {
            logger().error("unable to open {} file", csg_file);
            return EXIT_FAILURE;
        }
        std::string parse_error;
        if (!CSGTreeParser::parse(file, csg_tree, parse_error)) {
            logger().error("{}: {}", csg_file, parse_error);
            return EXIT_FAILURE;
        }
        file.close();

        CSGTreeParser::get_meshes(csg_tree, meshes, tree_with_ids);

        csg_Vs.resize(meshes.size());
        csg_Fs.resize(meshes.size());
        {
            geo::Mesh        tmp_mesh;
            std::vector<int> tmp_tags;
            for (int i = 0; i < meshes.size(); ++i) {
                if (!MeshIO::load_mesh(meshes[i], csg_Vs[i], csg_Fs[i], tmp_mesh, tmp_tags)) {
                    logger().error("unable to open {} file", meshes[i]);
                    return EXIT_FAILURE;
                }
            }
        }
        CSGTreeParser::merge(csg_Vs, csg_Fs, input_vertices, input_faces, sf_mesh, input_tags);

    }
    else {
        if (!MeshIO::load_mesh(
              params.input_path, input_vertices, input_faces, sf_mesh, input_tags)) {
            logger().error("Unable to load mesh at {}", params.input_path);
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return EXIT_FAILURE;
        }
        else if (input_vertices.empty() || input_faces.empty()) {
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return EXIT_FAILURE;
        }

        if (input_tags.size() != input_faces.size()) {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
    }
    AABBWrapper tree(sf_mesh);

    if (tetrahedralization(tree, input_vertices, input_faces, input_tags, mesh, skip_simplify) !=
        0) {
        MeshIO::write_mesh(output_mesh_name, mesh, false);
        return EXIT_FAILURE;
    }

    // Times the interior extraction below, recorded as wn_id. correct_tracked_surface_orientation
    // used to fall inside this window and is now inside tetrahedralization.
    timer.start();

    if (export_raw) {
        MatrixXs Vt;
        MatrixXi Ft;

        if (!csg_file.empty()) {
            int max_id = CSGTreeParser::get_max_id(tree_with_ids);

            for (int i = 0; i <= max_id; ++i) {
                get_tracked_surface(mesh, Vt, Ft, i);
                writeOBJ(out_prefix + "_" + std::to_string(i) + "_all.obj", Vt, Ft);
            }
        }
        else {
            get_tracked_surface(mesh, Vt, Ft);
            writeOBJ(out_prefix + "_all.obj", Vt, Ft);
        }
        MeshIO::write_mesh(out_prefix + "_all.msh", mesh, false);
    }

    MatrixXs Vt;
    MatrixXi Ft;
    get_tracked_surface(mesh, Vt, Ft);
    writeOBJ(out_prefix + "_tracked_surface.obj", Vt, Ft);

    if (!csg_file.empty())
        boolean_operation(mesh, tree_with_ids, csg_Vs, csg_Fs);
    else if (boolean_op >= 0)
        boolean_operation(mesh, boolean_op);
    else {
        if (params.smooth_open_boundary) {
            smooth_open_boundary(mesh, tree);
            for (auto& t : mesh.tets) {
                if (t.is_outside)
                    t.is_removed = true;
            }
        }
        else {
            if (!params.disable_filtering) {
                if (params.use_floodfill) {
                    filter_outside_floodfill(mesh);
                }
                else if (params.use_input_for_wn) {
                    filter_outside(mesh, input_vertices, input_faces);
                }
                else
                    filter_outside(mesh, Vt, Ft);
            }
        }
    }
    MatrixXd V_sf;
    MatrixXi F_sf;
    if (params.manifold_surface) {
        manifold_surface(mesh, V_sf, F_sf);
    }
    else {
        get_surface(mesh, V_sf, F_sf);
    }
    stats().record(StateInfo::wn_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy());
    logger().info("after winding number");
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    logger().info("winding number {}s", timer.getElapsedTimeInSec());
    logger().info("");

    std::vector<Scalar> colors;
    if (!nocolor) {
        colors.resize(mesh.tets.size(), -1);
        for (int i = 0; i < mesh.tets.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            colors[i] = mesh.tets[i].quality;
        }
    }
    MeshIO::write_mesh(output_mesh_name, mesh, false, colors, !nobinary, !csg_file.empty());
    writeOBJ(out_prefix + "_sf.obj", V_sf, F_sf);

    std::ofstream fout(params.log_path + "_" + params.postfix + ".csv");
    if (fout.good())
        fout << stats();
    fout.close();

    return EXIT_SUCCESS;
}
