// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <CLI/CLI.hpp>

#include <thread>

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

#include "AABBWrapper.h"
#include "FloatTetwild.h"
#include "LocalOperations.h"
#include "MeshImprovement.h"
#include "Statistics.h"
#include "CSGTreeParser.hpp"
#include "Mesh.hpp"
#include "MeshIO.hpp"

#include "Logger.hpp"

#include "Timer.h"
#include "ParallelFor.hpp"

#include "geo/geo_mesh.h"

using namespace floatTetWild;

namespace {

// Peak resident set size (physical memory use) in bytes, or zero if the value cannot be
// determined on this OS. Only the stats csv summary line below asks.
//
// From getCurrentAndPeakRSS by David Robert Nadeau, http://NadeauSoftware.com/, Creative Commons
// Attribution 3.0 Unported License, http://creativecommons.org/licenses/by/3.0/deed.en_US. Only
// getPeakRSS() is kept, and only for the platforms this project builds on. The original also had
// getCurrentRSS(), and branches for AIX and Solaris reading /proc/self/psinfo.
size_t getPeakRSS()
{
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;
#else
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif
#endif
}

// Write a mesh in an ascii obj file: V is #V by 3 vertex positions, F is #F by 3 indices into V.
// From libigl's igl/writeOBJ.h, Copyright (C) 2013 Alec Jacobson, MPL 2.0.
void writeOBJ(const std::string& path, const MatrixXd& V, const MatrixXi& F)
{
    std::ofstream s(path);
    if (!s.is_open()) {
        fprintf(stderr, "IOError: writeOBJ() could not open %s\n", path.c_str());
        return;
    }
    // Eigen wrote each block through IOFormat(FullPrecision, DontAlignCols, " ", "\n", prefix,
    // "", "", "\n"): a prefixed, space-separated row per line, and for an empty block nothing at
    // all but the trailing newline. FullPrecision resolved to digits10 significant digits, which
    // is what the vertices go out at; the indices are integers and ignore it.
    s.precision(std::numeric_limits<double>::digits10);
    if (V.rows() == 0) s << "\n";
    for (int i = 0; i < V.rows(); i++)
        s << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    if (F.rows() == 0) s << "\n";
    for (int i = 0; i < F.rows(); i++)
        s << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1 << "\n";
}

// Every recorded state as a row, then a summary row with id -1: the time the stages up to the
// winding number took together, the last state's counts and energies, and the peak memory.
void write_stats_csv(std::ostream& stream, const std::vector<StateInfo>& states)
{
    double time           = 0;
    int    cnt_uninserted = 0;
    for (const StateInfo& s : states) {
        stream << s.id << ", " << s.time << ", " << s.v_num << ", " << s.t_num << ", "
               << s.max_energy << ", " << s.avg_energy << ", " << s.cnt_fail_inserted_face
               << ", -1" << std::endl;
        if (s.cnt_fail_inserted_face >= 0)
            cnt_uninserted = s.cnt_fail_inserted_face;
        if (s.id < 6)
            time += s.time;
    }
    stream << -1 << ", " << time << ", " << states.back().v_num << ", " << states.back().t_num
           << ", " << states.back().max_energy << ", " << states.back().avg_energy << ", "
           << cnt_uninserted << ", " << getPeakRSS() / (1024 * 1024) << std::endl;
}

}  // namespace

int main(int argc, char** argv)
{
    bool skip_simplify = false;
    bool nobinary      = false;
    bool nocolor       = false;
    bool export_raw    = false;

    // Command-line state the mesher itself never reads. What the library consumes stays in
    // mesh.params.
    std::string input_path, output_path, tag_path, log_path, postfix;
    bool        is_quiet  = false;
    int         log_level = 2;  // Info
    bool do_smooth_open_boundary = false, do_manifold_surface = false, disable_filtering = false,
         use_floodfill = false, use_input_for_wn = false;

    Mesh        mesh;
    Parameters& params = mesh.params;

    CLI::App command_line {"float-tetwild"};
    command_line
      .add_option("-i,--input",
                  input_path,
                  "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")
      ->check(CLI::ExistingFile);
    command_line.add_option("-o,--output",
                            output_path,
                            "Output tetmesh OUTPUT in .msh format. (string, optional, default: "
                            "input_file+postfix+'.msh')");

    command_line.add_option("--tag", tag_path, "Tag input faces for Boolean operation.")
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

    command_line.add_option("--postfix", postfix, "(for debugging usage only)");
    command_line.add_option("--log", log_path, "Log info to given file.");
    command_line.add_option("--level", log_level, "Log level (0 = most verbose, 6 = off).");

    command_line.add_flag("-q,--is-quiet", is_quiet, "Mute console output. (optional)");
    command_line.add_flag("--skip-simplify", skip_simplify, "skip preprocessing.");
    command_line.add_flag("--no-binary", nobinary, "export meshes as ascii");
    command_line.add_flag("--no-color", nocolor, "don't export color");
    command_line.add_flag("--not-sort-input", params.not_sort_input, "(for debugging usage only)");

    command_line.add_flag(
      "--smooth-open-boundary", do_smooth_open_boundary, "Smooth the open boundary.");
    command_line.add_flag("--export-raw", export_raw, "Export raw output.");
    command_line.add_flag(
      "--manifold-surface", do_manifold_surface, "Force the output to be manifold.");
    command_line.add_flag("--coarsen", params.coarsen, "Coarsen the output as much as possible.");
    command_line.add_option("--csg", csg_file, "File containing a csg tree.")
      ->check(CLI::ExistingFile);

    command_line.add_flag(
      "--use-old-energy", floatTetWild::use_old_energy, "(for debugging usage only)");

    command_line.add_flag(
      "--disable-filtering", disable_filtering, "Disable filtering out outside elements.");
    command_line.add_flag(
      "--use-floodfill", use_floodfill, "Use flood-fill to extract interior volume.");
    command_line.add_flag(
      "--use-input-for-wn", use_input_for_wn, "Use input surface for winding number.");

    unsigned int max_threads = std::numeric_limits<unsigned int>::max();
    command_line.add_option("--max-threads", max_threads, "Maximum number of threads used");

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return command_line.exit(e);
    }

    const unsigned int num_threads =
      std::min(max_threads, std::max(1u, std::thread::hardware_concurrency()));
    floatTetWild::set_num_threads(num_threads);

    // The log file is opt-in: without --log, log_path is empty and logging goes to the console
    // only. It cannot default to output_path, which init would truncate — that path is the
    // output mesh when -o names one, and the *input* mesh when -o is omitted.
    Logger::init(!is_quiet, log_path);
    log_level = std::max(0, std::min(6, log_level));
    logger().set_level(log_level);
    logger().info("threads {}", num_threads);

    if (output_path.empty())
        output_path = input_path;

    // Everything but the main mesh is named after the output path plus the postfix.
    const std::string out_prefix = output_path + "_" + postfix;

    // An output path that already names a mesh file is taken as is.
    const auto ends_with = [&](const std::string& suffix) {
        return output_path.size() > suffix.size() &&
               output_path.compare(output_path.size() - suffix.size(),
                                          suffix.size(),
                                          suffix) == 0;
    };
    const std::string output_mesh_name =
      ends_with("msh") || ends_with("mesh") ? output_path : out_prefix + ".msh";

    std::vector<Vector3>  input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int>      input_tags;

    if (!tag_path.empty()) {
        std::string   line;
        std::ifstream fin(tag_path);
        if (fin.is_open()) {
            while (getline(fin, line)) {
                input_tags.push_back(std::stoi(line));
            }
            fin.close();
        }
    }

    Timer                              timer;
    geo::Mesh                          sf_mesh;
    CSGTree                            csg_tree;
    std::vector<std::vector<Vector3>>  csg_Vs;
    std::vector<std::vector<Vector3i>> csg_Fs;
    if (!csg_file.empty()) {
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

        const std::vector<std::string> meshes = CSGTreeParser::assign_mesh_ids(csg_tree);

        csg_Vs.resize(meshes.size());
        csg_Fs.resize(meshes.size());
        {
            geo::Mesh        tmp_mesh;
            std::vector<int> tmp_tags;
            for (int i = 0; i < meshes.size(); ++i) {
                if (!load_mesh(meshes[i], csg_Vs[i], csg_Fs[i], tmp_mesh, tmp_tags)) {
                    logger().error("unable to open {} file", meshes[i]);
                    return EXIT_FAILURE;
                }
            }
        }
        CSGTreeParser::merge(csg_Vs, csg_Fs, input_vertices, input_faces, sf_mesh, input_tags);
    }
    else {
        if (!load_mesh(input_path, input_vertices, input_faces, sf_mesh, input_tags))
            logger().error("Unable to load mesh at {}", input_path);
        if (input_vertices.empty() || input_faces.empty()) {
            write_mesh(output_mesh_name, mesh);
            return EXIT_FAILURE;
        }

        if (input_tags.size() != input_faces.size())
            input_tags.assign(input_faces.size(), 0);
    }
    AABBWrapper tree(sf_mesh);

    if (tetrahedralization(tree, input_vertices, input_faces, input_tags, mesh, skip_simplify) !=
        0) {
        write_mesh(output_mesh_name, mesh);
        return EXIT_FAILURE;
    }

    // Times the interior extraction below, recorded as wn_id. correct_tracked_surface_orientation
    // used to fall inside this window and is now inside tetrahedralization.
    timer.start();

    MatrixXd Vt;
    MatrixXi Ft;
    // Leaves the surface it wrote in Vt and Ft, which the winding number below reuses.
    const auto write_tracked_surface = [&](const std::string& name, int c_id = 0) {
        get_tracked_surface(mesh, Vt, Ft, c_id);
        writeOBJ(out_prefix + name, Vt, Ft);
    };

    if (export_raw) {
        if (!csg_file.empty()) {
            for (int i = 0; i < int(csg_Vs.size()); ++i)
                write_tracked_surface("_" + std::to_string(i) + "_all.obj", i);
        }
        else
            write_tracked_surface("_all.obj");
        write_mesh(out_prefix + "_all.msh", mesh);
    }

    write_tracked_surface("_tracked_surface.obj");

    if (!csg_file.empty())
        boolean_operation(mesh, csg_tree, csg_Vs, csg_Fs);
    else if (boolean_op >= 0)
        boolean_operation(mesh, boolean_op);
    else if (do_smooth_open_boundary) {
        smooth_open_boundary(mesh, tree);
        for (auto& t : mesh.tets) {
            if (t.is_outside)
                t.is_removed = true;
        }
    }
    else if (!disable_filtering) {
        if (use_floodfill)
            filter_outside_floodfill(mesh);
        else if (use_input_for_wn)
            filter_outside(mesh, input_vertices, input_faces);
        else
            filter_outside(mesh, Vt, Ft);
    }
    MatrixXd V_sf;
    MatrixXi F_sf;
    if (do_manifold_surface) {
        manifold_surface(mesh, V_sf, F_sf);
    }
    else {
        get_surface(mesh, V_sf, F_sf);
    }
    Scalar max_energy, avg_energy;
    mesh.get_max_avg_energy(max_energy, avg_energy);
    stats().push_back({StateInfo::wn_id,
                       timer.getElapsedTimeInSec(),
                       mesh.get_v_num(),
                       mesh.get_t_num(),
                       max_energy,
                       avg_energy});
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
    write_mesh(output_mesh_name, mesh, colors, !nobinary, !csg_file.empty());
    writeOBJ(out_prefix + "_sf.obj", V_sf, F_sf);

    // The stats csv sits next to the log file when --log was given, else next to the output.
    std::ofstream fout((log_path.empty() ? output_path : log_path) + "_" + postfix + ".csv");
    if (fout.good())
        write_stats_csv(fout, stats());
    fout.close();

    return EXIT_SUCCESS;
}
