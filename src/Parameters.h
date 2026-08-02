// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Types.hpp>


#include <floattetwild/Logger.hpp>

namespace floatTetWild {
    class Parameters {
    public:
        std::string log_path = "";
        std::string input_path = "";
        std::string output_path = "";
        std::string tag_path = "";
        std::string postfix = "";


        bool not_sort_input = false;
        bool is_quiet = false;
        int log_level = 2;  // Info

        bool smooth_open_boundary = false;
        bool manifold_surface = false;
        bool disable_filtering = false;
        bool use_floodfill = false;
        bool use_input_for_wn = false;
        bool coarsen = false;

        // Fractions of the bounding box diagonal, except ideal_edge_length_abs, which is an
        // absolute length and wins over ideal_edge_length_rel when it is set.
        Scalar box_scale = 1 / 15.0;
        Scalar eps_rel = 1e-3;
        Scalar ideal_edge_length_rel = 1 / 20.0;
        Scalar min_edge_len_rel = -1;
        Scalar ideal_edge_length_abs = 0.0;

        int max_its = 80;
        Scalar stop_energy = 10;

        int stage = 2;

        unsigned int num_threads = std::numeric_limits<unsigned int>::max();

        int stop_p = -1;

        Vector3 bbox_min;
        Vector3 bbox_max;
        Scalar bbox_diag_length;
        Scalar ideal_edge_length;
        Scalar eps_input;
        Scalar eps;
        Scalar eps_delta;
        Scalar eps_2;
        Scalar dd;

        Scalar split_threshold;
        Scalar collapse_threshold;
        Scalar split_threshold_2;
        Scalar collapse_threshold_2;

        Scalar eps_coplanar;
        Scalar eps_2_coplanar;
        Scalar eps_simplification;
        Scalar eps_2_simplification;
        Scalar dd_simplification;

        void init(Scalar bbox_diag_l) {
            if (stage > 5)
                stage = 5;

            bbox_diag_length = bbox_diag_l;

            if (ideal_edge_length_abs > 0.0) {
                ideal_edge_length = ideal_edge_length_abs;
                ideal_edge_length_rel = ideal_edge_length / bbox_diag_length;
            }
            else {
                ideal_edge_length = bbox_diag_length * ideal_edge_length_rel;
            }

            eps_input = bbox_diag_length * eps_rel;
            dd = eps_input;
            dd /= 1.5;

            // The sampled envelope check is conservative by the sampling error bound, the
            // circumradius dd/sqrt(3) of the triangle the samples sit on, so that budget comes
            // out of the epsilon the user asked for.
            double eps_usable = eps_input - dd / std::sqrt(3);
            eps_delta = eps_usable * 0.1;
            eps = eps_usable - eps_delta * (stage - 1);

            eps_2 = eps * eps;

            eps_coplanar = eps * 0.2;
            if (eps_coplanar > bbox_diag_length * 1e-6)
                eps_coplanar = bbox_diag_length * 1e-6;
            eps_2_coplanar = eps_coplanar * eps_coplanar;

            eps_simplification = eps * 0.8;
            eps_2_simplification = eps_simplification * eps_simplification;
            dd_simplification = dd / eps * eps_simplification;

            if (min_edge_len_rel < 0)
                min_edge_len_rel = eps_rel;

            split_threshold = ideal_edge_length * (4 / 3.0);
            collapse_threshold = ideal_edge_length * (4 / 5.0);
            split_threshold_2 = split_threshold * split_threshold;
            collapse_threshold_2 = collapse_threshold * collapse_threshold;

            logger().debug("bbox_diag_length = {}", bbox_diag_length);
            logger().debug("ideal_edge_length = {}", ideal_edge_length);

            logger().debug("stage = {}", stage);
            logger().debug("eps_input = {}", eps_input);
            logger().debug("eps = {}", eps);
            logger().debug("eps_simplification = {}", eps_simplification);
            logger().debug("eps_coplanar = {}", eps_coplanar);

            logger().debug("dd = {}", dd);
            logger().debug("dd_simplification = {}", dd_simplification);
        }
    };
}  // namespace floatTetWild
