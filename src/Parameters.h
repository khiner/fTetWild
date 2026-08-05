// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include "Types.hpp"

namespace floatTetWild {
    struct Parameters {
        bool not_sort_input = false;
        bool coarsen = false;

        // Fractions of the bounding box diagonal, except ideal_edge_length_abs, which is an
        // absolute length and wins over ideal_edge_length_rel when it is set.
        Scalar eps_rel = 1e-3;
        Scalar ideal_edge_length_rel = 1 / 20.0;
        Scalar min_edge_len_rel = -1;
        Scalar ideal_edge_length_abs = 0.0;

        int max_its = 80;
        Scalar stop_energy = 10;

        int stage = 2;

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

        Scalar split_threshold_2;
        Scalar collapse_threshold_2;

        Scalar eps_coplanar;
        Scalar eps_2_coplanar;
        Scalar eps_simplification;
        Scalar eps_2_simplification;
        Scalar dd_simplification;
    };
}  // namespace floatTetWild
