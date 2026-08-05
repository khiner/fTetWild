// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Types.hpp>

#include <vector>

namespace floatTetWild {

    // What a stage left behind when it finished: one row of the .csv written beside the output.
    // The fields are in the order they are written, which is also the order they are recorded in.
    struct StateInfo {
        enum OpId {
            init_id = 0,
            preprocessing_id = 1,
            tetrahedralization_id = 2,
            cutting_id = 3,
            optimization_id = 4,
            wn_id = 5,

            splitting_id = 6,
            collapsing_id = 7,
            swapping_id = 8,
            smoothing_id = 9
        };

        int id = 0;
        double time = 0;
        int v_num = 0;
        int t_num = 0;
        Scalar max_energy = 0;
        Scalar avg_energy = 0;
        int cnt_fail_inserted_face = -1;
    };

    // The one list every stage appends to, in the order the stages ran. Only the thread driving
    // the pipeline records, so there is nothing to lock.
    inline std::vector<StateInfo> &stats() {
        static std::vector<StateInfo> inst;
        return inst;
    }

} // namespace floatTetWild
