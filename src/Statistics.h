// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Types.hpp>

#include <cstddef>
#include <ostream>
#include <vector>

// Peak resident set size in bytes, from getRSS.c. Only the summary line below asks.
extern "C" size_t getPeakRSS();

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

    // Every state as a row, then a summary row with id -1: the time the stages up to the winding
    // number took together, the last state's counts and energies, and the peak memory.
    inline void write_stats_csv(std::ostream &stream, const std::vector<StateInfo> &states) {
        double time = 0;
        int cnt_uninserted = 0;
        for (const StateInfo &s : states) {
            stream << s.id << ", " << s.time << ", " << s.v_num << ", " << s.t_num << ", "
                   << s.max_energy << ", " << s.avg_energy << ", " << s.cnt_fail_inserted_face
                   << ", -1" << std::endl;
            if (s.cnt_fail_inserted_face >= 0)
                cnt_uninserted = s.cnt_fail_inserted_face;
            if (s.id < 6)
                time += s.time;
        }
        stream << -1 << ", " << time << ", " << states.back().v_num << ", "
               << states.back().t_num << ", " << states.back().max_energy << ", "
               << states.back().avg_energy << ", " << cnt_uninserted << ", "
               << getPeakRSS() / (1024 * 1024) << std::endl;
    }

} // namespace floatTetWild
