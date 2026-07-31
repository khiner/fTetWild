// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// This replaces libigl's igl::Timer, which reached for mach_absolute_time, QueryPerformanceCounter
// or gettimeofday behind three-way ifdefs. steady_clock is all three. The method names are the ones
// the call sites already use. Nothing here feeds the mesh: the times reach the log and the .csv
// beside the output, and neither is hashed.

#pragma once

#include <chrono>

namespace floatTetWild {

class Timer
{
  public:
    void start()
    {
        stopped_ = false;
        start_   = Clock::now();
    }

    void stop()
    {
        end_     = Clock::now();
        stopped_ = true;
    }

    // Seconds since start(), or between start() and stop() once stopped.
    double getElapsedTimeInSec() const
    {
        const Clock::time_point end = stopped_ ? end_ : Clock::now();
        return std::chrono::duration<double>(end - start_).count();
    }

    double getElapsedTime() const { return getElapsedTimeInSec(); }

  private:
    using Clock = std::chrono::steady_clock;

    bool              stopped_ = false;
    Clock::time_point start_{};
    Clock::time_point end_{};
};

}  // namespace floatTetWild
