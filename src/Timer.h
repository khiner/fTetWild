// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Replaces libigl's igl::Timer, which picked between mach_absolute_time, QueryPerformanceCounter and
// gettimeofday behind ifdefs. steady_clock is all three. The times reach the log and the .csv beside
// the output, neither of which is hashed.

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
