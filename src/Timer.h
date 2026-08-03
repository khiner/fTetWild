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
    void start() { start_ = Clock::now(); }

    // Seconds since start().
    double getElapsedTimeInSec() const
    {
        return std::chrono::duration<double>(Clock::now() - start_).count();
    }

  private:
    using Clock = std::chrono::steady_clock;

    Clock::time_point start_{};
};

}  // namespace floatTetWild
