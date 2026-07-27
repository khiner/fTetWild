// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_PARALLEL_FOR_H
#define FLOATTETWILD_PARALLEL_FOR_H

#include <floattetwild/default_num_threads.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <thread>
#include <vector>

namespace floatTetWild
{
  // Functional implementation of an open-mp style, parallel for loop with
  // accumulation.
  //
  // Inputs:
  //   loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
  //   prep_func  function handle taking n >= number of threads as only argument
  //   func  function handle taking iteration index i and thread id t as only
  //     arguments to compute inner block of for loop I.e.
  //     for(int i ...){ func(i,t); }
  //   accum_func  function handle taking thread index as only argument, to be
  //     called after all calls of func, e.g., for serial accumulation across
  //     all n (potential) threads, see n in description of prep_func.
  //   min_parallel  min size of loop_size such that parallel (non-serial)
  //     thread pooling should be attempted {0}
  // Returns true iff thread pool was invoked
  template<
    typename Index,
    typename PreFunctionType,
    typename FunctionType,
    typename AccumFunctionType>
  inline bool parallel_for(
    const Index loop_size,
    const PreFunctionType & prep_func,
    const FunctionType & func,
    const AccumFunctionType & accum_func,
    const size_t min_parallel=0)
  {
    assert(loop_size>=0);
    if(loop_size==0) return false;
    // Estimate number of threads in the pool
    // http://ideone.com/Z7zldb
#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
    const size_t nthreads = 1;
#else
    const size_t nthreads = floatTetWild::default_num_threads();
#endif
    if(loop_size<min_parallel || nthreads<=1)
    {
      // serial
      prep_func(1);
      for(Index i = 0;i<loop_size;i++) func(i,0);
      accum_func(0);
      return false;
    }else
    {
      // Size of a slice for the range functions
      Index slice =
        std::max(
          (Index)std::round((loop_size+1)/static_cast<double>(nthreads)),(Index)1);

      // [Helper] Inner loop
      const auto & range = [&func](const Index k1, const Index k2, const size_t t)
      {
        for(Index k = k1; k < k2; k++) func(k,t);
      };
      prep_func(nthreads);
      // Create pool and launch jobs
      std::vector<std::thread> pool;
      pool.reserve(nthreads);
      // Inner range extents
      Index i1 = 0;
      Index i2 = std::min(0 + slice, loop_size);
      {
        size_t t = 0;
        for (; t+1 < nthreads && i1 < loop_size; ++t)
        {
          pool.emplace_back(range, i1, i2, t);
          i1 = i2;
          i2 = std::min(i2 + slice, loop_size);
        }
        if (i1 < loop_size)
        {
          pool.emplace_back(range, i1, loop_size, t);
        }
      }
      // Wait for jobs to finish
      for (std::thread &t : pool) if (t.joinable()) t.join();
      // Accumulate across threads
      for(size_t t = 0;t<nthreads;t++)
      {
        accum_func(t);
      }
      return true;
    }
  }

  // \overload
  //
  // Inputs:
  //   func  function handle taking iteration index as only argument to compute
  //     inner block of for loop I.e. for(int i ...){ func(i); }
  template<typename Index, typename FunctionType >
  inline bool parallel_for(
    const Index loop_size,
    const FunctionType & func,
    const size_t min_parallel=0)
  {
    using namespace std;
    // no op preparation/accumulation
    const auto & no_op = [](const size_t /*n/t*/){};
    // two-parameter wrapper ignoring thread id
    const auto & wrapper = [&func](Index i,size_t /*t*/){ func(i); };
    return parallel_for(loop_size,no_op,wrapper,no_op,min_parallel);
  }
}

#endif
