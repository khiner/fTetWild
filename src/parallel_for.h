// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>          (parallel_for)
// Copyright (C) 2021 Jérémie Dumas <jeremie.dumas@ens-lyon.org>      (default_num_threads)
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// This is libigl's thread pool, not the mesher's. It stays because it is what
// FastWindingNumberForSoups.h was adapted onto, and swapping that vendored code over to
// ParallelFor.hpp would change how its BVH build is scheduled. Nothing else uses it.
#ifndef FLOATTETWILD_PARALLEL_FOR_H
#define FLOATTETWILD_PARALLEL_FOR_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <thread>
#include <vector>

namespace floatTetWild
{
  // Returns the default number of threads used by parallel_for. The value returned by the
  // first call to this function is cached. The following strategy is used to determine the
  // default number of threads:
  // 1. User-provided argument force_num_threads if != 0.
  // 2. Environment variable IGL_NUM_THREADS if > 0.
  // 3. Hardware concurrency if != 0.
  // 4. A fallback value of 8 is used otherwise.
  //
  // The environment variable keeps libigl's name so an existing IGL_NUM_THREADS in a user's
  // environment keeps working across this port.
  //
  // It is safe to call this method from multiple threads.
  inline unsigned int default_num_threads(unsigned int user_num_threads = 0)
  {
    // Thread-safe initialization using Meyers' singleton
    class MySingleton {
    public:
      static MySingleton &instance(unsigned int force_num_threads) {
        static MySingleton instance(force_num_threads);
        return instance;
      }

      unsigned int get_num_threads() const { return m_num_threads; }

    private:
      static const char* getenv_nowarning(const char* env_var)
      {
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
        return std::getenv(env_var);
#ifdef _MSC_VER
#pragma warning(pop)
#endif
      }

      MySingleton(unsigned int force_num_threads) {
        // User-defined default
        if (force_num_threads) {
          m_num_threads = force_num_threads;
          return;
        }
        // Set from env var
        if (const char *env_str = getenv_nowarning("IGL_NUM_THREADS")) {
          const int env_num_thread = atoi(env_str);
          if (env_num_thread > 0) {
            m_num_threads = static_cast<unsigned int>(env_num_thread);
            return;
          }
        }
        // Guess from hardware
        const unsigned int hw_num_threads = std::thread::hardware_concurrency();
        if (hw_num_threads) {
          m_num_threads = hw_num_threads;
          return;
        }
        // Fallback when std::thread::hardware_concurrency doesn't work
        m_num_threads = 8u;
      }

      unsigned int m_num_threads = 0;
    };

    return MySingleton::instance(user_num_threads).get_num_threads();
  }

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
