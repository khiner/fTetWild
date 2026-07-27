// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2021 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_DEFAULT_NUM_THREADS_H
#define FLOATTETWILD_DEFAULT_NUM_THREADS_H

#include <cstdlib>
#include <thread>

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
}

#endif
