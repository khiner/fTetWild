// The slice of geogram/basic/process.h that the vendored spatial sort uses.
//
// geogram runs these lambdas on its own ThreadGroup. Each one writes a distinct member and the
// ranges they sort are disjoint, so the result does not depend on how the work is spread; only
// the wall clock does.

#pragma once

#include "geo_basic.h"

#include <functional>
#include <thread>
#include <vector>

namespace floatTetWild {
namespace geo {

    namespace impl {
        inline void run_parallel(std::vector<std::function<void()>> tasks) {
            std::vector<std::thread> threads;
            threads.reserve(tasks.size() - 1);
            for(size_t i = 1; i < tasks.size(); ++i) {
                threads.emplace_back(tasks[i]);
            }
            tasks[0]();
            for(std::thread& thread : threads) {
                thread.join();
            }
        }
    }

    inline void parallel(
        std::function<void()> f1,
        std::function<void()> f2
    ) {
        impl::run_parallel({f1, f2});
    }

    inline void parallel(
        std::function<void()> f1,
        std::function<void()> f2,
        std::function<void()> f3,
        std::function<void()> f4
    ) {
        impl::run_parallel({f1, f2, f3, f4});
    }

    inline void parallel(
        std::function<void()> f1,
        std::function<void()> f2,
        std::function<void()> f3,
        std::function<void()> f4,
        std::function<void()> f5,
        std::function<void()> f6,
        std::function<void()> f7,
        std::function<void()> f8
    ) {
        impl::run_parallel({f1, f2, f3, f4, f5, f6, f7, f8});
    }
}
}
