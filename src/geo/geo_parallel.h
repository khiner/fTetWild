// The slice of geogram/basic/process.h that the vendored spatial sort and kd tree use.
//
// geogram runs these lambdas on its own ThreadGroup. Each one writes a distinct variable and the
// ranges they sort are disjoint, so the result does not depend on how the work is spread; only
// the wall clock does.

#pragma once

#include <functional>
#include <thread>
#include <vector>

namespace floatTetWild {
namespace geo {

    // Runs the given tasks at once. The spatial sort calls this with 2, 4 and 8 of them.
    template <typename... F>
    inline void parallel(F&&... fs) {
        std::vector<std::function<void()>> tasks = {
            std::function<void()>(std::forward<F>(fs))...};
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
}
