// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
//
// The slice of geogram/basic that the rest of the vendored code needs: numeric types, Sign,
// assertions and the index_t-sized vector. Renamed into floatTetWild::geo so a copy can coexist
// with anything else that still links geogram.

#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <thread>

namespace floatTetWild {
namespace geo {

    enum Sign {
        NEGATIVE = -1,
        ZERO = 0,
        POSITIVE = 1
    };

    inline Sign geo_sgn(double x) {
        return (x > 0) ? POSITIVE : (
            (x < 0) ? NEGATIVE : ZERO
        );
    }

    inline double geo_sqr(double x) {
        return x * x;
    }

    // The types for storing and manipulating indices, and the dummy index value.
    typedef uint32_t index_t;
    static const index_t NO_INDEX = index_t(-1);

    [[noreturn]] inline void geo_assertion_failed(
        const char* condition, const char* file, int line
    ) {
        std::fprintf(stderr, "Assertion failed: %s\n  at %s:%d\n", condition, file, line);
        std::abort();
    }

    // Runs the given tasks at once: the first on the calling thread, the rest each on their own.
    // The spatial sort calls this with 2, 4 and 8 of them, the kd tree with 2.
    //
    // geogram ran these lambdas on its own ThreadGroup. Each one writes a distinct variable and
    // the ranges they sort are disjoint, so the result does not depend on how the work is spread;
    // only the wall clock does.
    template <typename F0, typename... F>
    inline void parallel(F0&& f0, F&&... fs) {
        std::thread threads[] = {std::thread(std::forward<F>(fs))...};
        f0();
        for(std::thread& thread : threads) {
            thread.join();
        }
    }

}
}

#define geo_assert(x) {                                                       \
        if(!(x)) {                                                            \
            ::floatTetWild::geo::geo_assertion_failed(#x, __FILE__, __LINE__); \
        }                                                                     \
    }

