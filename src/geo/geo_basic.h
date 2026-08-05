// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
//
// The slice of geogram/basic that the rest of the vendored code needs: numeric types, Sign,
// assertions and the index_t-sized vector. Renamed into floatTetWild::geo so a copy can coexist
// with anything else that still links geogram.

#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

namespace floatTetWild {
namespace geo {

    enum Sign {
        NEGATIVE = -1,
        ZERO = 0,
        POSITIVE = 1
    };

    template <class T>
    inline Sign geo_sgn(const T& x) {
        return (x > 0) ? POSITIVE : (
            (x < 0) ? NEGATIVE : ZERO
        );
    }

    // A 32 bit integer between 0 and RAND_MAX. geogram keeps one default-seeded std::mt19937_64
    // for the whole process and never reseeds it, so the stream depends only on the order in
    // which callers draw from it. Delaunay3d is the only caller here, and it draws to pick a
    // starting tetrahedron for point location, so keeping the same stream keeps the same walk.
    inline int32_t random_int32() {
        static std::mt19937_64 engine;
        return std::uniform_int_distribution<int32_t>(0, RAND_MAX)(engine);
    }

    template <class T>
    inline T geo_sqr(T x) {
        return x * x;
    }

    // The types for storing and manipulating indices, and for coordinate indices, and the dummy
    // index value.
    typedef uint32_t index_t;
    typedef uint8_t coord_index_t;
    static const index_t NO_INDEX = index_t(-1);

    [[noreturn]] inline void geo_assertion_failed(
        const char* condition, const char* file, int line
    ) {
        std::cerr << "Assertion failed: " << condition << '\n'
                  << "  at " << file << ':' << line << std::endl;
        std::abort();
    }

}
}

#define geo_assert(x) {                                                       \
        if(!(x)) {                                                            \
            ::floatTetWild::geo::geo_assertion_failed(#x, __FILE__, __LINE__); \
        }                                                                     \
    }

#ifdef GEO_DEBUG
#define geo_debug_assert(x) geo_assert(x)
#else
#define geo_debug_assert(x)
#endif

namespace floatTetWild {
namespace geo {

    // std::vector with an index_t-typed size(), so that the vendored geogram code that assigns
    // size() to an index_t compiles without narrowing casts everywhere.
    template <class T>
    class vector : public ::std::vector<T> {
        typedef ::std::vector<T> baseclass;

    public:
        vector() : baseclass() {
        }

        explicit vector(index_t size) : baseclass(size) {
        }

        index_t size() const {
            return index_t(baseclass::size());
        }

        T& operator[] (index_t i) {
            geo_debug_assert(i < size());
            return baseclass::operator[] (i);
        }

        const T& operator[] (index_t i) const {
            geo_debug_assert(i < size());
            return baseclass::operator[] (i);
        }
    };

}
}
