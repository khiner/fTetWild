// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
//
// The slice of geogram/basic that the rest of the vendored code needs: numeric types, Sign,
// assertions and the index_t-sized vector. Renamed into floatTetWild::geo so a copy can coexist
// with anything else that still links geogram.

#pragma once

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <thread>
#include <vector>

namespace floatTetWild {
namespace geo {

    /**
     * \brief Integer constants that represent the sign of a value
     */
    enum Sign {
        NEGATIVE = -1,
        ZERO = 0,
        POSITIVE = 1
    };

    /**
     * \brief Gets the sign of a value
     */
    template <class T>
    inline Sign geo_sgn(const T& x) {
        return (x > 0) ? POSITIVE : (
            (x < 0) ? NEGATIVE : ZERO
        );
    }

    /**
     * \brief Returns a 32 bits integer between 0 and RAND_MAX.
     * \details geogram keeps one default-seeded std::mt19937_64 for the whole process and never
     *  reseeds it, so the stream depends only on the order in which callers draw from it.
     *  Delaunay3d is the only caller here, and it draws to pick a starting tetrahedron for point
     *  location, so keeping the same stream keeps the same walk.
     */
    inline int32_t random_int32() {
        static std::mt19937_64 engine;
        return std::uniform_int_distribution<int32_t>(0, RAND_MAX)(engine);
    }

    /**
     * \brief Gets the square value of a value
     */
    template <class T>
    inline T geo_sqr(T x) {
        return x * x;
    }

    /** \brief The type for storing and manipulating indices. */
    typedef uint32_t index_t;

    /** \brief The type for storing coordinate indices. */
    typedef uint8_t coord_index_t;

    /** \brief The dummy index value. */
    static const index_t NO_INDEX = index_t(-1);

    /**
     * \brief Prints an assertion failure and aborts.
     */
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

    /**
     * \brief std::vector with an index_t-typed size(), so that the vendored geogram code that
     *  assigns size() to an index_t compiles without narrowing casts everywhere.
     */
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

        // Null rather than unspecified when empty, which is what geogram's callers assume.
        T* data() {
            return size() == 0 ? nullptr : &(*this)[0];
        }

        const T* data() const {
            return size() == 0 ? nullptr : &(*this)[0];
        }
    };

}
}
