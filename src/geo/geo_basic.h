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

#define GEOGRAM_SPINLOCK_INIT ATOMIC_FLAG_INIT

// geogram/basic/common.h defines this per compiler.
#define geo_restrict __restrict__

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

    namespace Numeric {
        typedef void* pointer;
        typedef int8_t int8;
        typedef int16_t int16;
        typedef int32_t int32;
        typedef int64_t int64;
        typedef uint8_t uint8;
        typedef uint16_t uint16;
        typedef uint32_t uint32;
        typedef uint64_t uint64;
        typedef float float32;
        typedef double float64;

        inline float64 max_float64() {
            return std::numeric_limits<float64>::max();
        }

        inline float64 min_float64() {
            return -max_float64();
        }

        /**
         * \brief The random engine that random_int32() draws from.
         * \details geogram keeps one default-seeded std::mt19937_64 for the whole process and
         *  never reseeds it, so the stream depends only on the order in which callers draw from
         *  it. Delaunay3d is the only caller here, and it draws to pick a starting tetrahedron
         *  for point location, so keeping the same stream keeps the same walk.
         */
        inline std::mt19937_64& random_engine() {
            static std::mt19937_64 engine;
            return engine;
        }

        /**
         * \brief Returns a 32 bits integer between 0 and RAND_MAX.
         */
        inline int32 random_int32() {
            return std::uniform_int_distribution<int32>(0, RAND_MAX)(random_engine());
        }

        /**
         * \brief Place holder for optimizing internal number representation.
         * \details There are specializations for expansion_nt and for vecng.
         */
        template <class T> inline void optimize_number_representation(T& x) {
            (void)x;
        }

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

    inline index_t max_index_t() {
        return std::numeric_limits<index_t>::max();
    }

    /** \brief The type for storing and manipulating index differences. */
    typedef int32_t signed_index_t;

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

    [[noreturn]] inline void geo_should_not_have_reached(
        const char* file, int line
    ) {
        std::cerr << "Control should not have reached this point\n"
                  << "  at " << file << ':' << line << std::endl;
        std::abort();
    }

    /**
     * \brief Placeholder that suppresses unused-parameter warnings.
     */
    template <class T>
    inline void geo_argused(const T&) {
    }
}
}

#define geo_assert(x) {                                                       \
        if(!(x)) {                                                            \
            ::floatTetWild::geo::geo_assertion_failed(#x, __FILE__, __LINE__); \
        }                                                                     \
    }

#define geo_range_assert(x, min_val, max_val) \
    geo_assert((x) >= (min_val) && (x) <= (max_val))

#define geo_assert_not_reached \
    ::floatTetWild::geo::geo_should_not_have_reached(__FILE__, __LINE__)

#ifdef GEO_DEBUG
#define geo_debug_assert(x) geo_assert(x)
#define geo_debug_range_assert(x, min_val, max_val) \
    geo_range_assert(x, min_val, max_val)
#else
#define geo_debug_assert(x)
#define geo_debug_range_assert(x, min_val, max_val)
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

        explicit vector(index_t size, const T& val) : baseclass(size, val) {
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

        T& operator[] (signed_index_t i) {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        const T& operator[] (signed_index_t i) const {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        T* data() {
            return size() == 0 ? nullptr : &(*this)[0];
        }

        const T* data() const {
            return size() == 0 ? nullptr : &(*this)[0];
        }
    };

    /**
     * \brief A half-open range of indices, so that vendored code can write
     *  `for(index_t v: M.vertices)` and `for(index_t c: M.facets.corners(f))` as it does in
     *  geogram.
     */
    class index_range {
    public:
        class iterator {
        public:
            explicit iterator(index_t i) : i_(i) {
            }

            index_t operator* () const {
                return i_;
            }

            iterator& operator++ () {
                ++i_;
                return *this;
            }

            bool operator!= (const iterator& rhs) const {
                return i_ != rhs.i_;
            }

        private:
            index_t i_;
        };

        index_range(index_t begin, index_t end) : begin_(begin), end_(end) {
        }

        iterator begin() const {
            return iterator(begin_);
        }

        iterator end() const {
            return iterator(end_);
        }

    private:
        index_t begin_;
        index_t end_;
    };

    /**
     * \brief The slice of geogram's Memory:: that the vendored expansion allocator uses.
     */
    namespace Memory {
        typedef unsigned char byte;
        typedef byte* pointer;

        inline void copy(void* to, const void* from, size_t size) {
            ::memcpy(to, from, size);
        }
    }

    /**
     * \brief The slice of geogram's Process:: that the vendored expansion allocator uses.
     * \details geogram's spinlock is a std::atomic_flag on every platform we build for, so this
     *  is the same lock, not a substitute.
     */
    namespace Process {
        typedef std::atomic_flag spinlock;

        inline void acquire_spinlock(spinlock& x) {
            while(x.test_and_set(std::memory_order_acquire)) {
                // Spin.
            }
        }

        inline void release_spinlock(spinlock& x) {
            x.clear(std::memory_order_release);
        }

        /**
         * \brief How many threads can run at once.
         * \details Only used to decide whether a tree is built in parallel. The work is split
         *  into disjoint ranges either way, so this changes the wall clock and nothing else.
         */
        inline index_t maximum_concurrent_threads() {
            unsigned int n = std::thread::hardware_concurrency();
            return n == 0 ? 1 : index_t(n);
        }
    }

    template <>
    class vector<bool> : public ::std::vector<bool> {
        typedef ::std::vector<bool> baseclass;

    public:
        vector() : baseclass() {
        }

        explicit vector(index_t size) : baseclass(size) {
        }

        explicit vector(index_t size, bool val) : baseclass(size, val) {
        }

        index_t size() const {
            return index_t(baseclass::size());
        }
    };
}
}
