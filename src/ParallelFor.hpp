// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_PARALLELFOR_HPP
#define FLOATTETWILD_PARALLELFOR_HPP

#include <cstddef>
#include <functional>
#include <vector>

namespace floatTetWild {

// Number of threads the loops below may run on, the calling thread included. Set it once before any
// parallel work. 1 runs everything inline, through the same code a parallel run takes rather than a
// separate serial branch. Without a call, one thread per core.
//
// The loops below expect one driver thread, which is how the mesher runs them. A loop reached from
// inside another runs inline rather than nesting, but two threads starting a loop at the same time
// is not supported.
void set_num_threads(unsigned int n);

namespace detail {
// Splits [begin, end) into contiguous chunks and runs body(lo, hi, chunk) on each. Chunks are
// claimed dynamically, so which thread runs which chunk varies, but the split itself depends
// only on the range and the thread count, never on scheduling.
void parallel_ranges(size_t                                             begin,
                     size_t                                             end,
                     const std::function<void(size_t, size_t, size_t)>& body);

// How many chunks parallel_ranges makes for this range.
size_t chunk_count(size_t begin, size_t end);
}  // namespace detail

// Runs f(i) for every i in [begin, end). A range shorter than min_parallel runs inline, for the
// loops whose bodies are cheap enough that handing chunks to the pool costs more than the work.
template<typename F>
void parallel_for(size_t begin, size_t end, F&& f, size_t min_parallel = 0)
{
    if (end > begin && end - begin < min_parallel) {
        for (size_t i = begin; i < end; ++i)
            f(i);
        return;
    }
    detail::parallel_ranges(begin, end, [&f](size_t lo, size_t hi, size_t /*chunk*/) {
        for (size_t i = lo; i < hi; ++i)
            f(i);
    });
}

// Runs f(i, out) for every i in [begin, end), where f appends whatever index i contributes to out,
// and returns everything appended.
//
// Each chunk fills its own vector and the vectors are concatenated in chunk order, so the result is
// in index order however the threads interleaved. That is the order a single-threaded run produces,
// which is what the callers that feed a queue or a unique() rely on to stay reproducible.
template<typename T, typename F>
std::vector<T> parallel_collect(size_t begin, size_t end, F&& f)
{
    std::vector<std::vector<T>> per_chunk(detail::chunk_count(begin, end));
    detail::parallel_ranges(begin, end, [&](size_t lo, size_t hi, size_t chunk) {
        std::vector<T>& out = per_chunk[chunk];
        for (size_t i = lo; i < hi; ++i)
            f(i, out);
    });

    size_t total = 0;
    for (const std::vector<T>& chunk : per_chunk)
        total += chunk.size();

    std::vector<T> result;
    result.reserve(total);
    for (const std::vector<T>& chunk : per_chunk)
        result.insert(result.end(), chunk.begin(), chunk.end());
    return result;
}

}  // namespace floatTetWild

#endif  // FLOATTETWILD_PARALLELFOR_HPP
