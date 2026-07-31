// Vendored from geogram (https://github.com/BrunoLevy/geogram), Bruno Levy, INRIA.
// Original licence: BSD 3-clause, see LICENSE.geogram next to this file.
// Source: geogram/basic/algorithm.h
//
// The fixed seed used to be applied to geogram's own copy by
// cmake/patches/geogram_deterministic_shuffle.cmake. It is baked in here instead: geogram seeds a
// fresh std::mt19937 from std::random_device on every call, so the BRIO insertion order
// compute_BRIO_order picks would differ every run and every stage downstream of
// FloatTetDelaunay::tetrahedralize would diverge with it. The order is still drawn from the same
// distribution, just pinned to one draw.

#pragma once

#include "geo_basic.h"

#include <algorithm>
#include <random>

namespace floatTetWild {
namespace geo {

    /**
     * \brief Applies a random permutation to a sequence.
     * \details A drop-in replacement of std::random_shuffle(), that is deprecated since c++17.
     * \param[in] begin , end first position and one item past last position of the sequence to
     *  be randomly permuted
     */
    template <typename ITERATOR>
    inline void random_shuffle(
        const ITERATOR& begin, const ITERATOR& end
    ) {
        std::mt19937 urng(42);
        std::shuffle(begin, end, urng);
    }
}
}
