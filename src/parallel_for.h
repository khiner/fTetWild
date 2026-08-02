// igl::parallel_for's signature, over the mesher's thread pool. The signature is libigl's, MPL 2.0,
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>; the implementation is not. Keeping it
// is what lets the call sites in FastWindingNumberForSoups.h stay untouched.
//
// The four-argument form gives the body a slot id and prep_func a slot count. Slots are chunks here
// rather than threads, so the count differs from libigl's. Safe because every caller either writes
// disjointly by loop index or reduces exactly and order independently, over min/max box unions and
// integer span counts.
//
// One trap: BVH::computeFullBoundingBox never calls initBounds() on its per-slot boxes, and
// UT_Array zero-fills them because Box is POD, so the root box always contains the origin. That is
// upstream's bug and it stays -- it is invariant under slot count, and initBounds() would move the
// tree and change the mesh.
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_PARALLEL_FOR_H
#define FLOATTETWILD_PARALLEL_FOR_H

#include <floattetwild/ParallelFor.hpp>

#include <cassert>
#include <cstddef>

namespace floatTetWild
{
  // Runs func(i, slot) over [0, loop_size), after prep_func(nslots) and before one accum_func(slot)
  // per slot. Returns true iff it ran on more than one slot.
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

    const size_t n = size_t(loop_size);
    // 1 inside a parallel region, which is what runs a nested loop inline instead of re-entering
    // the pool. parallel_ranges() below asks again and gets the same answer.
    const size_t nslots = n < min_parallel ? 1 : detail::chunk_count(0, n);

    if(nslots <= 1)
    {
      prep_func(1);
      for(Index i = 0;i<loop_size;i++) func(i,0);
      accum_func(0);
      return false;
    }

    prep_func(nslots);
    detail::parallel_ranges(0, n, [&func](size_t lo, size_t hi, size_t slot)
    {
      for(size_t i = lo;i<hi;++i) func(Index(i),slot);
    });
    for(size_t slot = 0;slot<nslots;slot++)
    {
      accum_func(slot);
    }
    return true;
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
    // no op preparation/accumulation
    const auto & no_op = [](const size_t /*n/t*/){};
    // two-parameter wrapper ignoring the slot id
    const auto & wrapper = [&func](Index i,size_t /*t*/){ func(i); };
    return parallel_for(loop_size,no_op,wrapper,no_op,min_parallel);
  }
}

#endif
