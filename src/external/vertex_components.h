// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Only the overloads taking an adjacency matrix are vendored. libigl also has a
// face-list overload, which pulls in adjacency_matrix and is not reached here.
#ifndef FLOATTETWILD_VERTEX_COMPONENTS_H
#define FLOATTETWILD_VERTEX_COMPONENTS_H

#include <cassert>
#include <queue>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace floatTetWild
{
  // Compute connected components of a graph represented by a sparse adjacency matrix
  //
  // Inputs:
  //   A  n by n sparse adjacency matrix
  // Outputs:
  //   C  n list of component ids (starting with 0)
  //   counts  #components list of counts for each component
  template <typename DerivedA, typename DerivedC, typename Derivedcounts>
  inline void vertex_components(
    const Eigen::SparseCompressedBase<DerivedA> & A,
    Eigen::PlainObjectBase<DerivedC> & C,
    Eigen::PlainObjectBase<Derivedcounts> & counts)
  {
    using namespace Eigen;
    using namespace std;
    assert(A.rows() == A.cols() && "A should be square.");
    const size_t n = A.rows();
    Array<bool,Dynamic,1> seen = Array<bool,Dynamic,1>::Zero(n,1);
    C.resize(n,1);
    typename DerivedC::Scalar id = 0;
    vector<typename Derivedcounts::Scalar> vcounts;
    // breadth first search
    for(int k=0; k<A.outerSize(); ++k)
    {
      if(seen(k))
      {
        continue;
      }
      queue<int> Q;
      Q.push(k);
      vcounts.push_back(0);
      while(!Q.empty())
      {
        const int f = Q.front();
        Q.pop();
        if(seen(f))
        {
          continue;
        }
        seen(f) = true;
        C(f,0) = id;
        vcounts[id]++;
        // Iterate over inside
        for(typename DerivedA::InnerIterator it (A,f); it; ++it)
        {
          const int g = it.index();
          if(!seen(g) && it.value())
          {
            Q.push(g);
          }
        }
      }
      id++;
    }
    assert((size_t) id == vcounts.size());
    const size_t ncc = vcounts.size();
    assert((size_t)C.maxCoeff()+1 == ncc);
    counts.resize(ncc,1);
    for(size_t i = 0;i<ncc;i++)
    {
      counts(i) = vcounts[i];
    }
  }

  // \overload
  template <typename DerivedA, typename DerivedC>
  inline void vertex_components(
    const Eigen::SparseCompressedBase<DerivedA> & A,
    Eigen::PlainObjectBase<DerivedC> & C)
  {
    Eigen::VectorXi counts;
    return vertex_components(A,C,counts);
  }
}

#endif
