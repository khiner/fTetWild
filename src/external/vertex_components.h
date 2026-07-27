// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Only the overload taking an adjacency structure is vendored. libigl also has a face-list
// overload, which pulls in adjacency_matrix and is not reached here, and a variant returning the
// size of each component, which no caller read.
//
// libigl took a sparse adjacency matrix and walked each column, skipping stored zeros. The
// adjacency lists here are that column with the zeros already dropped, ascending, which is the
// order the sparse iterator produced.
#ifndef FLOATTETWILD_VERTEX_COMPONENTS_H
#define FLOATTETWILD_VERTEX_COMPONENTS_H

#include <floattetwild/Types.hpp>

#include <queue>
#include <vector>

namespace floatTetWild
{
  // Compute connected components of a graph
  //
  // Inputs:
  //   A  n lists of neighbours, each ascending
  // Outputs:
  //   C  n list of component ids (starting with 0)
  inline void vertex_components(
    const std::vector<std::vector<int> >& A,
    MatrixXi& C)
  {
    const int n = A.size();
    std::vector<bool> seen(n, false);
    C.resize(n,1);
    int id = 0;
    // breadth first search
    for(int k=0; k<n; ++k)
    {
      if(seen[k])
      {
        continue;
      }
      std::queue<int> Q;
      Q.push(k);
      while(!Q.empty())
      {
        const int f = Q.front();
        Q.pop();
        if(seen[f])
        {
          continue;
        }
        seen[f] = true;
        C(f,0) = id;
        for(const int g : A[f])
        {
          if(!seen[g])
          {
            Q.push(g);
          }
        }
      }
      id++;
    }
  }
}

#endif
