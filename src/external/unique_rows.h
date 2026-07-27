// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_UNIQUE_ROWS_H
#define FLOATTETWILD_UNIQUE_ROWS_H

#include <floattetwild/Types.hpp>
#include <floattetwild/sortrows.h>

#include <algorithm>
#include <vector>

namespace floatTetWild
{
  // Act like matlab's [C,IA,IC] = unique(X,'rows')
  //
  // Inputs:
  //  A  m by n matrix whose rows are to be unique'd
  // Outputs:
  //  C  #C by n matrix of unique rows of A
  //  IA  #C index vector so that C = A(IA,:);
  //  IC  #A index vector so that A = C(IC,:);
  template <typename T, typename U>
  inline void unique_rows(
    const MatrixX<T>& A,
    MatrixX<U>& C,
    MatrixXi& IA,
    MatrixXi& IC)
  {
    MatrixXi IM;
    MatrixX<T> sortA;
    sortrows(A,IM,sortA);

    const int num_rows = sortA.rows();
    const int num_cols = sortA.cols();
    std::vector<int> vIA(num_rows);
    for(int i=0;i<num_rows;i++)
    {
      vIA[i] = i;
    }

    const auto index_equal = [&sortA, num_cols](const int i, const int j) {
      for (int c=0; c<num_cols; c++) {
        if (sortA(i,c) != sortA(j,c))
          return false;
      }
      return true;
    };
    vIA.erase(std::unique(vIA.begin(),vIA.end(),index_equal),vIA.end());

    IC.resize(A.rows(),1);
    {
      int j = 0;
      for(int i = 0;i<num_rows;i++)
      {
        if(sortA.row(vIA[j]) != sortA.row(i))
        {
          j++;
        }
        IC(IM(i,0),0) = j;
      }
    }
    const int n_unique = vIA.size();
    C.resize(n_unique,A.cols());
    IA.resize(n_unique,1);
    // Reindex IA according to IM
    for(int i = 0;i<n_unique;i++)
    {
      IA(i,0) = IM(vIA[i],0);
      for(int j = 0;j<A.cols();j++) C(i,j) = U(A(IA(i,0),j));
    }
  }
}

#endif
