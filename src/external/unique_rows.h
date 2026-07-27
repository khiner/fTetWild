// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_UNIQUE_ROWS_H
#define FLOATTETWILD_UNIQUE_ROWS_H

#include <floattetwild/sortrows.h>

#include <algorithm>
#include <vector>
#include <Eigen/Core>

namespace floatTetWild
{
  // Act like matlab's [C,IA,IC] = unique(X,'rows')
  //
  // Inputs:
  //  A  m by n matrix whose entries are to unique'd according to rows
  // Outputs:
  //  C  #C vector of unique rows in A
  //  IA  #C index vector so that C = A(IA,:);
  //  IC  #A index vector so that A = C(IC,:);
  template <typename DerivedA, typename DerivedC, typename DerivedIA, typename DerivedIC>
  inline void unique_rows(
    const Eigen::DenseBase<DerivedA>& A,
    Eigen::PlainObjectBase<DerivedC>& C,
    Eigen::PlainObjectBase<DerivedIA>& IA,
    Eigen::PlainObjectBase<DerivedIC>& IC)
  {
    // IA and IC need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1
    static_assert(
      (DerivedIA::RowsAtCompileTime == 1 || DerivedIA::ColsAtCompileTime == 1) &&
      (DerivedIC::RowsAtCompileTime == 1 || DerivedIC::ColsAtCompileTime == 1),
      "IA and IC need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1");
    using namespace std;
    using namespace Eigen;
    VectorXi IM;
    Eigen::Matrix<typename DerivedA::Scalar, DerivedA::RowsAtCompileTime, DerivedA::ColsAtCompileTime> sortA;
    sortrows(A,true,sortA,IM);


    const int num_rows = sortA.rows();
    const int num_cols = sortA.cols();
    vector<int> vIA(num_rows);
    for(int i=0;i<num_rows;i++)
    {
      vIA[i] = i;
    }

    auto index_equal =
      //[&sortA, &num_cols]
      // using & so the warnings will shut up about &num_cols (which for some
      // templates is const at compile time and thus not required but for other
      // templates is not known until runtime and thus needed to be captured.
      [&]
      (const size_t i, const size_t j) {
      for (size_t c=0; c<num_cols; c++) {
        if (sortA(i,c) != sortA(j,c))
          return false;
      }
      return true;
    };
    vIA.erase(
      std::unique(
      vIA.begin(),
      vIA.end(),
      index_equal
      ),vIA.end());

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
    const int unique_rows = vIA.size();
    C.resize(unique_rows,A.cols());
    IA.resize(unique_rows,1);
    // Reindex IA according to IM
    for(int i = 0;i<unique_rows;i++)
    {
      IA(i,0) = IM(vIA[i],0);
      C.row(i) << A.row(IA(i,0));
    }
  }
}

#endif
