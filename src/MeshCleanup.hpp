// libigl's surface cleanup, in one unit. unique_rows stays here because Simplification and
// MeshCleanup instantiate it at different types; the rest lives in MeshCleanup.cpp.
//
// From libigl (https://github.com/libigl/libigl), MPL 2.0, the same licence as fTetWild:
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>  (round, sortrows,
//                                                             remove_duplicate_vertices)
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com>  (unique_rows)
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Types.hpp"

#include <algorithm>
#include <cstddef>
#include <numeric>
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
  template <typename T>
  inline void unique_rows(
    const MatrixX<T>& A,
    MatrixX<T>& C,
    MatrixXi& IA,
    MatrixXi& IC)
  {
    const int num_rows = A.rows();
    const int num_cols = A.cols();

    // matlab's sortrows(A) as a permutation, without building the sorted copy libigl did:
    // A(IM(i),:) is the i-th row in sorted order.
    MatrixXi IM(num_rows, 1);
    std::iota(IM.data(), IM.data() + num_rows, 0);
    const auto index_less_than = [&A, num_cols](int i, int j) {
      for (int c=0; c<num_cols; c++) {
        if (A(i, c) < A(j, c)) return true;
        else if (A(j,c) < A(i,c)) return false;
      }
      return false;
    };
    std::sort(IM.data(), IM.data()+IM.size(), index_less_than);

    std::vector<int> vIA(num_rows);
    std::iota(vIA.begin(), vIA.end(), 0);

    const auto index_equal = [&A, &IM, num_cols](const int i, const int j) {
      for (int c=0; c<num_cols; c++) {
        if (A(IM(i),c) != A(IM(j),c))
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
        if(!index_equal(vIA[j], i))
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
      for(int j = 0;j<A.cols();j++) C(i,j) = A(IA(i,0),j);
    }
  }

  // Remove duplicate vertices upto a uniqueness tolerance (epsilon), remapping the given faces
  // (F) --> (SF) so that SF index SV
  //
  // Inputs:
  //  V  #V by dim list of vertex positions
  //  F  #F by dim list of face indices into V
  //  epsilon  uniqueness tolerance used coordinate-wise: 1e0 --> integer
  //    match, 1e-1 --> match up to first decimal, ... , 0 --> exact match.
  // Outputs:
  //  SV  #SV by dim new list of vertex positions
  //  SF  #F by dim list of face indices into SV
  void remove_duplicate_vertices(
    const MatrixXd& V,
    const MatrixXi& F,
    double epsilon,
    MatrixXd& SV,
    MatrixXi& SF);

  // Consistently orient faces in orientable patches using BFS
  //
  // Inputs:
  //  F  #F by 3 list of faces
  // Outputs:
  //  FF  #F by 3 list of faces (OK if same as F)
  void bfs_orient(const MatrixXi &F, MatrixXi &FF);

}  // namespace floatTetWild
