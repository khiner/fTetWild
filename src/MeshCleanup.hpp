// Vendored from libigl (https://github.com/libigl/libigl), which carries its licence per file.
// The files gathered here were all MPL 2.0, the same licence as fTetWild.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>  (round, sortrows,
//                                                             remove_duplicate_vertices)
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com>  (unique_rows)
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// The surface cleanup libigl gave fTetWild, in one unit. Simplification and MeshImprovement
// instantiate the templates, so those stay here. The rest of the cluster -- sort2,
// vertex_components and orientable_patches -- has bfs_orient as its only caller and lives in
// MeshCleanup.cpp instead.
#pragma once

#include <floattetwild/Types.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace floatTetWild
{
namespace detail
{
  // Act like matlab's [Y,I] = sortrows(X). Only the ascending case is vendored: unique_rows,
  // the one caller, never asks for descending.
  //
  // Inputs:
  //  X  m by n matrix whose rows are to be sorted
  // Outputs:
  //  Y  m by n matrix whose rows are sorted (should not be same reference as X)
  //  IX  m list of indices so that Y = X(IX,:);
  template <typename T>
  inline void sortrows(const MatrixX<T>& X, MatrixXi& IX, MatrixX<T>& Y)
  {
    const int num_rows = X.rows();
    const int num_cols = X.cols();
    Y.resize(num_rows,num_cols);
    IX.resize(num_rows,1);
    for(int i = 0;i<num_rows;i++)
    {
      IX(i) = i;
    }
    const auto index_less_than = [&X, num_cols](int i, int j) {
      for (int c=0; c<num_cols; c++) {
        if (X.coeff(i, c) < X.coeff(j, c)) return true;
        else if (X.coeff(j,c) < X.coeff(i,c)) return false;
      }
      return false;
    };
    std::sort(IX.data(), IX.data()+IX.size(), index_less_than);
    for (int j=0; j<num_cols; j++) {
        for(int i = 0;i<num_rows;i++)
        {
            Y(i,j) = X(IX(i), j);
        }
    }
  }

  // Round a scalar value
  //
  // http://stackoverflow.com/a/485549
  template <typename T>
  inline T round(const T r)
  {
    return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5);
  }

  // Round a given matrix to nearest integers
  //
  // Inputs:
  //  X  m by n matrix of scalars
  // Outputs:
  //  Y  m by n matrix of rounded integers
  template <typename T>
  inline void round(const MatrixX<T>& X, MatrixX<T>& Y)
  {
    Y.resize(X.rows(),X.cols());
    for(int i = 0;i<X.rows();i++)
    {
      for(int j = 0;j<X.cols();j++)
      {
        Y(i,j) = floatTetWild::detail::round(X(i,j));
      }
    }
  }
}  // namespace detail

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
    detail::sortrows(A,IM,sortA);

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

  // Remove duplicate vertices upto a uniqueness tolerance (epsilon)
  //
  // Inputs:
  //  V  #V by dim list of vertex positions
  //  epsilon  uniqueness tolerance used coordinate-wise: 1e0 --> integer
  //    match, 1e-1 --> match up to first decimal, ... , 0 --> exact match.
  // Outputs:
  //  SV  #SV by dim new list of vertex positions
  //  SVI #SV by 1 list of indices so SV = V(SVI,:)
  //  SVJ #V by 1 list of indices so V = SV(SVJ,:)
  template <typename T, typename U>
  inline void remove_duplicate_vertices(
    const MatrixX<T>& V,
    const double epsilon,
    MatrixX<U>& SV,
    MatrixXi& SVI,
    MatrixXi& SVJ)
  {
    if(epsilon > 0)
    {
      // The rounded copy only feeds unique_rows: SV is gathered from V, not from it.
      MatrixX<T> scaled(V.rows(), V.cols());
      for(int i = 0;i<V.rows();i++)
        for(int j = 0;j<V.cols();j++)
          scaled(i,j) = V(i,j) / T(epsilon);
      MatrixX<T> rounded, unused;
      detail::round(scaled,rounded);
      unique_rows(rounded,unused,SVI,SVJ);
      SV.resize(SVI.rows(),V.cols());
      for(int i = 0;i<SVI.rows();i++)
        for(int j = 0;j<V.cols();j++)
          SV(i,j) = U(V(SVI(i),j));
    }else
    {
      unique_rows(V,SV,SVI,SVJ);
    }
  }

  // Wrapper that also remaps given faces (F) --> (SF) so that SF index SV
  //
  // Outputs:
  //  SF  #F by dim list of face indices into SV
  template <typename T, typename U>
  inline void remove_duplicate_vertices(
    const MatrixX<T>& V,
    const MatrixXi& F,
    const double epsilon,
    MatrixX<U>& SV,
    MatrixXi& SVI,
    MatrixXi& SVJ,
    MatrixXi& SF)
  {
    remove_duplicate_vertices(V,epsilon,SV,SVI,SVJ);
    SF.resize(F.rows(),F.cols());
    for(int f = 0;f<F.rows();f++)
    {
      for(int c = 0;c<F.cols();c++)
      {
        SF(f,c) = SVJ(F(f,c));
      }
    }
  }

  // Consistently orient faces in orientable patches using BFS
  //
  // Inputs:
  //  F  #F by 3 list of faces
  // Outputs:
  //  FF  #F by 3 list of faces (OK if same as F)
  //  C  #F list of component ids
  void bfs_orient(const MatrixXi &F, MatrixXi &FF, MatrixXi &C);

}  // namespace floatTetWild
