// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_REMOVE_DUPLICATE_VERTICES_H
#define FLOATTETWILD_REMOVE_DUPLICATE_VERTICES_H

#include <floattetwild/Types.hpp>
#include <floattetwild/round.h>
#include <floattetwild/unique_rows.h>

namespace floatTetWild
{
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
      round(scaled,rounded);
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
}

#endif
