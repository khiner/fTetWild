// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_REMOVE_DUPLICATE_VERTICES_H
#define FLOATTETWILD_REMOVE_DUPLICATE_VERTICES_H

#include <floattetwild/round.h>
#include <floattetwild/unique_rows.h>

#include <Eigen/Dense>

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
  template <
    typename DerivedV,
    typename DerivedSV,
    typename DerivedSVI,
    typename DerivedSVJ>
  inline void remove_duplicate_vertices(
    const Eigen::MatrixBase<DerivedV>& V,
    const double epsilon,
    Eigen::PlainObjectBase<DerivedSV>& SV,
    Eigen::PlainObjectBase<DerivedSVI>& SVI,
    Eigen::PlainObjectBase<DerivedSVJ>& SVJ)
  {
    static_assert(
      (DerivedSVI::RowsAtCompileTime == 1 || DerivedSVI::ColsAtCompileTime == 1) &&
      (DerivedSVJ::RowsAtCompileTime == 1 || DerivedSVJ::ColsAtCompileTime == 1),
      "SVI and SVJ need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1");
    if(epsilon > 0)
    {
      // The rounded copies only feed unique_rows, so their storage order is
      // immaterial: SV is gathered from V, not from them.
      Eigen::Matrix<typename DerivedV::Scalar,
                    DerivedV::RowsAtCompileTime,
                    DerivedV::ColsAtCompileTime,
                    DerivedV::Options> rV;
      round((V/(epsilon)).eval(),rV);
      Eigen::Matrix<typename DerivedV::Scalar,
                    Eigen::Dynamic,
                    DerivedV::ColsAtCompileTime,
                    DerivedV::Options> rSV;
      unique_rows(rV,rSV,SVI,SVJ);
      SV = V(SVI.derived(),Eigen::all);
    }else
    {
      unique_rows(V,SV,SVI,SVJ);
    }
  }

  // Wrapper that also remaps given faces (F) --> (SF) so that SF index SV
  //
  // Outputs:
  //  SF  #F by dim list of face indices into SV
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedSV,
    typename DerivedSVI,
    typename DerivedSVJ,
    typename DerivedSF>
  inline void remove_duplicate_vertices(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    const double epsilon,
    Eigen::PlainObjectBase<DerivedSV>& SV,
    Eigen::PlainObjectBase<DerivedSVI>& SVI,
    Eigen::PlainObjectBase<DerivedSVJ>& SVJ,
    Eigen::PlainObjectBase<DerivedSF>& SF)
  {
    // SVI and SVJ need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1
    static_assert(
      (DerivedSVI::RowsAtCompileTime == 1 || DerivedSVI::ColsAtCompileTime == 1) &&
      (DerivedSVJ::RowsAtCompileTime == 1 || DerivedSVJ::ColsAtCompileTime == 1),
      "SVI and SVJ need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1");
    using namespace Eigen;
    using namespace std;
    remove_duplicate_vertices(V,epsilon,SV,SVI,SVJ);
    SF.resizeLike(F);
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
