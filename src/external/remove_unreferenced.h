// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
//  remove_unreferenced.h
//  Preview3D
//
//  Created by Daniele Panozzo on 17/11/11.

#ifndef FLOATTETWILD_REMOVE_UNREFERENCED_H
#define FLOATTETWILD_REMOVE_UNREFERENCED_H

#include <algorithm>
#include <Eigen/Core>

namespace floatTetWild
{
  // Remove unreferenced vertices from V, updating F accordingly
  //
  // Inputs:
  //  n  number of vertices (possibly greater than F.maxCoeff()+1)
  //  F  #F by ss list of simplices (values of -1 are quietly skipped)
  // Outputs:
  //  I  #V by 1 list of indices such that: NF = IM(F) and NT = IM(T)
  //     and V(find(IM<=size(NV,1)),:) = NV
  //  J  #NV by 1 list, such that NV = V(J,:)
  template <
    typename DerivedF,
    typename DerivedI,
    typename DerivedJ>
  inline void remove_unreferenced(
    const size_t n,
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J)
  {
    // Mark referenced vertices
    typedef Eigen::Matrix<bool,Eigen::Dynamic,1> MatrixXb;
    MatrixXb mark = MatrixXb::Zero(n,1);
    for(int i=0; i<F.rows(); ++i)
    {
      for(int j=0; j<F.cols(); ++j)
      {
        if (F(i,j) != -1)
        {
          mark(F(i,j)) = 1;
        }
      }
    }

    // Sum the occupied cells
    int newsize = mark.count();

    I.resize(n,1);
    J.resize(newsize,1);

    // Do a pass on the marked vector and remove the unreferenced vertices
    int count = 0;
    for(int i=0;i<mark.rows();++i)
    {
      if (mark(i) == 1)
      {
        I(i) = count;
        J(count) = i;
        count++;
      }
      else
      {
        I(i) = -1;
      }
    }
  }

  // \overload
  //
  // Inputs:
  //  V  #V by dim list of mesh vertex positions
  // Outputs:
  //  NV  #NV by dim list of mesh vertex positions
  //  NF  #NF by ss list of simplices
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedNV,
    typename DerivedNF,
    typename DerivedI,
    typename DerivedJ>
  inline void remove_unreferenced(
    const Eigen::MatrixBase<DerivedV> &V,
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedNV> &NV,
    Eigen::PlainObjectBase<DerivedNF> &NF,
    Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedJ> &J)
  {
    using namespace std;
    const size_t n = V.rows();
    remove_unreferenced(n,F,I,J);
    NF = F;
    std::for_each(NF.data(),NF.data()+NF.size(),
      [&I](typename DerivedNF::Scalar & a){a=I(a);});
    NV = V(J.derived(),Eigen::all);
  }

  // \overload
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedNV,
    typename DerivedNF,
    typename DerivedI>
  inline void remove_unreferenced(
    const Eigen::MatrixBase<DerivedV> &V,
    const Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedNV> &NV,
    Eigen::PlainObjectBase<DerivedNF> &NF,
    Eigen::PlainObjectBase<DerivedI> &I)
  {
    Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,1> J;
    remove_unreferenced(V,F,NV,NF,I,J);
  }
}

#endif
