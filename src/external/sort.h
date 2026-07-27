// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Only sort2 is vendored. igl::sort's dim-based entry point dispatches on the
// inner dimension, and orientable_patches, the one caller here, always passes a
// Dynamic-by-2 matrix with dim == 2, so it lands on sort2 every time.
#ifndef FLOATTETWILD_SORT_H
#define FLOATTETWILD_SORT_H

#include <algorithm>
#include <Eigen/Core>

namespace floatTetWild
{
  // Special case of sort for a 2-element inner dimension
  //
  // Inputs:
  //  X  m by n matrix whose entries are to be sorted
  //  dim  dimension along which to sort: 1 sorts columns, 2 sorts rows
  //  ascending  sort ascending (true) or descending (false)
  // Outputs:
  //  Y  m by n matrix whose entries are sorted
  //  IX  m by n matrix of indices so that if dim = 1, Y(i,j) = X(IX(i,j),j),
  //    and if dim = 2, Y(i,j) = X(i,IX(i,j))
  template <typename DerivedX, typename DerivedY, typename DerivedIX>
  inline void sort2(
    const Eigen::DenseBase<DerivedX>& X,
    const int dim,
    const bool ascending,
    Eigen::PlainObjectBase<DerivedY>& Y,
    Eigen::PlainObjectBase<DerivedIX>& IX)
  {
    using namespace Eigen;
    using namespace std;
    typedef typename DerivedY::Scalar YScalar;
    Y = X.derived().template cast<YScalar>();


    // get number of columns (or rows)
    int num_outer = (dim == 1 ? X.cols() : X.rows() );
    // get number of rows (or columns)
    int num_inner = (dim == 1 ? X.rows() : X.cols() );
    assert(num_inner == 2);(void)num_inner;
    typedef typename DerivedIX::Scalar Index;
    IX.resizeLike(X);
    if(dim==1)
    {
      IX.row(0).setConstant(0);// = DerivedIX::Zero(1,IX.cols());
      IX.row(1).setConstant(1);// = DerivedIX::Ones (1,IX.cols());
    }else
    {
      IX.col(0).setConstant(0);// = DerivedIX::Zero(IX.rows(),1);
      IX.col(1).setConstant(1);// = DerivedIX::Ones (IX.rows(),1);
    }
    // loop over columns (or rows)
    for(int i = 0;i<num_outer;i++)
    {
      YScalar & a = (dim==1 ? Y(0,i) : Y(i,0));
      YScalar & b = (dim==1 ? Y(1,i) : Y(i,1));
      Index & ai = (dim==1 ? IX(0,i) : IX(i,0));
      Index & bi = (dim==1 ? IX(1,i) : IX(i,1));
      if((ascending && a>b) || (!ascending && a<b))
      {
        std::swap(a,b);
        std::swap(ai,bi);
      }
    }
  }
}

#endif
