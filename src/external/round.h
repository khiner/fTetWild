// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_ROUND_H
#define FLOATTETWILD_ROUND_H

#include <cmath>
#include <Eigen/Dense>

namespace floatTetWild
{
  // Round a scalar value
  //
  // http://stackoverflow.com/a/485549
  template <typename DerivedX>
  inline DerivedX round(const DerivedX r)
  {
    return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5);
  }

  // Round a given matrix to nearest integers
  //
  // Inputs:
  //  X  m by n matrix of scalars
  // Outputs:
  //  Y  m by n matrix of rounded integers
  template <typename DerivedX, typename DerivedY>
  inline void round(
    const Eigen::MatrixBase<DerivedX>& X,
    Eigen::PlainObjectBase<DerivedY>& Y)
  {
    Y.resizeLike(X);
    // loop over rows
    for(int i = 0;i<X.rows();i++)
    {
      // loop over cols
      for(int j = 0;j<X.cols();j++)
      {
        Y(i,j) = floatTetWild::round(X(i,j));
      }
    }
  }
}

#endif
