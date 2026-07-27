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
// Dynamic-by-2 matrix with dim == 2, so it lands on sort2 every time. Only that
// call survives here, so the dim argument is gone with it: rows are sorted.
#ifndef FLOATTETWILD_SORT_H
#define FLOATTETWILD_SORT_H

#include <floattetwild/Types.hpp>

#include <algorithm>
#include <cassert>

namespace floatTetWild
{
  // Sort each row of a #X by 2 matrix ascending. libigl also returned the permutation, which
  // orientable_patches, the one caller here, discarded.
  //
  // Inputs:
  //  X  m by 2 matrix whose rows are to be sorted
  // Outputs:
  //  Y  m by 2 matrix whose rows are sorted
  template <typename T>
  inline void sort2(const MatrixX<T>& X, MatrixX<T>& Y)
  {
    assert(X.cols() == 2);
    Y = X;
    for(int i = 0;i<X.rows();i++)
    {
      if(Y(i,0) > Y(i,1))
      {
        std::swap(Y(i,0), Y(i,1));
      }
    }
  }
}

#endif
