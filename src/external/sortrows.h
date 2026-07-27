// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Only the ascending case is vendored: unique_rows, the one caller here, never asks for descending.
#ifndef FLOATTETWILD_SORTROWS_H
#define FLOATTETWILD_SORTROWS_H

#include <floattetwild/Types.hpp>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace floatTetWild
{
  // Act like matlab's [Y,I] = sortrows(X)
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
}

#endif
