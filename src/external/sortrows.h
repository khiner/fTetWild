// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_SORTROWS_H
#define FLOATTETWILD_SORTROWS_H

#include <algorithm>
#include <Eigen/Core>

namespace floatTetWild
{
  // Act like matlab's [Y,I] = sortrows(X)
  //
  // Inputs:
  //  X  m by n matrix whose entries are to be sorted
  //  ascending  sort ascending (true, matlab default) or descending (false)
  // Outputs:
  //  Y  m by n matrix whose entries are sorted (should not be same reference as X)
  //  IX  m list of indices so that Y = X(IX,:);
  template <typename DerivedX, typename DerivedY, typename DerivedIX>
  inline void sortrows(
    const Eigen::DenseBase<DerivedX>& X,
    const bool ascending,
    Eigen::PlainObjectBase<DerivedY>& Y,
    Eigen::PlainObjectBase<DerivedIX>& IX)
  {
    // This is already 2x faster than matlab's builtin `sortrows`. I have tried
    // implementing a "multiple-pass" sort on each column, but see no performance
    // improvement.
    using namespace std;
    using namespace Eigen;
    // Resize output
    const size_t num_rows = X.rows();
    const size_t num_cols = X.cols();
    Y.resize(num_rows,num_cols);
    IX.resize(num_rows,1);
    for(int i = 0;i<num_rows;i++)
    {
      IX(i) = i;
    }
    if (ascending) {
      auto index_less_than = [&X, num_cols](size_t i, size_t j) {
        for (size_t c=0; c<num_cols; c++) {
          if (X.coeff(i, c) < X.coeff(j, c)) return true;
          else if (X.coeff(j,c) < X.coeff(i,c)) return false;
        }
        return false;
      };
        std::sort(
          IX.data(),
          IX.data()+IX.size(),
          index_less_than
          );
    } else {
      auto index_greater_than = [&X, num_cols](size_t i, size_t j) {
        for (size_t c=0; c<num_cols; c++) {
          if (X.coeff(i, c) > X.coeff(j, c)) return true;
          else if (X.coeff(j,c) > X.coeff(i,c)) return false;
        }
        return false;
      };
        std::sort(
          IX.data(),
          IX.data()+IX.size(),
          index_greater_than
          );
    }
    for (size_t j=0; j<num_cols; j++) {
        for(int i = 0;i<num_rows;i++)
        {
            Y(i,j) = X(IX(i), j);
        }
    }
  }
}

#endif
