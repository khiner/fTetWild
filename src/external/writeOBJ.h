// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_WRITEOBJ_H
#define FLOATTETWILD_WRITEOBJ_H

#include <floattetwild/Types.hpp>

#include <cassert>
#include <cstdio>
#include <limits>
#include <fstream>
#include <ostream>
#include <string>

namespace floatTetWild
{
  namespace writeobj_detail
  {
    // Eigen's IOFormat(FullPrecision, DontAlignCols, " ", "\n", prefix, "", "", "\n") over a
    // matrix: a prefixed, space-separated row per line, one trailing newline, and nothing at all
    // for an empty matrix but that newline. Coefficients go out through the stream's default
    // formatting, at digits10 significant digits for a floating point matrix, which is what
    // FullPrecision resolved to.
    template <typename T>
    inline void print_rows(std::ostream& s, const MatrixX<T>& m, const char* row_prefix)
    {
      if(m.size() == 0)
      {
        s << "\n";
        return;
      }
      const bool is_integer = std::numeric_limits<T>::is_integer;
      std::streamsize old_precision = 0;
      if(!is_integer)
        old_precision = s.precision(std::numeric_limits<T>::digits10);
      for(int i = 0;i<m.rows();i++)
      {
        s << row_prefix << m.coeff(i,0);
        for(int j = 1;j<m.cols();j++)
          s << " " << m.coeff(i,j);
        if(i < m.rows()-1)
          s << "\n";
      }
      s << "\n";
      if(!is_integer)
        s.precision(old_precision);
    }
  }

  // Write a mesh in an ascii obj file
  //
  // Inputs:
  //  str  path to outputfile
  //  V  #V by 3 mesh vertex positions
  //  F  #F by 3|4 mesh indices into V
  // Returns true on success, false on error
  template <typename T>
  inline bool writeOBJ(
    const std::string str,
    const MatrixX<T>& V,
    const MatrixX<int>& F)
  {
    assert(V.cols() == 3 && "V should have 3 columns");
    std::ofstream s(str);
    if(!s.is_open())
    {
      fprintf(stderr,"IOError: writeOBJ() could not open %s\n",str.c_str());
      return false;
    }
    MatrixX<int> F1(F.rows(),F.cols());
    for(int i = 0;i<F.rows();i++)
      for(int j = 0;j<F.cols();j++)
        F1(i,j) = F(i,j)+1;
    writeobj_detail::print_rows(s,V,"v ");
    writeobj_detail::print_rows(s,F1,"f ");
    return true;
  }
}

#endif
