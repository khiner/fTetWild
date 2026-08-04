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
#include <string>

namespace floatTetWild
{
  // Write a mesh in an ascii obj file
  //
  // Inputs:
  //  path  path to outputfile
  //  V  #V by 3 mesh vertex positions
  //  F  #F by 3 mesh indices into V
  // Returns true on success, false on error
  template <typename T>
  inline bool writeOBJ(
    const std::string& path,
    const MatrixX<T>& V,
    const MatrixXi& F)
  {
    assert(V.cols() == 3 && "V should have 3 columns");
    assert(F.cols() == 3 && "F should contain triangles");
    std::ofstream s(path);
    if(!s.is_open())
    {
      fprintf(stderr,"IOError: writeOBJ() could not open %s\n",path.c_str());
      return false;
    }
    // Eigen wrote each block through IOFormat(FullPrecision, DontAlignCols, " ", "\n", prefix,
    // "", "", "\n"): a prefixed, space-separated row per line, and for an empty block nothing at
    // all but the trailing newline. FullPrecision resolved to digits10 significant digits, which
    // is what the vertices go out at; the indices are integers and ignore it.
    s.precision(std::numeric_limits<T>::digits10);
    if(V.rows() == 0) s << "\n";
    for(int i = 0;i<V.rows();i++)
      s << "v " << V(i,0) << " " << V(i,1) << " " << V(i,2) << "\n";
    if(F.rows() == 0) s << "\n";
    for(int i = 0;i<F.rows();i++)
      s << "f " << F(i,0)+1 << " " << F(i,1)+1 << " " << F(i,2)+1 << "\n";
    return true;
  }
}

#endif
