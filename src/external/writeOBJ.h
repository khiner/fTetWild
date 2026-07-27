// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_WRITEOBJ_H
#define FLOATTETWILD_WRITEOBJ_H

#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <Eigen/Core>

namespace floatTetWild
{
  // Write a mesh in an ascii obj file
  //
  // Inputs:
  //  str  path to outputfile
  //  V  #V by 3 mesh vertex positions
  //  F  #F by 3|4 mesh indices into V
  // Returns true on success, false on error
  template <typename DerivedV, typename DerivedF>
  inline bool writeOBJ(
    const std::string str,
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F)
  {
    using namespace std;
    using namespace Eigen;
    assert(V.cols() == 3 && "V should have 3 columns");
    ofstream s(str);
    if(!s.is_open())
    {
      fprintf(stderr,"IOError: writeOBJ() could not open %s\n",str.c_str());
      return false;
    }
    s<<
      V.format(IOFormat(FullPrecision,DontAlignCols," ","\n","v ","","","\n"))<<
      (F.array()+1).format(IOFormat(FullPrecision,DontAlignCols," ","\n","f ","","","\n"));
    return true;
  }
}

#endif
