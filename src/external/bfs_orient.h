// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FLOATTETWILD_BFS_ORIENT_H
#define FLOATTETWILD_BFS_ORIENT_H

#include <floattetwild/Types.hpp>
namespace floatTetWild {
    // Consistently orient faces in orientable patches using BFS
    //
    // F = bfs_orient(F,V);
    //
    // Inputs:
    //  F  #F by 3 list of faces
    // Outputs:
    //  FF  #F by 3 list of faces (OK if same as F)
    //  C  #F list of component ids
    //
    //
    void bfs_orient(const MatrixXi &F, MatrixXi &FF, MatrixXi &C);
}
#endif
