// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Only the triangle-soup overloads fTetWild reaches are kept. libigl also has
// point-cloud and octree based overloads, which are not vendored.
#ifndef FLOATTETWILD_FAST_WINDING_NUMBER_H
#define FLOATTETWILD_FAST_WINDING_NUMBER_H

// Third-party header, so there is no source file to put compile options on. The first two lines
// let each compiler skip the other's warning names.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wgnu-anonymous-struct"
#pragma GCC diagnostic ignored "-Winvalid-offsetof"
#pragma GCC diagnostic ignored "-Wshadow"
#include <floattetwild/FastWindingNumberForSoups.h>
#pragma GCC diagnostic pop

#include <floattetwild/parallel_for.h>

#include <floattetwild/Types.hpp>

#include <cassert>
#include <cmath>
#include <vector>

namespace floatTetWild
{
  // igl::PI, which is M_PI when the platform defines it.
#ifdef M_PI
  constexpr double FastWindingNumberPI = M_PI;
#else
  constexpr double FastWindingNumberPI = 3.1415926535897932384626433832795;
#endif

  // Structure for caching precomputation for fast winding number for triangle soups
  struct FastWindingNumberBVH {
    FastWindingNumber::HDK_Sample::UT_SolidAngle<float,float> ut_solid_angle;
    // Need copies of these so they stay alive between calls.
    std::vector<FastWindingNumber::HDK_Sample::UT_Vector3T<float> > U;
    std::vector<int> F;
  };

  // Precomputation for computing approximate winding numbers of a triangle soup.
  //
  // Inputs:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of triangle mesh indices into rows of V
  //   order  Taylor series expansion order to use (e.g., 2)
  // Outputs:
  //   fwn_bvh  Precomputed bounding volume hierarchy
  inline void fast_winding_number(
    const MatrixXd & V,
    const MatrixXi & F,
    const int order,
    FastWindingNumberBVH & fwn_bvh)
  {
    assert(V.cols() == 3 && "V should be 3D");
    assert(F.cols() == 3 && "F should contain triangles");
    // Extra copies. Usuually this won't be the bottleneck.
    fwn_bvh.U.resize(V.rows());
    for(int i = 0;i<V.rows();i++)
    {
      for(int j = 0;j<3;j++)
      {
        fwn_bvh.U[i][j] = V(i,j);
      }
    }
    // Wouldn't need to copy if F is **RowMajor**
    fwn_bvh.F.resize(F.size());
    for(int f = 0;f<F.rows();f++)
    {
      for(int c = 0;c<F.cols();c++)
      {
        fwn_bvh.F[c+f*F.cols()] = F(f,c);
      }
    }
    fwn_bvh.ut_solid_angle.clear();
    fwn_bvh.ut_solid_angle.init(
       fwn_bvh.F.size()/3,
      &fwn_bvh.F[0],
       fwn_bvh.U.size(),
      &fwn_bvh.U[0],
      order);
  }

  // After precomputation, compute winding number at each of many points in a list.
  //
  // Inputs:
  //   fwn_bvh  Precomputed bounding volume hierarchy
  //   accuracy_scale  parameter controlling accuracy (e.g., 2)
  //   Q  #Q by 3 list of query positions
  // Outputs:
  //   W  #Q list of winding number values
  inline void fast_winding_number(
    const FastWindingNumberBVH & fwn_bvh,
    const float accuracy_scale,
    const MatrixXd & Q,
    MatrixXd & W)
  {
    assert(Q.cols() == 3 && "Q should be 3D");
    W.resize(Q.rows(),1);
    floatTetWild::parallel_for(Q.rows(),[&](int p)
    {
      FastWindingNumber::HDK_Sample::UT_Vector3T<float>Qp;
      Qp[0] = Q(p,0);
      Qp[1] = Q(p,1);
      Qp[2] = Q(p,2);
      W(p) = fwn_bvh.ut_solid_angle.computeSolidAngle(Qp,accuracy_scale) / (4.0*FastWindingNumberPI);
    },1000);
  }

  // Compute approximate winding number of a triangle soup mesh according to
  // "Fast Winding Numbers for Soups and Clouds" [Barill et al. 2018].
  //
  // Inputs:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of triangle mesh indices into rows of V
  //   Q  #Q by 3 list of query positions
  // Outputs:
  //   W  #Q list of winding number values
  inline void fast_winding_number(
    const MatrixXd & V,
    const MatrixXi & F,
    const MatrixXd & Q,
    MatrixXd & W)
  {
    FastWindingNumberBVH fwn_bvh;
    int order = 2;
    fast_winding_number(V,F,order,fwn_bvh);
    float accuracy_scale = 2;
    fast_winding_number(fwn_bvh,accuracy_scale,Q,W);
  }
}

#endif
