// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// libigl built the face-face adjacency as a sparse product uE2FT * uE2FT^T and left the
// non-manifold edges in as stored zeros for the callers to skip. The adjacency lists built here
// hold exactly the entries those callers kept: faces sharing at least one manifold edge, ascending,
// self excluded. A stored zero only ever meant "skip me", and a pair sharing both a manifold and a
// non-manifold edge stayed adjacent, which is why the manifold edges alone decide membership.
#ifndef FLOATTETWILD_ORIENTABLE_PATCHES_H
#define FLOATTETWILD_ORIENTABLE_PATCHES_H

#include <floattetwild/Types.hpp>
#include <floattetwild/sort.h>
#include <floattetwild/unique_rows.h>
#include <floattetwild/vertex_components.h>

#include <algorithm>
#include <cassert>
#include <vector>

namespace floatTetWild
{
  // Compute connected components of facets connected by manifold edges.
  //
  // Inputs:
  //   F  #F by 3 list of facets
  // Outputs:
  //   C  #F list of component ids
  //   A  #F lists of adjacent facets, each ascending
  inline void orientable_patches(
    const MatrixXi& F,
    MatrixXi& C,
    std::vector<std::vector<int> >& A)
  {
    assert(F.cols() == 3);
    const int nf = F.rows();

    // List of all "half"-edges: 3*#F by 2
    MatrixXi allE(nf*3,2);
    for(int f = 0;f<nf;f++)
    {
      allE(0*nf+f,0) = F(f,1); allE(0*nf+f,1) = F(f,2);
      allE(1*nf+f,0) = F(f,2); allE(1*nf+f,1) = F(f,0);
      allE(2*nf+f,0) = F(f,0); allE(2*nf+f,1) = F(f,1);
    }
    // Sort each row
    MatrixXi sortallE;
    sort2(allE,sortallE);
    //IC(i) tells us where to find sortallE(i,:) in uE:
    // so that sortallE(i,:) = uE(IC(i),:)
    MatrixXi uE, IA, IC;
    unique_rows(sortallE,uE,IA,IC);

    // Faces incident to each unique edge, each face listed once, ascending.
    std::vector<std::vector<int> > edge_faces(uE.rows());
    for(int e = 0;e<IC.rows();e++)
    {
      edge_faces[IC(e)].push_back(e%nf);
    }
    for(auto& faces : edge_faces)
    {
      std::sort(faces.begin(),faces.end());
      faces.erase(std::unique(faces.begin(),faces.end()),faces.end());
    }

    A.assign(nf,std::vector<int>());
    for(const auto& faces : edge_faces)
    {
      // Non-manifold edges join nothing.
      if(faces.size() > 2)
        continue;
      for(size_t i = 0;i<faces.size();i++)
        for(size_t j = 0;j<faces.size();j++)
          if(i != j)
            A[faces[i]].push_back(faces[j]);
    }
    for(auto& neighbours : A)
    {
      std::sort(neighbours.begin(),neighbours.end());
      neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
    }

    // graph connected components
    vertex_components(A,C);
  }
}

#endif
