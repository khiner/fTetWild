// From libigl (https://github.com/libigl/libigl), MPL 2.0, the same licence as fTetWild:
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>  (sort, orientable_patches, bfs_orient)
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>  (vertex_components)
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <floattetwild/MeshCleanup.hpp>

#include <algorithm>
#include <cassert>
#include <queue>
#include <vector>

namespace floatTetWild {
namespace {

// Sort each row of a #X by 2 matrix ascending. igl::sort dispatched on the inner dimension and
// returned the permutation too; orientable_patches, the one caller, always passes dim 2 and
// discards the permutation.
void sort2(const MatrixXi& X, MatrixXi& Y)
{
    assert(X.cols() == 2);
    Y = X;
    for (int i = 0; i < X.rows(); i++) {
        if (Y(i, 0) > Y(i, 1)) {
            std::swap(Y(i, 0), Y(i, 1));
        }
    }
}

// Compute connected components of a graph. libigl walked a sparse adjacency matrix column by
// column, skipping stored zeros; A is that column with the zeros already dropped, in the same
// ascending order the sparse iterator produced.
//
// Inputs:
//   A  n lists of neighbours, each ascending
// Outputs:
//   C  n list of component ids (starting with 0)
void vertex_components(const std::vector<std::vector<int>>& A, MatrixXi& C)
{
    const int n = A.size();
    std::vector<bool> seen(n, false);
    C.resize(n, 1);
    int id = 0;
    for (int k = 0; k < n; ++k) {
        if (seen[k]) {
            continue;
        }
        std::queue<int> Q;
        Q.push(k);
        while (!Q.empty()) {
            const int f = Q.front();
            Q.pop();
            if (seen[f]) {
                continue;
            }
            seen[f] = true;
            C(f, 0) = id;
            for (const int g : A[f]) {
                if (!seen[g]) {
                    Q.push(g);
                }
            }
        }
        id++;
    }
}

// Compute connected components of facets connected by manifold edges.
//
// libigl built face-face adjacency as a sparse uE2FT * uE2FT^T, leaving non-manifold edges in as
// stored zeros for callers to skip. A stored zero only ever meant "skip me", so A holds what those
// callers kept: faces sharing at least one manifold edge, ascending, self excluded.
//
// Inputs:
//   F  #F by 3 list of facets
// Outputs:
//   C  #F list of component ids
//   A  #F lists of adjacent facets, each ascending
void orientable_patches(const MatrixXi& F, MatrixXi& C, std::vector<std::vector<int>>& A)
{
    assert(F.cols() == 3);
    const int nf = F.rows();

    // List of all "half"-edges: 3*#F by 2
    MatrixXi allE(nf * 3, 2);
    for (int f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            allE(e * nf + f, 0) = F(f, (e + 1) % 3);
            allE(e * nf + f, 1) = F(f, (e + 2) % 3);
        }
    }
    MatrixXi sortallE;
    sort2(allE, sortallE);
    // IC(i) tells us where to find sortallE(i,:) in uE:
    // so that sortallE(i,:) = uE(IC(i),:)
    MatrixXi uE, IA, IC;
    unique_rows(sortallE, uE, IA, IC);

    // Faces incident to each unique edge, each face listed once, ascending.
    std::vector<std::vector<int>> edge_faces(uE.rows());
    for (int e = 0; e < IC.rows(); e++) {
        edge_faces[IC(e)].push_back(e % nf);
    }
    for (auto& faces : edge_faces) {
        std::sort(faces.begin(), faces.end());
        faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    }

    A.assign(nf, std::vector<int>());
    for (const auto& faces : edge_faces) {
        // Non-manifold edges join nothing.
        if (faces.size() > 2)
            continue;
        for (size_t i = 0; i < faces.size(); i++)
            for (size_t j = 0; j < faces.size(); j++)
                if (i != j)
                    A[faces[i]].push_back(faces[j]);
    }
    for (auto& neighbours : A) {
        std::sort(neighbours.begin(), neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
    }

    vertex_components(A, C);
}

}  // namespace

void bfs_orient(const MatrixXi &F, MatrixXi &FF, MatrixXi &C) {
    std::vector<std::vector<int>> A;
    orientable_patches(F, C, A);

    const int m = F.rows();
    const int num_cc = C.maxCoeff() + 1;
    std::vector<int> seen(m, 0);

    const int ES[3][2] = {{1, 2},
                          {2, 0},
                          {0, 1}};

    FF = F;

    // loop over patches. libigl ran this under `#pragma omp parallel for`, dropped with the rest of
    // the OpenMP in this tree.
    for (int c = 0; c < num_cc; c++) {
        std::queue<int> Q;
        int cnt = 0;
        for (int f = 0; f < FF.rows(); f++) {
            if (C(f) == c) {
                if (cnt == 0)
                    Q.push(f);
                cnt++;
            }
        }
        if (cnt < 5)
            continue;

        int cnt_inverted = 0;
        assert(!Q.empty());
        while (!Q.empty()) {
            const int f = Q.front();
            Q.pop();
            if (seen[f] > 0)
                continue;

            seen[f]++;
            for (const int n : A[f]) {
                for (int efi = 0; efi < 3; efi++) {
                    const Vector2i ef(FF(f, ES[efi][0]), FF(f, ES[efi][1]));
                    for (int eni = 0; eni < 3; eni++) {
                        const Vector2i en(FF(n, ES[eni][0]), FF(n, ES[eni][1]));
                        // Match (half-edges go same direction)
                        if (ef(0) == en(0) && ef(1) == en(1)) {
                            std::swap(FF(n, 0), FF(n, 2));
                            cnt_inverted++;
                        }
                    }
                }
                Q.push(n);
            }
        }
        if (cnt_inverted < cnt / 2)
            continue;

        for (int f = 0; f < FF.rows(); f++) {
            if (C(f) == c)
                std::swap(FF(f, 0), FF(f, 2));
        }
    }
}

}  // namespace floatTetWild
