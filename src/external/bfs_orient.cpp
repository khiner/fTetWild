// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include <floattetwild/bfs_orient.h>
#include <floattetwild/orientable_patches.h>
#include <queue>
#include <vector>

void floatTetWild::bfs_orient(const MatrixXi &F, MatrixXi &FF, MatrixXi &C) {
    std::vector<std::vector<int>> A;
    orientable_patches(F, C, A);

    // number of faces
    const int m = F.rows();
    // number of patches
    const int num_cc = C.maxCoeff() + 1;
    std::vector<int> seen(m, 0);

    // Edge sets
    const int ES[3][2] = {{1, 2},
                          {2, 0},
                          {0, 1}};

    if (((void *) &FF) != ((void *) &F))
        FF = F;

    // loop over patches
#pragma omp parallel for
    for (int c = 0; c < num_cc; c++) {
        std::queue<int> Q;
        // find first member of patch c
        int cnt = 0;
        for (int f = 0; f < FF.rows(); f++) {
            if (C(f) == c) {
                if (cnt == 0)
                    Q.push(f);
                cnt++;
//                break;
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
            // loop over neighbors of f
            for (const int n : A[f]) {
                // loop over edges of f
                for (int efi = 0; efi < 3; efi++) {
                    // efi'th edge of face f
                    const Vector2i ef(FF(f, ES[efi][0]), FF(f, ES[efi][1]));
                    // loop over edges of n
                    for (int eni = 0; eni < 3; eni++) {
                        // eni'th edge of face n
                        const Vector2i en(FF(n, ES[eni][0]), FF(n, ES[eni][1]));
                        // Match (half-edges go same direction)
                        if (ef(0) == en(0) && ef(1) == en(1)) {
                            // flip face n
                            std::swap(FF(n, 0), FF(n, 2));
                            cnt_inverted++;
                        }
                    }
                }
                // add neighbor to queue
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