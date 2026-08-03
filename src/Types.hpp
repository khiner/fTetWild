// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once


#include <floattetwild/Matrix.hpp>

#include <vector>


namespace floatTetWild {
#ifdef FLOAT_TETWILD_USE_FLOAT
    typedef float Scalar;
#define SCALAR_ZERO 1e-6
#define SCALAR_ZERO_2 1e-12
#define SCALAR_ZERO_3 1e-18
#else
    typedef double Scalar;
#define SCALAR_ZERO 1e-8
#define SCALAR_ZERO_2 1e-16
#define SCALAR_ZERO_3 1e-24
#endif

    // A dynamic vector is one of these with a single column. Eigen had a separate type for that
    // shape and these do not, so what would have been a VectorXd is a MatrixXd here.
    typedef MatrixX<Scalar> MatrixXs;
    typedef MatrixX<double> MatrixXd;
    typedef MatrixX<int>    MatrixXi;

    typedef Matrix33<Scalar> Matrix3;

    typedef Vector<Scalar, 3> Vector3;
    typedef Vector<Scalar, 2> Vector2;


    typedef Vector<int, 4> Vector4i;
    typedef Vector<int, 3> Vector3i;
    typedef Vector<int, 2> Vector2i;

    // A list of points or faces as a matrix with one row each, which is the form the surface
    // cleanup and the winding number take.
    template<typename T, int N>
    inline MatrixX<T> to_matrix(const std::vector<Vector<T, N>>& rows) {
        MatrixX<T> m(int(rows.size()), N);
        for (int i = 0; i < int(rows.size()); i++)
            m.row(i) = rows[i];
        return m;
    }
}
