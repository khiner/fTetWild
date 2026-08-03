// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_SIMPLIFICATION_H
#define FLOATTETWILD_SIMPLIFICATION_H

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Parameters.h>
#include <floattetwild/Types.hpp>

#include <vector>

namespace floatTetWild {
    void simplify(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, std::vector<int>& input_tags,
            const AABBWrapper& tree, const Parameters& params, bool skip_simplify = false);
}

#endif //FLOATTETWILD_SIMPLIFICATION_H
