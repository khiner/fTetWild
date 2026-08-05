// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_EDGESPLITTING_H
#define FLOATTETWILD_EDGESPLITTING_H

#include "Mesh.hpp"
#include "AABBWrapper.h"

namespace floatTetWild {
    void edge_splitting(Mesh &mesh, const AABBWrapper& tree);
}

#endif //FLOATTETWILD_EDGESPLITTING_H
