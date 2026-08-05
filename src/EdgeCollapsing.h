// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_EDGECOLLAPSING_H
#define FLOATTETWILD_EDGECOLLAPSING_H

#include "Mesh.hpp"
#include "AABBWrapper.h"

namespace floatTetWild {
    void edge_collapsing(Mesh& mesh, const AABBWrapper& tree);

    // tet_tss, when given, receives ts + 1 on every tet the collapse rewires, which is how
    // edge_collapsing tells whether anything moved near a previously failed edge.
    bool collapse_an_edge(Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree,
             std::vector<std::array<int, 2>>& new_edges, int ts, std::vector<int>* tet_tss);
}

#endif //FLOATTETWILD_EDGECOLLAPSING_H
