// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include "Types.hpp"

#include <array>
#include <vector>

// The generated cut tables, indexed by the 6-bit edge-cut pattern of a tet. An index whose
// pattern no real cut produces holds an empty list, which the insertion loop skips over.
namespace floatTetWild {
	namespace CutTable {
		const std::vector<std::vector<Vector4i>>& get_tet_confs(int idx);
		const std::vector<std::vector<Vector2i>>& get_diag_confs(int idx);
		const std::vector<std::vector<std::array<bool, 4>>>& get_surface_confs(int idx);
		const std::vector<std::vector<Vector4i>>& get_face_id_confs(int idx);
	}
}
