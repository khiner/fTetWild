// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Mesh.hpp>
#include <floattetwild/Types.hpp>
#include <floattetwild/geo_mesh.h>

namespace floatTetWild
{
	class MeshIO
	{
	public:
		static bool load_mesh(const std::string &path, std::vector<Vector3> &points, std::vector<Vector3i> &faces, geo::Mesh& input, std::vector<int> &flags);

		static void write_mesh(const std::string &path, const Mesh &mesh,
		        const bool only_interior = true, const std::vector<Scalar> &color = std::vector<Scalar>(), const bool binary = true, const bool separate_components = false);
	};
}
