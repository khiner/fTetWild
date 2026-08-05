// The fast winding number for triangle soups.
//
// Copyright (c) 2018 Side Effects Software Inc. MIT licence, reproduced in FastWindingNumber.cpp
// and in THIRD_PARTY.md.

#ifndef FLOATTETWILD_FASTWINDINGNUMBER_H
#define FLOATTETWILD_FASTWINDINGNUMBER_H

#include <vector>

namespace floatTetWild {
namespace FastWindingNumber {

// The signed solid angle each query point subtends against the triangle soup, approximated by
// "Fast Winding Numbers for Soups and Clouds" [Barill et al. 2018]: a four-wide BVH over the
// triangles, each node carrying a second order Taylor expansion of the angle its triangles
// subtend, descended only where the expansion is not accurate enough. Divide by 4 pi to get a
// winding number. points and queries hold three floats per point; triangles holds three indices
// into points per triangle. Returns one angle per query.
std::vector<float> solid_angles(const std::vector<float>& points,
                                const std::vector<int>&   triangles,
                                const std::vector<float>& queries);

}  // namespace FastWindingNumber
}  // namespace floatTetWild

#endif  // FLOATTETWILD_FASTWINDINGNUMBER_H
