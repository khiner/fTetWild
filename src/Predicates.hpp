#pragma once

#include <floattetwild/Types.hpp>

namespace floatTetWild {
namespace Predicates {

constexpr int ORI_POSITIVE = 1;
constexpr int ORI_ZERO     = 0;
constexpr int ORI_NEGATIVE = -1;
constexpr int ORI_UNKNOWN  = INT_MAX;

int    orient_3d(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4);
Scalar orient_3d_volume(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4);

int orient_2d(const Vector2& p1, const Vector2& p2, const Vector2& p3);

}  // namespace Predicates
}  // namespace floatTetWild
