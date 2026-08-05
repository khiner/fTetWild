#include "Predicates.hpp"

#include "geo/geo_predicates.h"

// The Shewchuk exact orient2d in predicates.c, which does not include these declarations itself.
extern "C" {
// Must be called once before orient2d.
void exactinit(void);

double orient2d(const double* pa, const double* pb, const double* pc);
}

namespace floatTetWild {
namespace Predicates {

int orient_3d(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4)
{
    // geogram's orient_3d takes the determinant of [p2-p1, p3-p1, p4-p1] where Shewchuk took
    // [p1-p4, p2-p4, p3-p4], so the two exact signs are each other's negation.
    return -geo::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
}

int orient_2d(const Vector2& p1, const Vector2& p2, const Vector2& p3)
{
    // Shewchuk's error bounds have to be computed once before any predicate runs.
    // libigl did this behind a Meyers singleton on every call, so do the same.
    static const bool initialized = []() { exactinit(); return true; }();
    (void) initialized;
    // The ORI_* constants are 1/0/-1, so the sign is the orientation.
    const Scalar r = orient2d(p1.data(), p2.data(), p3.data());
    return r > 0 ? ORI_POSITIVE : r < 0 ? ORI_NEGATIVE : ORI_ZERO;
}

}  // namespace Predicates
}  // namespace floatTetWild
