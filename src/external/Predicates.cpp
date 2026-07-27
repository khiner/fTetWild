#include <floattetwild/Predicates.hpp>

#include <floattetwild/predicates.h>

#include <geogram/delaunay/delaunay_3d.h>
namespace floatTetWild {
#define GEO_PREDICATES false

    namespace {
        // Shewchuk's error bounds have to be computed once before any predicate runs.
        // libigl did this behind a Meyers singleton on every call, so do the same.
        inline void init_predicates()
        {
            static const bool initialized = []() { exactinit(); return true; }();
            (void) initialized;
        }

        // libigl mapped the predicate's sign onto an enum and the callers below mapped
        // that back onto 1/0/-1. orient_3d_volume returns this normalised value
        // directly when it is not positive, so the normalisation has to stay.
        inline Scalar predicate_sign(Scalar r)
        {
            if (r > 0)
                return 1;
            else if (r < 0)
                return -1;
            else
                return 0;
        }
    }  // namespace

    const int Predicates::ORI_POSITIVE;
    const int Predicates::ORI_ZERO;
    const int Predicates::ORI_NEGATIVE;
    const int Predicates::ORI_UNKNOWN;

    int Predicates::orient_3d(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) {
#if GEO_PREDICATES
        const int result = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
        init_predicates();
        const Scalar result = predicate_sign(
          orient3d(const_cast<Scalar*>(p1.data()),
                   const_cast<Scalar*>(p2.data()),
                   const_cast<Scalar*>(p3.data()),
                   const_cast<Scalar*>(p4.data())));
#endif

        if (result > 0)
            return ORI_POSITIVE;
        else if (result < 0)
            return ORI_NEGATIVE;
        else
            return ORI_ZERO;
    }

    int Predicates::orient_3d_tolerance(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p) {
#if GEO_PREDICATES
        const int result = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p.data());
#else
        init_predicates();
        const Scalar result = predicate_sign(
          orient3d(const_cast<Scalar*>(p1.data()),
                   const_cast<Scalar*>(p2.data()),
                   const_cast<Scalar*>(p3.data()),
                   const_cast<Scalar*>(p.data())));
#endif

        if (result == 0)
            return ORI_ZERO;

        Vector3 n = ((p2 - p3).cross(p1 - p3)).normalized();
        Scalar d = std::abs(n.dot(p - p1));
        if (d <= SCALAR_ZERO)
            return Predicates::ORI_ZERO;

        if (result > 0)
            return ORI_POSITIVE;
        else
            return ORI_NEGATIVE;
    }

    Scalar Predicates::orient_3d_volume(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) {
#if GEO_PREDICATES
        const int ori = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
        init_predicates();
        const Scalar ori = predicate_sign(
          orient3d(const_cast<Scalar*>(p1.data()),
                   const_cast<Scalar*>(p2.data()),
                   const_cast<Scalar*>(p3.data()),
                   const_cast<Scalar*>(p4.data())));
#endif
        if (ori <= 0)
            return ori;
        else
            return (p1 - p4).dot((p2 - p4).cross(p3 - p4)) / 6;
    }

    int Predicates::orient_2d(const Vector2& p1, const Vector2& p2, const Vector2& p3) {
#if GEO_PREDICATES
        const int result = -GEO::PCK::orient_2d(p1.data(), p2.data(), p3.data());
#else
        init_predicates();
        const Scalar result = predicate_sign(
          orient2d(const_cast<Scalar*>(p1.data()),
                   const_cast<Scalar*>(p2.data()),
                   const_cast<Scalar*>(p3.data())));
#endif
        if (result > 0)
            return ORI_POSITIVE;
        else if (result < 0)
            return ORI_NEGATIVE;
        else
            return ORI_ZERO;
    }

} // namespace floatTetWild
