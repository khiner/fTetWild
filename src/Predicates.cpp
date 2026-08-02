#include <floattetwild/Predicates.hpp>

// The Shewchuk exact predicates in predicates.c, which defines REAL off the same
// FLOAT_TETWILD_USE_FLOAT switch Scalar follows and does not include these declarations itself.
extern "C" {
#ifdef FLOAT_TETWILD_USE_FLOAT
typedef float PredicatesReal;
#else
typedef double PredicatesReal;
#endif

// Must be called once before any of the predicates below.
void exactinit(void);

PredicatesReal orient2d(PredicatesReal* pa, PredicatesReal* pb, PredicatesReal* pc);
PredicatesReal orient3d(PredicatesReal* pa,
                        PredicatesReal* pb,
                        PredicatesReal* pc,
                        PredicatesReal* pd);
}

namespace floatTetWild {

    namespace {
        // Shewchuk's error bounds have to be computed once before any predicate runs.
        // libigl did this behind a Meyers singleton on every call, so do the same.
        inline void init_predicates()
        {
            static const bool initialized = []() { exactinit(); return true; }();
            (void) initialized;
        }

        // The ORI_* constants are 1/0/-1, so the sign is the orientation. orient_3d_volume
        // returns it as a Scalar when it is not positive, so the normalisation has to stay.
        inline int sign_of(Scalar r)
        {
            if (r > 0)
                return Predicates::ORI_POSITIVE;
            else if (r < 0)
                return Predicates::ORI_NEGATIVE;
            else
                return Predicates::ORI_ZERO;
        }

    }  // namespace

    const int Predicates::ORI_POSITIVE;
    const int Predicates::ORI_ZERO;
    const int Predicates::ORI_NEGATIVE;
    const int Predicates::ORI_UNKNOWN;

    int Predicates::orient_3d(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) {
        init_predicates();
        return sign_of(orient3d(const_cast<Scalar*>(p1.data()),
                                const_cast<Scalar*>(p2.data()),
                                const_cast<Scalar*>(p3.data()),
                                const_cast<Scalar*>(p4.data())));
    }

    Scalar Predicates::orient_3d_volume(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) {
        const int ori = orient_3d(p1, p2, p3, p4);

        if (ori <= 0)
            return ori;
        else
            return (p1 - p4).dot((p2 - p4).cross(p3 - p4)) / 6;
    }

    int Predicates::orient_2d(const Vector2& p1, const Vector2& p2, const Vector2& p3) {
        init_predicates();
        return sign_of(orient2d(const_cast<Scalar*>(p1.data()),
                                const_cast<Scalar*>(p2.data()),
                                const_cast<Scalar*>(p3.data())));
    }

} // namespace floatTetWild
