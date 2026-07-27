// Declarations for the Shewchuk exact predicates in predicates.c.
//
// REAL there is bound to fTetWild's Scalar by the same FLOAT_TETWILD_USE_FLOAT
// switch, so these signatures follow it.

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif
