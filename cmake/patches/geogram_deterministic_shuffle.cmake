# GEO::random_shuffle seeds a fresh std::mt19937 from std::random_device on every call, so the BRIO
# insertion order Delaunay::set_vertices picks differs every run and every stage downstream of
# FloatTetDelaunay::tetrahedralize diverges with it. Seed it fixed instead. Nothing else changes, so
# the order is still drawn from the same distribution, just pinned to one draw.
#
# Run as the geogram PATCH_COMMAND with GEOGRAM_SOURCE_DIR set.

set(algorithm_h "${GEOGRAM_SOURCE_DIR}/src/lib/geogram/basic/algorithm.h")

if(NOT EXISTS "${algorithm_h}")
    message(FATAL_ERROR "geogram patch: ${algorithm_h} does not exist")
endif()

file(READ "${algorithm_h}" contents)

string(FIND "${contents}" "FTETWILD_DETERMINISTIC_SHUFFLE" marker)
if(NOT marker EQUAL -1)
    return()
endif()

string(REPLACE
    "\tstd::random_device rng;\n\tstd::mt19937 urng(rng());"
    "\tstd::mt19937 urng(42); // FTETWILD_DETERMINISTIC_SHUFFLE, see cmake/patches"
    patched
    "${contents}")

if(patched STREQUAL contents)
    message(FATAL_ERROR
        "geogram patch: could not find the random_device seeding in GEO::random_shuffle. "
        "Check whether ${algorithm_h} still needs this, and update or drop the patch.")
endif()

file(WRITE "${algorithm_h}" "${patched}")
message(STATUS "geogram: seeded GEO::random_shuffle deterministically")
