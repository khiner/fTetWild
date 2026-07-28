# ##############################################################################
# Prepare dependencies
# ##############################################################################

# Use modern FetchContent for dependency management
include(FetchContent)

# Set a global preference for STATIC libraries over SHARED ones.
#set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)

# Set FetchContent properties for better performance
set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

# ##############################################################################
# Required libraries
# ##############################################################################

# CLI11
if(FLOAT_TETWILD_TOPLEVEL_PROJECT AND NOT TARGET CLI11::CLI11)
    FetchContent_Declare(
        cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11
        GIT_TAG v2.5.0
    )
    FetchContent_MakeAvailable(cli11)
endif()

# C++11 threads
find_package(Threads REQUIRED)

# winding number float_tetwild_download_windingnumber()
# set(windingnumber_SOURCES ${FLOAT_TETWILD_EXTERNAL}/windingnumber/SYS_Math.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/SYS_Types.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_Array.cpp
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_Array.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_ArrayImpl.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_BVH.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_BVHImpl.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_FixedVector.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_ParallelUtil.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_SmallArray.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_SolidAngle.cpp
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/UT_SolidAngle.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/VM_SIMD.h
# ${FLOAT_TETWILD_EXTERNAL}/windingnumber/VM_SSEFunc.h )
#
# add_library(fast_winding_number ${windingnumber_SOURCES})
# target_link_libraries(fast_winding_number PRIVATE tbb::tbb)
# target_compile_features(fast_winding_number PRIVATE ${CXX17_FEATURES})
# target_include_directories(fast_winding_number PUBLIC
# "${FLOAT_TETWILD_EXTERNAL}/")

