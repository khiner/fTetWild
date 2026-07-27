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

# Sanitizers
if(FLOAT_TETWILD_WITH_SANITIZERS)
    FetchContent_Declare(
        sanitizers-cmake
        GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
        GIT_TAG        6947cff3a9c9305eb9c16135dd81da3feb4bf87f
    )
    FetchContent_MakeAvailable(sanitizers-cmake)
    list(APPEND CMAKE_MODULE_PATH ${sanitizers-cmake_SOURCE_DIR}/cmake)
    find_package(Sanitizers)
endif()

# CLI11
if(FLOAT_TETWILD_TOPLEVEL_PROJECT AND NOT TARGET CLI11::CLI11)
    FetchContent_Declare(
        cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11
        GIT_TAG v2.5.0
    )
    FetchContent_MakeAvailable(cli11)
endif()

# Eigen
# libigl used to supply Eigen3::Eigen. Its recipe fetched tag 3.4.0 and wrapped
# the headers in an INTERFACE target rather than configuring Eigen's own CMake
# project, so keep both the version and that shape. SOURCE_SUBDIR points at a
# directory with no CMakeLists.txt, which skips Eigen's project, and dropping
# the recipe's FetchContent_Populate is what stops CMake 4.x rejecting it.
if(NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
        GIT_SHALLOW    TRUE
        SOURCE_SUBDIR  cmake
    )
    FetchContent_MakeAvailable(eigen)

    add_library(Eigen3_Eigen INTERFACE)
    add_library(Eigen3::Eigen ALIAS Eigen3_Eigen)
    target_include_directories(Eigen3_Eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
endif()



# C++11 threads
find_package(Threads REQUIRED)

# Json
if(NOT TARGET json)
    FetchContent_Declare(
        json
        GIT_REPOSITORY https://github.com/jdumas/json
        GIT_TAG        0901d33bf6e7dfe6f70fd9d142c8f5c6695c6c5b
    )
    FetchContent_MakeAvailable(json)
    
    # Create interface target if not provided by the library
    if(NOT TARGET json)
        add_library(json INTERFACE)
        target_include_directories(json SYSTEM
                                 INTERFACE ${json_SOURCE_DIR}/include)
    endif()
endif()

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

if(FLOAT_TETWILD_WITH_EXACT_ENVELOPE)
    FetchContent_Declare(
        exact_envelope
        GIT_REPOSITORY https://github.com/wangbolun300/fast-envelope
        GIT_TAG        520ee04b6c69a802db31d1fd3a3e6e382d10ef98
    )
    FetchContent_MakeAvailable(exact_envelope)
endif()
