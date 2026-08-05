# ##############################################################################
cmake_minimum_required(VERSION 3.20)
project(FloatTetwild)
# ##############################################################################

# Release unless asked for something else. With no build type CMake passes no optimisation flags
# at all, which is not just slower: the arithmetic differs, so script/check_corpus fails against a
# baseline recorded from an optimised build. An explicit -DCMAKE_BUILD_TYPE still wins, and
# multi-config generators pick their configuration at build time so they are left alone.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    # FORCE because project() has already created the cache entry, empty, and a plain set()
    # against an existing entry does nothing. The guard above is what keeps this from stomping
    # a build type the user asked for.
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Extra warnings
include(Warnings)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ##############################################################################
# Dependencies
# ##############################################################################

# The in-tree thread pool.
find_package(Threads REQUIRED)

# Argument parsing in main.cpp, so only the binary needs it.
include(FetchContent)
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(cli11)

# ##############################################################################
# FloatTetwild library
# ##############################################################################

# Sources are attached by src/CMakeLists.txt via target_sources().
add_library(${PROJECT_NAME})

# Headers are included relative to src/ (e.g. "geo/geo_mesh.h"), from the tests too.
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# Use C++14 (minimum for modern CMake practices)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

target_link_libraries(${PROJECT_NAME} PUBLIC Threads::Threads)

# ##############################################################################
# FloatTetwild binary
# ##############################################################################

add_executable(${PROJECT_NAME}_bin src/main.cpp)
target_compile_features(${PROJECT_NAME}_bin PUBLIC cxx_std_14)

target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME} CLI11::CLI11
                                               warnings::all)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
             PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)

# ##############################################################################
# Subdirectories
# ##############################################################################

# Other sources
add_subdirectory(src)

# Unit tests. include(CTest) already calls enable_testing().
include(CTest)
add_subdirectory(tests)
