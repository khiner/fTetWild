# Download Catch2 unit test framework using FetchContent
include(FetchContent)
FetchContent_Declare(
    Catch2
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.5.3.tar.gz
    URL_HASH MD5=1f51d817ce81d54b12e87d06e159305f
)
FetchContent_MakeAvailable(Catch2)

add_executable(unit_tests test_predicates.cpp test_table.cpp)

target_link_libraries(unit_tests PUBLIC FloatTetwild Catch2::Catch2WithMain warnings::all)

# Register tests
list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(unit_tests)
