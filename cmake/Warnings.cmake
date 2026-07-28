# ##############################################################################
# An INTERFACE target carrying the project's warning flags. Every flag is probed with
# check_cxx_compiler_flag first, so compiler-specific ones simply drop out where they are not
# understood.
# ##############################################################################

if(TARGET warnings::all)
    return()
endif()

set(MY_FLAGS
        -Wall
        -Wextra
        -pedantic
        -Wunused

        -Wno-long-long
        -Wpointer-arith
        -Wformat=2
        -Wuninitialized
        -Wcast-qual
        -Wmissing-noreturn
        -Wmissing-format-attribute
        -Wredundant-decls

        -Wno-unused-variable
        -Wunused-but-set-variable
        -Wno-unused-parameter
        -Wno-old-style-cast

        -Wshadow

        -Wstrict-null-sentinel
        -Woverloaded-virtual
        -Wsign-promo
        -Wstack-protector
        -Wstrict-aliasing
        -Wstrict-aliasing=2
        -Wswitch-default
        -Wswitch-enum
        -Wswitch-unreachable

        -Wcast-align
        -Wdisabled-optimization
        -Winvalid-pch
        -Wpacked
        -Wno-padded
        -Wstrict-overflow
        -Wstrict-overflow=2

        -Wctor-dtor-privacy
        -Wlogical-op
        -Wnoexcept

        -Wnon-virtual-dtor
        -Wdelete-non-virtual-dtor

        -Wno-sign-compare

        -Wno-error=stringop-overflow
        -Wnull-dereference
        -fdelete-null-pointer-checks
        -Wduplicated-cond
        -Wmisleading-indentation

        # Not warnings: these keep stack traces meaningful.
        -fno-omit-frame-pointer
        -fno-optimize-sibling-calls
)

# Flags above don't make sense for MSVC
if(MSVC)
    set(MY_FLAGS)
endif()

include(CheckCXXCompilerFlag)

add_library(warnings_all INTERFACE)
add_library(warnings::all ALIAS warnings_all)

foreach(FLAG IN ITEMS ${MY_FLAGS})
    string(REPLACE "=" "-" FLAG_VAR "${FLAG}")
    if(NOT DEFINED IS_SUPPORTED_${FLAG_VAR})
        check_cxx_compiler_flag("${FLAG}" IS_SUPPORTED_${FLAG_VAR})
    endif()
    if(IS_SUPPORTED_${FLAG_VAR})
        target_compile_options(warnings_all INTERFACE ${FLAG})
    endif()
endforeach()
