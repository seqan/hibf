# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.5...3.12)

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

message (STATUS "Finding HIBF (${HIBF_VERSION}) and checking requirements")

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (CheckCXXCompilerFlag)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (hibf_config_print text)
    message (STATUS "  ${text}")
endmacro ()

macro (hibf_config_error text)
    message (FATAL_ERROR "  ${text}")

endmacro ()

# ----------------------------------------------------------------------------
# Add submodules
# ----------------------------------------------------------------------------

set (HIBF_SUBMODULES_DIR
     "${HIBF_SOURCE_DIR}"
     CACHE STRING "Directory containing submodules.")
file (GLOB submodules ${HIBF_SUBMODULES_DIR}/submodules/*/include ${HIBF_SUBMODULES_DIR}/submodules/simde/simde
      ${HIBF_SUBMODULES_DIR}/simde/simde)
foreach (submodule ${submodules})
    if (IS_DIRECTORY ${submodule})
        hibf_config_print ("  …adding submodule include: ${submodule}")
        set (HIBF_DEPENDENCY_HEADER_PATHS ${submodule} ${HIBF_DEPENDENCY_HEADER_PATHS})
    endif ()
endforeach ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET 1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} ${HIBF_HEADER_PATH} ${HIBF_DEPENDENCY_HEADER_PATHS})
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Check supported compilers
# ----------------------------------------------------------------------------

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
    hibf_config_error ("GCC < 11 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17)
    hibf_config_error ("Clang < 17 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

# ----------------------------------------------------------------------------
# ccache
# ----------------------------------------------------------------------------

include ("${HIBF_SOURCE_DIR}/test/cmake/hibf_require_ccache.cmake")
hibf_require_ccache ()

# ----------------------------------------------------------------------------
# Require C++20
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})

set (CXXSTD_TEST_SOURCE
     "#if !defined (__cplusplus) || (__cplusplus < 202002)
      #error NOCXX20
      #endif
      int main() {}")

set (HIBF_FEATURE_CPP20_FLAG_BUILTIN "")
set (HIBF_FEATURE_CPP20_FLAG_STD20 "-std=c++20")
set (HIBF_FEATURE_CPP20_FLAG_STD2a "-std=c++2a")

set (HIBF_CPP20_FLAG "")

foreach (_FLAG BUILTIN STD20 STD2a)
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_SAVE} ${HIBF_FEATURE_CPP20_FLAG_${_FLAG}}")

    check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CPP20_FLAG_${_FLAG})

    if (CPP20_FLAG_${_FLAG})
        set (HIBF_CPP20_FLAG ${_FLAG})
        break ()
    endif ()
endforeach ()

set (CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})

if (HIBF_CPP20_FLAG STREQUAL "BUILTIN")
    hibf_config_print ("C++ Standard-20 support:    builtin")
elseif (HIBF_CPP20_FLAG)
    set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} ${HIBF_FEATURE_CPP20_FLAG_${HIBF_CPP20_FLAG}}")
    hibf_config_print ("C++ Standard-20 support:    via ${HIBF_FEATURE_CPP20_FLAG_${HIBF_CPP20_FLAG}}")
else ()
    hibf_config_error ("HIBF requires C++20, but your compiler does not support it.")
endif ()

# ----------------------------------------------------------------------------
# Required: OpenMP
# ----------------------------------------------------------------------------

check_cxx_compiler_flag ("-fopenmp" HIBF_HAS_OPENMP)
if (HIBF_HAS_OPENMP)
    set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} -fopenmp")
    hibf_config_print ("OpenMP Support:             via -fopenmp")
else ()
    hibf_config_error ("HIBF requires OpenMP, but your compiler does not support it.")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" HIBF_HAS_OPENMP_SIMD)
if (HIBF_HAS_OPENMP_SIMD)
    set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} -fopenmp-simd")
    set (HIBF_DEFINITIONS ${HIBF_DEFINITIONS} "-DSIMDE_ENABLE_OPENMP")
    hibf_config_print ("SIMD-OpenMP Support:        via -fopenmp-simd")
else ()
    hibf_config_print ("SIMD-OpenMP Support:        not found")
endif ()

check_cxx_compiler_flag ("-Wno-psabi" HIBF_SUPPRESS_GCC4_ABI)
if (HIBF_SUPPRESS_GCC4_ABI)
    set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} -Wno-psabi")
    hibf_config_print ("Suppressing GCC 4 warnings: via -Wno-psabi")
endif ()

# ----------------------------------------------------------------------------
# Optimizations
# ----------------------------------------------------------------------------

if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR "${CMAKE_BUILD_TYPE}" MATCHES "Coverage")
    set (HIBF_IS_DEBUG TRUE)
else ()
    set (HIBF_IS_DEBUG FALSE)
endif ()

option (HIBF_NATIVE_BUILD "Optimize build for current architecture." ON)
if (HIBF_IS_DEBUG)
    hibf_config_print ("Optimize build:             disabled")
elseif (HIBF_NATIVE_BUILD)
    set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} -march=native")
    hibf_config_print ("Optimize build:             via -march=native")
else ()
    check_cxx_compiler_flag ("-mpopcnt" HIBF_HAS_POPCNT)
    if (HIBF_HAS_POPCNT)
        set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} -mpopcnt")
        hibf_config_print ("Optimize build:             via -mpopcnt")
    else ()
        hibf_config_print ("Optimize build:             disabled")
    endif ()
endif ()

option (HIBF_LTO_BUILD "Use Link Time Optimisation." ON)
if (HIBF_IS_DEBUG OR NOT HIBF_LTO_BUILD)
    hibf_config_print ("Link Time Optimisation:     disabled")
else ()
    # CMake's check_ipo_supported uses hardcoded lto flags
    # macOS GCC supports -flto-auto, but not the hardcoded flag "-fno-fat-lto-objects"
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")

        set (HIBF_LTO_FLAGS "-flto=auto -ffat-lto-objects")
    else ()
        set (HIBF_LTO_FLAGS "-flto=auto")
    endif ()

    set (LTO_CMAKE_SOURCE
         "cmake_minimum_required(VERSION ${CMAKE_VERSION})\nproject(lto-test LANGUAGES CXX)
          cmake_policy(SET CMP0069 NEW)\nadd_library(foo foo.cpp)\nadd_executable(boo main.cpp)
          target_link_libraries(boo PUBLIC foo)")
    set (LTO_FOO_CPP "int foo(){return 0x42;}")
    set (LTO_MAIN_CPP "int foo();int main(){return foo();}")
    set (testdir "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/lto-test")
    set (bindir "${testdir}/bin")
    set (srcdir "${testdir}/src")
    file (MAKE_DIRECTORY "${bindir}")
    file (MAKE_DIRECTORY "${srcdir}")
    file (WRITE "${srcdir}/foo.cpp" "${LTO_FOO_CPP}")
    file (WRITE "${srcdir}/main.cpp" "${LTO_MAIN_CPP}")
    file (WRITE "${srcdir}/CMakeLists.txt" "${LTO_CMAKE_SOURCE}")
    try_compile (HIBF_HAS_LTO "${bindir}"
                 "${srcdir}" "lto-test"
                 CMAKE_FLAGS "-DCMAKE_VERBOSE_MAKEFILE=ON" "-DCMAKE_CXX_FLAGS:STRING=${HIBF_LTO_FLAGS}"
                 OUTPUT_VARIABLE output)
    if (HIBF_HAS_LTO)
        hibf_config_print ("Link Time Optimisation:     enabled")
        set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} ${HIBF_LTO_FLAGS}")
    else ()
        hibf_config_print ("Link Time Optimisation:     not available")
    endif ()
endif ()

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

set (THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads QUIET)

if (Threads_FOUND)
    if ("${CMAKE_THREAD_LIBS_INIT}" STREQUAL "")
        hibf_config_print ("Thread support:             builtin")
    else ()
        set (HIBF_LIBRARIES ${HIBF_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
        hibf_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    hibf_config_print ("Thread support:             not found")
endif ()

# ----------------------------------------------------------------------------
# xxHash dependency
# ----------------------------------------------------------------------------
find_path (HIBF_XXHASH_DIR
           NAMES xxhash.h
           HINTS "${HIBF_HEADER_PATH}/hibf/contrib/xxhash")

if (HIBF_XXHASH_DIR)
    hibf_config_print ("Required dependency:        xxHash found")
    set (HIBF_DEFINITIONS ${HIBF_DEFINITIONS} "-DXXH_INLINE_ALL")
else ()
    hibf_config_error ("Required dependency:        xxHash not found")
endif ()

unset (HIBF_XXHASH_DIR)

# ----------------------------------------------------------------------------
# robin-hood dependency
# ----------------------------------------------------------------------------
find_path (HIBF_ROBIN_HOOD_DIR
           NAMES robin_hood.hpp
           HINTS "${HIBF_HEADER_PATH}/hibf/contrib")

if (HIBF_ROBIN_HOOD_DIR)
    hibf_config_print ("Required dependency:        robin-hood found")
else ()
    hibf_config_error ("Required dependency:        robin-hood not found")
endif ()

unset (HIBF_ROBIN_HOOD_DIR)

# ----------------------------------------------------------------------------
# System dependencies
# ----------------------------------------------------------------------------

# librt
if ((${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "kFreeBSD")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "GNU"))
    set (HIBF_LIBRARIES ${HIBF_LIBRARIES} rt)
endif ()

# libexecinfo -- implicit
check_include_file_cxx (execinfo.h _HIBF_HAVE_EXECINFO)
mark_as_advanced (_HIBF_HAVE_EXECINFO)
if (_HIBF_HAVE_EXECINFO)
    hibf_config_print ("Optional dependency:        libexecinfo found")
    if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") OR (${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD"))
        set (HIBF_LIBRARIES ${HIBF_LIBRARIES} execinfo elf)
    endif ()
else ()
    hibf_config_print ("Optional dependency:        libexecinfo not found")
endif ()

# ----------------------------------------------------------------------------
# Perform compilability test of platform.hpp (tests some requirements)
# ----------------------------------------------------------------------------

set (CXXSTD_TEST_SOURCE "#include <hibf/platform.hpp>
                         int main() {}")

# using try_compile instead of check_cxx_source_compiles to capture output in case of failure
file (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "${CXXSTD_TEST_SOURCE}\n")

try_compile (HIBF_PLATFORM_TEST #
             ${CMAKE_BINARY_DIR}
             ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
             CMAKE_FLAGS "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_CXX_FLAGS} ${HIBF_CXX_FLAGS}"
                         "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_INCLUDE_PATH};${HIBF_HEADER_PATH};${HIBF_DEPENDENCY_HEADER_PATHS}"
             COMPILE_DEFINITIONS ${HIBF_DEFINITIONS}
             LINK_LIBRARIES ${HIBF_LIBRARIES}
             OUTPUT_VARIABLE HIBF_PLATFORM_TEST_OUTPUT)

if (HIBF_PLATFORM_TEST)
    hibf_config_print ("HIBF platform.hpp build:    passed")
else ()
    hibf_config_error ("HIBF platform.hpp build:    failed\n${HIBF_PLATFORM_TEST_OUTPUT}")
endif ()

separate_arguments (HIBF_CXX_FLAGS UNIX_COMMAND "${HIBF_CXX_FLAGS}")
