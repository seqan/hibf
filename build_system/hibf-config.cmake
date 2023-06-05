# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------
#
# This CMake module will try to find SeqAn and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package (HIBF [REQUIRED] ...)
#
# Since this makes a difference for CMAKE, pay attention to the case
# ("HIBF", "HIBF" and "hibf" are all valid, but other names not).
#
# SeqAn has the following platform requirements:
#
#   C++20
#   pthread
#
# SeqAn requires the following libraries:
#
#   SDSL      -- the succinct data structure library
#
# SeqAn has the following optional dependencies:
#
#   ZLIB      -- zlib compression library
#   BZip2     -- libbz2 compression library
#   Cereal    -- Serialisation library
#
# If you don't wish for these to be detected (and used), you may define HIBF_NO_ZLIB,
# HIBF_NO_BZIP2, and HIBF_NO_CEREAL respectively.
#
# If you wish to require the presence of ZLIB or BZip2, just check for the module before
# finding HIBF, e.g. "find_package (ZLIB REQUIRED)" and "find_package (BZip2 REQUIRED)".
# If you wish to require the presence of CEREAL, you may define HIBF_CEREAL.
#
# Once the search has been performed, the following variables will be set.
#
#   HIBF_FOUND            -- Indicate whether SeqAn was found and requirements met.
#
#   HIBF_VERSION          -- The version as string, e.g. "3.0.0"
#   HIBF_VERSION_MAJOR    -- e.g. 3
#   HIBF_VERSION_MINOR    -- e.g. 0
#   HIBF_VERSION_PATCH    -- e.g. 0
#
#   HIBF_INCLUDE_DIRS     -- to be passed to include_directories ()
#   HIBF_LIBRARIES        -- to be passed to target_link_libraries ()
#   HIBF_DEFINITIONS      -- to be passed to add_definitions ()
#   HIBF_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   hibf::hibf          -- interface target where
#                                  target_link_libraries(target hibf::hibf)
#                              automatically sets
#                                  target_include_directories(target $HIBF_INCLUDE_DIRS),
#                                  target_link_libraries(target $HIBF_LIBRARIES),
#                                  target_compile_definitions(target $HIBF_DEFINITIONS) and
#                                  target_compile_options(target $HIBF_CXX_FLAGS)
#                              for a target.
#
#   [IMPORTED]: https://cmake.org/cmake/help/v3.10/prop_tgt/IMPORTED.html#prop_tgt:IMPORTED
#
# ============================================================================

cmake_minimum_required (VERSION 3.4...3.12)

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

# make output globally quiet if required by find_package, this effects cmake functions like `check_*`
set (CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set (CMAKE_REQUIRED_QUIET ${${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY})

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Finding HIBF and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (CheckCXXCompilerFlag)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (hibf_config_print text)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (hibf_config_error text)
    if (${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
        message (FATAL_ERROR ${text})
    else ()
        if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message (WARNING ${text})
        endif ()
        return ()
    endif ()
endmacro ()

# ----------------------------------------------------------------------------
# Find HIBF include path
# ----------------------------------------------------------------------------

# Note that hibf-config.cmake can be standalone and thus HIBF_CLONE_DIR might be empty.
# * `HIBF_CLONE_DIR` was already found in hibf-config-version.cmake
# * `HIBF_INCLUDE_DIR` was already found in hibf-config-version.cmake
find_path (HIBF_SUBMODULES_DIR
           NAMES submodules
           HINTS "${HIBF_CLONE_DIR}" "${HIBF_INCLUDE_DIR}/hibf")

if (HIBF_INCLUDE_DIR)
    hibf_config_print ("HIBF include dir found:     ${HIBF_INCLUDE_DIR}")
else ()
    hibf_config_error ("HIBF include directory could not be found (HIBF_INCLUDE_DIR: '${HIBF_INCLUDE_DIR}')")
endif ()

find_path (HIBF_SOURCE_DIR
           NAMES hierarchical_interleaved_bloom_filter.cpp
           HINTS "${HIBF_CLONE_DIR}/src")

if (HIBF_SOURCE_DIR)
    hibf_config_print ("HIBF source dir found:      ${HIBF_SOURCE_DIR}")
else ()
    hibf_config_error ("HIBF source directory could not be found (HIBF_SOURCE_DIR: '${HIBF_SOURCE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Detect if we are a clone of repository and if yes auto-add submodules
# ----------------------------------------------------------------------------

if (HIBF_CLONE_DIR)
    hibf_config_print ("Detected as running from a repository checkout…")
endif ()

if (HIBF_SUBMODULES_DIR)
    file (GLOB submodules ${HIBF_SUBMODULES_DIR}/submodules/*/include ${HIBF_SUBMODULES_DIR}/submodules/simde/simde)
    foreach (submodule ${submodules})
        if (IS_DIRECTORY ${submodule})
            hibf_config_print ("  …adding submodule include: ${submodule}")
            set (HIBF_DEPENDENCY_INCLUDE_DIRS ${submodule} ${HIBF_DEPENDENCY_INCLUDE_DIRS})
        endif ()
    endforeach ()
endif ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET 1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} ${HIBF_INCLUDE_DIR} ${HIBF_DEPENDENCY_INCLUDE_DIRS})
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Check supported compilers
# ----------------------------------------------------------------------------

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message (FATAL_ERROR "GCC < 10 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

option (HIBF_DISABLE_COMPILER_CHECK "Skips the check for supported compilers." OFF)

if (NOT HIBF_DISABLE_COMPILER_CHECK)
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message (FATAL_ERROR "Only GCC is supported. "
                             "The detected compiler version is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. "
                             "You can disable this error by passing -DHIBF_DISABLE_COMPILER_CHECK=ON to CMake.")
    endif ()
else ()
    set (HIBF_DEFINITIONS ${HIBF_DEFINITIONS} "-DHIBF_DISABLE_COMPILER_CHECK")
endif ()

# ----------------------------------------------------------------------------
# Require C++20
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})

set (CXXSTD_TEST_SOURCE
     "#if !defined (__cplusplus) || (__cplusplus < 201709)
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
    set (HIBF_CXX_FLAGS "${HIBF_CXX_FLAGS} -fopenmp-simd -DSIMDE_ENABLE_OPENMP")
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

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

set (THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads QUIET)

if (Threads_FOUND)
    set (HIBF_LIBRARIES ${HIBF_LIBRARIES} Threads::Threads)
    if ("${CMAKE_THREAD_LIBS_INIT}" STREQUAL "")
        hibf_config_print ("Thread support:             builtin.")
    else ()
        hibf_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    hibf_config_print ("Thread support:             not found.")
endif ()

# ----------------------------------------------------------------------------
# xxHash dependency
# ----------------------------------------------------------------------------
find_path (HIBF_XXHASH_DIR
           NAMES xxhash.h
           HINTS "${HIBF_INCLUDE_DIR}/hibf/contrib/xxhash")

if (HIBF_XXHASH_DIR)
    hibf_config_print ("Required dependency:        xxHash found.")
    set (HIBF_DEFINITIONS ${HIBF_DEFINITIONS} "-DXXH_INLINE_ALL")
else ()
    hibf_config_error ("Required dependency xxHash not found.")
endif ()

unset (HIBF_XXHASH_DIR)

# ----------------------------------------------------------------------------
# robin-hood dependency
# ----------------------------------------------------------------------------
find_path (HIBF_ROBIN_HOOD_DIR
           NAMES robin_hood.hpp
           HINTS "${HIBF_INCLUDE_DIR}/hibf/contrib")

if (HIBF_ROBIN_HOOD_DIR)
    hibf_config_print ("Required dependency:        robin-hood found.")
else ()
    hibf_config_error ("Required dependency robin-hood not found.")
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
    hibf_config_print ("Optional dependency:        libexecinfo found.")
    if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") OR (${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD"))
        set (HIBF_LIBRARIES ${HIBF_LIBRARIES} execinfo elf)
    endif ()
else ()
    hibf_config_print ("Optional dependency:        libexecinfo not found.")
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
                         "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_INCLUDE_PATH};${HIBF_INCLUDE_DIR};${HIBF_DEPENDENCY_INCLUDE_DIRS}"
             COMPILE_DEFINITIONS ${HIBF_DEFINITIONS}
             LINK_LIBRARIES ${HIBF_LIBRARIES}
             OUTPUT_VARIABLE HIBF_PLATFORM_TEST_OUTPUT)

if (HIBF_PLATFORM_TEST)
    hibf_config_print ("HIBF platform.hpp build:    passed.")
else ()
    hibf_config_error ("HIBF platform.hpp build:    failed!\n\
                        ${HIBF_PLATFORM_TEST_OUTPUT}")
endif ()

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS HIBF_INCLUDE_DIR)

# Set HIBF_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
# This needs to be done, because `find_package(HIBF)` might be called in any case-sensitive way and we want to
# guarantee that HIBF_* are always set.
foreach (package_var
         FOUND
         DIR
         ROOT
         CONFIG
         VERSION
         VERSION_MAJOR
         VERSION_MINOR
         VERSION_PATCH
         VERSION_TWEAK
         VERSION_COUNT)
    set (HIBF_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
endforeach ()

# propagate HIBF_INCLUDE_DIR into HIBF_INCLUDE_DIRS
set (HIBF_INCLUDE_DIRS ${HIBF_INCLUDE_DIR} ${HIBF_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

if (HIBF_FOUND AND NOT TARGET hibf::hibf)
    separate_arguments (HIBF_CXX_FLAGS_LIST UNIX_COMMAND "${HIBF_CXX_FLAGS}")

    add_library (hibf_hibf STATIC
                 ${HIBF_SOURCE_DIR}/hierarchical_interleaved_bloom_filter.cpp
                 ${HIBF_SOURCE_DIR}/detail/layout/simple_binning.cpp
                 ${HIBF_SOURCE_DIR}/detail/layout/execute.cpp
                 ${HIBF_SOURCE_DIR}/detail/layout/output.cpp
                 ${HIBF_SOURCE_DIR}/detail/layout/compute_fp_correction.cpp
                 ${HIBF_SOURCE_DIR}/detail/layout/hierarchical_binning.cpp
                 ${HIBF_SOURCE_DIR}/detail/sketch/toolbox.cpp
                 ${HIBF_SOURCE_DIR}/detail/sketch/hyperloglog.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/initialise_build_tree.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/insert_into_ibf.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/parse_chopper_pack_header.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/compute_kmers.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/read_chopper_pack_file.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/update_header_node_data.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/parse_chopper_pack_line.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/construct_ibf.cpp
                 ${HIBF_SOURCE_DIR}/detail/build/update_content_node_data.cpp)
    target_compile_definitions (hibf_hibf PUBLIC ${HIBF_DEFINITIONS})
    target_compile_options (hibf_hibf PUBLIC ${HIBF_CXX_FLAGS_LIST})
    target_link_options (hibf_hibf PUBLIC ${HIBF_CXX_FLAGS_LIST})
    target_link_libraries (hibf_hibf PUBLIC "${HIBF_LIBRARIES}")
    # include hibf/include/ as -I, because hibf should never produce warnings.
    target_include_directories (hibf_hibf PUBLIC "${HIBF_INCLUDE_DIR}")
    # include everything except hibf/include/ as -isystem, i.e.
    # a system header which suppresses warnings of external libraries.
    target_include_directories (hibf_hibf SYSTEM PUBLIC "${HIBF_DEPENDENCY_INCLUDE_DIRS}")
    add_library (hibf::hibf ALIAS hibf_hibf)
endif ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (HIBF_FIND_DEBUG)
    message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
    message ("")
    message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
    message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
    message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
    message ("  HIBF_INCLUDE_DIR          ${HIBF_INCLUDE_DIR}")
    message ("")
    message ("  ${CMAKE_FIND_PACKAGE_NAME}_FOUND                ${${CMAKE_FIND_PACKAGE_NAME}_FOUND}")
    message ("  HIBF_HAS_ZLIB             ${ZLIB_FOUND}")
    message ("  HIBF_HAS_BZIP2            ${BZIP2_FOUND}")
    message ("")
    message ("  HIBF_INCLUDE_DIRS         ${HIBF_INCLUDE_DIRS}")
    message ("  HIBF_LIBRARIES            ${HIBF_LIBRARIES}")
    message ("  HIBF_DEFINITIONS          ${HIBF_DEFINITIONS}")
    message ("  HIBF_CXX_FLAGS            ${HIBF_CXX_FLAGS}")
    message ("")
    message ("  HIBF_VERSION              ${HIBF_VERSION}")
    message ("  HIBF_VERSION_MAJOR        ${HIBF_VERSION_MAJOR}")
    message ("  HIBF_VERSION_MINOR        ${HIBF_VERSION_MINOR}")
    message ("  HIBF_VERSION_PATCH        ${HIBF_VERSION_PATCH}")
endif ()
