# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library_template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------
#
# This CMake module will try to find SeqAn and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package (Library_Template [REQUIRED] ...)
#
# Since this makes a difference for CMAKE, pay attention to the case
# ("Library_Template", "LIBRARY_TEMPLATE" and "library_template" are all valid, but other names not).
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
# If you don't wish for these to be detected (and used), you may define LIBRARY_TEMPLATE_NO_ZLIB,
# LIBRARY_TEMPLATE_NO_BZIP2, and LIBRARY_TEMPLATE_NO_CEREAL respectively.
#
# If you wish to require the presence of ZLIB or BZip2, just check for the module before
# finding Library_Template, e.g. "find_package (ZLIB REQUIRED)" and "find_package (BZip2 REQUIRED)".
# If you wish to require the presence of CEREAL, you may define LIBRARY_TEMPLATE_CEREAL.
#
# Once the search has been performed, the following variables will be set.
#
#   LIBRARY_TEMPLATE_FOUND            -- Indicate whether SeqAn was found and requirements met.
#
#   LIBRARY_TEMPLATE_VERSION          -- The version as string, e.g. "3.0.0"
#   LIBRARY_TEMPLATE_VERSION_MAJOR    -- e.g. 3
#   LIBRARY_TEMPLATE_VERSION_MINOR    -- e.g. 0
#   LIBRARY_TEMPLATE_VERSION_PATCH    -- e.g. 0
#
#   LIBRARY_TEMPLATE_INCLUDE_DIRS     -- to be passed to include_directories ()
#   LIBRARY_TEMPLATE_LIBRARIES        -- to be passed to target_link_libraries ()
#   LIBRARY_TEMPLATE_DEFINITIONS      -- to be passed to add_definitions ()
#   LIBRARY_TEMPLATE_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   library_template::library_template          -- interface target where
#                                  target_link_libraries(target library_template::library_template)
#                              automatically sets
#                                  target_include_directories(target $LIBRARY_TEMPLATE_INCLUDE_DIRS),
#                                  target_link_libraries(target $LIBRARY_TEMPLATE_LIBRARIES),
#                                  target_compile_definitions(target $LIBRARY_TEMPLATE_DEFINITIONS) and
#                                  target_compile_options(target $LIBRARY_TEMPLATE_CXX_FLAGS)
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
    message (STATUS "${ColourBold}Finding Library_Template and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (library_template_config_print text)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (library_template_config_error text)
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
# Find Library_Template include path
# ----------------------------------------------------------------------------

# Note that library_template-config.cmake can be standalone and thus LIBRARY_TEMPLATE_CLONE_DIR might be empty.
# * `LIBRARY_TEMPLATE_CLONE_DIR` was already found in library_template-config-version.cmake
# * `LIBRARY_TEMPLATE_INCLUDE_DIR` was already found in library_template-config-version.cmake
find_path (LIBRARY_TEMPLATE_SUBMODULES_DIR
           NAMES submodules
           HINTS "${LIBRARY_TEMPLATE_CLONE_DIR}" "${LIBRARY_TEMPLATE_INCLUDE_DIR}/library_template")

if (LIBRARY_TEMPLATE_INCLUDE_DIR)
    library_template_config_print ("Library_Template include dir found:   ${LIBRARY_TEMPLATE_INCLUDE_DIR}")
else ()
    library_template_config_error (
        "Library_Template include directory could not be found (LIBRARY_TEMPLATE_INCLUDE_DIR: '${LIBRARY_TEMPLATE_INCLUDE_DIR}')"
    )
endif ()

# ----------------------------------------------------------------------------
# Detect if we are a clone of repository and if yes auto-add submodules
# ----------------------------------------------------------------------------

if (LIBRARY_TEMPLATE_CLONE_DIR)
    library_template_config_print ("Detected as running from a repository checkout…")
endif ()

if (LIBRARY_TEMPLATE_SUBMODULES_DIR)
    file (GLOB submodules ${LIBRARY_TEMPLATE_SUBMODULES_DIR}/submodules/*/include)
    foreach (submodule ${submodules})
        if (IS_DIRECTORY ${submodule})
            library_template_config_print ("  …adding submodule include:  ${submodule}")
            set (LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS ${submodule} ${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS})
        endif ()
    endforeach ()
endif ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET 1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} ${LIBRARY_TEMPLATE_INCLUDE_DIR}
                             ${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS})
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Force-deactivate optional dependencies
# ----------------------------------------------------------------------------

# These two are "opt-in", because detected by CMake
# If you want to force-require these, just do find_package (zlib REQUIRED) before find_package (library_template)
option (LIBRARY_TEMPLATE_NO_ZLIB "Don't use ZLIB, even if present." OFF)
option (LIBRARY_TEMPLATE_NO_BZIP2 "Don't use BZip2, even if present." OFF)

# ----------------------------------------------------------------------------
# Check supported compilers
# ----------------------------------------------------------------------------

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message (FATAL_ERROR "GCC < 10 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

option (LIBRARY_TEMPLATE_DISABLE_COMPILER_CHECK "Skips the check for supported compilers." OFF)

if (NOT LIBRARY_TEMPLATE_DISABLE_COMPILER_CHECK)
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message (FATAL_ERROR "Only GCC is supported. "
                             "The detected compiler version is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. "
                             "You can disable this error by passing -DLIBRARY_TEMPLATE_DISABLE_COMPILER_CHECK=ON to CMake."
        )
    endif ()
else ()
    set (LIBRARY_TEMPLATE_DEFINITIONS ${LIBRARY_TEMPLATE_DEFINITIONS} "-DLIBRARY_TEMPLATE_DISABLE_COMPILER_CHECK")
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

set (LIBRARY_TEMPLATE_FEATURE_CPP20_FLAG_BUILTIN "")
set (LIBRARY_TEMPLATE_FEATURE_CPP20_FLAG_STD20 "-std=c++20")
set (LIBRARY_TEMPLATE_FEATURE_CPP20_FLAG_STD2a "-std=c++2a")

set (LIBRARY_TEMPLATE_CPP20_FLAG "")

foreach (_FLAG BUILTIN STD20 STD2a)
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_SAVE} ${LIBRARY_TEMPLATE_FEATURE_CPP20_FLAG_${_FLAG}}")

    check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CPP20_FLAG_${_FLAG})

    if (CPP20_FLAG_${_FLAG})
        set (LIBRARY_TEMPLATE_CPP20_FLAG ${_FLAG})
        break ()
    endif ()
endforeach ()

set (CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})

if (LIBRARY_TEMPLATE_CPP20_FLAG STREQUAL "BUILTIN")
    library_template_config_print ("C++ Standard-20 support:    builtin")
elseif (LIBRARY_TEMPLATE_CPP20_FLAG)
    set (LIBRARY_TEMPLATE_CXX_FLAGS
         "${LIBRARY_TEMPLATE_CXX_FLAGS} ${LIBRARY_TEMPLATE_FEATURE_CPP20_FLAG_${LIBRARY_TEMPLATE_CPP20_FLAG}}")
    library_template_config_print (
        "C++ Standard-20 support:    via ${LIBRARY_TEMPLATE_FEATURE_CPP20_FLAG_${LIBRARY_TEMPLATE_CPP20_FLAG}}")
else ()
    library_template_config_error ("Library_Template requires C++20, but your compiler does not support it.")
endif ()

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

set (THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads QUIET)

if (Threads_FOUND)
    set (LIBRARY_TEMPLATE_LIBRARIES ${LIBRARY_TEMPLATE_LIBRARIES} Threads::Threads)
    if ("${CMAKE_THREAD_LIBS_INIT}" STREQUAL "")
        library_template_config_print ("Thread support:             builtin.")
    else ()
        library_template_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    library_template_config_print ("Thread support:             not found.")
endif ()

# ----------------------------------------------------------------------------
# ZLIB dependency
# ----------------------------------------------------------------------------

if (NOT LIBRARY_TEMPLATE_NO_ZLIB)
    find_package (ZLIB QUIET)
endif ()

if (ZLIB_FOUND)
    set (LIBRARY_TEMPLATE_LIBRARIES ${LIBRARY_TEMPLATE_LIBRARIES} ${ZLIB_LIBRARIES})
    set (LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS ${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
    set (LIBRARY_TEMPLATE_DEFINITIONS ${LIBRARY_TEMPLATE_DEFINITIONS} "-DLIBRARY_TEMPLATE_HAS_ZLIB=1")
    library_template_config_print ("Optional dependency:        ZLIB-${ZLIB_VERSION_STRING} found.")
else ()
    library_template_config_print ("Optional dependency:        ZLIB not found.")
endif ()

# ----------------------------------------------------------------------------
# BZip2 dependency
# ----------------------------------------------------------------------------

if (NOT LIBRARY_TEMPLATE_NO_BZIP2)
    find_package (BZip2 QUIET)
endif ()

if (NOT ZLIB_FOUND AND BZIP2_FOUND)
    # NOTE (marehr): iostream_bzip2 uses the type `uInt`, which is defined by
    # `zlib`. Therefore, `bzip2` will cause a ton of errors without `zlib`.
    message (AUTHOR_WARNING "Disabling BZip2 [which was successfully found], "
                            "because ZLIB was not found. BZip2 depends on ZLIB.")
    unset (BZIP2_FOUND)
endif ()

if (BZIP2_FOUND)
    set (LIBRARY_TEMPLATE_LIBRARIES ${LIBRARY_TEMPLATE_LIBRARIES} ${BZIP2_LIBRARIES})
    set (LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS ${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS} ${BZIP2_INCLUDE_DIRS})
    set (LIBRARY_TEMPLATE_DEFINITIONS ${LIBRARY_TEMPLATE_DEFINITIONS} "-DLIBRARY_TEMPLATE_HAS_BZIP2=1")
    library_template_config_print ("Optional dependency:        BZip2-${BZIP2_VERSION_STRING} found.")
else ()
    library_template_config_print ("Optional dependency:        BZip2 not found.")
endif ()

# ----------------------------------------------------------------------------
# System dependencies
# ----------------------------------------------------------------------------

# librt
if ((${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "kFreeBSD")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "GNU"))
    set (LIBRARY_TEMPLATE_LIBRARIES ${LIBRARY_TEMPLATE_LIBRARIES} rt)
endif ()

# libexecinfo -- implicit
check_include_file_cxx (execinfo.h _LIBRARY_TEMPLATE_HAVE_EXECINFO)
mark_as_advanced (_LIBRARY_TEMPLATE_HAVE_EXECINFO)
if (_LIBRARY_TEMPLATE_HAVE_EXECINFO)
    library_template_config_print ("Optional dependency:        libexecinfo found.")
    if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") OR (${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD"))
        set (LIBRARY_TEMPLATE_LIBRARIES ${LIBRARY_TEMPLATE_LIBRARIES} execinfo elf)
    endif ()
else ()
    library_template_config_print ("Optional dependency:        libexecinfo not found.")
endif ()

# ----------------------------------------------------------------------------
# Perform compilability test of platform.hpp (tests some requirements)
# ----------------------------------------------------------------------------

set (CXXSTD_TEST_SOURCE "#include <library_template/platform.hpp>
                         int main() {}")

# using try_compile instead of check_cxx_source_compiles to capture output in case of failure
file (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "${CXXSTD_TEST_SOURCE}\n")

try_compile (LIBRARY_TEMPLATE_PLATFORM_TEST #
             ${CMAKE_BINARY_DIR}
             ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
             CMAKE_FLAGS "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_CXX_FLAGS} ${LIBRARY_TEMPLATE_CXX_FLAGS}"
                         "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_INCLUDE_PATH};${LIBRARY_TEMPLATE_INCLUDE_DIR};${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS}"
             COMPILE_DEFINITIONS ${LIBRARY_TEMPLATE_DEFINITIONS}
             LINK_LIBRARIES ${LIBRARY_TEMPLATE_LIBRARIES}
             OUTPUT_VARIABLE LIBRARY_TEMPLATE_PLATFORM_TEST_OUTPUT)

if (LIBRARY_TEMPLATE_PLATFORM_TEST)
    library_template_config_print ("Library_Template platform.hpp build:  passed.")
else ()
    library_template_config_error ("Library_Template platform.hpp build:  failed!\n\
                        ${LIBRARY_TEMPLATE_PLATFORM_TEST_OUTPUT}")
endif ()

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS LIBRARY_TEMPLATE_INCLUDE_DIR)

# Set LIBRARY_TEMPLATE_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
# This needs to be done, because `find_package(Library_Template)` might be called in any case-sensitive way and we want to
# guarantee that LIBRARY_TEMPLATE_* are always set.
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
    set (LIBRARY_TEMPLATE_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
endforeach ()

# propagate LIBRARY_TEMPLATE_INCLUDE_DIR into LIBRARY_TEMPLATE_INCLUDE_DIRS
set (LIBRARY_TEMPLATE_INCLUDE_DIRS ${LIBRARY_TEMPLATE_INCLUDE_DIR} ${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

if (LIBRARY_TEMPLATE_FOUND AND NOT TARGET library_template::library_template)
    separate_arguments (LIBRARY_TEMPLATE_CXX_FLAGS_LIST UNIX_COMMAND "${LIBRARY_TEMPLATE_CXX_FLAGS}")

    add_library (library_template_library_template INTERFACE)
    target_compile_definitions (library_template_library_template INTERFACE ${LIBRARY_TEMPLATE_DEFINITIONS})
    target_compile_options (library_template_library_template INTERFACE ${LIBRARY_TEMPLATE_CXX_FLAGS_LIST})
    target_link_libraries (library_template_library_template INTERFACE "${LIBRARY_TEMPLATE_LIBRARIES}")
    # include library_template/include/ as -I, because library_template should never produce warnings.
    target_include_directories (library_template_library_template INTERFACE "${LIBRARY_TEMPLATE_INCLUDE_DIR}")
    # include everything except library_template/include/ as -isystem, i.e.
    # a system header which suppresses warnings of external libraries.
    target_include_directories (library_template_library_template SYSTEM
                                INTERFACE "${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS}")
    add_library (library_template::library_template ALIAS library_template_library_template)
endif ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (LIBRARY_TEMPLATE_FIND_DEBUG)
    message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
    message ("")
    message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
    message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
    message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
    message ("  LIBRARY_TEMPLATE_INCLUDE_DIR          ${LIBRARY_TEMPLATE_INCLUDE_DIR}")
    message ("")
    message ("  ${CMAKE_FIND_PACKAGE_NAME}_FOUND                ${${CMAKE_FIND_PACKAGE_NAME}_FOUND}")
    message ("  LIBRARY_TEMPLATE_HAS_ZLIB             ${ZLIB_FOUND}")
    message ("  LIBRARY_TEMPLATE_HAS_BZIP2            ${BZIP2_FOUND}")
    message ("")
    message ("  LIBRARY_TEMPLATE_INCLUDE_DIRS         ${LIBRARY_TEMPLATE_INCLUDE_DIRS}")
    message ("  LIBRARY_TEMPLATE_LIBRARIES            ${LIBRARY_TEMPLATE_LIBRARIES}")
    message ("  LIBRARY_TEMPLATE_DEFINITIONS          ${LIBRARY_TEMPLATE_DEFINITIONS}")
    message ("  LIBRARY_TEMPLATE_CXX_FLAGS            ${LIBRARY_TEMPLATE_CXX_FLAGS}")
    message ("")
    message ("  LIBRARY_TEMPLATE_VERSION              ${LIBRARY_TEMPLATE_VERSION}")
    message ("  LIBRARY_TEMPLATE_VERSION_MAJOR        ${LIBRARY_TEMPLATE_VERSION_MAJOR}")
    message ("  LIBRARY_TEMPLATE_VERSION_MINOR        ${LIBRARY_TEMPLATE_VERSION_MINOR}")
    message ("  LIBRARY_TEMPLATE_VERSION_PATCH        ${LIBRARY_TEMPLATE_VERSION_PATCH}")
endif ()
