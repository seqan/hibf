# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

message (STATUS "Finding HIBF (${HIBF_VERSION}) and checking requirements")

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

set (CPM_INDENT "  CMake Package Manager CPM: ")
include (${CMAKE_CURRENT_LIST_DIR}/CPM.cmake)
CPMUsePackageLock (${CMAKE_CURRENT_LIST_DIR}/package-lock.cmake)

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (CheckCXXSourceRuns)
include (CheckCXXCompilerFlag)
include (CMakeDependentOption)
include (CheckIPOSupported)

# ----------------------------------------------------------------------------
# Find or add dependencies
# ----------------------------------------------------------------------------

CPMGetPackage (cereal)

if (cereal_ADDED)
    set_target_properties (cereal PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                                             $<TARGET_PROPERTY:cereal,INTERFACE_INCLUDE_DIRECTORIES>)
endif ()

include (${CMAKE_CURRENT_LIST_DIR}/simde.cmake)

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
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET 1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} ${HIBF_HEADER_PATH})
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Check supported compilers
# ----------------------------------------------------------------------------

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12)
    hibf_config_error ("GCC < 12 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
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
# Require C++23
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
set (HIBF_CXX_FLAGS "")

set (CXXSTD_TEST_SOURCE
     "#if !defined (__cplusplus) || (__cplusplus < 202100)
      #error NOCXX23
      #endif
      int main() {}")

set (HIBF_FEATURE_CPP23_FLAG_BUILTIN "")
set (HIBF_FEATURE_CPP23_FLAG_STD23 "-std=c++23")

set (HIBF_CPP23_FLAG "")

foreach (_FLAG BUILTIN STD23)
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_SAVE} ${HIBF_FEATURE_CPP23_FLAG_${_FLAG}}")

    check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CPP23_FLAG_${_FLAG})

    if (CPP23_FLAG_${_FLAG})
        set (HIBF_CPP23_FLAG ${_FLAG})
        break ()
    endif ()
endforeach ()

set (CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})

if (HIBF_CPP23_FLAG STREQUAL "BUILTIN")
    hibf_config_print ("C++ Standard-23 support:    builtin")
elseif (HIBF_CPP23_FLAG)
    list (APPEND HIBF_CXX_FLAGS "${HIBF_FEATURE_CPP23_FLAG_${HIBF_CPP23_FLAG}}")
    hibf_config_print ("C++ Standard-23 support:    via ${HIBF_FEATURE_CPP23_FLAG_${HIBF_CPP23_FLAG}}")
else ()
    hibf_config_error ("HIBF requires C++23, but your compiler does not support it.")
endif ()

# ----------------------------------------------------------------------------
# Required: 128 bit unsigned integer extension
# ----------------------------------------------------------------------------

set (HIBF_128BIT_TEST_SOURCE
     "#if !defined (__SIZEOF_INT128__)
      #error NO128BITS
      #endif
      int main() {}")

check_cxx_source_compiles ("${HIBF_128BIT_TEST_SOURCE}" HIBF_128BIT_SUPPORT)

if (HIBF_128BIT_SUPPORT)
    hibf_config_print ("128 bit support:            builtin")
else ()
    hibf_config_error ("HIBF requires 128 bit integers, but your compiler does not support it.")
endif ()

# ----------------------------------------------------------------------------
# Required: OpenMP
# ----------------------------------------------------------------------------

find_package (OpenMP QUIET COMPONENTS CXX)

if (OpenMP_FOUND)
    set (HIBF_DEFINITIONS ${HIBF_DEFINITIONS} "-DSIMDE_ENABLE_OPENMP")
    hibf_config_print ("OpenMP support:             via ${OpenMP_CXX_FLAGS}")
else ()
    hibf_config_error ("HIBF requires OpenMP, but your compiler does not support it.")
endif ()

check_cxx_compiler_flag ("-Wno-psabi" HIBF_SUPPRESS_GCC4_ABI)
if (HIBF_SUPPRESS_GCC4_ABI)
    list (APPEND HIBF_CXX_FLAGS "-Wno-psabi")
    hibf_config_print ("Suppressing GCC 4 warnings: via -Wno-psabi")
endif ()

# ----------------------------------------------------------------------------
# Optional: seqan::hibf::bit_vector::resize_for_overwrite
# ----------------------------------------------------------------------------

set (HIBF_UNINITIALISED_RESIZE_TEST_SOURCE
     "#include <cstddef>
      #include <cstdint>
      #include <vector>

      struct my_vector : public std::vector<uint64_t>
      {
          using base_t = std::vector<uint64_t>;

      #ifdef _LIBCPP_VERSION
          inline void resize_for_overwrite(size_t const size)
          {
              struct fake_vector
              {
                  using allocator_t = typename base_t::allocator_type;
                  using pointer = typename std::allocator_traits<allocator_t>::pointer;

                  pointer begin;
                  pointer end;
                  std::__compressed_pair<pointer, allocator_t> end_cap;
              };

              static_assert(sizeof(fake_vector) == sizeof(base_t));
              static_assert(alignof(fake_vector) == alignof(base_t));

              if (size > base_t::capacity())
                  base_t::reserve(size);

              fake_vector & vec = reinterpret_cast<fake_vector &>(*this);
              vec.end = vec.begin + size;
          }
      #else
          inline void resize_for_overwrite(size_t const size)
          {
              if (size > base_t::capacity())
                  base_t::reserve(size);

              this->_M_impl._M_finish = this->_M_impl._M_start + size;
          }
      #endif
      };

      int main()
      {
          my_vector vec{};
          vec.resize_for_overwrite(10u);
          return vec.size() != 10u;
      }
     ")

check_cxx_source_runs ("${HIBF_UNINITIALISED_RESIZE_TEST_SOURCE}" HIBF_UNINITIALISED_RESIZE_SUPPORT)

if (HIBF_UNINITIALISED_RESIZE_SUPPORT)
    hibf_config_print ("Unitialised resize support: enabled")
    set (HIBF_DEFINITIONS ${HIBF_DEFINITIONS} "-DHIBF_UNINITIALISED_RESIZE")
else ()
    hibf_config_print ("Unitialised resize support: disabled")
endif ()

# ----------------------------------------------------------------------------
# Optimizations
# ----------------------------------------------------------------------------

if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR "${CMAKE_BUILD_TYPE}" MATCHES "Coverage")
    set (HIBF_IS_DEBUG TRUE)
else ()
    set (HIBF_IS_DEBUG FALSE)
endif ()

check_cxx_compiler_flag ("-march=native" HIBF_HAS_MARCH_NATIVE)
cmake_dependent_option (HIBF_NATIVE_BUILD "Optimize build for current architecture." ON "HIBF_HAS_MARCH_NATIVE" OFF)

if (HIBF_NATIVE_BUILD)
    list (APPEND HIBF_CXX_FLAGS "-march=native")
    hibf_config_print ("Optimize build:             via -march=native")
else ()
    check_cxx_compiler_flag ("-mpopcnt" HIBF_HAS_POPCNT)
    if (HIBF_HAS_POPCNT)
        list (APPEND HIBF_CXX_FLAGS "-mpopcnt")
        hibf_config_print ("Optimize build:             via -mpopcnt")
    else ()
        hibf_config_print ("Optimize build:             disabled")
    endif ()
endif ()

check_ipo_supported (RESULT HIBF_HAS_LTO OUTPUT HIBF_HAS_LTO_OUTPUT)
cmake_dependent_option (HIBF_LTO_BUILD "Use Link Time Optimisation." ON "HIBF_HAS_LTO" OFF)
cmake_dependent_option (HIBF_DEV_CHECK_LTO "LTO check." ON
                        "HIBF_LTO_BUILD;NOT HIBF_IS_TOP_LEVEL;NOT CMAKE_INTERPROCEDURAL_OPTIMIZATION" OFF)

if (HIBF_DEV_CHECK_LTO)
    message (FATAL_ERROR " HIBF heavily benefits from Link Time Optimisation (LTO).\n"
                         " Your compiler supports LTO, but is has not been enabled for your project.\n \n"
                         " Add the following at the beginning of your CMakeLists.txt:\n"
                         " ```\n"
                         " include (CheckIPOSupported)\n"
                         " check_ipo_supported (RESULT result OUTPUT output)\n"
                         " if (result)\n"
                         "     set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)\n"
                         " endif ()\n"
                         " ```"
                         " \n \n"
                         " For development purposes, this error can be disabled via `HIBF_DEV_CHECK_LTO=OFF`.")
elseif (HIBF_LTO_BUILD)
    set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    hibf_config_print ("Link Time Optimisation:     enabled")
else ()
    set (HIBF_LTO_LOG "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/hibf.lto.log")
    file (WRITE "${HIBF_LTO_LOG}" "${HIBF_HAS_LTO_OUTPUT}")
    hibf_config_print ("Link Time Optimisation:     not supported")
    hibf_config_print ("  See ${HIBF_LTO_LOG}")
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
        list (APPEND HIBF_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
        hibf_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    hibf_config_print ("Thread support:             not found")
endif ()

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
    list (APPEND HIBF_LIBRARIES "rt")
endif ()

# libexecinfo -- implicit
check_include_file_cxx (execinfo.h _HIBF_HAVE_EXECINFO)
mark_as_advanced (_HIBF_HAVE_EXECINFO)
if (_HIBF_HAVE_EXECINFO)
    hibf_config_print ("Optional dependency:        libexecinfo found")
    if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") OR (${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD"))
        list (APPEND HIBF_LIBRARIES "execinfo")
        list (APPEND HIBF_LIBRARIES "elf")
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
list (JOIN HIBF_CXX_FLAGS " " HIBF_TRY_COMPILE_CXX_FLAGS)
# cmake-format: off
try_compile (HIBF_PLATFORM_TEST
             SOURCE_FROM_CONTENT src.cxx "#include <hibf/platform.hpp>\nint main() {}"
             CMAKE_FLAGS "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_CXX_FLAGS} ${HIBF_TRY_COMPILE_CXX_FLAGS}"
                         "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_INCLUDE_PATH};${HIBF_HEADER_PATH}"
             COMPILE_DEFINITIONS ${HIBF_DEFINITIONS}
             LINK_LIBRARIES ${HIBF_LIBRARIES}
             OUTPUT_VARIABLE HIBF_PLATFORM_TEST_OUTPUT)
# cmake-format: on
if (HIBF_PLATFORM_TEST)
    hibf_config_print ("HIBF platform.hpp build:    passed")
else ()
    set (HIBF_PLATFORM_TEST_LOG "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/hibf.platform.compile.log")
    file (WRITE "${HIBF_PLATFORM_TEST_LOG}" "${HIBF_PLATFORM_TEST_OUTPUT}")
    hibf_config_error ("HIBF platform.hpp build:    failed\n    See ${HIBF_PLATFORM_TEST_LOG}")
endif ()
