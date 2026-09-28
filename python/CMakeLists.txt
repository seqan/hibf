# SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Python bindings via nanobind. Included from the top-level CMakeLists.txt if HIBF_PYTHON_BINDINGS is ON.
# Usually built via `pip install .` (see pyproject.toml).

cmake_minimum_required (VERSION 3.20...3.31)

find_package (Python 3.10 REQUIRED
              COMPONENTS Interpreter Development.Module
              OPTIONAL_COMPONENTS Development.SABIModule)

execute_process (COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                 OUTPUT_STRIP_TRAILING_WHITESPACE
                 OUTPUT_VARIABLE nanobind_ROOT)
find_package (nanobind CONFIG REQUIRED)

# nanobind <= 3.1.0 uses stderr without including <cstdio>. libstdc++ includes it transitively, libc++ does not.
# Fixed upstream in wjakob/nanobind@942e75d. Remove once the minimum required nanobind version contains the fix.
add_compile_options (-include cstdio)

# The static library is linked into a shared Python extension.
set_target_properties (hibf PROPERTIES POSITION_INDEPENDENT_CODE ON)

# NOMINSIZE: The hot insertion path of hibf::insert_iterator is inlined into the bindings; do not optimise for size.
nanobind_add_module (_hibf STABLE_ABI NB_STATIC NOMINSIZE src/hibf_python.cpp)
target_link_libraries (_hibf PRIVATE seqan::hibf)

# cmake-format: off
nanobind_add_stub (_hibf_stub
                   MODULE _hibf
                   OUTPUT _hibf.pyi
                   PYTHON_PATH $<TARGET_FILE_DIR:_hibf>
                   DEPENDS _hibf
                   MARKER_FILE py.typed)
# cmake-format: on

install (TARGETS _hibf LIBRARY DESTINATION hibf)
install (FILES "${CMAKE_CURRENT_BINARY_DIR}/_hibf.pyi" "${CMAKE_CURRENT_BINARY_DIR}/py.typed" DESTINATION hibf)
