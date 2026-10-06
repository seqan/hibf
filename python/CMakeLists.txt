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

# Licences of the third-party code compiled into the module. Files installed into SKBUILD_METADATA_DIR end up in the
# .dist-info directory of the wheel. Dependencies found on the system instead of via CPM are not covered.
if (DEFINED SKBUILD_METADATA_DIR)
    set (HIBF_THIRD_PARTY_LICENSES "${SKBUILD_METADATA_DIR}/licenses/third_party")
    set (HIBF_GENERATED_LICENSES "${CMAKE_CURRENT_BINARY_DIR}/licenses")

    execute_process (COMMAND "${Python_EXECUTABLE}" -c
                             "import importlib.metadata as m; d = m.distribution('nanobind'); \
print(next(d.locate_file(f) for f in d.files if f.name == 'LICENSE'))"
                     OUTPUT_VARIABLE HIBF_NANOBIND_LICENSE
                     OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

    foreach (license IN ITEMS "nanobind|${HIBF_NANOBIND_LICENSE}" "cereal|${cereal_SOURCE_DIR}/LICENSE"
                              "simde|${simde_SOURCE_DIR}/COPYING")
        string (REPLACE "|" ";" license "${license}")
        list (GET license 0 name)
        list (GET license 1 file)
        if (EXISTS "${file}" AND NOT IS_DIRECTORY "${file}")
            install (FILES "${file}" DESTINATION "${HIBF_THIRD_PARTY_LICENSES}" RENAME "${name}.txt")
        else ()
            message (WARNING "Licence of ${name} not found. It will be missing from the wheel.")
        endif ()
    endforeach ()

    # robin_hood.hpp contains its licence in a comment.
    file (READ "${HIBF_SOURCE_DIR}/include/hibf/contrib/robin_hood.hpp" robin_hood_header LIMIT 4096)
    string (REGEX MATCH "// Copyright[^\n]*\n(//[^\n]*\n)*" robin_hood_license "${robin_hood_header}")
    string (REGEX REPLACE "(^|\n)// ?" "\\1" robin_hood_license "${robin_hood_license}")
    file (WRITE "${HIBF_GENERATED_LICENSES}/robin_hood.txt" "${robin_hood_license}")

    # The HyperLogLog implementation states its copyright and licence via SPDX.
    file (READ "${HIBF_SOURCE_DIR}/LICENSES/MIT.txt" hyperloglog_license)
    string (REPLACE "<year> <copyright holders>" "2013 Hideaki Ohno <hide.o.j55{at}gmail.com>" hyperloglog_license
                    "${hyperloglog_license}")
    file (WRITE "${HIBF_GENERATED_LICENSES}/hyperloglog.txt" "${hyperloglog_license}")

    install (FILES "${HIBF_GENERATED_LICENSES}/robin_hood.txt" "${HIBF_GENERATED_LICENSES}/hyperloglog.txt"
             DESTINATION "${HIBF_THIRD_PARTY_LICENSES}")
endif ()
