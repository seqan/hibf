# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library-template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# This file describes where and which parts of Library_Template should be installed to.

cmake_minimum_required (VERSION 3.14)

include (GNUInstallDirs)

# install documentation files in /share/doc
install (FILES "${LIBRARY_TEMPLATE_CLONE_DIR}/CHANGELOG.md" #
               "${LIBRARY_TEMPLATE_CLONE_DIR}/CODE_OF_CONDUCT.md" #
               "${LIBRARY_TEMPLATE_CLONE_DIR}/CONTRIBUTING.md" #
               "${LIBRARY_TEMPLATE_CLONE_DIR}/LICENSE.md" #
               "${LIBRARY_TEMPLATE_CLONE_DIR}/README.md"
         TYPE DOC)

# install cmake files in /share/cmake
install (FILES "${LIBRARY_TEMPLATE_CLONE_DIR}/build_system/library_template-config.cmake"
               "${LIBRARY_TEMPLATE_CLONE_DIR}/build_system/library_template-config-version.cmake"
         DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/library_template")

# install library_template header files in /include/library_template
install (DIRECTORY "${LIBRARY_TEMPLATE_INCLUDE_DIR}/library_template" TYPE INCLUDE)

# install submodule header files, e.g. all external dependencies in /home/user/library_template/submodules/*,
# in /include/library_template/submodules/<submodule>/include
foreach (submodule_dir ${LIBRARY_TEMPLATE_DEPENDENCY_INCLUDE_DIRS})
    # e.g. submodule_dir: (1) /home/user/library_template/submodules/sdsl-lite/include or (2) /usr/include
    # strip /home/user/library_template/submodules/ and /include part.
    file (RELATIVE_PATH submodule "${LIBRARY_TEMPLATE_SUBMODULES_DIR}/submodules" "${submodule_dir}/..")
    # submodule is either a single module name, like sdsl-lite or a relative path to a folder ../../../usr
    # skip relative folders and only keep submodules that reside in the submodules folder
    if (NOT submodule MATCHES "^\\.\\.") # skip relative folders
        install (DIRECTORY "${submodule_dir}"
                 DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/library_template/submodules/${submodule}")
    endif ()
endforeach ()
