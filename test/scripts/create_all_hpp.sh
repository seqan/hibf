#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Usage: create_all_hpp.sh
# Will create all.hpp in all subdirectories

# Directories for which no all.hpp shall be created.
# The paths are relative to the root directory.
# Exclusions also apply to subdirectories.
EXCLUSIONS=(
    include/hibf/contrib
    include/hibf/cereal
)

# The directory of this script.
SCRIPT_DIR=$(dirname $(readlink -f $0))
# The root directory of the project.
HIBF_DIR=$(realpath "${SCRIPT_DIR}/../..")
# The current year.
CURRENT_YEAR=$(date "+%Y")

if [[ $# -ne 0 ]]; then
    echo "Usage: create_all_hpp.sh"
    echo "This command does not take any arguments."
    exit 1
fi

if [[ ! -d ${HIBF_DIR} ]]; then
    echo "The directory ${HIBF_DIR} does not exist."
    exit 1
fi

if [[ ! -f ${HIBF_DIR}/include/hibf/version.hpp ]]; then
    echo "Cannot find ${HIBF_DIR}/include/hibf/version.hpp."
    exit 1
fi

# Check whether a directory is excluded.
is_excluded()
{
    directory="$1"
    for entry in "${EXCLUSIONS[@]}"; do
        if [ ! -z $(echo "${directory}" | grep $entry) ]; then
            return 0
        fi
    done
    return 1
}

# The expected license header.
print_header()
{
cat << EOF
// SPDX-FileCopyrightText: 2006-${CURRENT_YEAR}, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-${CURRENT_YEAR}, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
EOF
}

# Generate the include list for the all.hpp file.
print_includes()
{
    directory="$1"
    echo -ne "\n#pragma once\n\n"
    hasIncludes=0

    for path in $(find ${directory}/* -maxdepth 0 -not -name "all.hpp"); do
        include_path="$(realpath --relative-to ${HIBF_DIR}/include "${path}")"
        # Include all headers in the directory.
        if [ -f "${path}" ] && [[ "${path}" == *.hpp ]]; then
            echo "#include <${include_path}>"
        fi
        # Include all `all.hpp` files in the subdirectories, unless excluded.
        if [ -d "${path}" ] && ! is_excluded "${path}"; then
            echo "#include <${include_path}/all.hpp>"
        fi
    done
}

# Iterate through all folders
for directory in $(find ${HIBF_DIR}/include/hibf -type d); do
    if is_excluded "${directory}"; then
        continue
    fi

    # create all.hpp
    if [ ! -f "${directory}/all.hpp" ]; then
        echo "Missing all.hpp file in ${directory}"
        (
            print_header
            print_includes "${directory}"

        ) > "${directory}/all.hpp"
    else
        # Check whether license text has changed.
        HEADER_LINES=$(print_header | wc -l)
        licence_changes=$(head -n ${HEADER_LINES} "${directory}/all.hpp" | diff /dev/stdin <( print_header) | wc -c)

        # Check whether include list has changed.
        INCLUDE_LINES=$(print_includes "${directory}" | wc -l)
        include_changes=$(tail -n ${INCLUDE_LINES} "${directory}/all.hpp" | diff /dev/stdin <( print_includes "${directory}" ) | wc -c)

        if [ "$licence_changes" -gt 0 ] || [ "$include_changes" -gt 0 ]; then
            echo "License text or includes changed, updating ${directory}/all.hpp"
            tmp_file=$(mktemp)
            # Create new license text and include list, but keep doxygen documentation.
            cp "${directory}/all.hpp" "${tmp_file}"
            (
                PRAGMA_LINE=$(cat "${tmp_file}" | grep -n "#pragma once" | cut -f 1 -d ':')
                print_header
                tail ${tmp_file} -n +$(expr 1 + ${HEADER_LINES}) | grep -v "#include" | grep -v "#pragma once"
                print_includes "${directory}"
            ) | cat -s | sed -e :a -e '/^\n*$/{$d;N;};/\n$/ba' >  "${directory}/all.hpp"
            rm $tmp_file
        fi
    fi
done
