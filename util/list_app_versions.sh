#!/bin/bash
# ----------------------------------------------------------------------------
#                  SeqAn - The Library for Sequence Analysis
# ----------------------------------------------------------------------------
# Copyright (c) 2006-2026, Knut Reinert, FU Berlin
# All rights reserved.
#
# License: BSD 3-clause

set -Eeuo pipefail

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
SEQAN_ROOT="$(realpath "${SCRIPT_DIR}/..")"

while IFS= read -r -d '' file; do
    # Extract version
    version=$(grep -oP 'set \(SEQAN_APP_VERSION "\K[^"]+' "${file}" || true)

    if [[ -z "${version}" ]]; then
        continue
    fi

    # Extract all install TARGETS commands (handling multi-line)
    # Read entire file and extract targets between "install (TARGETS" and "DESTINATION"
    targets=$(awk '
        /install \(TARGETS/ {
            in_install = 1
            line = $0
        }
        in_install {
            if ($0 !~ /install \(TARGETS/) {
                line = line " " $0
            }
            if (/DESTINATION/) {
                # Extract content between TARGETS and DESTINATION
                match(line, /install \(TARGETS[ \t\n]+([^D]+)DESTINATION/, arr)
                if (arr[1] != "") {
                    gsub(/[ \t\n]+/, " ", arr[1])
                    gsub(/^ +| +$/, "", arr[1])
                    print arr[1]
                }
                in_install = 0
                line = ""
            }
        }
    ' "${file}")

    # Output each target with version
    if [[ -n "${targets}" ]]; then
        while IFS= read -r target_line; do
            if [[ -n "${target_line}" ]]; then
                # Split by spaces to get individual target names
                for target in ${target_line}; do
                    echo "${target}	${version}"
                done
            fi
        done <<< "${targets}"
    fi
done < <(find "${SEQAN_ROOT}/apps/" -type f -name "CMakeLists.txt" -print0) | sort
