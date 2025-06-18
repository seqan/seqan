#!/bin/bash

set -Eeuo pipefail

usage="\
SYNOPSIS
    adjust_copyright_years.sh [new_year=$(date +%Y)]

DESCRIPTION
    Updates the copyright year of files that are formatted in a certain way.

EXAMPLES
    ./util/bin/adjust_copyright_years.sh
    ./util/bin/adjust_copyright_years.sh 2026 # Overwrite default year (current year)
"

if [[ $# -gt 1 ]]; then
    echo -e "$usage"
    exit 1
fi

YEAR="${1:-$(date +%Y)}"

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
SEQAN_ROOT="$(realpath "${SCRIPT_DIR}/../..")"

find "${SEQAN_ROOT}" \
    \( \
        -not -path '*/\.git/*' -and \
        -not -path '*/build/*' -and \
        -not -path '*/util/py_lib/seqan/dox/*' \
    \) \
    -and \
    \( \
        -name "*.bat" -or \
        -name "*.cmake" -or \
        -name "*.cpp" -or \
        -name "*.h" -or \
        -name "*.hpp" -or \
        -name "*.md" -or \
        -name "*.rst" -or \
        -name "*.sh" -or \
        -name "*.txt" -or \
        -name "*.yml" -or \
        -name "COPYRIGHT" -or \
        -name "LICENSE" \
    \) \
    -type f \
    -exec sed -i -e "s/20\(..\)-20../20\1-${YEAR}/g" {} \;
