#!/bin/bash
# ----------------------------------------------------------------------------
#                  SeqAn - The Library for Sequence Analysis
# ----------------------------------------------------------------------------
# Copyright (c) 2006-2025, Knut Reinert, FU Berlin
# All rights reserved.
#
# License: BSD 3-clause

set -Eeuo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <seqan_root>"
    exit 1
fi

directory=$(readlink -m "$1")

if [[ ! -d "${directory}" ]]; then
    echo "Not a valid directory: ${directory}"
    exit 1
fi

if [[ ! -f "${directory}/apps/CMakeLists.txt" ]]; then
    echo "This does not seem to be the SeqAn root directory: ${directory}"
    exit 1
fi

while IFS= read -r -d '' file; do
    # shellcheck disable=SC2016
    sed -i -E 's@set \(SEQAN_APP_VERSION "([0-9]+)\.([0-9]+)\.([0-9]+)"\)@echo "set (SEQAN_APP_VERSION \\"\1.\2.$((\3 + 1))\\")"@ge' "${file}"
done < <(find "${directory}/apps/" -type f -name "CMakeLists.txt" -print0)
