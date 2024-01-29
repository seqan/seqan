#!/usr/bin/env bash
# ==========================================================================
#                 SeqAn - The Library for Sequence Analysis
# ==========================================================================
# Copyright (c) 2006-2024, Knut Reinert, FU Berlin
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Knut Reinert or the FU Berlin nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
# ==========================================================================
#
# Usage: link_check.sh <SeqAn root directory>
# Will output the status of links in the repository.
#
# Of main interest are broken links and those with a "Link STATUS" message.
# Some URLs may not be properly matched by the regex.
#
# The general workflow is to first run the script and then check the non-working links by searching the occurrence
# within the codebase and verifying that they are indeed broken.

do_check ()
{
    RESPONSE=$(curl --http2 -Is -A 'Mozilla/5.0' $1) # HTTP2 is the default.
    if ! [[ "$RESPONSE" =~ ^HTTP.* ]]; then # If this does not work,
        RESPONSE=$(curl --http1.1 -Is -A 'Mozilla/5.0' $1) # fall back to HTTP1.1.
    fi

    HEADER=($(echo $RESPONSE | head -1)); # May look like: HTTP/2 200
    STATUS=${HEADER[1]}
    case "$STATUS" in
        200) echo "Link OK         :" $1;;
        301) echo "Link PERM MOVED :" $1;;
        302) echo "Link TEMP MOVED :" $1;;
        404) echo "Link BROKE      :" $1;;
        429) sleep 5; do_check $1;;
        *)   echo "Link STATUS" $STATUS ":" $1;;
    esac
}

if [[ $# -ne 1 ]]; then
    echo "Usage: link_check.sh <SeqAn root directory>"
    exit 1
fi

if [[ ! -d $1 ]]; then
    echo "The directory $1 does not exist."
    exit 1
fi

if [[ ! -f $1/include/seqan/version.h ]]; then
    echo "The directory $1 does not seem to be the SeqAn root directory."
    echo "Cannot find $1/include/seqan/version.h."
    exit 1
fi

# -o: print only the matching part
# -h: do not print filenames
# -r: recursive
# --binary-files=without-match: do not print matches in binary files
# --exclude-dir: exclude directories, must be directory names, not paths (util/py_lib wouldn't work)
for URL in $(grep -ohr --binary-files=without-match --exclude-dir={.git,.vscode,py_lib,build,jquery-bbq} "https*://[a-zA-Z0-9./#+?=_%:-]*[a-zA-Z0-9/#+?=_%:-]" $1 | sort | uniq)
do
  do_check $URL
done
