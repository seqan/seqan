#!/bin/bash

# Sorting of GFF files.
#
# USAGE: sort_gff.sh [OPTIONS] -i IN.roi -o OUT.roi
#
# Options:
#   -r     reverse orientation
#   -p     sort by position (ref, start, end) -- default

# The parameters that we will pass to sort.
SORT_POS_ARGS="-k 1,1 -k 4,4n -k 5,5n"
SORT_POS_ARGS_REV="-k 1,1r -k 4,4nr -k 5,5nr"

# The arguments will go here.
SORT_BY=beginPos
REVERSE=

# Parse option values.
while getopts "pc:i:o:n:r" opt; do
    case $opt in
        i)
            IN_FILE=$OPTARG
            ;;
        o)
            OUT_FILE=$OPTARG
            ;;
        r)
            REVERSE=r
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

# Check that -i or -o are given.
if [[ "$IN_FILE" == "" || "$OUT_FILE" == "" ]]; then
    echo "Missing option -i or -o" >&2
    exit 1
fi

# Setup sort args.
case $SORT_BY in
    beginPos)
        if [[ "$REVERSE" == "r" ]]; then
            SORT_ARGS=${SORT_POS_ARGS_REV}
        else
            SORT_ARGS=${SORT_POS_ARGS}
        fi
        ;;
esac

# Execute sorting.
#echo "OUT_FILE=${OUT_FILE}" 2>&2
#echo "SORT_ARGS=${SORT_ARGS}" 1>&2
(
    export LC_ALL=C
    #echo "sort ${SORT_ARGS} <(grep -v '^#' ${IN_FILE});" 1>&2
    sort ${SORT_ARGS} <(grep -v '^#' ${IN_FILE});
) > ${OUT_FILE}

exit $?
