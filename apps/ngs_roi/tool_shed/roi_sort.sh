#!/bin/bash

# Sorting of ROI files.
#
# USAGE: sort_roi.sh [OPTIONS] -i IN.roi -o OUT.roi
#
# Options:
#   -r     reverse orientation
#   -p     sort by position (ref, start, end) -- default
#   -c COL sort by column COL
#   -n COL sort by column COL, compare as number

# Maximal number of header lines.
MAX_HEADER=1000

# The parameters that we will pass to sort.
SORT_POS_ARGS="-k 1,1 -k 2,2n -k 3,3n"
SORT_POS_ARGS_REV="-k 1,1r -k 2,2nr -k 3,3nr"
SORT_COL_ARGS="-k"

# The arguments will go here.
SORT_BY=beginPos
SORT_COL=0
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
        p)
            SORT_BY=beginPos
            ;;
        c)
            SORT_BY=c
            SORT_COL=$OPTARG
            ;;
        n)
            SORT_BY=n
            SORT_COL=$OPTARG
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
    c)
        SORT_ARGS="-k ${SORT_COL},${SORT_COL}${REVERSE}"
        ;;
    n)
        SORT_ARGS="-k ${SORT_COL},${SORT_COL}g${REVERSE}"
        ;;
    beginPos)
        if [[ "$REVERSE" == "r" ]]; then
            SORT_ARGS=${SORT_POS_ARGS_REV}
        else
            SORT_ARGS=${SORT_POS_ARGS}
        fi
        ;;
esac

# Execute sorting.
#echo "SORT_ARGS=${SORT_ARGS}" 1>&2
(
    export LC_ALL=C
    head -n ${MAX_HEADER} ${IN_FILE} | grep '^#';
    #echo "sort ${SORT_ARGS} <(grep -v '^#' ${IN_FILE});" 1>&2
    sort ${SORT_ARGS} <(grep -v '^#' ${IN_FILE});
) > ${OUT_FILE}

exit $?
