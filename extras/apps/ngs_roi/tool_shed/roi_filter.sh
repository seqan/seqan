#!/bin/bash

# Sorting of ROI files.
#
# USAGE: sort_filter.sh [-n] -c COL -o OP -v VALUE <IN.roi >OUT.roi
#
# Read lines from IN.roi, write header lines.  Otherwise, check that the
# given column is in relation OP to the value and only write if this is
# fulfilled.  Use -n to compare as number.

# Maximal number of header lines.
MAX_HEADER=1000

# The arguments will go here.
COL=
OP=
VALUE=
TXT=1

# Parse option values.
while getopts "c:o:v:n" opt; do
    case $opt in
        c )
            COL=$OPTARG
            ;;
        o )
            case $OPTARG in
                eq)
                    OP="=="
                    ;;
                geq)
                    OP=">="
                    ;;
                leq)
                    OP="<="
                    ;;
                gt)
                    OP=">"
                    ;;
                lt)
                    OP="<"
                    ;;
                neq)
                    OP="!="
                    ;;
                * )
                    echo "Invalid operator $OPTARG" >&2
                    exit 1
                    ;;
            esac
            ;;
        v )
            VALUE=$OPTARG
            ;;
        n )
            TXT=0
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

shift $(($OPTIND - 1))

# Check that -i or -o are given.
if [ "x$COL" = "x" ] || [ "x$OP" = "x" ] || [ "x$VALUE" = "x" ] ; then
    echo "Missing option -c, -o, or -v" >&2
    echo "Col: $COL" >&2
    echo "OP: $OP" >&2
    echo "VALUE: $VALUE" >&2
    exit 1
fi

if [[ "x$TXT" == "x1" ]]; then
    VALUE="\"$VALUE\""
fi

# Execute filtering.
#echo awk "/^#/ { print \$0 } !/^#/ { if (\$$COL $OP $VALUE) print \$0 }" >&2
gawk "/^#/ { print \$0 } !/^#/ { if (\$$COL $OP $VALUE) print \$0 }"

exit $?
