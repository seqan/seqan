#!/bin/bash

# Plot ROIs 9 per page to PDF with links to IGV.
#
# USAGE: roi_plot_9.sh -i IN.roi -o OUT.pdf
#
# Maximal number of header lines.
MAX_HEADER=1000

# The parameters that we will pass to sort.
# The arguments will go here.
SORT_BY=beginPos
SORT_COL=0
REVERSE=

# Current directory.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

helpme () {
   echo "host: ${HOST}"
   echo "Plot ROIs 9 per page to PDF with links to IGV."
   echo ""
   echo "USAGE: roi_plot_9.sh -i IN.roi -o OUT.pdf"
   echo ""
}

# Parse option values.
while getopts "i:o:" opt; do
    case $opt in
        i)
            IN_FILE=$OPTARG
            ;;
        o)
            OUT_FILE=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            helpme
            exit 1
            ;;
    esac
done

# Check that -i or -o are given.
if [[ "$IN_FILE" == "" || "$OUT_FILE" == "" ]]; then
    echo "Missing option -i or -o" >&2
    exit 1
fi

USED_INFILE=${IN_FILE}
if [ "x${IN_FILE##*.}" == "xgz" ]; then
   zcat $IN_FILE > ${IN_FILE}.tmp.roi
   USED_INFILE=${IN_FILE}.tmp.roi
fi

   gawk -v fileName=${IN_FILE}.tmp.ps \
        -f ${DIR}/plot.awk \
        ${USED_INFILE} |gnuplot  2> /dev/null
   gawk -v roiFile=${USED_INFILE} \
        -f ${DIR}/ps2pswLinks.gawk \
        ${IN_FILE}.tmp.ps > ${IN_FILE}.tmp.ln.ps
   ps2pdf ${IN_FILE}.tmp.ln.ps ${OUT_FILE}
#   rm ${IN_FILE}.tmp.ps
#   rm ${IN_FILE}.tmp.ln.ps
if [ "x${IN_FILE##*.}" == "xgz" ]; then
   echo rm $USED_INFILE
fi


