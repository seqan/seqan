#!/bin/sh
#
# Output generation script for breakpoint_calculator

BREAKPOINT_CALCULATOR=../../../../build/debug_linux/extras/apps/breakpoint_calculator/breakpoint_calculator

# ============================================================
# First Section
# ============================================================

${BREAKPOINT_CALCULATOR} -d2 alignment.xmfa > d2_xmfa.stdout
${BREAKPOINT_CALCULATOR} -d2 -d alignment.maf > d2_maf.stdout
${BREAKPOINT_CALCULATOR} -d3 -f xmfa alignment.xmfa > d3_xmfa.stdout
${BREAKPOINT_CALCULATOR} -d3 -d alignment.maf > d3_maf.stdout
