# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/CMakeLists.txt
#
# CMakeLists.txt file for the core apps.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps CXX)

# Register all apps in core.  The apps will be built in SEQAN_RELEASE mode or
# the apps called ${app} in APP:${app} mode.
seqan_register_apps ()
