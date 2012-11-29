# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/CMakeLists.txt
#
# CMakeLists.txt file for the extras apps.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps)

# Register all apps in extras.  The apps will be built in SEQAN_RELEASE mode or
# the apps called ${app} in APP:${app} mode.
seqan_register_apps ()
