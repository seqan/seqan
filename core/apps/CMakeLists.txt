# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/CMakeLists.txt
#
# CMakeLists.txt file for the core apps.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_apps)
message (STATUS "Configuring core/apps")

# Register all apps in core.  The apps will be built in SEQAN_RELEASE mode or
# the apps called ${app} in APP:${app} mode.
seqan_register_apps ()
