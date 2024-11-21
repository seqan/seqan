# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/CMakeLists.txt
#
# CMakeLists.txt file for the core apps.
# ===========================================================================

cmake_minimum_required (VERSION 3.12)
project (seqan_apps CXX)

# Register all apps in core.  The apps will be built in SEQAN_RELEASE mode or
# the apps called ${app} in APP:${app} mode.

# Strip the build path from app binaries. These are included by macros such as FAIL().
# Used for packaging, e.g., Debian.
if (DEFINED SEQAN_STRIP_APPS_BUILD_PATH)
    set (SEQAN_CXX_FLAGS_BACKUP ${SEQAN_CXX_FLAGS})
    set (SEQAN_CXX_FLAGS ${SEQAN_CXX_FLAGS}  -ffile-prefix-map=${SEQAN_STRIP_APPS_BUILD_PATH}=.)
    seqan_register_apps ()
    set (SEQAN_CXX_FLAGS ${SEQAN_CXX_FLAGS_BACKUP})
    unset (SEQAN_CXX_FLAGS_BACKUP)
else ()
    seqan_register_apps ()
endif ()
