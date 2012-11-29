# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/extras/CMakeLists.txt
#
# CMakeLists.txt file for the extras demos.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_demos)
message (STATUS "Configuring extras/demos")

# Recursively register all demos in core and prefix their target/binary names
# with "demo_".  The demos will be built in DEVELOP mode, copied in
# SEQAN_RELEASE mode and left alone otherwise.
seqan_register_demos (demo_)
