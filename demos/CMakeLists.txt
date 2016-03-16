# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /demos/CMakeLists.txt
#
# CMakeLists.txt file for the core demos.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_demos CXX)

# Recursively register all demos in core and prefix their target/binary names
# with "demo_".  The demos will be built in DEVELOP mode, copied in
# SEQAN_RELEASE mode and left alone otherwise.
seqan_register_demos (demo_)
