cmake_minimum_required (VERSION 2.6)
project (seqan_demos_tutorials CXX)

# Supress unused parameter warnings for these demos.
if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

seqan_add_all_executables(demos_tutorials tutorial_)

seqan_add_all_subdirectories()
