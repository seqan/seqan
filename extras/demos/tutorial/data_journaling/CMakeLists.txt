cmake_minimum_required (VERSION 2.6)
project (seqan_extras_demos_tutorials_data_journaling)

# Supress unused parameter warnings for these demos.
if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

seqan_setup_apps(extras_demos_tutorials_data_journaling)
seqan_add_all_executables(extras_demos_tutorials_data_journaling tutorial_data_journaling_)