cmake_minimum_required (VERSION 2.8)
project (seqan_extras_demos_tutorials_annotation_store)

if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

seqan_setup_apps(core_extras_tutorials_annotation_store)
seqan_add_all_executables(core_extras_tutorials_annotation_store tutorial_)
