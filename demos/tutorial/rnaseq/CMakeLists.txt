cmake_minimum_required (VERSION 2.6)
project (seqan_demos_tutorials_rnaseq)

if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

if (OpenMP_FLAG_DETECTED)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OpenMP_FLAG_DETECTED)

seqan_setup_apps(demos_tutorials_rnaseq)
seqan_add_all_executables(demos_tutorials_rnaseq tutorial_)
