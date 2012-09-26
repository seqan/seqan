cmake_minimum_required (VERSION 2.6)
project (seqan_core_demos_tutorials_rnaseq)

if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

if (OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

seqan_setup_apps(core_demos_tutorials_rnaseq)
seqan_add_all_executables(core_demos_tutorials_rnaseq tutorial_)
