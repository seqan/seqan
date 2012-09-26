cmake_minimum_required (VERSION 2.6)
project (seqan_extras_demos)

seqan_setup_apps(extra_demos)

if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

seqan_add_all_executables(extras_demos)

seqan_add_all_subdirectories()
