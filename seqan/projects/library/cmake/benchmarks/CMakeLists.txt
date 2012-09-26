cmake_minimum_required (VERSION 2.6)

################################################################################
# For the benchmark apps, use assertions in debug mode, no assertions in
# release mode and never go into testing mode.
################################################################################

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=0")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")

################################################################################
# Include the CMakeLists.txt for each benchmark's app directory.
################################################################################

get_filename_component (SEQAN_BENCHMARK_ABS ${SEQAN_LIBRARY}/../benchmarks ABSOLUTE)
file (GLOB SEQAN_BENCHMARKS ${SEQAN_BENCHMARK_ABS}/[A-z]*)

foreach (BENCHMARK_DIR ${SEQAN_BENCHMARKS})
  get_filename_component(BENCHMARK_NAME ${BENCHMARK_DIR} NAME)
  if (EXISTS ${BENCHMARK_DIR}/apps/CMakeLists.txt)
    add_subdirectory(${BENCHMARK_DIR}/apps ${BENCHMARK_NAME})
  endif (EXISTS ${BENCHMARK_DIR}/apps/CMakeLists.txt)
endforeach (BENCHMARK_DIR ${SEQAN_BENCHMARK})
