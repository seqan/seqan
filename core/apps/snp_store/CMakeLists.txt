cmake_minimum_required (VERSION 2.6)
project (core_apps_snp_store)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../extras/include)

if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  seqan_add_executable(snp_store snp_store.cpp snp_store.h)
  install(TARGETS snp_store
          DESTINATION bin)
else (Boost_FOUND)
  message(STATUS "Boost not found, not building snp_store.")
endif (Boost_FOUND)

if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_snp_store COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)

