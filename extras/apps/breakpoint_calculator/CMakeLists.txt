cmake_minimum_required (VERSION 2.6)
project (sandbox_bkehr_apps_breakpoint_calculator)

# The FindLemon.cmake file is in the current directory.  Thus, add the latter to
# the CMAKE_MODULE_PATH.
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(LEMON QUIET)

if (LEMON_FOUND)
  # Add app tests if Python interpreter could be found.
  if(PYTHONINTERP_FOUND)
    add_test(NAME app_test_breakpoint_calculator COMMAND ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
      ${CMAKE_BINARY_DIR})
  endif(PYTHONINTERP_FOUND)

  # add lemon include path
  include_directories(${LEMON_INCLUDE_DIRS})

  # Update the list of file names below if you add source files to your application.
  seqan_add_executable(breakpoint_calculator breakpoint_calculator.cpp
                                             breakpoint_calculator.h
                                             parse_alignment.h
                                             breakpoint_counts.h)
else (LEMON_FOUND)
  message(STATUS "LEMON library not found, not building breakpoint_calculator.")
endif (LEMON_FOUND)
