cmake_minimum_required (VERSION 2.6)
project (core_apps_pair_align)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(pair_align pair_align.cpp)

# Install pair_align in /bin directory
install(TARGETS pair_align
        DESTINATION bin)

# Run app tests if python was found
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_pair_align COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  #message(STATUS "add_test(NAME app_test_pair_align COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif(PYTHONINTERP_FOUND)