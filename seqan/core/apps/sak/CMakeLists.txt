cmake_minimum_required (VERSION 2.6)
project (core_apps_sak)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(sak sak.cpp)

# Install sak in bin/ directory.
install(TARGETS sak
        DESTINATION bin)

# Run app tests if python was found.
if (PYTHONINTERP_FOUND)
  add_test(NAME app_test_sak COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  #message(STATUS "add_test(NAME app_test_sak COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif (PYTHONINTERP_FOUND)
