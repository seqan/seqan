cmake_minimum_required (VERSION 2.6)
project (core_apps_razers2)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(razers2 razers.cpp
                             razers.h
                             param_tabs.h
                             param_tabs.cpp
                             param_tabs.inc
                             outputFormat.h
                             paramChooser.h
                             razers_matepairs.h
                             readSimulator.h)

# Install razers2 in /bin directory
install(TARGETS razers2
        DESTINATION bin)

# Run app tests if python was found
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_razers2 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  #message(STATUS "add_test(NAME app_test_razers2 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif(PYTHONINTERP_FOUND)
