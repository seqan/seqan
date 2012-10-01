cmake_minimum_required (VERSION 2.6)
project (core_apps_seqan_tcoffee)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(seqan_tcoffee seqan_tcoffee.cpp)

# Install seqan_tcoffee in /bin directory
install(TARGETS seqan_tcoffee
        DESTINATION bin)

# Run app test if python was found 
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_seqan_tcoffee COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  #message(STATUS  "add_test(NAME app_test_seqan_tcoffee COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif(PYTHONINTERP_FOUND)
