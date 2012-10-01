cmake_minimum_required (VERSION 2.6)
project (core_apps_stellar)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(stellar stellar.cpp
			     stellar.h
                             stellar_extension.h
			     stellar_output.h
			     stellar_types.h)

install(TARGETS stellar
        DESTINATION bin)

# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_stellar COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)

