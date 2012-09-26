cmake_minimum_required (VERSION 2.6)
project (core_apps_splazers)


# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_splazers COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)



# Update the list of file names below if you add source files to your application.
seqan_add_executable(splazers splazers.cpp
			      outputFormat.h
			      paramChooser.h
			      razers.h
			      razers_matepairs.h
			      razers_parallel.h
			      razers_spliced.h
			      readSimulator.h)

install(TARGETS splazers
        DESTINATION bin)
