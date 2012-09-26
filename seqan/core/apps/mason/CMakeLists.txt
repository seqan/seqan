cmake_minimum_required (VERSION 2.6)
project (core_apps_mason)

if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(mason mason.cpp
			   mason.h
			   simulate_454.h
			   simulate_454_base_calling.h
			   simulate_illumina.h
			   simulate_sanger.h
			   store_config.h
			   util.h)

# Install mason in /bin directory
install(TARGETS mason
        DESTINATION bin)

# Run app tests if python was found
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_mason COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  #message(STATUS "add_test(NAME app_test_mason COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif(PYTHONINTERP_FOUND)

