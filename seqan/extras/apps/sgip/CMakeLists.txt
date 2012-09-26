cmake_minimum_required (VERSION 2.6)
project (extras_apps_sgip)

# Do not build SGIP on MSVC <10.
set(BUILD_SGIP TRUE)
if (MSVC AND NOT MSVC10)
    set(BUILD_SGIP FALSE)
endif (MSVC AND NOT MSVC10)

if (BUILD_SGIP)
	seqan_add_executable(sgip sgip.cpp sgip.h sgip_base.h sgip_output.h)
endif (BUILD_SGIP)

# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_sgip COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)
