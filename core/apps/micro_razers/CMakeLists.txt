cmake_minimum_required (VERSION 2.6)
project (core_apps_micro_razers)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(micro_razers micro_razers.cpp outputFormat.h ../razers/paramChooser.h ../razers/param_tabs.h ../razers/param_tabs.cpp ../razers/param_tabs.inc
)

# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_micro_razers COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)

install(TARGETS micro_razers
        DESTINATION bin)



