cmake_minimum_required (VERSION 2.6)
project (core_apps_tree_recon)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(tree_recon tree_recon.cpp)

# Install tree_recon in bin/ directory.
install(TARGETS tree_recon
        DESTINATION bin)

# Run app tests if python was found.
if (PYTHONINTERP_FOUND)
  add_test(NAME app_test_tree_recon COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  #message(STATUS "add_test(NAME app_test_tree_recon COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif (PYTHONINTERP_FOUND)
