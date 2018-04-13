# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /util/py_lib/CMakeLists.txt
#
# CMakeLists.txt file for Python stuff.
# ===========================================================================
# Kicks off tests using Python nosetests.
# ===========================================================================

# Look for Python and stop if it could not be found.
# require python 2.7, not python3
set(PythonInterp_FIND_VERSION 2.7)
set(PythonInterp_FIND_VERSION_MAJOR 2)
set(PythonInterp_FIND_VERSION_MINOR 7)
set(PythonInterp_FIND_VERSION_COUNT 2)
find_package (PythonInterp)

if (NOT PYTHONINTERP_FOUND)
  message (STATUS "  Python not found, cannot test py_lib.")
  return ()
endif (NOT PYTHONINTERP_FOUND)

# Look for nosetests and stop if it cannot be found
execute_process(COMMAND ${PYTHON_EXECUTABLE} -m "nose"
                RESULT_VARIABLE PYTHON_NOSETESTS_NOT_FOUND
                OUTPUT_VARIABLE _IGNORED
                ERROR_VARIABLE _IGNORED)

if (PYTHON_NOSETESTS_NOT_FOUND)
  message (STATUS "  Python nosetests ('import nose' failed), cannot add tests for seqan.dox")
  return ()
endif (PYTHON_NOSETESTS_NOT_FOUND)

# Adding test for dox.
message (STATUS "  adding nosetests for seqan.dox")
add_test (NAME test_py_lib_dox
          COMMAND nosetests
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/seqan/dox)
