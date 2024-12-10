# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /util/py_lib/CMakeLists.txt
#
# CMakeLists.txt file for Python stuff.
# ===========================================================================
# Kicks off tests using Python nose2.
# ===========================================================================

# Look for Python and stop if it could not be found.
if (NOT Python3_Interpreter_FOUND)
  find_package (Python3)
endif ()

if (NOT Python3_Interpreter_FOUND)
  message (STATUS "  Skipping py_lib tests")
  return ()
endif (NOT Python3_Interpreter_FOUND)

# Look for nose2 and stop if it cannot be found
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import nose2"
                RESULT_VARIABLE PYTHON_NOSE2_NOT_FOUND
                OUTPUT_VARIABLE _IGNORED
                ERROR_VARIABLE _IGNORED)

if (PYTHON_NOSE2_NOT_FOUND)
  message (STATUS "Could NOT find Python3 module nose2")
  message (STATUS "  Skipping py_lib tests")
  return ()
endif (PYTHON_NOSE2_NOT_FOUND)

# Adding test for dox.
message (STATUS "Found Python3 module nose2")
add_test (NAME test_py_lib_dox
          COMMAND nose2
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/seqan/dox)
