# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /util/py_lib/CMakeLists.txt
#
# CMakeLists.txt file for Python stuff.
# ===========================================================================
# Kicks off tests using Python nosetests.
# ===========================================================================

find_package (PythonInterp)

if (PYTHONINTERP_FOUND)
  # Adding test for dox.
  message (STATUS "  nosetests for seqan.dox")
  add_test (NAME test_py_lib_dox
            COMMAND nosetests
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/seqan/dox)
endif (PYTHONINTERP_FOUND)
