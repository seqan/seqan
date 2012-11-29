# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/CMakeLists.txt
#
# CMakeLists.txt file for the sandbox folder.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox)

# Add all subdirectories with a CMakeLists.txt in them.

file (GLOB ENTRIES
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)
foreach (ENTRY ${ENTRIES})
  if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
    if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
      add_subdirectory (${ENTRY})
    endif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
  endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
endforeach (ENTRY ${ENTRIES})
