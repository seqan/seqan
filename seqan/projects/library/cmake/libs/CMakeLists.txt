cmake_minimum_required (VERSION 2.6)

################################################################################
# Include the CMakeLists.txt for each lib's directory.
################################################################################

get_filename_component (SEQAN_LIB_ABS ${SEQAN_LIBRARY}/lib ABSOLUTE)
file (GLOB SEQAN_LIBRARIES ${SEQAN_LIB_ABS}/[A-z]*)

foreach (LIBRARY_DIR ${SEQAN_LIBRARIES})
  get_filename_component(LIBRARY_NAME ${LIBRARY_DIR} NAME)
  if (EXISTS ${LIBRARY_DIR}/CMakeLists.txt)
    add_subdirectory(${LIBRARY_DIR} ${LIBRARY_NAME})
  endif (EXISTS ${LIBRARY_DIR}/CMakeLists.txt)
endforeach (LIBRARY_DIR ${SEQAN_LIBRARIES})
