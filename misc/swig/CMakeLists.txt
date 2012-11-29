cmake_minimum_required (VERSION 2.6)
project (SwigBinding)
#################################################################################
## get svn version
#################################################################################

find_package (Subversion)
if (Subversion_FOUND)
  Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
  set (SUBVERSION_REVISION ${ER_WC_REVISION})
endif (Subversion_FOUND)


#################################################################################
## Set Path Variables
#################################################################################
SET (SEQAN_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../projects/library")

if (NOT SEQAN_LIBRARY)
  get_filename_component (SEQAN_LIBRARY_ABS ../../projects/library ABSOLUTE)
  set (SEQAN_LIBRARY ${SEQAN_LIBRARY_ABS} CACHE PATH
       "Path to SeqAn's parent folder (contains apps,demos,seqan)")
endif (NOT SEQAN_LIBRARY)


set (SEQAN_TARGET "Seqan")
################################################################################
# Define SeqAn Library Related Variables
#
# We define some variables that contain the base path of the library,
# a list of all headers in SeqAn, the modules and generated forwards.
################################################################################
get_filename_component (SEQAN_BASE_ABS ${SEQAN_LIBRARY}/seqan ABSOLUTE)

# A sorted list of all SeqAn headers goes into SEQAN_HEADERS.  We
# exclude the generated headers.
file (GLOB SEQAN_HEADERS_TMP ${SEQAN_BASE_ABS}/[A-z]*/[A-z]*.h)
file (GLOB SEQAN_SUPER_HEADERS ${SEQAN_BASE_ABS}/[A-z]*.h)
foreach (HEADER ${SEQAN_HEADERS_TMP})
  if (NOT ${HEADER} MATCHES ".*generated.*")
    list (APPEND SEQAN_HEADERS ${HEADER})
  endif (NOT ${HEADER} MATCHES ".*generated.*")
endforeach (HEADER ${SEQAN_HEADERS})
list (SORT SEQAN_HEADERS)

# A sorted list of all SeqAn modules goes into SEQAN_MODULES.
file (GLOB SEQAN_BASE_CONTENTS
      RELATIVE ${SEQAN_BASE_ABS}
      ${SEQAN_BASE_ABS}/[A-z]*)
foreach (ENTRY ${SEQAN_BASE_CONTENTS})
  if (IS_DIRECTORY ${SEQAN_BASE_ABS}/${ENTRY})
    list (APPEND SEQAN_MODULES ${ENTRY})
  endif (IS_DIRECTORY ${SEQAN_BASE_ABS}/${ENTRY})
endforeach (ENTRY ${SEQAN_BASE_CONTENTS})
list (SORT SEQAN_MODULES)
list (REMOVE_DUPLICATES SEQAN_MODULES)
set (SEQAN_TARGET "SeqAn")
################################################################################
# CMake Bug Workarounds
################################################################################

# For Non-GUI generators there is a bug in cmake: The SOURCE command in
# add_custom_target is not recognized there
set (NONGUI_GENERATORS "Unix Makefiles" "MinGW Makefiles")
list (FIND NONGUI_GENERATORS ${CMAKE_GENERATOR} FOUND)
if (FOUND EQUAL -1)
  set (SEQAN_GUI_SOURCES SOURCES ${SEQAN_HEADERS} ${SEQAN_SUPER_HEADERS})
endif (FOUND EQUAL -1)

################################################################################
# Forwards Generation For g++ Version >= 4.0
################################################################################
if (CMAKE_COMPILER_IS_GNUCXX)
  # Build a list of generated forwards headers.  Goes into SEQAN_FORWARDS.
  foreach (MODULE ${SEQAN_MODULES})
    list (APPEND SEQAN_FORWARDS
          ${SEQAN_BASE_ABS}/${MODULE}/${MODULE}_generated_forwards.h)
  endforeach (MODULE ${SEQAN_MODULES})

  # Now tell CMake that the forward headers can be generated with
  # build_forwards.py
  add_custom_command (
    OUTPUT ${SEQAN_FORWARDS}
    COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_LIBRARY}/../../misc/build_forwards.py
            ${SEQAN_BASE_ABS} all
    DEPENDS ${SEQAN_HEADERS})
endif (CMAKE_COMPILER_IS_GNUCXX)

################################################################################
# SeqAn Pseudo Target
#
# We create a target ${SEQAN_TARGET} that contains all SeqAn headers,
# sources and the generated forwards.
################################################################################

# Note that SEQAN_FORWARDS is not defined when not compiling with g++
# which implies no forwards generation for non-g++ compilers.
add_custom_target (${SEQAN_TARGET} ALL
  DEPENDS ${SEQAN_HEADERS}
          ${SEQAN_FORWARDS}
          ${SEQAN_GUI_SOURCES}
)
################################################################################
# SWIG Build
#
# This is an example how to compile the interface swig file with cmake
# 
################################################################################

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
#find Python
FIND_PACKAGE(PythonLibs)

################################################################################
# Set Path Variables
################################################################################
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${SEQAN_LIBRARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

################################################################################
# Set SWIG Flags -> use "" for the default behaivour. This depends on the 
# target language 
################################################################################

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(alignment.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(pattern_matching.i PROPERTIES CPLUSPLUS ON)
################################################################################
# indicate the files needed to be wrapped: target_name  target_language interface_file wrapper_file
################################################################################

SWIG_ADD_MODULE(alignment python alignment.i alignment.h)
SWIG_LINK_LIBRARIES(alignment ${PYTHON_LIBRARIES})
IF (WIN32)
  SET_TARGET_PROPERTIES ( _alignment PROPERTIES SUFFIX ".pyd")
ENDIF (WIN32)

INSTALL ( TARGETS _alignment LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

INSTALL ( FILES ${CMAKE_BINARY_DIR}/alignment.py DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

INSTALL(FILES example.py  DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin )

#SWIG_ADD_MODULE(pattern_matching python pattern_matching.i pattern_matching.h)
#SWIG_LINK_LIBRARIES(pattern_matching ${PYTHON_LIBRARIES})

#IF (WIN32)
##  SET_TARGET_PROPERTIES ( _pattern_matching PROPERTIES SUFFIX ".pyd")
#ENDIF (WIN32)

SET(CPACK_SET_DESTDIR "ON")
IF (WIN32)
    SET(CMAKE_INSTALL_PREFIX .)
ENDIF (WIN32)

SET(CPACK_GENERATOR "ZIP")
set (CPACK_PACKAGE_VERSION_PATCH "${SUBVERSION_REVISION}")
IF (WIN32)
    SET(CMAKE_INSTALL_PREFIX .)
    SET(CPACK_GENERATOR "ZIP")
ENDIF (WIN32)

include (InstallRequiredSystemLibraries)
INCLUDE (CPack) 
################################################################################
# to generate the package call  cpack -C CPackConfig.cmake 
################################################################################


