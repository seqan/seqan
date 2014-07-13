# Taken from
# https://bitbucket.org/alpar/lemon-project-template/src/368136e07d23/cmake/FindLEMON.cmake

SET(LEMON_ROOT_DIR "" CACHE PATH "LEMON root directory")

FIND_PATH(LEMON_INCLUDE_DIR
  lemon/core.h
  HINTS ${LEMON_ROOT_DIR}/include
        $ENV{LEMON_ROOT_DIR}/include
        "C:/Program Files/LEMON/include"
        
)
FIND_LIBRARY(LEMON_LIBRARY
  NAMES lemon emon
  HINTS ${LEMON_ROOT_DIR}/lib
        $ENV{LEMON_ROOT_DIR}/lib
        "C:/Program Files/LEMON/lib" 
        
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LEMON DEFAULT_MSG LEMON_LIBRARY LEMON_INCLUDE_DIR)

IF(LEMON_FOUND)
  SET(LEMON_INCLUDE_DIRS ${LEMON_INCLUDE_DIR})
  SET(LEMON_LIBRARIES ${LEMON_LIBRARY})
ENDIF(LEMON_FOUND)

MARK_AS_ADVANCED(LEMON_LIBRARY LEMON_INCLUDE_DIR)
