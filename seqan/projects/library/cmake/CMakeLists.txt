cmake_minimum_required (VERSION 2.6)
project (Seqan)

################################################################################
# Load Package Finder Modules
################################################################################

# Add FindTBB directory to CMake's module path
set(CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindTBB/"
    "${CMAKE_CURRENT_SOURCE_DIR}/CMake/cuda"
    "${CMAKE_CURRENT_SOURCE_DIR}/CMake/"
    ${CMAKE_MODULE_PATH})
find_package (PythonInterp)

# Locate BZIP2.
find_package(BZip2)
if (BZIP2_FOUND)
  include_directories(${BZIP_INCLUDE_DIRS})
  add_definitions(-DSEQAN_HAS_BZIP2=1)
endif (BZIP2_FOUND)

# Hacky conditional inclusion of samtools.
find_package(ZLIB QUIET)
if (ZLIB_FOUND)
  add_definitions(-DSEQAN_HAS_ZLIB=1)
  set(C_CXX_FLAGS "-DSEQAN_HAS_SAMTOOLS=1")
  set(SAMTOOLS_CXX_FLAGS "-DSEQAN_HAS_SAMTOOLS=1")
  set(SAMTOOLS_LIBRARIES "bam")
  include_directories(${ZLIB_INCLUDE_DIRS})
else (ZLIB_FOUND)
  message("WARNING: zlib not found!")
  set(SAMTOOLS_FOUND "0")
endif (ZLIB_FOUND)

include_directories(${CMAKE_INCLUDE_PATH})

################################################################################
# Set Path Variables
################################################################################

# Set SEQAN_LIBRARY if it has not been set yet (e.g. by using ccmake).
if (NOT SEQAN_LIBRARY)
  get_filename_component (SEQAN_LIBRARY_ABS .. ABSOLUTE)
  set (SEQAN_LIBRARY ${SEQAN_LIBRARY_ABS} CACHE PATH
       "Path to SeqAn's parent folder (contains apps,demos,seqan)")
endif (NOT SEQAN_LIBRARY)

# Path to and not including the "projects" directory.
get_filename_component (SEQAN_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../.. ABSOLUTE)
get_filename_component (SEQAN_BINARY_ROOT ${CMAKE_CURRENT_BINARY_DIR}/../../.. ABSOLUTE)

set (SEQAN_TARGET "Seqan")

################################################################################
# Include paths for vendor copies of libraries.
################################################################################

INCLUDE_DIRECTORIES(${SEQAN_LIBRARY}/lib/samtools)

################################################################################
# Try to find execinfo.h
################################################################################

include(CheckIncludeFiles)
check_include_files(execinfo.h HAVE_EXECINFO)
if (HAVE_EXECINFO)
  add_definitions(-DSEQAN_HAS_EXECINFO=1)
endif (HAVE_EXECINFO)

################################################################################
# Define SeqAn Library Related Variables
#
# We define some variables that contain the base path of the library,
# a list of all headers in SeqAn, the modules and generated forwards.
################################################################################

# Absolute path to "seqan" directory goes into SEQAN_BASE_ABS.
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
add_custom_target (${SEQAN_TARGET}
  DEPENDS ${SEQAN_HEADERS}
          ${SEQAN_FORWARDS}
          ${SEQAN_GUI_SOURCES}
)

# Group library headers into modules.
foreach (SEQAN_HEADER ${SEQAN_HEADERS})
  file (RELATIVE_PATH SEQAN_HEADER_REL ${SEQAN_BASE_ABS} ${SEQAN_HEADER})
  get_filename_component (SEQAN_MODULE ${SEQAN_HEADER_REL} PATH)
  source_group (${SEQAN_MODULE} FILES ${SEQAN_HEADER})
endforeach (SEQAN_HEADER)

################################################################################
# Set Compiler Flags
################################################################################

# Make sure the compiler can find include files from the SeqAn library.
include_directories (${SEQAN_LIBRARY})

# We enable maximal warnings for the GCC.
if (CMAKE_COMPILER_IS_GNUCXX)
  # For the GCC, enable warnings.
  set (CMAKE_CXX_WARNING_LEVEL 4)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wno-long-long -fstrict-aliasing -Wstrict-aliasing")

  # Determine GCC version.
  message("Determining GCC version.")
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
               ARGS --version
               OUTPUT_VARIABLE _GCC_VERSION)
  STRING(REGEX REPLACE ".* ([0-9])\\.([0-9])\\.([0-9]) .*" "\\1\\2\\3"
         _GCC_VERSION ${_GCC_VERSION})
  message("  GCC version is ${_GCC_VERSION}")

  # Add -Wno-longlong if the GCC version is < 4.0.0.  Add -pedantic flag
  # but disable warnings for variadic macros with GCC >= 4.0.0.  Earlier
  # versions warn because of anonymous variadic macros in pedantic mode
  # but do not have a flag to disable these warnings.
  if (400 GREATER _GCC_VERSION)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
  else (400 GREATER _GCC_VERSION)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wno-variadic-macros")
  endif (400 GREATER _GCC_VERSION)

  # Force GCC to keep the frame pointer when debugging is enabled.
  # This is mainly important for 64 bit but does not get into the way
  # on 32 bit either at minimal performance impact.
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELDEBUG "${CMAKE_CXX_FLAGS_RELEASE} -g")# -fno-omit-frame-pointer")
  
  # Pass CXX flags to flags.
  #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_CXX_FLAGS_=\"${CMAKE_CXX_FLAGS}\"")
endif (CMAKE_COMPILER_IS_GNUCXX)

if (MSVC)
  # Warning level 3 for MSVC is disabled for now to see how much really bad warnings there are.
  #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3  /wd4996 -D_CRT_SECURE_NO_WARNINGS")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W2  /wd4996 -D_CRT_SECURE_NO_WARNINGS")
  # Disable Microsoft C++ language extensions.
  # TODO(holtgrew): Re-disable again, Microsoft's header do not compile with this option! Look for a workaround.
  #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Za")
endif (MSVC)

################################################################################
# Set Default Build Type
################################################################################

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

################################################################################
# Define Convenience Macros
################################################################################

macro(SEQAN_ADD_EXECUTABLE seqan_target)

	# create a new target
	add_executable(${seqan_target} ${ARGN})

	# add dependency to SeqAn library
	add_dependencies (${seqan_target} ${SEQAN_TARGET})

	# optionally link against runtime-library
	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
		target_link_libraries (${seqan_target} rt)
	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

endmacro(SEQAN_ADD_EXECUTABLE seqan_target)

################################################################################
# Include Sub-Projects
################################################################################

add_subdirectory (libs)
add_subdirectory (apps)
add_subdirectory (demos)
add_subdirectory (tests)
add_subdirectory (benchmarks)
