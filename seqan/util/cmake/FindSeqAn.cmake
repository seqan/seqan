# ===========================================================================
# FindSeqAn.cmake -- Utilities for building SeqAn applications.
# ===========================================================================

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------

# Use c++11 standard for compilation
option (SEQAN_C++11_STANDARD "Use the c++11 standard for compilation." OFF)
if (SEQAN_C++11_STANDARD)
  if (APPLE)
    set (CMAKE_XCODE_ATTRIBUTE_GCC_VERSION com.apple.compilers.llvm.clang.macports CACHE STRING "" FORCE)
  endif (APPLE)
  if (UNIX)
    list (APPEND CMAKE_CXX_FLAGS -std=c++0x)
  endif (UNIX)
endif (SEQAN_C++11_STANDARD)

# ---------------------------------------------------------------------------
# Macro seqan_workshop_instrumentation_cmake ()
# ---------------------------------------------------------------------------

macro(seqan_instrumentation_cmake)
	if(SEQAN_INSTRUMENTATION)
	    message(STATUS "Prepare SeqAn Usability Analyzer data collection...")
	    if(CMAKE_HOST_WIN32 AND NOT PYTHONINTERP_FOUND)
	        execute_process(COMMAND ${CMAKE_SOURCE_DIR}/misc/seqan_instrumentation/py2exe/dist/seqan_instrumentation.exe cmake ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})
	        add_custom_target(seqan_instrumentation_build ${CMAKE_SOURCE_DIR}/misc/seqan_instrumentation/py2exe/dist/seqan_instrumentation.exe build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
	                          COMMENT "Build Instrumentation...")
	    else(CMAKE_HOST_WIN32 AND NOT PYTHONINTERP_FOUND)
	        execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/misc/seqan_instrumentation/bin/seqan_instrumentation.py cmake ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})
	        add_custom_target(seqan_instrumentation_build ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/misc/seqan_instrumentation/bin/seqan_instrumentation.py build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
	                          COMMENT "Build Instrumentation...")
	    endif(CMAKE_HOST_WIN32 AND NOT PYTHONINTERP_FOUND)
	else(SEQAN_INSTRUMENTATION)
		#message(STATUS "SeqAn Usability Analyzer data collection deactivated")
	endif(SEQAN_INSTRUMENTATION)
endmacro (seqan_instrumentation_cmake)


# ---------------------------------------------------------------------------
# Macro seqan_workshop_instrumentation_target (TARGET)
# ---------------------------------------------------------------------------

macro (seqan_instrumentation_target TARGET)
	if(SEQAN_INSTRUMENTATION)
	    add_dependencies(${TARGET} seqan_instrumentation_build)

	    if(CMAKE_HOST_WIN32 AND NOT PYTHONINTERP_FOUND)
	      add_custom_command(TARGET ${TARGET}
	                         PRE_BUILD
	                         COMMAND ${CMAKE_SOURCE_DIR}/misc/seqan_instrumentation/py2exe/dist/seqan_instrumentation.exe pre_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
	                         COMMENT "Pre Build Instrumentation...")
	      add_custom_command(TARGET ${TARGET}
	                         POST_BUILD
	                         COMMAND ${CMAKE_SOURCE_DIR}/misc/seqan_instrumentation/py2exe/dist/seqan_instrumentation.exe post_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
	                         COMMENT "Post Build Instrumentation...")
	    else(CMAKE_HOST_WIN32 AND NOT PYTHONINTERP_FOUND)
	      add_custom_command(TARGET ${TARGET}
	                         PRE_BUILD
	                         COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/misc/seqan_instrumentation/bin/seqan_instrumentation.py pre_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
	                         COMMENT "Pre Build Instrumentation...")
	      add_custom_command(TARGET ${TARGET}
	                         POST_BUILD
	                         COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/misc/seqan_instrumentation/bin/seqan_instrumentation.py post_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
	                         COMMENT "Post Build Instrumentation...")
	    endif(CMAKE_HOST_WIN32 AND NOT PYTHONINTERP_FOUND)
	endif(SEQAN_INSTRUMENTATION)
endmacro (seqan_instrumentation_target TARGET)


# ---------------------------------------------------------------------------
# Function seqan_get_version()
#
# Sets the variables SEQAN_VERSION, SEQAN_VERSION_MAJOR, SEQAN_VERSION_MINOR,
# SEQAN_VERSION_PATCH, determined from seqan/version.h
# ---------------------------------------------------------------------------

macro (seqan_get_version)
  try_run(_SEQAN_RUN_RESULT
          _SEQAN_COMPILE_RESULT
          ${CMAKE_BINARY_DIR}/CMakeFiles/SeqAnVersion
          ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/SeqAnVersion.cpp
          CMAKE_FLAGS -DINCLUDE_DIRECTORIES:STRING=${SEQAN_INCLUDE_DIR_FOR_SeqAnCore}
          COMPILE_OUTPUT_VARIABLE _COMPILE_OUTPUT
          RUN_OUTPUT_VARIABLE _RUN_OUTPUT)
  if (NOT _RUN_OUTPUT)
	message("")
	message("ERROR: Could not determine SeqAn version.")
	message("COMPILE OUTPUT:")
	message(${_COMPILE_OUTPUT})
  endif (NOT _RUN_OUTPUT)
  string(REGEX REPLACE ".*SEQAN_VERSION_MAJOR:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_MAJOR ${_RUN_OUTPUT})
  string(REGEX REPLACE ".*SEQAN_VERSION_MINOR:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_MINOR ${_RUN_OUTPUT})
  string(REGEX REPLACE ".*SEQAN_VERSION_PATCH:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_PATCH ${_RUN_OUTPUT})
  string(REGEX REPLACE ".*SEQAN_VERSION_PRE_RELEASE:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_PRE_RELEASE ${_RUN_OUTPUT})
    set(SEQAN_VERSION "${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}")
  if (SEQAN_VERSION_PRE_RELEASE STREQUAL 1)
    set(SEQAN_VERSION "pre${SEQAN_VERSION}")
  endif (SEQAN_VERSION_PRE_RELEASE STREQUAL 1)
endmacro (seqan_get_version)

# ---------------------------------------------------------------------------
# Macro seqan_setup_global ()
# ---------------------------------------------------------------------------

# Global setup for SeqAn.
#
# This consists of setting the warning level and suppressing some warnings.

macro (seqan_setup_global)
    # This is used for calling anything in util.
    set (SEQAN_ROOT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL
         "Used to get the directory to util, for example." FORCE)

    # -----------------------------------------------------------------------
    # Check whether we compile with CLANG.
    # -----------------------------------------------------------------------
    set (COMPILER_IS_CLANG FALSE)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set (COMPILER_IS_CLANG TRUE)
    endif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

    # -----------------------------------------------------------------------
    # Fix CMAKE_COMPILER_IS_GNUCXX for MinGW.
    # -----------------------------------------------------------------------

    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
      set (CMAKE_COMPILER_IS_GNUCXX TRUE)
    endif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")

    # -----------------------------------------------------------------------
    # GCC Setup
    # -----------------------------------------------------------------------
    if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
        # For the GCC, enable warnings.
        set (CMAKE_CXX_WARNING_LEVEL 4)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wno-long-long -fstrict-aliasing -Wstrict-aliasing")
        add_definitions (-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64)

        # Determine GCC version.
        # message("Determining GCC version.")
        EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
                     ARGS --version
                     OUTPUT_VARIABLE __GCC_VERSION)
        STRING(REGEX REPLACE ".*([0-9])\\.([0-9])\\.([0-9]).*" "\\1\\2\\3"
               _GCC_VERSION ${__GCC_VERSION})
        # message("  GCC version is ${_GCC_VERSION}")

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
        set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
        set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g -fno-omit-frame-pointer")
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fomit-frame-pointer")

        # Pass CXX flags to flags.
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_CXX_FLAGS_=\"${CMAKE_CXX_FLAGS}\"")
    endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

    # -----------------------------------------------------------------------
    # Windows Setup
    # -----------------------------------------------------------------------
    if (WIN32)
        # Always set NOMINMAX such that <Windows.h> does not define min/max
        # as macros.
        add_definitions(-DNOMINMAX)
    endif (WIN32)

    # -----------------------------------------------------------------------
    # Visual Studio Setup
    # -----------------------------------------------------------------------
    if (MSVC)
        # Warning level 3 for MSVC is disabled for now to see how much really bad warnings there are.
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W2")
    endif (MSVC)

    # -----------------------------------------------------------------------
    # Instrumentation
    # -----------------------------------------------------------------------
    seqan_instrumentation_cmake()
endmacro (seqan_setup_global)

# ---------------------------------------------------------------------------
# Macro seqan_setup_includes ()
# ---------------------------------------------------------------------------

# Setup an "include" directory, i.e. add it to the include directories via
# add_includes() and register forward building for everything inside a
# subdirectory "seqan".  Also registers it in variable SEQAN_LIBRARY_TARGETS.

function (seqan_setup_includes REL_PATH TARGET_NAME)
    set (PATH ${CMAKE_CURRENT_SOURCE_DIR}/${REL_PATH})
    set (PATH_BIN ${CMAKE_CURRENT_BINARY_DIR}/${REL_PATH})
    file (GLOB HEADERS ${PATH}/seqan/[A-z]*/[A-z]*.h)
    file (GLOB SUPER_HEADERS ${PATH}/seqan/[A-z]*.h)
	
    set (SEQAN_INCLUDE_DIR_FOR_${TARGET_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${REL_PATH} CACHE INTERNAL "asdf" FORCE)
    # message("SEQAN_INCLUDE_DIR_FOR_${TARGET_NAME} <- ${SEQAN_INCLUDE_DIR_FOR_${TARGET_NAME}}")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${REL_PATH})

    set (SEQAN_LIBRARY_TARGETS ${TARGET_NAME} CACHE INTERNAL "asdf" FORCE)

    # Sort headers.
    if (HEADERS)
        list (SORT HEADERS)
    endif (HEADERS)
	
    # ---------------------------------------------------------------------------
    # GUI Setup Stuff.
    # ---------------------------------------------------------------------------
    get_filename_component (BASE_ABS ${PATH}/seqan ABSOLUTE)
    
    # CMake bug workaround: For Non-GUI generators there is a bug in cmake.
    # The SOURCE command in add_custom_target is not recognized there.
    set (NONGUI_GENERATORS "Unix Makefiles" "MinGW Makefiles")
    list (FIND NONGUI_GENERATORS ${CMAKE_GENERATOR} FOUND)
    if (FOUND EQUAL -1)
        set (GUI_SOURCES SOURCES ${HEADERS} ${SUPER_HEADERS})
    endif (FOUND EQUAL -1)

    # Add SeqAn Pseudo Target for GUIs.
    #
    # This target contains all headers, forwards and umbrella headers.
    add_custom_target(
        ${TARGET_NAME}
        DEPENDS ${HEADERS}
                ${GUI_SOURCES}
    )
    
    # Group library headers into modules.  The CMake documentation says this
    # is mostly (only?) used for Visual Studio Projects.
    foreach (HEADER ${HEADERS})
        file (RELATIVE_PATH HEADER_REL ${BASE_ABS} ${HEADER})
        get_filename_component (MODULE ${HEADER_REL} PATH)
        source_group (${MODULE} FILES ${HEADER})
        # message("source_group(${MODULE} FILES ${HEADER})")    
    endforeach (HEADER ${HEADERS})

    # -----------------------------------------------------------------------
    # Installation
    # -----------------------------------------------------------------------
    foreach (HEADER ${HEADERS} ${SUPER_HEADERS})
        string(REPLACE ${CMAKE_CURRENT_BINARY_DIR}/${REL_PATH}/seqan "" NEW_PATH ${HEADER})
        string(REPLACE ${BASE_ABS} "" NEW_PATH ${NEW_PATH})
        string(REPLACE "//" "/" NEW_PATH ${NEW_PATH})
        install(FILES ${HEADER}
                RENAME seqan${NEW_PATH}
                DESTINATION include)
        #message("install(FILES ${HEADER} RENAME seqan${NEW_PATH} DESTINATION include COMPONENT dev)")
    endforeach()
endfunction (seqan_setup_includes)

# ---------------------------------------------------------------------------
# Macro seqan_make_seqan_available ()
# ---------------------------------------------------------------------------

# Register a seqan extension that has been previously defined available.
# TODO(holtgrew): Is order really important at all?

function (seqan_make_seqan_available TARGET_NAME)
   set (SEQAN_LIBRARY_TARGETS ${SEQAN_LIBRARY_TARGETS} ${TARGET_NAME} CACHE INTERNAL "asdf" FORCE)
   foreach (x ${SEQAN_INCLUDE_DIR_FOR_${TARGET_NAME}})
     # message(include_directories(${x}))
     include_directories(${x})
   endforeach (x SEQAN_INCLUDE_DIR_FOR_${TARGET_NAME})
endfunction (seqan_make_seqan_available TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_setup_tests ()
# ---------------------------------------------------------------------------

# Switch to testing mode.  This function should be called in the
# CMakeLists.txt in the tests directories before including subdirectories.
#
# The following will happen:
#
# * A target with the name TEST_TARGET will be created and dependencies
#   of this target to subsequently targets added by seqan_add_test will
#   be added.
# * Setting definitions SEQAN_ENABLE_DEBUG=1 and SEQAN_ENABLE_TESTING=1.
# * If the ${MODEL} variable is NightlyCoverage OR ExperimentalCoverage,
#   and the compiler is GCC C++ then symbols for test coverate are added.

macro (seqan_setup_tests TEST_TARGET)
    # Setup flags for tests.
    add_definitions(-DSEQAN_ENABLE_DEBUG=1)
    add_definitions(-DSEQAN_ENABLE_TESTING=1)
    
    # Add a target for the tests.
    add_custom_target(${TEST_TARGET})
    # Create a CMake variable for storing the current test target.
    set (SEQAN_CURRENT_TEST_TARGET ${TEST_TARGET} CACHE INTERNAL
         "Test target, communicated to seqan_add_test_executable" FORCE)
    # message (STATUS "TEST_TARGET " ${TEST_TARGET})
    # message (STATUS "SEQAN_CURRENT_TEST_TARGET " ${SEQAN_CURRENT_TEST_TARGET})
    
    # Conditionally enable coverage mode.
    if (MODEL STREQUAL "NightlyCoverage")
        if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
            add_definitions(-DSEQAN_ENABLE_CHECKPOINTS=0)
        endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
    endif (MODEL STREQUAL "NightlyCoverage")
    if (MODEL STREQUAL "ExperimentalCoverage")
        if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
            add_definitions(-DSEQAN_ENABLE_CHECKPOINTS=0)
        endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
    endif (MODEL STREQUAL "ExperimentalCoverage")
endmacro (seqan_setup_tests)

# ---------------------------------------------------------------------------
# Macro seqan_setup_apps ()
# ---------------------------------------------------------------------------

# Initialize "apps" area.  This function should be called in the
# CMakeLists.txt in the apps directories before including subdirectories.
#
# The following will happen:
#
# * A target with the name APP_TARGET will be created and dependencies
#   of this target to subsequently targets added by seqan_add_test will
#   be added.
# * Setup the correct flags for Debug, Release and RelWithDebInfo mode.

macro (seqan_setup_apps APP_TARGET)
    # Set flags for SeqAn.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=0")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")

    # Add a target for the tests.
    add_custom_target(${APP_TARGET})
    # Create a CMake variable for storing the current test target.
    set (SEQAN_CURRENT_APP_TARGET ${APP_TARGET} CACHE INTERNAL
         "App target, communicated to seqan_add_executable" FORCE)
endmacro (seqan_setup_apps)

# ---------------------------------------------------------------------------
# Macro seqan_setup_demos ()
# ---------------------------------------------------------------------------

# Initialize "demos" area.  This function should be called in the
# CMakeLists.txt in the demos directories before including subdirectories.
#
# The following will happen:
#
# * A target with the name DEMO_TARGET will be created and dependencies of
#   this target to subsequently targets added by seqan_add_executable will be
#   added.
# * Setup the correct flags for Debug, Release and RelWithDebInfo mode.

macro (seqan_setup_demos DEMO_TARGET)
    # Set flags for SeqAn.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=0")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")

    # Add a target for the tests.
    add_custom_target(${DEMO_TARGET})
    # Create a CMake variable for storing the current test target.
    set (SEQAN_CURRENT_APP_TARGET ${DEMO_TARGET} CACHE INTERNAL
         "App target, communicated to seqan_add_executable" FORCE)
endmacro (seqan_setup_demos)

# ---------------------------------------------------------------------------
# Macro seqan_find_dependencies ()
# ---------------------------------------------------------------------------

# Try to find all dependencies using the find() function.
#
# Currently the libraries SeqAn can use to extend its functionality are:
#
#  * zlib
#  * bzlib
#  * OpenMP
#  * CUDA

macro (seqan_find_dependencies)
  # Register external include directory.

  find_package (ZLIB QUIET)
  if (ZLIB_FOUND)
    add_definitions(-DSEQAN_HAS_ZLIB=1)
    include_directories(${ZLIB_INCLUDE_DIRS})
  endif (ZLIB_FOUND)

  find_package (BZip2 QUIET)
  if (BZIP2_FOUND)
    include_directories(${BZIP_INCLUDE_DIRS})
    add_definitions(-DSEQAN_HAS_BZIP2=1)
  endif (BZIP2_FOUND)

  # search OpenMP flags for only non-clang compilers
  if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    find_package (OpenMP QUIET)
  endif (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  find_package (CUDA QUIET)
  find_package (Boost)
  if (Boost_FOUND)
    #    include_directories(${Boost_INCLUDE_DIRS})
  endif (Boost_FOUND)

  if (CUDA_FOUND)
    add_definitions(-DSEQAN_HAS_CUDA=1)
  endif (CUDA_FOUND)

  include(CheckIncludeFiles)
  check_include_files(execinfo.h HAVE_EXECINFO)
  if (HAVE_EXECINFO)
    add_definitions(-DSEQAN_HAS_EXECINFO=1)
  endif (HAVE_EXECINFO)
endmacro (seqan_find_dependencies)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_subdirectories ()
# ---------------------------------------------------------------------------

# This macro calls add_subdirectory() for all subdirectories below
# ${SEQAN_CURRENT_SOURCE_DIR} if they contain a CMakeLists.txt.
#
# Example:
#
#   seqan_add_all_subdirectories(seqan_tests)

macro (seqan_add_all_subdirectories)
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)
    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            # message(STATUS "Going into ${ENTRY}")
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                add_subdirectory(${ENTRY})
            endif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
        endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
    endforeach (ENTRY ${ENTRIES})
endmacro (seqan_add_all_subdirectories)

# ---------------------------------------------------------------------------
# Macro seqan_project_just_include_all_subdirs ()
# ---------------------------------------------------------------------------

# Call this macro as the only command in a CMakeLists.txt file to include
# all subdirectories.
#
# Args:
#   PROJECT_NAME  String, the project name of this subdirectory.
#
# Example:
#
#   # Only line in file tests/CMakeLists.txt
#   seqan_project_just_include_all_subdirs(seqan_tests)

macro (seqan_project_just_include_all_subdirs PROJECT_NAME)
    cmake_minimum_required (VERSION 2.6)
    project(${PROJECT_NAME})
    seqan_add_all_subdirectories()
endmacro (seqan_project_just_include_all_subdirs)

# ---------------------------------------------------------------------------
# Macro seqan_add_executable (TARGET_NAME source1.cpp source2.[cpp|h] ...)
#       seqan_add_executable (TARGET_NAME)
# ---------------------------------------------------------------------------

# Create a SeqAn executable from the given source files.  If no such files
# are given then all files in the current directory will be used.

macro (seqan_add_executable TARGET_NAME)
    # TODO(holtgrew): Use all files in directory if ${ARGC} == 0
    # message(STATUS "add_executable (${ARGV})")
    add_executable (${ARGV})

    add_dependencies(${SEQAN_CURRENT_APP_TARGET} ${ARGV0})
    
    # Link against librt on Linux.
  	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    	target_link_libraries (${TARGET_NAME} rt)
  	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  	
  	# Link against stdc++ on Mac OS X (clang seems not to do this automatically)
  	if (APPLE)
  		target_link_libraries (${TARGET_NAME} stdc++)
  	endif (APPLE)
  	
  	# Dependencies on all registered seqan extensions.
  	add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})
    # message(STATUS "add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})")

    # -----------------------------------------------------------------------
    # Instrumentation
    # -----------------------------------------------------------------------
    seqan_instrumentation_target(${ARGV0})

  	# Link against zlib and bzlib if found.
  	if (ZLIB_FOUND)
  	    include_directories (${ZLIB_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${ZLIB_LIBRARIES})
  	endif (ZLIB_FOUND)
  	if (BZIP2_FOUND)
  	    include_directories (${BZIP2_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${BZIP2_LIBRARIES})
  	endif (BZIP2_FOUND)
endmacro (seqan_add_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_add_test_executable (TARGET_NAME source1.cpp source2.[cpp|h])
#       seqan_add_test_executable (TARGET_NAME)
# ---------------------------------------------------------------------------

# Create a SeqAn executable from the given source files.  If no such files
# are given then all files in the current directory will be used.
#
# Also the test will be registered with add_test and depend on the current
# test target, i.e. from closest seqan_setup_tests(TEST_TARGET_NAME) call.

macro (seqan_add_test_executable TARGET_NAME)
    # TODO(holtgrew): Use all files in directory if ${ARGC} == 0
    # message(STATUS "add_executable (${ARGV})")
    add_executable (${ARGV})
    add_test(test_${ARGV0} ${ARGV0})
                     
    add_dependencies(${SEQAN_CURRENT_TEST_TARGET} ${ARGV0})
    
    # Link against librt on Linux.
  	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  		  target_link_libraries (${ARGV0} rt)
  	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

  	# Dependencies on all registered seqan extensions.
    # message("add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})")
  	add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})

    # -----------------------------------------------------------------------
    # Instrumentation
    # -----------------------------------------------------------------------
    seqan_instrumentation_target(${ARGV0})

  	# Link against zlib and bzlib if found.
  	if (ZLIB_FOUND)
  	    include_directories (${ZLIB_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${ZLIB_LIBRARIES})
  	endif (ZLIB_FOUND)
  	if (BZIP2_FOUND)
  	    include_directories (${BZIP2_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${BZIP2_LIBRARIES})
  	endif (BZIP2_FOUND)
endmacro (seqan_add_test_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Function seqan_add_cuda_executable (TARGET_NAME source1.cu source2.[cu|cpp|h])
#          seqan_add_cuda_executable (TARGET_NAME)
# ---------------------------------------------------------------------------

function (seqan_add_cuda_executable TARGET_NAME)
    if (CUDA_FOUND)
        # TODO(holtgrew): Use all files in directory if ${ARGC} == 0

        # -------------------------------------------------------------------
        # Set CUDA variables
        # -------------------------------------------------------------------
        set (CUDA_PROPAGATE_HOST_FLAGS OFF)
        set (CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
        #set (CUDA_NVCC_FLAGS "-pg")
        list (APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS "cu")
        cuda_include_directories(${CUDA_CUT_INCLUDE_DIR})
        #string (REGEX REPLACE "\\-(W( |all)|pedantic)" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})
        string (REGEX REPLACE "\\-pedantic" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})

        # -------------------------------------------------------------------
        # Go on normally as in apps.
        # -------------------------------------------------------------------
        # message(STATUS "cuda_add_executable (${ARGV})")
        cuda_add_executable (${ARGV})

        add_dependencies(${SEQAN_CURRENT_APP_TARGET} ${ARGV0})
    
        # Link against librt on Linux.
      	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        	target_link_libraries (${TARGET_NAME} rt)
      	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  	
      	# Dependencies on all registered seqan extensions.
      	add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})

        # -----------------------------------------------------------------------
        # Instrumentation
        # -----------------------------------------------------------------------
	seqan_instrumentation_target(${ARGV0})
	
      	# Link against zlib and bzlib if found.
      	if (ZLIB_FOUND)
      	    include_directories (${ZLIB_INCLUDE_DIR})
      	    target_link_libraries (${TARGET_NAME} ${ZLIB_LIBRARIES})
      	endif (ZLIB_FOUND)
      	if (BZIP2_FOUND)
      	    include_directories (${BZIP2_INCLUDE_DIR})
      	    target_link_libraries (${TARGET_NAME} ${BZIP2_LIBRARIES})
      	endif (BZIP2_FOUND)
    endif (CUDA_FOUND)
endfunction (seqan_add_cuda_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Function seqan_add_cuda_test_executable (TARGET_NAME source1.cu source2.[cu|cpp|h])
#          seqan_add_cuda_test_executable (TARGET_NAME)
# ---------------------------------------------------------------------------

# Create a SeqAn CUDA executable from the given source files.  If no such
# files are given then all files in the current directory will be used.
#
# Also the test will be registered with add_test and depend on the current
# test target, i.e. from closest seqan_setup_tests(TEST_TARGET_NAME) call.
#
# Only adds target if CUDA_FOUND is set to true.

function (seqan_add_cuda_test_executable TARGET_NAME)
    if (CUDA_FOUND)
        # TODO(holtgrew): Use all files in directory if ${ARGC} == 0

        # -------------------------------------------------------------------
        # Set CUDA variables
        # -------------------------------------------------------------------
        set (CUDA_PROPAGATE_HOST_FLAGS OFF)
        set (CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
        #set (CUDA_NVCC_FLAGS "-pg")
        list (APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS "cu")
        cuda_include_directories(${CUDA_CUT_INCLUDE_DIR})
        #string (REGEX REPLACE "\\-(W( |all)|pedantic)" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})
        string (REGEX REPLACE "\\-pedantic" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})

        # -------------------------------------------------------------------
        # Go on normally as in apps.
        # -------------------------------------------------------------------
        # message(STATUS "cuda_add_executable (${ARGV})")
        cuda_add_executable (${ARGV})
        add_test(test_${ARGV0} ${ARGV0})

        add_dependencies(${SEQAN_CURRENT_TEST_TARGET} ${ARGV0})
    
        # Link against librt on Linux.
      	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        	target_link_libraries (${TARGET_NAME} rt)
      	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  	
      	# Dependencies on all registered seqan extensions.
      	add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})

        # -----------------------------------------------------------------------
        # Instrumentation
        # -----------------------------------------------------------------------
	seqan_instrumentation_target(${ARGV0})
	
      	# Link against zlib and bzlib if found.
      	if (ZLIB_FOUND)
      	    include_directories (${ZLIB_INCLUDE_DIR})
      	    target_link_libraries (${TARGET_NAME} ${ZLIB_LIBRARIES})
      	endif (ZLIB_FOUND)
      	if (BZIP2_FOUND)
      	    include_directories (${BZIP2_INCLUDE_DIR})
      	    target_link_libraries (${TARGET_NAME} ${BZIP2_LIBRARIES})
      	endif (BZIP2_FOUND)
    endif (CUDA_FOUND)
endfunction (seqan_add_cuda_test_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_executables ([TARGET])
# ---------------------------------------------------------------------------

# This macro calls seqan_add_executable() for all all .cpp files in the
# current directory.  If the optional TARGET parameter is given, a target
# with this name is created and the executable targets will depend on this
# target.
#
# Example:
#
#   seqan_add_all_executables()
#   seqan_add_all_executables(depend_on_this)
#   seqan_add_all_executables(depend_on_this prefix)

macro (seqan_add_all_executables)
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cpp)
    if (${ARGC} GREATER 1)
        set (PREFIX ${ARGV1})
    endif (${ARGC} GREATER 1)
    if (${ARGC} GREATER 0)
        if (TARGET ${ARGV0})  # Add target only if it does not exist yet.
        else (TARGET ${ARGV0})
            add_custom_target(${ARGV0})
        endif (TARGET ${ARGV0})
    endif (${ARGC} GREATER 0)
    foreach (ENTRY ${ENTRIES})
        get_filename_component(BIN_NAME ${ENTRY} NAME_WE)
        seqan_add_executable(${PREFIX}${BIN_NAME} ${ENTRY})
        if (${ARGC} GREATER 0)
            add_dependencies(${ARGV0} ${PREFIX}${BIN_NAME})
        endif (${ARGC} GREATER 0)
    endforeach (ENTRY ${ENTRIES})
endmacro (seqan_add_all_executables)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_cuda_executables ([TARGET])
# ---------------------------------------------------------------------------

# This macro calls seqan_add_cuda_executable() for all all .cu files in the
# current directory.  If the optional TARGET parameter is given, a target
# with this name is created and the executable targets will depend on this
# target.
#
# Only adds target if CUDA_FOUND is set to true.
#
# Example:
#
#   seqan_add_all_cuda_executables()
#   seqan_add_all_cuda_executables(depend_on_this)

macro (seqan_add_all_cuda_executables)
    if (CUDA_FOUND)
        file (GLOB ENTRIES
              RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
              ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cu)
        if (${ARGC} GREATER 0)
            if (TARGET ${ARGV0})  # Add target only if it does not exist yet.
            else (TARGET ${ARGV0})
                add_custom_target(${ARGV0})
            endif (TARGET ${ARGV0})
        endif (${ARGC} GREATER 0)
        foreach (ENTRY ${ENTRIES})
            get_filename_component(BIN_NAME ${ENTRY} NAME_WE)
            # message("seqan_add_cuda_executable(${BIN_NAME} ${ENTRY})")
            seqan_add_cuda_executable(${BIN_NAME} ${ENTRY})
            if (${ARGC} GREATER 0)
                add_dependencies(${ARGV0} ${BIN_NAME})
            endif (${ARGC} GREATER 0)
        endforeach (ENTRY ${ENTRIES})
    endif (CUDA_FOUND)
endmacro (seqan_add_all_cuda_executables)
