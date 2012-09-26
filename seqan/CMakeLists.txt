project (seqan)

# ===========================================================================
# We require CMake 2.8.2.  Since it is available in the currently stable
# Debian, it should be available everywhere else, too.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)

# ===========================================================================
# Iff the file INST or SUA exists in the source root
# turn on SeqAn instrumentation for experimental data collection.
# ===========================================================================

if(EXISTS ${CMAKE_SOURCE_DIR}/INST OR EXISTS ${CMAKE_SOURCE_DIR}/SUA OR EXISTS ${CMAKE_SOURCE_DIR}/../SUA)
  set(SEQAN_INSTRUMENTATION YES)
else(EXISTS ${CMAKE_SOURCE_DIR}/INST OR EXISTS ${CMAKE_SOURCE_DIR}/SUA OR EXISTS ${CMAKE_SOURCE_DIR}/../SUA)
  set(SEQAN_INSTRUMENTATION NO)
endif(EXISTS ${CMAKE_SOURCE_DIR}/INST OR EXISTS ${CMAKE_SOURCE_DIR}/SUA OR EXISTS ${CMAKE_SOURCE_DIR}/../SUA)


# ===========================================================================
# Warn People About the Perils Of In-Source Builds
# ===========================================================================

if (EXISTS "$CMAKE_BINARY_DIR/CMakeLists.txt")
  if (EXISTS "$CMAKE_BINARY_DIR/I_WANT_IN_SOURCE_BUILDS")
    set(IN_SOURCE_BUILD YES)
  else (EXISTS "$CMAKE_BINARY_DIR/I_WANT_IN_SOURCE_BUILDS")
    message("CAUTION: In-source builds are not allowed by default!")
    message("")
    message("(If you want to use Eclipse and have read the 'Getting Started,")
    message("with Eclipse Tutorial' then see the instructions below.)")
    message("")
    message("Instead of doing in source builds, go to into the directory build")
    message("and call CMake from there.")
    message("We recommend you to create an additional directory level, this")
    message("allows you to have multiple build configurations (e.g. Debug,")
    message("Release) or project types.")
    message("")
    message("  cd build")
    message("  mkdir Debug")
    message("  cd Debug")
    message("  cmake ../.. -DCMAKE_BUILD_TYPE=Debug")
    message("")
    message("NOTE: CMake has created some files in this directory that you have")
    message("to remove before executing the above commands.")
    message("")
    message("  [Linux, Mac Os X]  rm -rf CMakeCache.txt CMakeFiles")
    message("  [Windows]          del CMakeCache.txt ; rmdir /S CMakeFiles")
    message("")
    message("If you want to use Eclipse then you might want to have in-source")
    message("builds because of limitations of Eclipse CDT4.  In this case")
    message("Create a file in this directory called I_WANT_IN_SOURCE_BUILDS:")
    message("")
    message("  [Linux, Mac Os X]  touch I_WANT_IN_SOURCE_BUILDS")
    message("  [Windows]          type nul >> I_WANT_IN_SOURCE_BUILDS")
    message("")
    message(FATAL_ERROR "Please follow the instructions above.")
  endif (EXISTS "$CMAKE_BINARY_DIR/I_WANT_IN_SOURCE_BUILDS")
endif (EXISTS "$CMAKE_BINARY_DIR/CMakeLists.txt")

# ===========================================================================
# Add contrib to CMAKE_FIND_ROOT_PATH.
# ===========================================================================

if (WIN32)
    # The contrib version we look for.
	set (_SEQAN_CONTRIB_VERSION D20111031)
	set (_SEQAN_CONTRIB_DIR "seqan-contrib-${_SEQAN_CONTRIB_VERSION}")

	# Try to figure out where the user installed the contrib.  We expect
	# it to be either in C:\, or one of the Program Files dirs.
	#
	# First, look into Program Files on 64 bit.
	if (DEFINED ENV{ProgramW6432})
		if (IS_DIRECTORY "$ENV{ProgramW6432}/${_SEQAN_CONTRIB_DIR}-x64")
			set (SEQAN_CONTRIB_BASE "$ENV{ProgramW6432}/${_SEQAN_CONTRIB_DIR}-x64")
		endif (IS_DIRECTORY "$ENV{ProgramW6432}/${_SEQAN_CONTRIB_DIR}-x64")
	endif (DEFINED ENV{ProgramW6432})
	# Try out Program Files for 32bit Windows.
	if (NOT DEFINED SEQAN_CONTRIB_BASE)
		if (IS_DIRECTORY "$ENV{ProgramFiles}/${_SEQAN_CONTRIB_DIR}-x86")
			set (SEQAN_CONTRIB_BASE "$ENV{ProgramFiles}/${_SEQAN_CONTRIB_DIR}-x86")
		endif (IS_DIRECTORY "$ENV{ProgramFiles}/${_SEQAN_CONTRIB_DIR}-x86")
	endif (NOT DEFINED SEQAN_CONTRIB_BASE)
	# Try out on C:/.
	if (NOT DEFINED SEQAN_CONTRIB_BASE)
		if (IS_DIRECTORY "C:/${_SEQAN_CONTRIB_DIR}-x86")
			set (SEQAN_CONTRIB_BASE "C:/${_SEQAN_CONTRIB_DIR}-x86")
		elseif (IS_DIRECTORY "C:/${_SEQAN_CONTRIB_DIR}-x64")
			set (SEQAN_CONTRIB_BASE "C:/${_SEQAN_CONTRIB_DIR}-x64")
		endif (IS_DIRECTORY "C:/${_SEQAN_CONTRIB_DIR}-x86")
	endif (NOT DEFINED SEQAN_CONTRIB_BASE)

	# Debug help.
    #if (NOT DEFINED SEQAN_CONTRIB_BASE)
    #	message("SEQAN_CONTRIB_BASE is undefined!")
    #else (NOT DEFINED SEQAN_CONTRIB_BASE)
    #	message("SEQAN_CONTRIB_BASE is ${SEQAN_CONTRIB_BASE}")
    #endif (NOT DEFINED SEQAN_CONTRIB_BASE)

	# Try to figure out the generator.
	if (IS_DIRECTORY ${SEQAN_CONTRIB_BASE})
		if (CMAKE_GENERATOR MATCHES "^Visual Studio .*")
			string (REGEX REPLACE "^Visual Studio ([0-9]+).*$" "\\1" SEQAN_CONTRIB_VARIANT ${CMAKE_GENERATOR})
			set (SEQAN_CONTRIB_VARIANT "vs${SEQAN_CONTRIB_VARIANT}")
		elseif (CMAKE_GENERATOR STREQUAL "MinGW Makefiles")
			set (SEQAN_CONTRIB_VARIANT mingw)
		endif (CMAKE_GENERATOR MATCHES "^Visual Studio .*")

        #message(STATUS "SEQAN_CONTRIB_BASE    is ${SEQAN_CONTRIB_BASE}")
        #message(STATUS "SEQAN_CONTRIB_VARIANT is ${SEQAN_CONTRIB_VARIANT}")

		# Compose contrib path.
		set(SEQAN_CONTRIB_PATH "${SEQAN_CONTRIB_BASE}/${SEQAN_CONTRIB_VARIANT}")

		# Extend CMAKE_FIND_ROOT_PATH.
		if (IS_DIRECTORY ${SEQAN_CONTRIB_PATH})
			set (CMAKE_FIND_ROOT_PATH ${SEQAN_CONTRIB_PATH} ${CMAKE_FIND_ROOT_PATH})
		endif (IS_DIRECTORY ${SEQAN_CONTRIB_PATH})
	endif (IS_DIRECTORY ${SEQAN_CONTRIB_BASE})

    #message(STATUS "CMAKE_FIND_ROOT_PATH is \"${CMAKE_FIND_ROOT_PATH}\".")
endif (WIN32)

# ===========================================================================
# Setup Variables
# ===========================================================================

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Debug CACHE STRING
         "Choose the type of build, options are: Debug Release RelDebug ..."
         FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#message (STATUS "CMAKE_BUILD_TYPE == ${CMAKE_BUILD_TYPE}")

# ===========================================================================
# Setup Modules and Find Packages.
# ===========================================================================

# Note that the find modules in util/make take precedence
# over the built-in ones from CMake.  This is used in
# FindBZip2, for example which we slightly modify.
set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/FindTBB
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/cuda
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake
    ${CMAKE_MODULE_PATH})

#message (STATUS "CMAKE_MODULE_PATH is ${CMAKE_MODULE_PATH}")

include(CTest)

find_package (PythonInterp)
find_package (SeqAn)

# ===========================================================================
# Prerequisites
# ===========================================================================

seqan_find_dependencies()  # Find SeqAn optional dependencies.
seqan_setup_global()

# ===========================================================================
# Include Subdirectories
# ===========================================================================

add_subdirectory(core)
add_subdirectory(docs)
add_subdirectory(extras)
add_subdirectory(sandbox)

# ===========================================================================
# Include CPack
# ===========================================================================

include(${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/package.cmake)

# ===========================================================================
# Check Optional Libraries, Print Usage
# ===========================================================================

include(FeatureSummary)
