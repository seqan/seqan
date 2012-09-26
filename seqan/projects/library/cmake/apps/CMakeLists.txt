cmake_minimum_required (VERSION 2.6)
project (Seqan_Apps)

################################################################################
# For the apps, use assertions in debug mode, no assertions in
# release mode and never go into testing mode.
################################################################################

# Use CMake's package location functionality to find the directories of
# TBB and Boost.
find_package (TBB QUIET)
find_package (CUDA)
find_package (OpenMP QUIET)
find_package (Boost 1.35.0)  # >= 1.35 for math/distributions, however broken?
find_package (TR1)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${SAMTOOLS_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=0")
set (CMAKE_CXX_FLAGS_RELDEBUG "${CMAKE_CXX_FLAGS_RELDEBUG} -DSEQAN_ENABLE_DEBUG=0")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")

# We need Boost >= 1.35.0.  However, the FindBoost module seems to be
# broken in this respect, so we will set Boost_FOUND to FALSE if the version
# is not large enough on our own.
if (Boost_FOUND)
  if (Boost_VERSION LESS 103500 )
    set (Boost_FOUND FALSE)
  endif (Boost_VERSION LESS 103500 )
endif (Boost_FOUND)

################################################################################
# Set CUDA variables
################################################################################

set (CUDA_PROPAGATE_HOST_FLAGS OFF)
set (CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
#set (CUDA_NVCC_FLAGS "-pg")
list (APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS "cu")
cuda_include_directories(${CUDA_CUT_INCLUDE_DIR})
#string (REGEX REPLACE "\\-(W( |all)|pedantic)" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})
string (REGEX REPLACE "\\-pedantic" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})

################################################################################
# Set path variables
################################################################################

# If Boost is not found then Boost_INCLUDE_DIRS will be set
# to "NOTFOUND" which will make CMake complain on Windows.
if (Boost_FOUND)
  set (BOOST_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
endif (Boost_FOUND)

# Make sure the compiler can find include files from the SeqAn, Boost and TBB library.
include_directories (${SEQAN_LIBRARY}
                     ${TBB_INCLUDE_DIRS}
                     ${BOOST_INCLUDE_DIRS})
link_directories (${SEQAN_LIBRARY}/lib
                  ${TBB_LIBRARY_DIRS})

# Make sure the compiler can find object files from the TBB library.
##link_directories (${SEQAN_LIBRARY}/lib)
# ^-- currently not working because cmake appends Debug|Release to ../lib/

################################################################################
# Applications
################################################################################

# TODO(holtgrew): Disable apps if dependencies are missing.  This would mean removing the automatic generation of targets but this invites errors anyway.

# pseudo-target to build all apps
add_custom_target (apps)

# Find targets automatically
get_filename_component (SEQAN_APPS_ABS ${SEQAN_LIBRARY}/apps ABSOLUTE)
file (GLOB SEQAN_APPS ${SEQAN_APPS_ABS}/[A-z]*)

file (GLOB INSEGT_HEADERS ${SEQAN_LIBRARY}/apps/insegt/[A-z]*.h)
file (GLOB RAZERS_HEADERS ${SEQAN_LIBRARY}/apps/razers/[A-z]*.h)
file (GLOB RAZERS2_HEADERS ${SEQAN_LIBRARY}/apps/razers2/[A-z]*.h)
file (GLOB RAZERS3_HEADERS ${SEQAN_LIBRARY}/apps/razers3/[A-z]*.h)
file (GLOB SEQAN_TCOFFEE_HEADERS ${SEQAN_LIBRARY}/apps/seqan_tcoffee/[A-z]*.h)

# Manually add predefined targets
list (APPEND SEQAN_APPS ${SEQAN_LIBRARY}/apps/paramChooser)
list (APPEND SEQAN_APPS ${SEQAN_LIBRARY}/apps/transrectify)
list (APPEND SEQAN_APPS ${SEQAN_LIBRARY}/apps/transsplice)
list (APPEND SEQAN_APPS ${SEQAN_LIBRARY}/apps/transform)
list (APPEND SEQAN_APPS ${SEQAN_LIBRARY}/apps/chopsuey)
set  (APP_transrectify_SOURCES ${SEQAN_LIBRARY}/apps/transquant/transrectify.cpp)
set  (APP_transsplice_SOURCES ${SEQAN_LIBRARY}/apps/transquant/transsplice.cpp)
set  (APP_transform_SOURCES ${SEQAN_LIBRARY}/apps/transquant/transform.cpp)
set  (APP_chopsuey_SOURCES ${SEQAN_LIBRARY}/apps/transquant/chopsuey.cpp)

# Manually set predefined sets of source files for some targets
set (APP_razers_SOURCES ${SEQAN_LIBRARY}/apps/razers/razers.cpp ${RAZERS_HEADERS})
set (APP_razers2_SOURCES ${SEQAN_LIBRARY}/apps/razers2/razers.cpp ${RAZERS2_HEADERS})
set (APP_micro_razers_SOURCES ${SEQAN_LIBRARY}/apps/micro_razers/micro_razers.cpp ${RAZERS_HEADERS})
set (APP_paramChooser_SOURCES ${SEQAN_LIBRARY}/apps/razers/paramChooser.cpp ${RAZERS_HEADERS})
set (APP_insegt_SOURCES ${SEQAN_LIBRARY}/apps/insegt/insegt.cpp ${INSEGT_HEADERS})
if (Boost_FOUND)
  set (APP_transquant_SOURCES ${SEQAN_LIBRARY}/apps/transquant/transquant.cpp)
else (Boost_FOUND)
  message ("WARNING Boost not found, not building transquant!")
  # Disable building of transquant.
  set (APP_transquant_NOBUILD TRUE)
  list (REMOVE_ITEM SEQAN_APPS ${SEQAN_LIBRARY}/apps/transquant)
endif (Boost_FOUND)
if (OpenMP_CXX_FLAGS)
  if (Boost_FOUND OR TR1_SHARED_PTR_FOUND)
    if (NOT MINGW)
      set (APP_razers3_SOURCES ${SEQAN_LIBRARY}/apps/razers3/razers.cpp ${RAZERS3_HEADERS})
    else (NOT MINGW)
      message ("RazerS 3 is not built on MinGW since __sync_val_compare_and_swap is missing.")
      set (APP_razers3_NOBUILD true)
      list (REMOVE_ITEM SEQAN_APPS ${SEQAN_LIBRARY}/apps/razers3)
    endif (NOT MINGW)
  else (Boost_FOUND OR TR1_SHARED_PTR_FOUND)
    # TODO(holtgrew): Is there a better way to specify complex dependencies?
    message ("WARNING neither TR1 nor Boost found, not building razers3!")
    # Disable building of RazerS 3.
    set (APP_razers3_NOBUILD true)
    list (REMOVE_ITEM SEQAN_APPS ${SEQAN_LIBRARY}/apps/razers3)
  endif (Boost_FOUND OR TR1_SHARED_PTR_FOUND)
else (OpenMP_CXX_FLAGS)
  message ("WARNING OpenMP not found, not building razers3!")
  # Disable building of RazerS 3.
  set (APP_razers3_NOBUILD true)
  list (REMOVE_ITEM SEQAN_APPS ${SEQAN_LIBRARY}/apps/razers3)
endif (OpenMP_CXX_FLAGS)

################################################################################
# Auto-generate headers.
################################################################################

foreach (APP_DIR ${SEQAN_APPS})
  if (EXISTS ${APP_DIR}/config.h.in)
    configure_file (${APP_DIR}/config.h.in
                    ${APP_DIR}/config.h)
  endif (EXISTS ${APP_DIR}/config.h.in)
endforeach (APP_DIR ${SEQAN_APPS})

################################################################################
# Build applications.
################################################################################

# Make sure the compiler can find object files from the TBB library.
link_directories (${SEQAN_LIBRARY}/lib)

foreach (APP_DIR ${SEQAN_APPS})
	file (RELATIVE_PATH APP_DIR_REL ${SEQAN_APPS_ABS} ${APP_DIR})

	# if there is no predefined set of files for this target
	# find the source files automatically
	if (NOT APP_${APP_DIR_REL}_SOURCES)
		file (GLOB APP_${APP_DIR_REL}_SOURCES ${APP_DIR}/[A-z]*.cpp ${APP_DIR}/[A-z]*.h ${APP_DIR}/[A-z]*.cu)
	endif (NOT APP_${APP_DIR_REL}_SOURCES)

	foreach (SRC ${APP_${APP_DIR_REL}_SOURCES})
		if (NOT EXISTS ${SRC})
#			message ("WARNING ${SRC} not found. Disabling ${APP_DIR_REL}.")
			set (APP_${APP_DIR_REL}_NOBUILD TRUE)
		endif (NOT EXISTS ${SRC})
	endforeach (SRC)

	# auto-detect CUDA targets
	foreach (SRC ${APP_${APP_DIR_REL}_SOURCES})
		get_filename_component(SRC_EXT ${SRC} EXT)
		if (${SRC_EXT} STREQUAL ".cu")
#			message ("Found CUDA target: ${APP_DIR_REL}")
			set (APP_${APP_DIR_REL}_CUDA TRUE)
			if (NOT CUDA_FOUND)
				set (APP_${APP_DIR_REL}_NOBUILD TRUE)
			endif (NOT CUDA_FOUND)
		endif (${SRC_EXT} STREQUAL ".cu")
	endforeach (SRC)

	# add target, dependencies and linked libraries
	set (SAVED_FLAGS ${CMAKE_CXX_FLAGS})
	if (APP_${APP_DIR_REL}_SOURCES)
		if (APP_${APP_DIR_REL}_NOBUILD)
			message ("Not building ${APP_DIR_REL} -- it is deactivated.")
		else (APP_${APP_DIR_REL}_NOBUILD)
			list (SORT APP_${APP_DIR_REL}_SOURCES)
			if (APP_${APP_DIR_REL}_CUDA)
			  message ("Found CUDA target: ${APP_DIR_REL}")
			  set (CMAKE_CXX_FLAGS ${CUDA_CXX_FLAGS})
			  cuda_add_executable(${APP_DIR_REL} ${APP_${APP_DIR_REL}_SOURCES})
			  target_link_libraries(${APP_DIR_REL} ${CUDA_CUT_LIBRARY})
			else (APP_${APP_DIR_REL}_CUDA)
			  add_executable (${APP_DIR_REL} ${APP_${APP_DIR_REL}_SOURCES})
			endif (APP_${APP_DIR_REL}_CUDA)
			add_dependencies (${APP_DIR_REL} ${SEQAN_TARGET})
			add_dependencies (apps ${APP_DIR_REL})
            # Link against libz everywhere, used for reading BAM, for example.
            if (ZLIB_FOUND)
			  target_link_libraries (${APP_DIR_REL} z)
            endif (ZLIB_FOUND)
            # Link against librt on Linux for real-time clock.
			if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
				target_link_libraries (${APP_DIR_REL} rt)
			endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
			set (CMAKE_CXX_FLAGS ${SAVED_FLAGS})
		endif (APP_${APP_DIR_REL}_NOBUILD)
	endif (APP_${APP_DIR_REL}_SOURCES)
endforeach (APP_DIR)
