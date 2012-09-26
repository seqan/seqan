cmake_minimum_required (VERSION 2.6)
project (Seqan_Demos)

################################################################################
# For the demos, use assertions in debug mode, no assertions in
# release mode and never go into testing mode.
################################################################################

set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=0")

################################################################################
# Set path variables
################################################################################

include_directories (${SEQAN_LIBRARY})
link_directories (${SEQAN_LIBRARY}/lib)

################################################################################
# Demos
################################################################################

add_custom_target(demos)

file(GLOB SEQAN_DEMOS ${SEQAN_LIBRARY}/demos/[A-z]*.cpp)
foreach (DEMOFILE ${SEQAN_DEMOS})
	get_filename_component (DEMO ${DEMOFILE} NAME_WE)
	seqan_add_executable (${DEMO} ${DEMOFILE})
    if (ZLIB_FOUND)
      target_link_libraries(${DEMO} ${ZLIB_LIBRARIES})
    endif (ZLIB_FOUND)
    if (BZIP2_FOUND)
      target_link_libraries(${DEMO} ${BZIP2_LIBRARIES})
    endif (BZIP2_FOUND)
	add_dependencies (demos ${DEMO})
endforeach (DEMOFILE)

################################################################################
# Tutorials
################################################################################

add_custom_target(tutorials)

file(GLOB SEQAN_DEMOS ${SEQAN_LIBRARY}/demos/tutorial/[A-z]*.cpp
                      ${SEQAN_LIBRARY}/demos/tutorial/*/[A-z]*.cpp)
foreach (DEMOFILE ${SEQAN_DEMOS})
	get_filename_component (DEMO ${DEMOFILE} NAME_WE)
	seqan_add_executable (tutorial_${DEMO} ${DEMOFILE})
	add_dependencies (tutorials tutorial_${DEMO})
endforeach (DEMOFILE)

################################################################################
# HowTos
################################################################################

add_custom_target(howtos)

file(GLOB SEQAN_DEMOS ${SEQAN_LIBRARY}/demos/howto/[A-z]*.cpp
                      ${SEQAN_LIBRARY}/demos/howto/*/[A-z]*.cpp)
foreach (DEMOFILE ${SEQAN_DEMOS})
	get_filename_component (DEMO ${DEMOFILE} NAME_WE)
	seqan_add_executable (howto_${DEMO} ${DEMOFILE})
	add_dependencies (howtos howto_${DEMO})
endforeach (DEMOFILE)

