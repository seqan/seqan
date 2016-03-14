# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/pair_align/lib/CMakeLists.txt
#
# CMakeLists.txt file for pair_align_lib.
# ===========================================================================

# ---------------------------------------------------------------------------
# Configure static libraries for different alignment configurations.
# ---------------------------------------------------------------------------

# Add library for lcs and local option.
add_library(pair_align_lcs OBJECT pair_align_lcs.cpp)
add_library(pair_align_local OBJECT pair_align_local.cpp)

# Add lcs and local target names to the lib target list.
set(LIB_TARGETS "$<TARGET_OBJECTS:pair_align_lcs>")
list(APPEND LIB_TARGETS "$<TARGET_OBJECTS:pair_align_local>")

# Add library for global options. 
# For every combination of AlignConfig flags there is a separate library, which 
# is build independently and thus saving resources during compilation.
set(FLAG_VALUES "0" "1")
foreach(GAPS_T ${FLAG_VALUES})
	foreach(GAPS_L ${FLAG_VALUES})
		foreach(GAPS_R ${FLAG_VALUES})
			foreach(GAPS_B ${FLAG_VALUES})
				# Set target name composed of the four gap configurations.
				set(TARGET "pair_align_global_${GAPS_T}${GAPS_L}${GAPS_R}${GAPS_B}")
				
				# Build library from source -> Note: STATIC could be replaced with OBJECT only building 
				# the object files not linking them into a lirary. However this feature
				# is shipped with cmake 2.8.8, which is not supported at the moment.
                add_library(${TARGET} OBJECT pair_align_global.cpp)

				# Add properties for the four free-gap options.
				set_target_properties (${TARGET} PROPERTIES COMPILE_DEFINITIONS "SUFFIX_GAP_TOP=${GAPS_T}")
				set_property (TARGET ${TARGET} APPEND PROPERTY COMPILE_DEFINITIONS "SUFFIX_GAP_LEFT=${GAPS_L}")
				set_property (TARGET ${TARGET} APPEND PROPERTY COMPILE_DEFINITIONS "SUFFIX_GAP_RIGHT=${GAPS_R}")
				set_property (TARGET ${TARGET} APPEND PROPERTY COMPILE_DEFINITIONS "SUFFIX_GAP_BOTTOM=${GAPS_B}")
				
				# Add target name to lib targets.
                list(APPEND LIB_TARGETS "$<TARGET_OBJECTS:${TARGET}>")
			endforeach(GAPS_B)
		endforeach(GAPS_R)
	endforeach(GAPS_L)
endforeach(GAPS_T)

# ---------------------------------------------------------------------------
# Build static library linking all the object libraries.
# ---------------------------------------------------------------------------

# This part can be used instead for cmake version greater or equal to 2.8.8. 
# This requires to update the script above as well.
add_library(pair_align_lib pair_align_lib.cpp pair_align_lib.h ${LIB_TARGETS})
target_link_libraries(pair_align_lib ${SEQAN_LIBRARIES})
