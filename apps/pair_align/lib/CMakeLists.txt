# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/pair_align/lib/CMakeLists.txt
#
# CMakeLists.txt file for pair_align_lib.
# ===========================================================================

# ---------------------------------------------------------------------------
# Generate all AlignConfig configurations.
# ---------------------------------------------------------------------------

set(FLAG_VALUES "f" "t")
foreach(GAPS_T ${FLAG_VALUES})
	foreach(GAPS_L ${FLAG_VALUES})
		foreach(GAPS_R ${FLAG_VALUES})
			foreach(GAPS_B ${FLAG_VALUES})
				list(APPEND ALIGN_CONFIG_FLAGS "${GAPS_T}_${GAPS_L}_${GAPS_R}_${GAPS_B}")
			endforeach(GAPS_B)
		endforeach(GAPS_R)
	endforeach(GAPS_L)
endforeach(GAPS_T)

# ---------------------------------------------------------------------------
# Generate target names for global alignment with the align_config as 
# the suffix of the target name.
# ---------------------------------------------------------------------------

foreach(CONFIG ${ALIGN_CONFIG_FLAGS})
	list(APPEND GLOBAL_ALIGN_TARGETS "pair_align_global_${CONFIG}")
endforeach(CONFIG)


# ---------------------------------------------------------------------------
# Add static libraries for all pair_align configurations, which are later
# linked to by the pair_align_lib.
#
# Note with cmake 2.8.8 and higher we could use the OBJECT feature of
# add_library, which only builds the object files without linking them into
# libraries.
# ---------------------------------------------------------------------------

# Add library for lcs option.
add_library(pair_align_lcs STATIC pair_align_lcs.cpp)
# Add library for local option.
add_library(pair_align_local STATIC pair_align_local.cpp)
# Add library for all global options.
foreach(TARGET ${GLOBAL_ALIGN_TARGETS})	
	add_library(${TARGET} STATIC pair_align_global.cpp)
endforeach(TARGET)

# ---------------------------------------------------------------------------
# Add COMPILE_DEFINITIONS in order to generate all pairAlignGlobal_<config>
# interfaces using preprocessor macros in pair_align_global.cpp.
# ---------------------------------------------------------------------------

# Extract number of different configurations
list(LENGTH GLOBAL_ALIGN_TARGETS NUM_TARGETS)
math(EXPR MAX_RANGE "${NUM_TARGETS} - 1")

# Enumerate all possible configurations and add compile definitions accordingly.
foreach(INDEX RANGE ${MAX_RANGE})
	list(GET GLOBAL_ALIGN_TARGETS ${INDEX} TARGET)
	list(GET ALIGN_CONFIG_FLAGS ${INDEX} CONFIG)
	
	# Form a list out of the suffix string.	
	string(REPLACE "f" "0" CONFIG_1 ${CONFIG})
	string(REPLACE "t" "1" CONFIG_2 ${CONFIG_1})
	string(REPLACE "_" ";" CONFIG_LIST ${CONFIG_2})
	list(GET CONFIG_LIST 0 GAP_T)
	list(GET CONFIG_LIST 1 GAP_L)
	list(GET CONFIG_LIST 2 GAP_R)
	list(GET CONFIG_LIST 3 GAP_B)

	# Add properties for the four free-gap options.
	set_target_properties (${TARGET} PROPERTIES COMPILE_DEFINITIONS "SUFFIX_GAP_TOP=${GAP_T}")
	set_property (TARGET ${TARGET} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "SUFFIX_GAP_LEFT=${GAP_L}")
	set_property (TARGET ${TARGET} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "SUFFIX_GAP_RIGHT=${GAP_R}")
	set_property (TARGET ${TARGET} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "SUFFIX_GAP_BOTTOM=${GAP_B}")
endforeach(INDEX)

# ---------------------------------------------------------------------------
# Set references to the corresponding object libraries.
# ---------------------------------------------------------------------------

# Add lcs and local reference.
set(LIB_TARGETS pair_align_lcs)
list(APPEND LIB_TARGETS pair_align_local)

# Add reference for all global configurations.
foreach(TARGET ${GLOBAL_ALIGN_TARGETS})
	list(APPEND LIB_TARGETS ${TARGET})
endforeach(TARGET)

## The commented parts can be used instead if minimum cmake version is 2.8.8 
## or higher. Then the OBJECT feature is available and here we would link 
## against the object files generated in the steps before.
##set(OBJECT_REF "$<TARGET_OBJECTS:pair_align_lcs>")
##list(APPEND OBJECT_REF "$<TARGET_OBJECTS:pair_align_local>")

## Add reference for all global configurations.
##foreach(TARGET ${GLOBAL_ALIGN_TARGETS})
##	list(APPEND OBJECT_REF "$<TARGET_OBJECTS:${TARGET}>")
##endforeach(TARGET)

# ---------------------------------------------------------------------------
# Build static library linking all the object libraries.
# ---------------------------------------------------------------------------

add_library(pair_align_lib pair_align_lib.cpp pair_align_lib.h)
target_link_libraries(pair_align_lib ${LIB_TARGETS} ${SEQAN_LIBRARIES})

## Can be enabled instead for cmake version greater or equal to 2.8.8.
##add_library(pair_align_lib pair_align_lib.cpp pair_align_lib.h ${OBJECT_REF})
##target_link_libraries(pair_align_lib ${SEQAN_LIBRARIES})
