cmake_minimum_required (VERSION 2.6)
project (core_apps_rabema)

# We need the module bam_io from extras, so add the necessary include
# directory here.
include_directories(../../../extras/include)


# Preparation of SAM files, i.e. filling SEQ field of each record.
seqan_add_executable(rabema_prepare_sam prepare_sam.cpp sorting.h)

if (ZLIB_FOUND)
    # Building of GSI file from "perfect" SAM file and reference sequence.
    seqan_add_executable(rabema_build_gold_standard rabema_build_gold_standard.cpp
                         curve_smoothing.h
                         find_hamming_simple_ext.h
                         find_myers_ukkonen_ext.h
                         find_myers_ukkonen_reads.h
                         io_gsi.h
                         rabema_stats.h
                         ref_id_mapping.h
                         sorting.h)

    # Comparison of arbitrary read mapper SAM against GSI file.
    seqan_add_executable(rabema_evaluate rabema_evaluate.cpp
                         curve_smoothing.h
                         find_hamming_simple_ext.h
                         find_myers_ukkonen_ext.h
                         find_myers_ukkonen_reads.h
                         io_gsi.h
                         rabema_stats.h
                         ref_id_mapping.h
                         sorting.h)

    # Install rabema in /bin directory
    install(TARGETS rabema_prepare_sam rabema_build_gold_standard rabema_evaluate
            DESTINATION bin)
else (ZLIB_FOUND)
    message(STATUS "INFO: Rabema not build because zlib was not found.")
endif (ZLIB_FOUND)

# Small tool for helping in debugging.
seqan_add_executable(do_search do_search.cpp find_myers_ukkonen_ext.h)


# Run app tests if python was found
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_rabema COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  #message(STATUS "add_test(NAME app_test_rabema COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif(PYTHONINTERP_FOUND)

