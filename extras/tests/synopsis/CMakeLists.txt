cmake_minimum_required (VERSION 2.6)
project (extras_tests_synopsis)

if (TR1_UNORDERED_MAP_FOUND)
    # We need a configuration header file.  First, configure this file, then
    # add the target directory to the inclue path.
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
                    ${CMAKE_CURRENT_BINARY_DIR}/config.h)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})

    # Update the list of file names below if you add source files to your test.
    seqan_add_test_executable(test_synopsis test_synopsis.cpp test_synopsis_counter_buckets.h test_synopsis_hot_list.h test_synopsis_histogram_set.h)
else (TR1_UNORDERED_MAP_FOUND)
    message(STATUS "Not building test_synopsis: Could not find unordered_map implementation.")
endif (TR1_UNORDERED_MAP_FOUND)
