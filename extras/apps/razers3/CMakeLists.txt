cmake_minimum_required (VERSION 2.6)
project (sandbox_weese_apps_razers3)

#if (OPENMP_FOUND)
  # We need a configuration header file.  First, configure this file, then
  # add the target directory to the inclue path.
  #configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  #                ${CMAKE_CURRENT_BINARY_DIR}/config.h)
  #include_directories(${CMAKE_CURRENT_BINARY_DIR})

  # Enable OpenMP.
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fno-exceptions")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# disable RazerS3 for MinGW and VS8 and earlier
if (NOT MINGW AND NOT(MSVC AND MSVC_VERSION LESS 1500))
  # Add executable for RazerS 3.
  seqan_add_executable (razers3 razers.cpp
                                razers.h
                                job_queue.h
                                outputFormat.h
                                parallel_job_queue.h
                                parallel_misc.h
                                parallel_store.h
                                paramChooser.h
                                param_tabs.cpp
                                param_tabs.inc
                                param_tabs.h
                                profile_timeline.h
                                razers_match_filter.h
                                razers_paired_match_filter.h
                                razers_matepairs.h
                                razers_matepairs_parallel.h
                                razers_parallel.h
                                razers_window.h
                                readSimulator.h)
else (NOT MINGW AND NOT(MSVC AND MSVC_VERSION LESS 1500))
	message (STATUS "Not building RazerS 3: MinGW does not provide __sync_val_compare_and_swap.")
endif (NOT MINGW AND NOT(MSVC AND MSVC_VERSION LESS 1500))

#  seqan_add_executable (compPHSens compPHSens.cpp razers.h)
  seqan_add_executable (simulate_reads simulate_reads.cpp)
  seqan_add_executable (quality2prob quality2prob.cpp)

  #elseif (OPENMP_FOUND)
  #  message(STATUS "Not building RazerS3, since OpenMP could not be found.")
  #endif (OPENMP_FOUND)
