# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/anise/CMakeLists.txt
#
# CMakeLists.txt file for anise.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_anise)
message (STATUS "Configuring extras/apps/anise")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for C++11.
find_package (CXX11)

if (NOT CXX11_FOUND)
    message (STATUS "  C++11 not found, not building ANISE.")
	return ()
endif ()

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# Search for LEMON libraries.
find_package (LEMON QUIET)

if (NOT LEMON_FOUND)
  message (STATUS "  LEMON library not found: Not building anise.")
  return ()
endif (NOT LEMON_FOUND)

# Search for Boost >= 1.51.0.
find_package (Boost 1.51.0 REQUIRED)

# If the BOOST library could not be found then do not configure BASIL.
if (NOT Boost_FOUND)
  message (STATUS "  BOOST library not found: Not building basil.")
  return ()
endif (NOT Boost_FOUND)

set (CXX_FLAGS "${CXX_FLAGS} -fno-omit-frame-pointers")

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories for SeqAn and other dependencies.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM ${LEMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})

# Also add include directories for shared stuff.
include_directories (shared .)

# Add CXX flags found by find_package  .
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} -pthread")

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS} -DSEQAN_NO_GLOBAL_EXCEPTION_HANDLER=1)

# Library with the progress indicator.
add_library (progress_indicator
             shared/progress_indicator.h
             shared/progress_indicator.cpp)

# Assembler module used in ANISE.
add_library (anise_asm
             asm/ag_to_frag_store.h
             asm/best_overlap_graph.cpp
             asm/best_overlap_graph.h
             asm/bog_unitigger.h
             asm/bog_unitigger.cpp
             asm/consensus_builder.cpp
             asm/consensus_builder.h
             asm/contig_graph.cpp
             asm/contig_graph.h
             asm/mira_unitigger.h
             asm/mira_unitigger.cpp
             asm/options.h
             asm/overlap_candidate_generator.cpp
             asm/overlap_store.cpp
             asm/overlapper.cpp
             asm/overlapper.h
             asm/store_utils.cpp
             asm/store_utils.h
             asm/unitigger.cpp
             asm/unitigger.h
             asm/unitigger_impl.h
             asm/unitigger_impl.cpp)

# Module used for scaffolding in ANISE assembler.
add_library (anise_scaffolder
             scaffolder/gpm.h
             scaffolder/gpm_options.h
             scaffolder/link_bundler.cpp
             scaffolder/mate_link.cpp
             scaffolder/mate_link.h
             scaffolder/path_merger.cpp
             scaffolder/scaffolding_result.h
             scaffolder/scaffolding_result.cpp
             scaffolder/transitive_reducer.cpp
             scaffolder/utils.h
             scaffolder/overlap_resolution.h
             scaffolder/overlap_resolution.cpp
             )
target_link_libraries (anise_scaffolder anise_asm)

# Library with repeat separation.
add_library (anise_rep_sep
             rep_sep/clique.cpp
             rep_sep/clique.h
             rep_sep/clique_enum.cpp
             rep_sep/clique_enum.h
             rep_sep/clique_set.cpp
             rep_sep/clique_set.h
             rep_sep/cluster_linking.cpp
             rep_sep/cluster_linking.h
             rep_sep/feature_map.cpp
             rep_sep/feature_map.h
             rep_sep/feature_vector.cpp
             rep_sep/feature_vector.h
             rep_sep/local_variation_store.cpp
             rep_sep/local_variation_store.h
             rep_sep/merge_read_set.cpp
             rep_sep/merge_read_set.h
             rep_sep/pair_based_sep.cpp
             rep_sep/pair_based_sep.h
             rep_sep/profile_computation.cpp
             rep_sep/profile_computation.h
             rep_sep/read.cpp
             rep_sep/read.cpp
             rep_sep/read.h
             rep_sep/read_correction.cpp
             rep_sep/read_correction.h
             rep_sep/read_cover.cpp
             rep_sep/read_cover.h
             rep_sep/read_set.cpp
             rep_sep/read_set.h
             rep_sep/rep_sep.cpp
             rep_sep/rep_sep.h
             rep_sep/separating_columns.cpp
             rep_sep/separating_columns.h
             rep_sep/separation_sites.cpp
             rep_sep/separation_sites.h
             rep_sep/shared_typedefs.h
             rep_sep/utils.h
             rep_sep/weighted_feature_vector.cpp
             rep_sep/weighted_feature_vector.h)
target_link_libraries (anise_rep_sep anise_asm)

# Library with repeat resolution (reconstruct global haplotypes of repeats from
# local haplotype reconstruction).
add_library (anise_rep_solv
             rep_solv/contig_graph.h
             rep_solv/contig_graph.cpp
             rep_solv/options.h
             rep_solv/dead_branch.h
             rep_solv/dead_branch.cpp
             rep_solv/directed_tree.h
             rep_solv/directed_tree.cpp
             rep_solv/path_enumeration.h
             rep_solv/path_enumeration.cpp
             rep_solv/scaffolding.h
             rep_solv/scaffolding.cpp
             rep_solv/graph_expansion.h
             rep_solv/graph_expansion.cpp)
target_link_libraries (anise_rep_solv anise_asm anise_rep_sep anise_scaffolder)

# Evaluation app.
add_executable (lba eval/lba.cpp )
target_link_libraries (lba progress_indicator ${SEQAN_LIBRARIES})

# ANISE application.
add_executable (anise
                anise/anise_app.cpp
                anise/anise_app.h
                anise/anise_options.cpp
                anise/anise_options.h
                anise/app_state.cpp
                anise/app_state.h
                anise/assembler.cpp
                anise/assembler.h
                anise/assembly_substep.cpp
                anise/assembly_substep.h
                anise/bam_tag_names.cpp
                anise/bam_tag_names.h
                anise/file_name_tokens.cpp
                anise/file_name_tokens.h
                anise/finishing_step.cpp
                anise/finishing_step.h
                anise/library_info.cpp
                anise/library_info.h
                anise/parallel_utils.cpp
                anise/parallel_utils.h
                anise/preparation_step.cpp
                anise/preparation_step.h
                anise/read_mapping_impl.h
                anise/read_mapping_impl_all.cpp
                anise/read_mapping_impl_all_best.cpp
                anise/read_mapping_impl_any_best.cpp
                anise/read_mapping_substep.cpp
                anise/read_mapping_substep.h
                anise/realignment.cpp
                anise/realignment.h
                anise/site_data.cpp
                anise/site_data.h
                anise/site_data_updater.cpp
                anise/site_data_updater.h
                anise/site_state.cpp
                anise/site_state.h
                anise/store_utils.cpp
                anise/store_utils.h
                anise/temporary_file_manager.cpp
                anise/temporary_file_manager.h
                anise/time_log.cpp
                anise/time_log.h
                shared/bam_io_utils.h
                shared/bam_io_utils.cpp
                shared/streaming_exception.h
                shared/utils.h)
target_link_libraries (anise progress_indicator anise_rep_solv anise_asm anise_rep_sep anise_scaffolder
                             ${SEQAN_LIBRARIES})

add_executable (basil
                basil/bam_filter_pipeline.cpp
                basil/bam_filter_pipeline.h
                basil/basil_app.cpp
                basil/basil_app.h
                basil/basil_stream.cpp
                basil/clipping_clustering.cpp
                basil/clipping_clustering.h
                basil/cluster_matching.cpp
                basil/cluster_matching.h
                basil/filter_clean_pairs.cpp
                basil/filter_clean_pairs.h
                basil/filter_discordant.cpp
                basil/filter_discordant.h
                basil/filter_high_coverage.cpp
                basil/filter_high_coverage.h
                basil/filter_realignment.cpp
                basil/filter_realignment.h
                basil/filter_remove_clipping.cpp
                basil/filter_remove_clipping.h
                basil/filter_reorder.cpp
                basil/filter_reorder.h
                basil/filter_shared_types.h
                basil/thread_safe_queue.h
                basil/oea_clustering.cpp
                basil/oea_clustering.h
                basil/utils.h)
target_link_libraries (basil progress_indicator ${SEQAN_LIBRARIES})

add_executable (oea_extractor
                shared/orphan_extractor.h
                basil/oea_extractor.cpp
                basil/bam_filter_pipeline.cpp
                basil/bam_filter_pipeline.h
                basil/filter_clean_pairs.cpp
                basil/filter_clean_pairs.h
                basil/filter_discordant.cpp
                basil/filter_discordant.h
                basil/filter_high_coverage.cpp
                basil/filter_high_coverage.h
                basil/filter_realignment.cpp
                basil/filter_realignment.h
                basil/filter_remove_clipping.cpp
                basil/filter_remove_clipping.h
                basil/filter_reorder.cpp
                basil/filter_reorder.h
                basil/filter_shared_types.h
                basil/thread_safe_queue.h
                basil/utils.h)
target_link_libraries (oea_extractor progress_indicator ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (anise)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install anise in ${PREFIX}/bin directory
install (TARGETS anise
                 basil
                 lba
         DESTINATION bin)
install (PROGRAMS scripts/filter_basil.py
                  scripts/extract_spanning.awk
                  scripts/best_blat.py
         DESTINATION bin)

# Install example files
install (FILES example/left.fq.gz
               example/right.fq.gz
               example/ref.fa
               example/ins.fa
         DESTINATION example)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/anise for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (anise)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:anise")
  set (CPACK_PACKAGE_NAME "anise")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "anise")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (anise "anise")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:anise")

