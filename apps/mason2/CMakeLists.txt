# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /mason2/CMakeLists.txt
#
# CMakeLists.txt file for mason2.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_mason2 CXX)
message (STATUS "Configuring apps/mason2")

set (SEQAN_APP_VERSION "2.0.9")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
    find_package (ZLIB)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

if (NOT ZLIB_FOUND)
  message (STATUS "  ZLIB not found: MASON2 will not be built.")
  return ()
endif (NOT ZLIB_FOUND)

# Enable C++11 in GCC and/or CLANG.
#if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
#    add_definitions ("-std=c++0x")
#endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# We define a library for the reusable parts of Mason.
add_library (mason_sim
             external_split_merge.h
             external_split_merge.cpp
             fragment_generation.h
             genomic_variants.h
             genomic_variants.cpp
             individual_genome.h
             mason_options.cpp
             mason_options.h
             mason_types.h
             methylation_levels.cpp
             methylation_levels.h
             omp_wrapper.h
             sequencing.h
             simulate_454.cpp
             simulate_base.cpp
             simulate_genome.cpp
             simulate_genome.h
             simulate_illumina.cpp
             simulate_sanger.cpp
             variation_size_tsv.h
             vcf_materialization.cpp
             vcf_materialization.h)
target_link_libraries (mason_sim ${SEQAN_LIBRARIES})


# Simulate random genomes.
add_executable (mason_genome mason_genome.cpp)

# Simulate methylation levels for gnome.
add_executable (mason_methylation mason_methylation.cpp)

# Simulate sequencing of fragments.
add_executable (mason_frag_sequencing mason_frag_sequencing.cpp)

# Simulation of structural variants into FASTA sequences.
add_executable (mason_variator mason_variator.cpp)

# The materializer reads in FASTA and VCF and writes out FASTA.
add_executable (mason_materializer mason_materializer.cpp)

# The end-to-end read simulator includes materialization and fragment sampling.
add_executable (mason_simulator mason_simulator.cpp)

# Simulation splicing given a reference, GFF/GTF file, and possibly VCF for variants.
add_executable (mason_splicing mason_splicing.cpp)

# Tests for some of the library functionality.
add_executable (mason_tests mason_tests.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (mason_frag_sequencing ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_genome ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_materializer ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_methylation ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_simulator ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_splicing ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_tests ${SEQAN_LIBRARIES} mason_sim)
target_link_libraries (mason_variator ${SEQAN_LIBRARIES} mason_sim)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include mason executables in CTD structure.
set (SEQAN_CTD_EXECUTABLES
    ${SEQAN_CTD_EXECUTABLES}
    mason_genome
    mason_frag_sequencing
    mason_variator
    mason_materializer
    mason_simulator
    CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (mason2)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install Mason binaries in ${PREFIX}/bin directory
install (TARGETS mason_frag_sequencing
                 mason_genome
                 mason_materializer
                 mason_methylation
                 mason_simulator
                 mason_splicing
                 mason_variator
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/mason2 for SeqAn release builds.
install (FILES LICENSE
               README
               README.mason_frag_sequencing
               README.mason_genome
               README.mason_materializer
               README.mason_methylation
               README.mason_simulator
               README.mason_splicing
               README.mason_variator
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/adeno_virus.fa
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (mason2)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:mason2")
  set (CPACK_PACKAGE_NAME "mason2")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Versatile package for the simulation of nucleic sequence.")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (mason2 "mason2")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:mason2")
