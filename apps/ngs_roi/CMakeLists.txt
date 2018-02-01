# ===========================================================================
#                  NGS: Regions of Interest Analysis
# ===========================================================================
# File: /sandbox/jagla/apps/ngs_roi/CMakeLists.txt
#
# CMakeLists.txt file for ngs_roi.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_ngs_roi CXX)
message (STATUS "Configuring apps/ngs_roi")

set (SEQAN_APP_VERSION "0.2.12")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.

if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (ZLIB)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

if (NOT ZLIB_FOUND)
  message (STATUS "  ZLIB not found: ngs_roi will not be built.")
  return ()
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable (roi_plot_thumbnails roi_plot_thumbnails.cpp png_canvas.h)
target_link_libraries (roi_plot_thumbnails ${SEQAN_LIBRARIES})

# bam2roi converts BAM file into ROI files.
add_executable (bam2roi bam2roi.cpp roi_builder.h roi_builder.cpp record_ext.h)
target_link_libraries (bam2roi ${SEQAN_LIBRARIES})

# roi_feature_projection allows the projection of ROIs to features (BED/GFF
# intervals or GFF exons).
add_executable (roi_feature_projection
                roi_feature_projection.cpp
                project_spliced.cpp
                project_spliced.h
                project_interval.cpp
                project_interval.h)
target_link_libraries (roi_feature_projection ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (ngs_roi)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install binaries into ${PREFIX}/bin directory.
install (TARGETS bam2roi
                 roi_feature_projection
                 roi_plot_thumbnails
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install helper scripts into ${PREFIX}/bin directory.
install (FILES # Scripts for sorting.
                 tool_shed/bed_sort.sh
                 tool_shed/gff_sort.sh
                 tool_shed/roi_sort.sh
                 # PDF plot script helper.
                 tool_shed/roi_plot_9.sh
                 tool_shed/plot.awk
                 tool_shed/ps2pswLinks.gawk
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/ngs_roi for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# Example files.
install (FILES example/dmel.bed
               example/dmel.gtf
               example/example.bam
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# R package for Linux and Mac.
install (FILES R/ngsroi_0.1.tar.gz
               R/ngsroi_0.1.zip
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/R)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (ngs_roi)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:ngs_roi")
  set (CPACK_PACKAGE_NAME "ngs_roi")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "ngs_roi")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Bernd Jagla <bernd.jagla@pasteur.fr>")
  set (CPACK_PACKAGE_VENDOR "Bernd Jagla <bernd.jagla@pasteur.fr>")

  seqan_configure_cpack_app (ngs_roi "ngs_roi")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:ngs_roi")

