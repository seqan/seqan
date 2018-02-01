# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/sak/CMakeLists.txt
#
# CMakeLists.txt file for sak.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_sak CXX)
message (STATUS "Configuring apps/sak")

set (SEQAN_APP_VERSION "0.4.8")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn CONFIG REQUIRED)
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

# Update the list of file names below if you add source files to your application.
add_executable (sak sak.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (sak ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Generate text and man page information with sak_doc target.
add_custom_command (OUTPUT README.sak.txt
                    COMMAND sak --export-help=txt > README.sak.txt
                    DEPENDS sak)
add_custom_command (OUTPUT sak.1
                    COMMAND sak --export-help=man > sak.1
                    DEPENDS sak)
add_custom_target (sak_doc DEPENDS README.sak.txt sak.1)

# Build documentation on installation.
install (CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build ${CMAKE_BINARY_DIR} --target sak_doc)")

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (sak)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install sak in ${PREFIX}/bin directory
install (TARGETS sak
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/sak for SeqAn release builds.
install (FILES LICENSE
               README
               ${CMAKE_CURRENT_BINARY_DIR}/README.sak.txt
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/sak.1
         DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (sak)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------
	
# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} sak CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:sak")
  set (CPACK_PACKAGE_NAME "sak")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "sak - Swiss Army Knife")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (sak "sak")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:sak")

