# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
# Copyright (c) 2006-2012, Knut Reinert, FU Berlin
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Knut Reinert or the FU Berlin nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
# ============================================================================
# Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
# ============================================================================
# CMake code for generation of CTD structure.
#
# All executables in the list SEQAN_CTD_EXECUTABLES will be included in the
# resulting CTD structure.
#
# The following environment variables configure the output:
#
#   WORKFLOW_PLUGIN_DIR -- Output directory for CTD structure.  Defaults to
#                          ${CMAKE_BINARY_DIR}/workflow_plugin_dir
# ============================================================================

# ============================================================================
# Dependency Check
# ============================================================================

# If Java cannot be found, we disable CTD support.

find_package (Java)
if (NOT Java_JAR_EXECUTABLE)
  message (STATUS "jar binary not found, disabling CTD support.")
  return ()
endif ()

# ============================================================================
# Variable Setup
# ============================================================================

# Get path that all binaries are placed in.  With MSVC, we have to extend that
# path with the configuration name.
set (SEQAN_BIN_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
if (MSVC)
  set (SEQAN_BIN_DIR "${SEQAN_BIN_DIR}/$(ConfigurationName)")
endif ()

# The user-definable setting for the output plugin dir.
set (WORKFLOW_PLUGIN_DIR ${CMAKE_BINARY_DIR}/workflow_plugin_dir CACHE PATH
     "Directory containing the generated plugin-sources for the SeqAn Workflow/KNIME package")
# Shortcut to "descriptors" below target directory.
set (CTD_PATH ${WORKFLOW_PLUGIN_DIR}/descriptors)
# Shortcut to "payload" below target directory.
set (PAYLOAD_PATH ${WORKFLOW_PLUGIN_DIR}/payload)

# We will create the contents of the payload directory temporarily within the
# output directory.
set (PAYLOAD_TMP_PATH ${WORKFLOW_PLUGIN_DIR}/payload.tmp)
set (PAYLOAD_TMP_BIN_PATH ${PAYLOAD_TMP_PATH}/bin)

# ============================================================================
# Targets for preparing payload
# ============================================================================

# ----------------------------------------------------------------------------
# Targets for creating payload "bin" dir
# ----------------------------------------------------------------------------

# Target for creating the CTD output directories.  Triggers building of CTD
# executables through dependency.
add_custom_target (prepare_ctd_payload_tmp_bin
                   # These two commands are equivalent to rm -rf.
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_TMP_PATH}
                   COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_TMP_PATH}
                   # Create payload directory.
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_TMP_PATH}
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_TMP_BIN_PATH}
                   # Make sure that all binaries that we want to have CTDs for are built.
                   DEPENDS ${SEQAN_CTD_EXECUTABLES})

# For each white-listed executable in SEQAN_CTD_EXECUTABLES, get its path and
# (1) copy it into the temporary payload directory's bin subdirectory, and (2)
# create a target target_ctd_ctds_${EXECUTABLE} that the CTD creation target
# below will depend on (through the list PREPARE_CTD_CTDS_TARGETS).
foreach (_BINARY ${SEQAN_CTD_EXECUTABLES})
  # Get platform-dependent executable path.
  set (_BINARY_PATH "${SEQAN_BIN_DIR}/${_BINARY}")
  if (WIN32)
    set (_BINARY_PATH "${_BINARY_PATH}.exe)")
  endif ()
  add_custom_command (TARGET prepare_ctd_payload_tmp_bin POST_BUILD
                      COMMAND ${CMAKE_COMMAND} -E copy ${_BINARY_PATH} ${PAYLOAD_TMP_BIN_PATH})
  add_custom_target (target_ctd_ctds_${_BINARY}
                     COMMAND ${_BINARY_PATH} --write-ctd ${CTD_PATH}/${_BINARY}.ctd)
  list (APPEND PREPARE_CTD_CTDS_TARGETS target_ctd_ctds_${_BINARY})
endforeach ()

# NOTE: When adding lib, share, etc. directories to payload, add targets
# following the pattern for the binary files and register as dependency
# everywhere below.

# ----------------------------------------------------------------------------
# Targets for creating payload archive
# ----------------------------------------------------------------------------

# Get system name and word size for the payload archive name.

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set (SEQAN_PLATFORM "lnx")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set (SEQAN_PLATFORM "win")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set (SEQAN_PLATFORM "mac")
else ()
  message (FATAL_ERROR "Unknown platform name ${SEQAN_PLATFORM}")
endif ()

if (SEQAN_SYSTEM_PROCESSOR MATCHES ".*64.*")
  set (SEQAN_SYSTEM_WORDSIZE "64")
else ()
  set (SEQAN_SYSTEM_WORDSIZE "32")
endif ()

# We use the jar command for creating the payload files in a platform independent manner.

set (_ZIP_NAME "binaries_${SEQAN_PLATFORM}_${SEQAN_SYSTEM_WORDSIZE}.zip")
set (_ZIP_PATH "${PAYLOAD_PATH}")
add_custom_target (prepare_ctd_payload
                   # rm -rf
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${_ZIP_PATH}
                   COMMAND ${CMAKE_COMMAND} -E remove_directory ${_ZIP_PATH}
                   # mkdir
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${_ZIP_PATH}
                   # compress
                   COMMAND ${Java_JAR_EXECUTABLE} cfvM ${_ZIP_PATH}/${_ZIP_NAME} -C ${PAYLOAD_TMP_PATH} .
                   # remove temporary files
                   COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_TMP_PATH}
                   DEPENDS prepare_ctd_payload_tmp_bin)

# ============================================================================
# Targets for creating Eclipse plugin files.
# ============================================================================

# ----------------------------------------------------------------------------
# Copy static files (LICENSE etc.)
# ----------------------------------------------------------------------------

add_custom_target (prepare_ctd_static_files
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${WORKFLOW_PLUGIN_DIR}
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/util/cmake/ctd/COPYRIGHT ${WORKFLOW_PLUGIN_DIR}
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/util/cmake/ctd/DESCRIPTION ${WORKFLOW_PLUGIN_DIR}
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/util/cmake/ctd/LICENSE ${WORKFLOW_PLUGIN_DIR}
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/util/cmake/ctd/mimetypes.xml ${WORKFLOW_PLUGIN_DIR}/descriptors)

add_custom_target (prepare_ctd_icons
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${WORKFLOW_PLUGIN_DIR}/icons
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/util/cmake/ctd/icons ${WORKFLOW_PLUGIN_DIR}/icons)

# ----------------------------------------------------------------------------
# Configure plugin.properties.
# ----------------------------------------------------------------------------

# If possible, get latest change date from SeqAn SVN.
find_package(Subversion)
if (Subversion_FOUND)
  file (TO_CMAKE_PATH "${CMAKE_SOURCE_DIR}" _SEQAN_SOURCE_DIR)
  Subversion_WC_INFO (${_SEQAN_SOURCE_DIR} SEQAN)
  string(REGEX REPLACE "^([0-9]+)-([0-9]+)-([0-9]+) ([0-9]+):([0-9]+).*"
    "\\1\\2\\3\\4\\5" SEQAN_LAST_CHANGE_DATE "${SEQAN_WC_LAST_CHANGED_DATE}")
  message (STATUS "SEQAN_LAST_CHANGE_DATE ${SEQAN_LAST_CHANGE_DATE}")
  message (STATUS "SEQAN_VERSION_STRING ${SEQAN_VERSION_STRING}")
  set (CF_SEQAN_VERSION ${SEQAN_VERSION_STRING}.${SEQAN_LAST_CHANGE_DATE})
else ()
  set (CF_SEQAN_VERSION "${SEQAN_VERSION_STRING}")
endif ()

# Configure the file.
add_custom_target (prepare_ctd_properties
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CTD_PATH}
                   COMMAND ${CMAKE_COMMAND} -DSEQAN_SOURCE_DIR=${CMAKE_SOURCE_DIR} -DWORKFLOW_PLUGIN_DIR=${WORKFLOW_PLUGIN_DIR}
                                            -DCF_SEQAN_VERSION=${CF_SEQAN_VERSION}
                                            -P ${CMAKE_SOURCE_DIR}/util/cmake/ctd/configure_profile_properties.cmake)

# ----------------------------------------------------------------------------
# Create *.ctd files.
# ----------------------------------------------------------------------------

add_custom_target (prepare_ctd_ctds
                   DEPENDS ${PREPARE_CTD_CTDS_TARGETS})

# ============================================================================
# Top-level target
# ============================================================================

# Depends on all targets for preparing the workflow plugin.

add_custom_target (prepare_workflow_plugin
                   DEPENDS prepare_ctd_payload prepare_ctd_static_files
                           prepare_ctd_properties prepare_ctd_ctds prepare_ctd_icons)
