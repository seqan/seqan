# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
# Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

# Optionally find java for zip support
find_package (Java QUIET)

# Create the payload binary ZIP file.
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set (SEQAN_PLATFORM "lnx")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set (SEQAN_PLATFORM "win")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set (SEQAN_PLATFORM "mac")
else ()
  message (STATUS "Unsupported platform ${CMAKE_SYSTEM_NAME}, disabling CTD support.")
  return()
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
set (PAYLOAD_TMP_PATH ${CMAKE_BINARY_DIR}/CMakeFiles/payload.tmp)
set (PAYLOAD_TMP_BIN_PATH ${PAYLOAD_TMP_PATH}/bin)

# ============================================================================
# Creating directory structure.
# ============================================================================

# Create directory: workflow_plugin_dir
add_custom_command (OUTPUT ${WORKFLOW_PLUGIN_DIR}
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${WORKFLOW_PLUGIN_DIR})
# Create directory: workflow_plugin_dir/icons
add_custom_command (OUTPUT ${WORKFLOW_PLUGIN_DIR}/icons
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${WORKFLOW_PLUGIN_DIR}/icons
                    DEPENDS ${WORKFLOW_PLUGIN_DIR})
# Create directory: workflow_plugin_dir/descriptors
add_custom_command (OUTPUT ${CTD_PATH}
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${CTD_PATH}
                    DEPENDS ${WORKFLOW_PLUGIN_DIR})
# Create directory: workflow_plugin_dir/payload
add_custom_command (OUTPUT ${PAYLOAD_PATH}
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_PATH}
                    DEPENDS ${WORKFLOW_PLUGIN_DIR})

# Create directory: workflow_plugin_dir/payload.tmp
add_custom_command (OUTPUT ${PAYLOAD_TMP_PATH}
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_TMP_PATH}
                    DEPENDS ${WORKFLOW_PLUGIN_DIR})
# Create directory: workflow_plugin_dir/payload.tmp/bin
add_custom_command (OUTPUT ${PAYLOAD_TMP_BIN_PATH}
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_TMP_BIN_PATH}
                    DEPENDS ${PAYLOAD_TMP_PATH})

# ============================================================================
# Creating payload data.
# ============================================================================

# Binaries.
foreach (_BINARY ${SEQAN_CTD_EXECUTABLES})
  set (_TARGET ${_BINARY})
  set (_BINARY_PATH "${SEQAN_BIN_DIR}/${_BINARY}")
  if (WIN32)
    set (_BINARY "${_BINARY}.exe")
    set (_BINARY_PATH "${_BINARY_PATH}.exe")
  endif ()

  list (APPEND TMP_PAYLOAD_FILES "${PAYLOAD_TMP_BIN_PATH}/${_BINARY}")
  add_custom_command (OUTPUT ${PAYLOAD_TMP_BIN_PATH}/${_BINARY}
                      COMMAND ${CMAKE_COMMAND} -E copy "${_BINARY_PATH}" "${PAYLOAD_TMP_BIN_PATH}/${_BINARY}"
                      DEPENDS ${_TARGET}
                              ${PAYLOAD_TMP_BIN_PATH})
endforeach ()

# binaries.ini file.
add_custom_command (OUTPUT "${PAYLOAD_TMP_PATH}/binaries.ini"
                    COMMAND ${CMAKE_COMMAND} -E touch ${PAYLOAD_TMP_PATH}/binaries.ini
                    DEPENDS ${PAYLOAD_TMP_PATH})

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set (SEQAN_SYSTEM_WORDSIZE "64")
else ()
  set (SEQAN_SYSTEM_WORDSIZE "32")
endif ()

set (_ZIP_NAME "binaries_${SEQAN_PLATFORM}_${SEQAN_SYSTEM_WORDSIZE}.zip")
set (_ZIP_PATH "${PAYLOAD_PATH}")

# ----------------------------------------------------------------------------
# different ways to zip
# ----------------------------------------------------------------------------

if (NOT CMAKE_VERSION VERSION_LESS "3.3") # internal since 3.3
    add_custom_command (OUTPUT ${_ZIP_PATH}/${_ZIP_NAME}
                        COMMAND  ${CMAKE_COMMAND} -E tar "cfv" ${_ZIP_PATH}/${_ZIP_NAME} --format=zip .
                        WORKING_DIRECTORY ${PAYLOAD_TMP_PATH}
                        DEPENDS ${PAYLOAD_PATH}
                                ${TMP_PAYLOAD_FILES}
                                ${PAYLOAD_TMP_PATH}/binaries.ini)
elseif (Java_JAR_EXECUTABLE) # use java
    add_custom_command (OUTPUT ${_ZIP_PATH}/${_ZIP_NAME}
                        COMMAND ${Java_JAR_EXECUTABLE} cfvM ${_ZIP_PATH}/${_ZIP_NAME} -C ${PAYLOAD_TMP_PATH} .
                        DEPENDS ${PAYLOAD_PATH}
                                ${TMP_PAYLOAD_FILES}
                                ${PAYLOAD_TMP_PATH}/binaries.ini)
elseif (NOT CMAKE_SYSTEM_NAME STREQUAL "Windows") # use unix zip
    add_custom_command (OUTPUT ${_ZIP_PATH}/${_ZIP_NAME}
                        COMMAND zip -r ${_ZIP_PATH}/${_ZIP_NAME} .
                        WORKING_DIRECTORY ${PAYLOAD_TMP_PATH}
                        DEPENDS ${PAYLOAD_PATH}
                                ${TMP_PAYLOAD_FILES}
                                ${PAYLOAD_TMP_PATH}/binaries.ini)
else ()
    message (STATUS "Either update cmake to >= 3.3 or install java to create CTDs!")
    return ()
endif ()

# ============================================================================
# CTDs and other descriptors contents.
# ============================================================================

# descriptors/mimetypes.xml
# Note: The mimetypes.xml file is deprecated but we keep it here for backward compatibilty.
add_custom_command (OUTPUT ${CTD_PATH}/mimetypes.xml
                    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/util/cmake/ctd/mimetypes.xml"
                                                     "${CTD_PATH}/mimetypes.xml"
                    DEPENDS ${CTD_PATH}
                            ${CMAKE_SOURCE_DIR}/util/cmake/ctd/mimetypes.xml)
list (APPEND DESCRIPTOR_FILES ${CTD_PATH}/mimetypes.xml)

# descriptors/mime.types
add_custom_command (OUTPUT ${CTD_PATH}/mime.types
                    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/util/cmake/ctd/mime.types"
                                                     "${CTD_PATH}/mime.types"
                    DEPENDS ${CTD_PATH}
                            ${CMAKE_SOURCE_DIR}/util/cmake/ctd/mime.types)
list (APPEND DESCRIPTOR_FILES ${CTD_PATH}/mime.types)

# *.ctd
foreach (_BINARY ${SEQAN_CTD_EXECUTABLES})
  set (_BINARY_PATH "${SEQAN_BIN_DIR}/${_BINARY}")
  if (WIN32)
    set (_BINARY_PATH "${_BINARY_PATH}.exe")
  endif ()

  add_custom_command (OUTPUT ${CTD_PATH}/${_BINARY}.ctd
                      COMMAND ${_BINARY_PATH} --write-ctd "${CTD_PATH}/${_BINARY}.ctd"
                      DEPENDS ${CTD_PATH}
                              ${_BINARY})
  list (APPEND DESCRIPTOR_FILES ${CTD_PATH}/${_BINARY}.ctd)
endforeach ()

# ============================================================================
# Eclipse Plugin Files.
# ============================================================================

# plugin.properties
add_custom_command (OUTPUT ${WORKFLOW_PLUGIN_DIR}/plugin.properties
                    COMMAND ${CMAKE_COMMAND} "-DSEQAN_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
                                             "-DWORKFLOW_PLUGIN_DIR=${WORKFLOW_PLUGIN_DIR}"
                                             "-DSEQAN_VERSION_STRING=${SEQAN_VERSION_STRING}"
                                             "-DSEQAN_DATE=${SEQAN_DATE}"
                                             "-DCTD_PLUGIN_PACKAGE=${CTD_PLUGIN_PACKAGE}"
                                             "-DCTD_PLUGIN_NAME=${CTD_PLUGIN_NAME}"
                                             -P "${CMAKE_SOURCE_DIR}/util/cmake/ctd/configure_profile_properties.cmake"
                    DEPENDS ${WORKFLOW_PLUGIN_DIR}
                            ${CMAKE_SOURCE_DIR}/util/cmake/ctd/plugin.properties.in)

# ============================================================================
# Static Files.
# ============================================================================

# Static files in plugin root.
foreach (_FILE COPYRIGHT DESCRIPTION LICENSE)
  add_custom_command (OUTPUT ${WORKFLOW_PLUGIN_DIR}/${_FILE}
                      COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/util/cmake/ctd/${_FILE}"
                                                       "${WORKFLOW_PLUGIN_DIR}/${_FILE}"
                      DEPENDS ${WORKFLOW_PLUGIN_DIR}
                              ${CMAKE_SOURCE_DIR}/util/cmake/ctd/${_FILE})
  list (APPEND STATIC_FILES ${WORKFLOW_PLUGIN_DIR}/${_FILE})
endforeach ()

# Icon files.
foreach (_FILE category.png splash.png)
  add_custom_command (OUTPUT ${WORKFLOW_PLUGIN_DIR}/icons/${_FILE}
                      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/util/cmake/ctd/icons/${_FILE}
                                                       ${WORKFLOW_PLUGIN_DIR}/icons/${_FILE}
                      DEPENDS ${WORKFLOW_PLUGIN_DIR}/icons
                              ${CMAKE_SOURCE_DIR}/util/cmake/ctd/icons/${FILE})
  list (APPEND ICON_FILES ${WORKFLOW_PLUGIN_DIR}/icons/${_FILE})
endforeach ()

# ============================================================================
# Master target.
# ============================================================================

add_custom_target (prepare_workflow_plugin
                   DEPENDS ${DESCRIPTOR_FILES}
                           ${STATIC_FILES}
                           ${ICON_FILES}
                           ${WORKFLOW_PLUGIN_DIR}/plugin.properties
                           ${_ZIP_PATH}/${_ZIP_NAME})
