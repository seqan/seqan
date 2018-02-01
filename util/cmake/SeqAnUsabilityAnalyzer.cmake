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
# Define the macros required for using the instrumentation for the SeqAn
# Usability Analyzer.
# ============================================================================

# ---------------------------------------------------------------------------
# Macro seqan_setup_sua ()
#
# Setup target for the SUA data collection initialization.
# ---------------------------------------------------------------------------

function (seqan_setup_sua)
  if (NOT SEQAN_INSTRUMENTATION)
    return ()  # Do not process further when instrumentation is disabled.
  endif (NOT SEQAN_INSTRUMENTATION)

  message(STATUS "Prepare SeqAn Usability Analyzer data collection...")
#  if (CMAKE_HOST_WIN32 AND NOT PYTHON_EXECUTABLE)
#    # Use EXE for instrumentation with bundled Python runtime.
#    execute_process (COMMAND ${CMAKE_SOURCE_DIR}/misc/seqan_instrumentation/py2exe/dist/seqan_instrumentation.exe cmake ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})
#    add_custom_target (seqan_sua_target ${CMAKE_SOURCE_DIR}/misc/seqan_instrumentation/py2exe/dist/seqan_instrumentation.exe build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
#                       COMMENT "Build Instrumentation...")
#  else (CMAKE_HOST_WIN32 AND NOT PYTHON_EXECUTABLE)
    # Use system's Python runtime.
    execute_process (COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/misc/seqan_instrumentation/bin/seqan_instrumentation.py cmake ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})
    add_custom_target (seqan_sua_target ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/misc/seqan_instrumentation/bin/seqan_instrumentation.py build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
                       COMMENT "Build Instrumentation...")
#  endif (CMAKE_HOST_WIN32 AND NOT PYTHON_EXECUTABLE)
endfunction (seqan_setup_sua)

# ---------------------------------------------------------------------------
# Macro seqan_add_sua_dependency (TARGET)
#
# Add hooks into the SUA data collection before and after building TARGET.
# ---------------------------------------------------------------------------

function (seqan_add_sua_dependency TARGET)
  if (NOT SEQAN_INSTRUMENTATION)
    return ()  # Do not add if instrumentation is disabled.
  endif (NOT SEQAN_INSTRUMENTATION)

  # Add dependency on the instrumentation target.
  add_dependencies (${TARGET} seqan_sua_target)

  # Add hooks before and after building.
#  if (CMAKE_HOST_WIN32 AND NOT PYTHON_EXECUTABLE)
#    set (_INST ${CMAKE_SOURCE_DIR}/misc/seqan_instrumentation/py2exe/dist/seqan_instrumentation.exe)
#    add_custom_command (TARGET ${TARGET}
#                        PRE_BUILD
#                        COMMAND ${_INST} pre_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
#                        COMMENT "Pre Build Instrumentation...")
#    add_custom_command (TARGET ${TARGET}
#                        POST_BUILD
#                        COMMAND ${_INST} post_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
#                        COMMENT "Post Build Instrumentation...")
#  else (CMAKE_HOST_WIN32 AND NOT PYTHON_EXECUTABLE)
    set (_INST ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/misc/seqan_instrumentation/bin/seqan_instrumentation.py)
    add_custom_command (TARGET ${TARGET}
                        PRE_BUILD
                        COMMAND ${_INST} pre_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
                        COMMENT "Pre Build Instrumentation...")
    add_custom_command (TARGET ${TARGET}
                        POST_BUILD
                        COMMAND ${_INST} post_build ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${TARGET}
                        COMMENT "Post Build Instrumentation...")
#  endif (CMAKE_HOST_WIN32 AND NOT PYTHON_EXECUTABLE)
endfunction (seqan_add_sua_dependency TARGET)
