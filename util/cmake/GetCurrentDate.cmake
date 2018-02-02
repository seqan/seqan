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
# Get current date from compiled C++ program.
#
# The C++ program GetCurrentDate.cpp belongs together with this module.
#
# The resulting date will be written to the variables CURRENT_YEAR,
# CURRENT_MONTH, and CURRENT_DAY.
# ============================================================================

if (NOT CMAKE_CURRENT_LIST_DIR)  # CMAKE_CURRENT_LIST_DIR only from cmake 2.8.3.
  get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
endif (NOT CMAKE_CURRENT_LIST_DIR)

try_run(_GET_CURRENT_DATE_RUN_RESULT
        _GET_CURRENT_DATE_COMPILE_RESULT
        ${CMAKE_BINARY_DIR}/CMakeFiles/GetCurrentDate
        ${CMAKE_CURRENT_LIST_DIR}/GetCurrentDate.cpp
        COMPILE_OUTPUT_VARIABLE _GET_CURRENT_DATE_COMPILE_OUTPUT
        RUN_OUTPUT_VARIABLE _GET_CURRENT_DATE_RUN_OUTPUT)
if (NOT _GET_CURRENT_DATE_RUN_OUTPUT)
  message(FATAL_ERROR "Could not determine current date!.")
endif ()

string(REGEX REPLACE ".*YEAR ([0-9a-zA-Z]+).*" "\\1" CURRENT_YEAR ${_GET_CURRENT_DATE_RUN_OUTPUT})
string(REGEX REPLACE ".*MONTH ([0-9a-zA-Z]+).*" "\\1" CURRENT_MONTH ${_GET_CURRENT_DATE_RUN_OUTPUT})
string(REGEX REPLACE ".*DAY ([0-9a-zA-Z]+).*" "\\1" CURRENT_DAY ${_GET_CURRENT_DATE_RUN_OUTPUT})
