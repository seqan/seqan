# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
# Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
# Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
# ============================================================================
#
#.rst:
# SeqanSimdUtility
# ----------------
#
# Utility to create binaries for different SIMD implementations and target
# architectures. Disables architectures if the compiler can't handle it or bugs
# in the code generation are known.
#
# .. code-block:: cmake
#
#   add_executable (
#     test_simd_vector
#     test_simd_vector.cpp
#     test_simd_vector.h)
#   target_link_libraries (test_simd_vector ${SEQAN_LIBRARIES})
#   add_simd_platform_tests(test_simd_vector)
#
# Could build the following binaries and add tests for them
#
# ::
#
#   Built target test_simd_vector
#   Built target test_simd_vector_sse4
#   Built target test_simd_vector_avx2
#   Built target test_simd_vector_umesimd
#   Built target test_simd_vector_umesimd_sse4
#   Built target test_simd_vector_umesimd_avx2
#   Built target test_simd_vector_umesimd_avx512_knl
#
# Functions
# +++++++++++
#
# ::
#
#   add_simd_platform_tests(target)
#       - builds target on different simd implementations and architectures
#       - adds targets which use the seqan::simd default implementation
#       - adds targets which use the UME::SIMD library instead of the seqan
#         default implementation if UMESIMD_FOUND is set
#       - adds tests for different architectures if the Intel® Software
#         Development Emulator could be found, i.e. if SDE_FOUND is set
#
#         NOTE: You don't have to add_test the target (e.g. test_simd_vector)
#         itself, because it will be done by this function.
#
#   add_simd_executables(target blacklist) [not required to use explicitly]
#       - adds executables for sse4, avx2, etc. as `target`_sse4, `target`_avx2,
#         etc. (Will clone the SOURCE, FLAGS and other properties from `target`)
#       - blacklist: list of architectures that shouldn't be added (e.g. sse4,
#         avx2, avx512_knl, avx512)
#
#   add_simd_tests(target blacklist) [not required to use explicitly]
#       - adds tests for executables (targets) specified by add_simd_executables
#       - blacklist: list of architectures which executables shouldn't be tested
#         (e.g. sse4, avx2, avx512_knl, avx512)
#

find_package (SDE)
find_package (Umesimd)

macro(transfer_target_property property source_target target_target)
    get_target_property(_property_value ${source_target} ${property})

    # message(STATUS "${source_target}: ${property} == ${_property_value}")

    if(_property_value)
        # message(STATUS "${target_target}: set ${property} = ${_property_value}")
        set_target_properties(${target_target} PROPERTIES ${property} "${_property_value}")
    endif()
endmacro()

macro(clone_target source_target target_target)
    get_target_property(_SOURCES ${source_target} SOURCES)

    # message(STATUS "transfer properties from ${source_target} to ${target_target}")
    # message(STATUS "${source_target}: SOURCES == ${_SOURCES}")
    # message(STATUS "${target_target}: set SOURCES = ${_SOURCES}")
    add_executable(${target_target} ${_SOURCES})

    # https://cmake.org/cmake/help/v3.4/manual/cmake-properties.7.html#properties-on-targets
    transfer_target_property(SOURCE_DIR ${source_target} ${target_target})
    transfer_target_property(COMPILE_DEFINITIONS ${source_target} ${target_target})
    transfer_target_property(COMPILE_FEATURES ${source_target} ${target_target})
    transfer_target_property(COMPILE_FLAGS ${source_target} ${target_target})
    transfer_target_property(COMPILE_OPTIONS ${source_target} ${target_target})
    transfer_target_property(LINK_LIBRARIES ${source_target} ${target_target})
    transfer_target_property(LINK_FLAGS ${source_target} ${target_target})
    transfer_target_property(CXX_EXTENSIONS ${source_target} ${target_target})
    transfer_target_property(CXX_STANDARD ${source_target} ${target_target})
    transfer_target_property(CXX_STANDARD_REQUIRED ${source_target} ${target_target})
    transfer_target_property(FOLDER ${source_target} ${target_target})
    transfer_target_property(INCLUDE_DIRECTORIES ${source_target} ${target_target})
endmacro(clone_target)

macro(add_simd_executables target blacklist)
    if (NOT (";${blacklist};" MATCHES ";sse4;"))
        clone_target(${target} "${target}_sse4")
        target_compile_options("${target}_sse4" PRIVATE -msse4)
    endif()

    if (NOT (";${blacklist};" MATCHES ";avx2;"))
        clone_target(${target} "${target}_avx2")
        target_compile_options("${target}_avx2" PRIVATE -mavx2)
    endif()

    if (NOT (";${blacklist};" MATCHES ";avx512_knl;"))
        clone_target(${target} "${target}_avx512_knl")
        target_compile_options("${target}_avx512_knl" PRIVATE -mavx512f -mavx512cd -mavx512er -mavx512pf)
    endif()

    if (NOT (";${blacklist};" MATCHES ";avx512;"))
        clone_target(${target} "${target}_avx512")
        target_compile_options("${target}_avx512" PRIVATE -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi)
    endif()
endmacro(add_simd_executables)

macro(add_simd_tests target blacklist)
    add_test(NAME "test_${target}" COMMAND $<TARGET_FILE:${target}>)

    if (NOT SDE_FOUND)
        message (STATUS "Intel Software Development Emulator not found, not building platform emulated tests.")
        return()
    endif ()

    if (TARGET "${target}_sse4" AND NOT (";${blacklist};" MATCHES ";sse4;"))
        add_test(NAME "test_${target}_sse4" COMMAND ${SDE_EXECUTABLE} -snb -- $<TARGET_FILE:${target}>_sse4)
    endif()

    if (TARGET "${target}_avx2" AND NOT (";${blacklist};" MATCHES ";avx2;"))
        add_test(NAME "test_${target}_avx2" COMMAND ${SDE_EXECUTABLE} -hsw -- $<TARGET_FILE:${target}>_avx2)
    endif()

    if (TARGET "${target}_avx512_knl" AND NOT (";${blacklist};" MATCHES ";avx512_knl;"))
        add_test(NAME "test_${target}_avx512_knl" COMMAND  ${SDE_EXECUTABLE} -knl -- $<TARGET_FILE:${target}>_avx512_knl)
    endif()

    if (TARGET "${target}_avx512" AND NOT (";${blacklist};" MATCHES ";avx512;"))
        add_test(NAME "test_${target}_avx512" COMMAND  ${SDE_EXECUTABLE} -skx -- $<TARGET_FILE:${target}>_avx512)
    endif()
endmacro(add_simd_tests)

macro(add_simd_platform_tests target)
    if (UMESIMD_FOUND)
        clone_target("${target}" "${target}_umesimd")
        target_compile_definitions("${target}_umesimd" PUBLIC SEQAN_UMESIMD_ENABLED=1)
    endif()

    # seqan-simd doesn't support avx512, but will fallback to avx2
    set(seqansimd_compile_blacklist "")
    set(seqansimd_test_blacklist "")

    # ume-simd has some problems with clang
    set(umesimd_compile_blacklist "")
    set(umesimd_test_blacklist "")

    # clang <= 3.9.x produces executables using invalid instructions for avx512
    if (COMPILER_CLANG AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0)
        set(umesimd_test_blacklist "avx512;avx512_knl")
        set(seqansimd_test_blacklist "${umesimd_test_blacklist}")
    endif()

    add_simd_executables("${target}" "${seqansimd_compile_blacklist}")
    add_simd_tests("${target}" "${seqansimd_test_blacklist}")

    if (UMESIMD_FOUND)
        add_simd_executables("${target}_umesimd" "${umesimd_compile_blacklist}")
        add_simd_tests("${target}_umesimd" "${umesimd_test_blacklist}")
    endif()
endmacro(add_simd_platform_tests)
