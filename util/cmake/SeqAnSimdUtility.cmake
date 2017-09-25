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
#         avx2, avx512_knl, avx512_skx, avx512_cnl)
#
#   add_simd_tests(target blacklist) [not required to use explicitly]
#       - adds tests for executables (targets) specified by add_simd_executables
#       - blacklist: list of architectures which executables shouldn't be tested
#         (e.g. sse4, avx2, avx512_knl, avx512_skx, avx512_cnl)
#

find_package (SeqAn CONFIG REQUIRED)
find_package (SDE)
find_package (Umesimd)

set(SEQAN_SIMD_SUPPORTED_EXTENSIONS "sse4;avx2;avx512_knl;avx512_skx;avx512_cnl")

if (COMPILER_MSVC)
    set(SEQAN_SIMD_SSE4_OPTIONS /arch:AVX)
    set(SEQAN_SIMD_AVX2_OPTIONS /arch:AVX2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS "")
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS "")
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS "")
elseif (COMPILER_WINTEL)
    set(SEQAN_SIMD_SSE4_OPTIONS /QxSSE4.2)
    set(SEQAN_SIMD_AVX2_OPTIONS /QxCORE-AVX2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS /QxMIC-AVX512)
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS /QxCORE-AVX512)
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS /QxCORE-AVX512)
elseif (COMPILER_LINTEL)
    set(SEQAN_SIMD_SSE4_OPTIONS -xSSE4.2)
    set(SEQAN_SIMD_AVX2_OPTIONS -xCORE-AVX2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS -xMIC-AVX512)
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS -xCORE-AVX512)
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS -xCORE-AVX512)
else()
    set(SEQAN_SIMD_SSE4_OPTIONS -msse4)
    set(SEQAN_SIMD_AVX2_OPTIONS -mavx2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS -mavx512f -mavx512cd -mavx512er -mavx512pf)
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl)
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS "${SEQAN_SIMD_AVX512_SKX_OPTIONS}" -mavx512ifma -mavx512vbmi)
endif()

set(SEQAN_SIMD_AVX2_SOURCE
"#include <immintrin.h>

int main() {
  _mm256_add_epi32(__m256i{}, __m256i{});
  return 0;
}")

set(SEQAN_SIMD_AVX512_KNL_SOURCE
"#include <cstdint>
#include <immintrin.h>

int main() {
  alignas(64) int64_t a[]{0,1,2,3,4,5,6,7};
  _mm512_mask_load_epi64(__m512i{}, 0, a);
  return 0;
}")

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
        target_compile_options("${target}_sse4" PRIVATE "${SEQAN_SIMD_SSE4_OPTIONS}")
    endif()

    if (NOT (";${blacklist};" MATCHES ";avx2;"))
        clone_target(${target} "${target}_avx2")
        target_compile_options("${target}_avx2" PRIVATE "${SEQAN_SIMD_AVX2_OPTIONS}")
    endif()

    if (NOT (";${blacklist};" MATCHES ";avx512_knl;"))
        clone_target(${target} "${target}_avx512_knl")
        target_compile_options("${target}_avx512_knl" PRIVATE "${SEQAN_SIMD_AVX512_KNL_OPTIONS}")

        if (COMPILER_MSVC)
            message(STATUS "avx512_knl not supported on msvc")
        endif()
    endif()

    if (NOT (";${blacklist};" MATCHES ";avx512_skx;"))
        clone_target(${target} "${target}_avx512_skx")
        target_compile_options("${target}_avx512_skx" PRIVATE "${SEQAN_SIMD_AVX512_SKX_OPTIONS}")

        if (COMPILER_MSVC)
            message(STATUS "avx512_skx not supported on msvc")
        endif()
    endif()

    if (NOT (";${blacklist};" MATCHES ";avx512_cnl;"))
        clone_target(${target} "${target}_avx512_cnl")
        target_compile_options("${target}_avx512_cnl" PRIVATE "${SEQAN_SIMD_AVX512_CNL_OPTIONS}")

        if (COMPILER_MSVC)
            message(STATUS "avx512_cnl not supported on msvc")
        endif()
    endif()
endmacro(add_simd_executables)

macro(add_simd_tests target blacklist)
    add_test(NAME "test_${target}" COMMAND $<TARGET_FILE:${target}>)

    if (NOT SDE_FOUND)
        message (STATUS "Intel Software Development Emulator not found, not building platform emulated tests.")
    else ()
        if (TARGET "${target}_sse4" AND NOT (";${blacklist};" MATCHES ";sse4;"))
            add_test(NAME "test_${target}_sse4" COMMAND ${SDE_EXECUTABLE} -snb -- $<TARGET_FILE:${target}_sse4>)
        endif()

        if (TARGET "${target}_avx2" AND NOT (";${blacklist};" MATCHES ";avx2;"))
            add_test(NAME "test_${target}_avx2" COMMAND ${SDE_EXECUTABLE} -hsw -- $<TARGET_FILE:${target}_avx2>)
        endif()

        if (TARGET "${target}_avx512_knl" AND NOT (";${blacklist};" MATCHES ";avx512_knl;"))
            add_test(NAME "test_${target}_avx512_knl" COMMAND  ${SDE_EXECUTABLE} -knl -- $<TARGET_FILE:${target}_avx512_knl>)
        endif()

        if (TARGET "${target}_avx512_skx" AND NOT (";${blacklist};" MATCHES ";avx512_skx;"))
            add_test(NAME "test_${target}_avx512_skx" COMMAND  ${SDE_EXECUTABLE} -skx -- $<TARGET_FILE:${target}_avx512_skx>)
        endif()

        if (TARGET "${target}_avx512_cnl" AND NOT (";${blacklist};" MATCHES ";avx512_cnl;"))
            add_test(NAME "test_${target}_avx512_cnl" COMMAND  ${SDE_EXECUTABLE} -cnl -- $<TARGET_FILE:${target}_avx512_cnl>)
        endif()
    endif()
endmacro(add_simd_tests)

macro(add_simd_platform_tests target)
    if (UMESIMD_FOUND)
        clone_target("${target}" "${target}_umesimd")
        target_include_directories("${target}_umesimd" PUBLIC "${UMESIMD_INCLUDE_DIR}")
        target_compile_definitions("${target}_umesimd" PUBLIC SEQAN_UMESIMD_ENABLED=1)
    endif()

    # We don't disable AVX512 even though seqan-simd doesn't support AVX512 (it
    # will fallback to AVX2 in source code), but it replaces some intrinsics
    # with ones introduced in AVX512.
    set(seqansimd_compile_blacklist "")
    set(seqansimd_test_blacklist "")

    set(umesimd_compile_blacklist "")
    set(umesimd_test_blacklist "")

    if (COMPILER_GCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
        ## seqan-simd:
        # avx2: gcc4.9 has a problem to differentiate `__vector(32) signed char`
        #       from `__vector(16) signed char` (some conflicts in the abi).
        #       Theoretically, one could build avx2, but one wouldn't be allowed
        #       to mix it with sse4 in the complete program. This can be fixed
        #       by adding `-fabi-version=6` to the compiler flags, but since we
        #       don't know which consequences this has, we disable it.
        #
        ## ume-simd:
        # avx512_knl/avx512_skx/avx512_cnl:
        #           gcc4.9 has no complete avx512_knl intrinsics support, which
        #           are used in umesimd
        #   error: ‘__mmask32’ does not name a type
        #   error: ‘_mm512_castsi128_si512’ was not declared in this scope
        set(umesimd_compile_blacklist "avx512_knl;avx512_skx;avx512_cnl")
        set(seqansimd_compile_blacklist "avx2;avx512_knl;avx512_skx;avx512_cnl")
    endif()

    if (COMPILER_CLANG AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0.1)
            message(STATUS "Clang 4.0.0 crashes for AVX512_skx (seqan-simd only), see https://llvm.org/bugs/show_bug.cgi?id=31731")
            set(umesimd_compile_blacklist "")
            set(seqansimd_compile_blacklist "avx512_skx;avx512_cnl")
        endif()

        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.0.0)
            message(AUTHOR_WARNING "Clang >=4.0.1 reevaluate if AVX512_skx (seqan-simd only) binaries are working.")
            set(umesimd_compile_blacklist "")
            set(seqansimd_compile_blacklist "")
        endif()
    endif()

    # Build the executables, but don't execute them, because clang <= 3.9.x
    # produces executables which contain invalid instructions for AVX512.
    if (COMPILER_CLANG AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0)

        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9.1)
            message(STATUS "Clang 3.9.0 produces executables which contain invalid instructions for AVX512_skx and AVX512_knl")
            set(umesimd_test_blacklist "avx512_knl;avx512_skx;avx512_cnl")
            set(seqansimd_test_blacklist "${umesimd_test_blacklist}")
        endif()

        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9.2)
            message(STATUS "Clang 3.9.1 produces executables which contain invalid instructions for AVX512_skx (seqan-simd only), see https://llvm.org/bugs/show_bug.cgi?id=31731")
            set(umesimd_test_blacklist "")
            set(seqansimd_test_blacklist "avx512_skx;avx512_cnl")
        endif()

        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3.9.1)
            message(AUTHOR_WARNING "Clang >=3.9.2 reevaluate if AVX512_skx (seqan-simd only) binaries are working.")
            set(umesimd_test_blacklist "")
            set(seqansimd_test_blacklist "")
        endif()

    endif()

    if (COMPILER_CLANG AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9)
        ## all-simd wrappers
        # avx512_cnl: clang compiler <= 3.8.x can't handle avx512 compiler flags
        #   error : unknown argument: '-mavx512ifma'
        #   error : unknown argument: '-mavx512vbmi'
        #
        ## ume-simd:
        # avx512_knl: clang compiler <= 3.8.x has no complete avx512_knl
        #             intrinsics support, which are used in umesimd
        #   error : use of undeclared identifier '_mm512_castsi128_si512'
        #   error : use of undeclared identifier '_mm512_mask_mov_epi32'
        #
        ## seqan-simd:
        # avx512_knl: clang compiler 3.6.x crashes on avx512_knl, thus disabling
        #             it also for 3.7.x and 3.8.x, where avx512_knl works
        set(umesimd_compile_blacklist "avx512_knl;avx512_skx;avx512_cnl")
        set(seqansimd_compile_blacklist "avx512_knl;avx512_skx;avx512_cnl")

        ## seqan-simd:
        # avx2: clang 3.7.x fails test cases:
        #       SimdVectorTestCommon_ShuffleConstant1 type parameter (un-)signed char __vector(32) FAILED
        #       SimdVectorTestCommon_ShuffleConstant1 type parameter (un-)signed short __vector(16) FAILED
        if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6) AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8))
            set(seqansimd_compile_blacklist "avx2;${seqansimd_compile_blacklist}")
        endif()

        if (STDLIB_VS)
            ## ume-simd:
            # avx2: Also, it can't handle avx2 on windows, because some
            #       intrinsics are not implemented yet:
            #   %16 = call <8 x i32> @llvm.x86.avx2.gather.d.d.256(<8 x i32> %7, i8* %9, <8 x i32> %13, <8 x i32> %15, i8 4), !dbg !7025
            set(umesimd_compile_blacklist "avx2;${umesimd_compile_blacklist}")

            ## seqan-simd: Similar, simd-seqan doesn't work, even though
            #              it should theoretically be possible, but some
            #              intrinsic are not yet implemented.
            # sse4: fatal error C1001: An internal error has occurred in the
            #       compiler.
            #   clang!DllGetC2Telemetry()+0x1d7c28
            #   clang!LLVM_IR_InvokeCompilerPassW()+0xd127
            #   clang!DllGetC2Telemetry()+0xfa13e
            #   clang!crt_at_quick_exit()+0x104
            #   clang!BaseThreadInitThunk()+0x12
            # avx2: 'mm256_unpackhi_epi128': Intrinsic not yet implemented:
            #   %8 = call <4 x i64> @llvm.x86.avx2.vperm2i128(<4 x i64> %5, <4 x i64> %7, i8 49), !dbg !5915
            #   fatal error C1001: An internal error has occurred in the compiler.
            set(seqansimd_compile_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
        endif()
    endif()

    if (COMPILER_MSVC)
        # msvc 2015/2017 only supports /arch:AVX and /arch:AVX2, thus don't
        # compile AVX512
        set(umesimd_compile_blacklist "avx512_knl;avx512_skx;avx512_cnl")

        # Don't compile simd-seqan, because it uses a compiler extensions which
        # is not supported by msvc
        set(seqansimd_compile_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
    endif()

    if (COMPILER_WINTEL)
        # Don't compile simd-seqan, because wintel compiler doesn't support
        # vector extension even though lintel does.
        set(seqansimd_compile_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
    endif()

    if (COMPILER_LINTEL AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0.2)
        ## seqan-simd: intel compiler crashes on <= 17.0.1 at least in the test
        #              case test_align_simd, others may work (e.g.
        #              test_simd_vector)
        #   ": internal error: ** The compiler has encountered an unexpected
        #   problem. ** Segmentation violation signal raised. ** Access
        #   violation or stack overflow. Please contact Intel Support for
        #   assistance.
        if(NOT ("${target}" MATCHES "test_simd_vector"))
            set(seqansimd_compile_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
        endif()
    endif()

    if(COMPILER_LINTEL AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16.0.2)
        ## seqan-simd: intel compiler crashes on <= 16.0.1
        #   ": internal error: ** The compiler has encountered an unexpected
        #   problem. ** Segmentation violation signal raised. ** Access
        #   violation or stack overflow. Please contact Intel Support for
        #   assistance.
        set(seqansimd_compile_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
    endif()

    if(COMPILER_GCC AND DEFINED ENV{TRAVIS} AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0.0)
        # Disable avx2 on travis, because sometimes gcc 5.x and 6.x crashes
        set(seqansimd_compile_blacklist "${seqansimd_compile_blacklist};avx2")
        set(umesimd_compile_blacklist "${umesimd_compile_blacklist};avx2")
        message(STATUS "Disabling avx2 on travis, because gcc<=6.x crashes occasionally")
    endif()

    # This checks if a simple avx2 program builds
    set(CMAKE_REQUIRED_FLAGS "${SEQAN_SIMD_AVX2_OPTIONS}")
    check_cxx_source_compiles("${SEQAN_SIMD_AVX2_SOURCE}" AVX2_SOURCE_COMPILES)
    unset(CMAKE_REQUIRED_FLAGS)
    if (NOT AVX2_SOURCE_COMPILES)
        set(seqansimd_compile_blacklist "${seqansimd_compile_blacklist};avx2")
        set(umesimd_compile_blacklist "${umesimd_compile_blacklist};avx2")
        message(STATUS "Your compiler/linker doesn't support avx2")
    endif()

    # This check is primarily for travis, which has an old linker
    set(CMAKE_REQUIRED_FLAGS "${SEQAN_SIMD_AVX512_KNL_OPTIONS}")
    check_cxx_source_compiles("${SEQAN_SIMD_AVX512_KNL_SOURCE}" KNL_SOURCE_COMPILES)
    unset(CMAKE_REQUIRED_FLAGS)
    if (NOT KNL_SOURCE_COMPILES)
        set(seqansimd_compile_blacklist "${seqansimd_compile_blacklist};avx512_knl;avx512_skx;avx512_cnl")
        set(umesimd_compile_blacklist "${umesimd_compile_blacklist};avx512_knl;avx512_skx;avx512_cnl")
        message(STATUS "Your compiler/linker doesn't support avx512")
    endif()

    # don't use simd on 32bit architectures at all
    if (CMAKE_SIZEOF_VOID_P EQUAL 4)
        set(seqansimd_compile_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
        set(umesimd_compile_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
        message(STATUS "SIMD acceleration is only available on 64bit systems")
    endif()

    add_simd_executables("${target}" "${seqansimd_compile_blacklist}")
    add_simd_tests("${target}" "${seqansimd_test_blacklist}")

    if (UMESIMD_FOUND)
        add_simd_executables("${target}_umesimd" "${umesimd_compile_blacklist}")
        add_simd_tests("${target}_umesimd" "${umesimd_test_blacklist}")
    endif()
endmacro(add_simd_platform_tests)
