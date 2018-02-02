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

find_package (SDE)
find_package (Umesimd)

include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)

set(SEQAN_SIMD_UTILITY_VERBOSE OFF)
set(SEQAN_SIMD_SUPPORTED_EXTENSIONS "sse4;avx2;avx512_knl;avx512_skx;avx512_cnl")

if (COMPILER_MSVC)
    set(SEQAN_SIMD_SSE4_FLAGS "/arch:AVX")
    set(SEQAN_SIMD_AVX2_FLAGS "/arch:AVX2")
    set(SEQAN_SIMD_AVX512_KNL_FLAGS "")
    set(SEQAN_SIMD_AVX512_SKX_FLAGS "")
    set(SEQAN_SIMD_AVX512_CNL_FLAGS "")
elseif (COMPILER_WINTEL)
    set(SEQAN_SIMD_SSE4_FLAGS "/QxSSE4.2")
    set(SEQAN_SIMD_AVX2_FLAGS "/QxCORE-AVX2")
    set(SEQAN_SIMD_AVX512_KNL_FLAGS "/QxMIC-AVX512")
    set(SEQAN_SIMD_AVX512_SKX_FLAGS "/QxCORE-AVX512")
    set(SEQAN_SIMD_AVX512_CNL_FLAGS "/QxCORE-AVX512")
elseif (COMPILER_LINTEL)
    set(SEQAN_SIMD_SSE4_FLAGS "-xSSE4.2")
    set(SEQAN_SIMD_AVX2_FLAGS "-xCORE-AVX2")
    set(SEQAN_SIMD_AVX512_KNL_FLAGS "-xMIC-AVX512")
    set(SEQAN_SIMD_AVX512_SKX_FLAGS "-xCORE-AVX512")
    set(SEQAN_SIMD_AVX512_CNL_FLAGS "-xCORE-AVX512")
else()
    set(SEQAN_SIMD_SSE4_FLAGS "-msse4")
    set(SEQAN_SIMD_AVX2_FLAGS "-mavx2")
    set(SEQAN_SIMD_AVX512_KNL_FLAGS "-mavx512f -mavx512cd -mavx512er -mavx512pf")
    set(SEQAN_SIMD_AVX512_SKX_FLAGS "-mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl")
    set(SEQAN_SIMD_AVX512_CNL_FLAGS "${SEQAN_SIMD_AVX512_SKX_FLAGS} -mavx512ifma -mavx512vbmi")
endif()

# offer flags (string) as list for functions like target_compile_options
string (REPLACE " " ";" SEQAN_SIMD_SSE4_OPTIONS "${SEQAN_SIMD_SSE4_FLAGS}")
string (REPLACE " " ";" SEQAN_SIMD_AVX2_OPTIONS "${SEQAN_SIMD_AVX2_FLAGS}")
string (REPLACE " " ";" SEQAN_SIMD_AVX512_KNL_OPTIONS "${SEQAN_SIMD_AVX512_KNL_FLAGS}")
string (REPLACE " " ";" SEQAN_SIMD_AVX512_SKX_OPTIONS "${SEQAN_SIMD_AVX512_SKX_FLAGS}")
string (REPLACE " " ";" SEQAN_SIMD_AVX512_CNL_OPTIONS "${SEQAN_SIMD_AVX512_CNL_FLAGS}")

set(SEQAN_SIMD_SSE4_SDE_OPTIONS "-snb")
set(SEQAN_SIMD_AVX2_SDE_OPTIONS "-hsw")
set(SEQAN_SIMD_AVX512_KNL_SDE_OPTIONS "-knl")
set(SEQAN_SIMD_AVX512_SKX_SDE_OPTIONS "-skx")
set(SEQAN_SIMD_AVX512_CNL_SDE_OPTIONS "-cnl")

set(SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS)
set(SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD)
set(SEQAN_SIMD_HOST_CPU_SUPPORTS)

set(SEQAN_SIMD_SSE4_SOURCE
"#include <cstdint>
#include <immintrin.h>
#include <iostream>
#include <random>

int main() {
  std::random_device r;
  std::default_random_engine e(r());
  std::uniform_int_distribution<int32_t> d(1, 10);

  alignas(16) int32_t t[]{0,0,0,0};
  volatile auto a = _mm_set_epi32(d(e),d(e),d(e),d(e));
  volatile auto b = _mm_set_epi32(d(e),d(e),d(e),d(e));
  volatile auto z = _mm_max_epi32(a, b);
  _mm_store_si128((__m128i*)t, z);
  std::cout << \"(\" << t[0] << \", \" << t[1] << \", \" << t[2] << \", \" << t[3] << \")\" << std::endl;
  return 0;
}")

set(SEQAN_SIMD_AVX2_SOURCE
"#include <cstdint>
#include <immintrin.h>
#include <iostream>
#include <random>

int main() {
  std::random_device r;
  std::default_random_engine e(r());
  std::uniform_int_distribution<int32_t> d(1, 10);

  alignas(32) int64_t t[]{0,0,0,0};
  volatile auto a = _mm256_set_epi64x(d(e),d(e),d(e),d(e));
  volatile auto b = _mm256_set_epi64x(d(e),d(e),d(e),d(e));
  volatile auto z = _mm256_add_epi64(a, b);
  _mm256_store_si256((__m256i*)t, z);
  std::cout << \"(\" << t[0] << \", \" << t[1] << \", \" << t[2] << \", \" << t[3] << \")\" << std::endl;
  return 0;
}")

set(SEQAN_SIMD_AVX512_KNL_SOURCE
"#include <cstdint>
#include <immintrin.h>
#include <iostream>
#include <random>

int main() {
  std::random_device r;
  std::default_random_engine e(r());
  std::uniform_int_distribution<int64_t> d(1, 10);

  { // Some avx2 bug on AppleClang when compiling with avx512_skx or higher
    using TSimdVector = int8_t __attribute__ ((__vector_size__(32)));

    TSimdVector a{0u}, b{0u};

    for (auto i = 0; i < 32; ++i)
    {
        a[i] = (i - 1) * 3;
        b[i] = 32 - i;
    }
    TSimdVector c = _mm256_cmpeq_epi8(a, b);

    for (auto i = 0; i < 32; ++i)
    {
        if (c[i] == (a[i] == b[i]))
            std::cout << \"Failed!\" << std::endl;
    }
  }

  {
    alignas(64) uint64_t s[]{9,9,9,9,9,9,9,9};
    alignas(64) uint64_t t[]{0,0,0,0,0,0,0,0};

    // gcc 4.9 does not know _mm512_cmpgt_epu64_mask
    volatile auto a = _mm512_setr_epi64(d(e),d(e),d(e),d(e),d(e),d(e),d(e),d(e));
    volatile auto m = _mm512_cmpgt_epu64_mask(a, _mm512_set1_epi64(4)); // m = a > 4
    volatile auto z = _mm512_mask_load_epi64(a, m, s); // (a > 4) ? s : a
    _mm512_store_epi64(t, z);

    std::cout << \"(\" << t[0] << \", \" << t[1] << \", \" << t[2] << \", \" << t[3] << \", ...)\" << std::endl;
  }
  return 0;
}")

set(SEQAN_SIMD_AVX512_SKX_SOURCE "${SEQAN_SIMD_AVX512_KNL_SOURCE}")
set(SEQAN_SIMD_AVX512_CNL_SOURCE "${SEQAN_SIMD_AVX512_KNL_SOURCE}")

set(SEQAN_SIMD_SEQANSIMD_SOURCE
"#include <cstdint>
#include <iostream>
using int32x4_t = int32_t __attribute__ ((__vector_size__(4 * sizeof(int32_t)))); // SSE4 = 128bit
using int32x8_t = int32_t __attribute__ ((__vector_size__(8 * sizeof(int32_t)))); // AVX2 = 256bit

// gcc 4.9 bug (-fabi-version=6 (or =0) avoids this error with a change in mangling)
template <typename vector_t>
struct LENGTH;

template <>
struct LENGTH<int32x4_t>{
    static const std::size_t VALUE = 4;
};

template <>
struct LENGTH<int32x8_t>{
    static const std::size_t VALUE = 8;
};

// icc 16.0.0, 16.0.1 bug
template <typename TSimdVector, typename TValue>
void assign(TSimdVector & a, int index, TValue value) { a[index] = value; }

// icc >= 16.0.0 & <17.0.2 bug
struct simd_t { int32x4_t simd_vector; };

namespace ns {
template <typename T>
T & value(T * me) { return *me; }
}

template <typename value_t>
void destruct(value_t * p) { }

template <typename value_t>
struct string_t {
    value_t * value;

    template <typename T1>
    bool static assert(const T1 & value1, const T1 & value2) { return value1 <= value2; }

    string_t() : value(0) { assert(value, value); }
    ~string_t() { destruct(&ns::value(value)); }
};

template <typename value_t>
struct holder_t {
    value_t * value;

    holder_t() { value = new value_t; }
    ~holder_t() { delete value; }
};

struct matrix_t {  holder_t<string_t<simd_t>> cells; };

int main() {
  int32x4_t a{0,1,2,3}, b{4,3,2,1};
  int32x8_t x{0,1,2,3,4,5,6,7}, y{4,3,2,1,0,-1,-2,-3};
  auto c = a + b;
  auto z = x + y;

  // gcc 4.9 bug
  constexpr auto length1 = LENGTH<int32x4_t>::VALUE;
  constexpr auto length2 = LENGTH<int32x8_t>::VALUE;
  static_assert(length1 == 4u, \"\");
  static_assert(length2 == 8u, \"\");
  std::cout << \"length1: \" << length1 << std::endl;
  std::cout << \"length2: \" << length2 << std::endl;

  // icc 16.0.0, 16.0.1 bug
  assign(a, 0, 4);

  // icc >= 16.0.0 & <17.0.2 bug
  holder_t<matrix_t> matrix;
  return 0;
}")

set(SEQAN_SIMD_SEQANSIMD_AVX512_KNL_SOURCE
"#include <x86intrin.h>
#include <iostream>


int main() {
  // clang bug 4.0.0, https://bugs.llvm.org//show_bug.cgi?id=31731
  // -std=c++14 -mavx512bw -O3

  {
    using int8x32_t = signed char __attribute__ ((__vector_size__(32)));
    unsigned length = sizeof(int8x32_t) / sizeof(char);
    int8x32_t a{}, b{};
    for (auto i = 0u; i < length; ++i) { a[i] = i-1; b[i] = -i; }

    auto c = a < b;
    for(auto i = 0u; i < length; ++i)
      std::cout << (int)c[i] << std::endl;
  }

  // gcc 5.0 bug
  // -std=c++14 -mavx512f -O3
  {
    using int8x64_t = signed char __attribute__ ((__vector_size__(64)));
    unsigned length = sizeof(int8x64_t) / sizeof(char);
    int8x64_t a{}, b{};
    for (auto i = 0u; i < length; ++i) { a[i] = i-1; b[i] = -i; }

    auto c = a == b;
    for(auto i = 0u; i < length; ++i)
      std::cout << (int)c[i] << std::endl;
  }
  return 0;
}")

# list1 and list2
macro(list_intersect output list1 list2)
    set(${output})
    foreach(item ${list1})
        if (";${list2};" MATCHES ";${item};")
            list(APPEND ${output} "${item}")
        endif()
    endforeach()
endmacro()

# list1\list2 = list1 and not list2
macro(list_complement output list1 list2)
    set(${output})
    foreach(item ${list1})
        if (NOT(";${list2};" MATCHES ";${item};"))
            list(APPEND ${output} "${item}")
        endif()
    endforeach()
endmacro()

# Returns list of all simd extension prior to `simd_ext`
macro(simd_list_version_less output simd_ext)
  set(${output})

  foreach(item ${SEQAN_SIMD_SUPPORTED_EXTENSIONS})
      if (";${item};" MATCHES ";${simd_ext};")
          break()
      endif()
      list(APPEND ${output} "${item}")
  endforeach()
endmacro()

# Returns list of all simd extension after to `simd_ext`
macro(simd_list_version_greater output simd_ext)
  set(${output})
  simd_list_version_less(${output} "${simd_ext}")
  set(${output} ${${output}} "${simd_ext}")
  list_complement(${output} "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}" "${${output}}")
endmacro()

# simd specialized try-compile macro
macro(check_cxx_simd_source_runs SOURCE VAR)
    if(NOT DEFINED "${VAR}")
        # empty flags means that this simd extension is not supported by the compiler
        if (CMAKE_REQUIRED_FLAGS)
            # -O3 is needed to trigger a compiler bug on clang 4.0.0
            if (COMPILER_GCC OR COMPILER_CLANG OR COMPILER_LINTEL)
                set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -O3")
            endif()

            # CMAKE_REQUIRED_QUIET available since cmake 3.1.0
            set(CMAKE_REQUIRED_QUIET ON)
            check_cxx_source_runs("${SOURCE}" ${VAR})
            set(CMAKE_REQUIRED_QUIET OFF)
        endif()

        message(STATUS "Performing Test ${VAR}")
        if (NOT CMAKE_REQUIRED_FLAGS)
            set(${VAR} 0 CACHE INTERNAL "Test ${VAR}")
            message(STATUS "Performing Test ${VAR} - No flags to enable simd architecture")
        elseif (${VAR})
            message(STATUS "Performing Test ${VAR} - Success (Compiled & Host CPU support)")
        elseif (NOT ${VAR}_COMPILED)
            message(STATUS "Performing Test ${VAR} - Failed")
        elseif (NOT ${VAR})
            message(STATUS "Performing Test ${VAR} - Success (Compiled)")
        endif()
    endif()
endmacro()

# detect simd support and cache it
macro(detect_simd_support)
    # Executing simd on 32bit architectures has issues
    if (NOT DEFINED SEQAN_SIMD_DETECTED AND CMAKE_SIZEOF_VOID_P EQUAL 4)
        message(STATUS "SIMD acceleration is only available on 64bit systems")
        set(SEQAN_SIMD_DETECTED 1 CACHE INTERNAL "Detected SEQAN_SIMD_DETECTED")
        set(SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS "" CACHE INTERNAL "Test SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS")
        set(SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD "" CACHE INTERNAL "Test SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD")
        set(SEQAN_SIMD_HOST_CPU_SUPPORTS "" CACHE INTERNAL "Test SEQAN_SIMD_HOST_CPU_SUPPORTS")
    elseif (NOT DEFINED SEQAN_SIMD_DETECTED)
        # First try to compile AND run e.g. SEQAN_SIMD_SSE4_SOURCE, otherwise try to
        # only compile the program
        foreach(simd_ext ${SEQAN_SIMD_SUPPORTED_EXTENSIONS})
            string(TOUPPER ${simd_ext} SIMD_EXT)
            set(CMAKE_REQUIRED_FLAGS "${SEQAN_SIMD_${SIMD_EXT}_FLAGS}")
            check_cxx_simd_source_runs("${SEQAN_SIMD_${SIMD_EXT}_SOURCE}" ${SIMD_EXT}_SOURCE_RUNS)

            if (${SIMD_EXT}_SOURCE_RUNS)
                list(APPEND SEQAN_SIMD_HOST_CPU_SUPPORTS "${simd_ext}")
            endif()

            if (${SIMD_EXT}_SOURCE_RUNS_COMPILED)
                list(APPEND SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS "${simd_ext}")
            elseif (CMAKE_REQUIRED_FLAGS)
                message(STATUS "=> Abort simd detection of newer instruction sets")
                break()
            endif()
        endforeach()

        # test seqan simd
        set(CMAKE_REQUIRED_FLAGS "")
        check_cxx_source_compiles("${SEQAN_SIMD_SEQANSIMD_SOURCE}" SEQAN_SIMD_SEQANSIMD_SUPPORTED)

        # try-compile known compiler crashes/errors with seqan-simd and exclude them
        if (SEQAN_SIMD_SEQANSIMD_SUPPORTED)
            set(CMAKE_REQUIRED_FLAGS "${SEQAN_SIMD_AVX512_KNL_FLAGS}")
            check_cxx_simd_source_runs("${SEQAN_SIMD_SEQANSIMD_AVX512_KNL_SOURCE}" SEQANSIMD_AVX512_KNL_SOURCE_RUNS)

            if (SEQANSIMD_AVX512_KNL_SOURCE_RUNS_COMPILED)
                set(_SEQANSIMD_SUPPORTED "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}")
            else ()
                simd_list_version_less(_SEQANSIMD_SUPPORTED "avx512_knl")
            endif()
        endif()

        # seqan-simd also uses intrinsics, so exclude those that failed the intrinsic try-compile
        list_intersect(SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD "${_SEQANSIMD_SUPPORTED}" "${SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS}")

        # cache results
        set(SEQAN_SIMD_DETECTED 1 CACHE INTERNAL "Detected SEQAN_SIMD_DETECTED")
        set(SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS "${SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS}" CACHE INTERNAL "Test SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS")
        set(SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD "${SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD}" CACHE INTERNAL "Test SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD")
        set(SEQAN_SIMD_HOST_CPU_SUPPORTS "${SEQAN_SIMD_HOST_CPU_SUPPORTS}" CACHE INTERNAL "Test SEQAN_SIMD_HOST_CPU_SUPPORTS")

        if(SEQAN_SIMD_UTILITY_VERBOSE)
            message(STATUS "SEQAN_SIMD_SEQANSIMD_SUPPORTED: ${SEQAN_SIMD_SEQANSIMD_SUPPORTED}")
            message(STATUS "SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS: ${SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS}")
            message(STATUS "SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD: ${SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD}")
            message(STATUS "SEQAN_SIMD_HOST_CPU_SUPPORTS: ${SEQAN_SIMD_HOST_CPU_SUPPORTS}")
        endif()

        unset(CMAKE_REQUIRED_FLAGS)
        unset(SIMD_EXT)
        unset(_SEQANSIMD_SUPPORTED)
    endif()
endmacro()

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
    set(properies "SOURCE_DIR;COMPILE_DEFINITIONS;COMPILE_FEATURES;COMPILE_FLAGS"
                "COMPILE_OPTIONS;LINK_LIBRARIES;LINK_FLAGS;CXX_EXTENSIONS"
                "CXX_STANDARD;CXX_STANDARD_REQUIRED;FOLDER;INCLUDE_DIRECTORIES")

    foreach(property ${properies})
        transfer_target_property(${property} ${source_target} ${target_target})
    endforeach()
endmacro(clone_target)

macro(add_simd_executables target blacklist)
    foreach(simd_ext ${SEQAN_SIMD_SUPPORTED_EXTENSIONS})
        string(TOUPPER ${simd_ext} SIMD_EXT)
        if (NOT (";${blacklist};" MATCHES ";${simd_ext};"))
            # i.e. clone_target(${target} "${target}_avx2")
            clone_target(${target} "${target}_${simd_ext}")
            # i.e. target_compile_options("${target}_avx2" PRIVATE "${SEQAN_SIMD_AVX2_OPTIONS}")
            target_compile_options("${target}_${simd_ext}" PRIVATE "${SEQAN_SIMD_${SIMD_EXT}_OPTIONS}")

            # empty FLAGS means no support for this simd architecture
            if (NOT SEQAN_SIMD_${SIMD_EXT}_OPTIONS)
                message(STATUS "${simd_ext} not supported on ${CMAKE_CXX_COMPILER_ID} (no flags)")
            endif()
        endif()
        unset(SIMD_EXT)
    endforeach()
endmacro(add_simd_executables)

macro(add_simd_tests target blacklist)
    add_test(NAME "test_${target}" COMMAND $<TARGET_FILE:${target}>)

    # simd extensions supported by the host CPU
    foreach(simd_ext ${SEQAN_SIMD_HOST_CPU_SUPPORTS})
        string(TOUPPER ${simd_ext} SIMD_EXT)
        if (TARGET "${target}_${simd_ext}" AND NOT (";${blacklist};" MATCHES ";${simd_ext};"))
            # expands as
            # add_test(NAME "test_test_simd_vector_sse4" COMMAND /usr/bin/sde64 -snb -- $<TARGET_FILE:test_simd_vector_sse4>)
            add_test(NAME "test_${target}_${simd_ext}_host" COMMAND $<TARGET_FILE:${target}_${simd_ext}>)
        endif()
    endforeach()

    if (SDE_FOUND)
        foreach(simd_ext ${SEQAN_SIMD_SUPPORTED_EXTENSIONS})
            string(TOUPPER ${simd_ext} SIMD_EXT)
            if (TARGET "${target}_${simd_ext}" AND NOT (";${blacklist};" MATCHES ";${simd_ext};"))
                # expands as
                # add_test(NAME "test_test_simd_vector_sse4" COMMAND /usr/bin/sde64 -snb -- $<TARGET_FILE:test_simd_vector_sse4>)
                add_test(NAME "test_${target}_${simd_ext}_sde" COMMAND ${SDE_EXECUTABLE} ${SEQAN_SIMD_${SIMD_EXT}_SDE_OPTIONS} -- $<TARGET_FILE:${target}_${simd_ext}>)
            endif()
        endforeach()
    else ()
        message (STATUS "Intel Software Development Emulator not found, not building platform emulated tests.")
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

    if (COMPILER_CLANG)
        # clang 4.x
        if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0.2) AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0) AND (";${SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD};" MATCHES ";avx512_skx;"))
            message(AUTHOR_WARNING "Clang 4.x; reevaluate if AVX512_skx (seqan-simd only) binaries are working. "
                    "An earlier version had an Internal Compiler Error (https://llvm.org/bugs/show_bug.cgi?id=31731), "
                    "which was fixed, the produced binaries might work now. (clang 5.0 is known to work)")
        elseif ((CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3.9) AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0.2))
            simd_list_version_greater(seqansimd_compile_blacklist avx512_knl)
            set(reason_for_disabled_test "Clang 4.0.x produces executables that fail the basic vector test `test_simd_vector`")
        # clang =3.9.0
        elseif (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9) AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9.1)
            # Build the executables, but don't execute them, because clang <= 3.9.x
            # produces executables which contain invalid instructions for AVX512.
            simd_list_version_greater(umesimd_test_blacklist avx2)
            simd_list_version_greater(seqansimd_test_blacklist avx2)
            set(reason_for_disabled_test "Clang 3.9.0 produces executables that contain invalid instructions for at least AVX512_knl and AVX512_skx")
        # clang =3.9.1
        elseif (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9.1) AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.9.2)
            simd_list_version_greater(seqansimd_test_blacklist avx512_knl)
            set(reason_for_disabled_test "Clang 3.9.1 produces executables that contain invalid instructions for at least AVX512_skx (seqan-simd only), see https://llvm.org/bugs/show_bug.cgi?id=31731")
        # clang >=3.9.2
        elseif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3.9.1 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0.0)
            message(AUTHOR_WARNING "Clang >=3.9.2 reevaluate if AVX512_skx (seqan-simd only) binaries are working.")
        # clang 3.7.x
        elseif(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.7) AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8))
            ## seqan-simd:
            # avx2: clang 3.7.x fails test cases:
            #       SimdVectorTestCommon_ShuffleConstant1 type parameter (un-)signed char __vector(32) FAILED
            #       SimdVectorTestCommon_ShuffleConstant1 type parameter (un-)signed short __vector(16) FAILED
            simd_list_version_greater(seqansimd_test_blacklist sse4)
            set(reason_for_disabled_test "Clang 3.7.x produces executables that fail the basic vector test `test_simd_vector`")
        endif()
    endif()

    if(COMPILER_GCC AND DEFINED ENV{TRAVIS} AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0.0)
        # Disable avx2,avx512_knl,... on travis, because sometimes gcc 5.x and 6.x crashes
        simd_list_version_greater(_travis_compile_blacklist sse4)
        set(seqansimd_compile_blacklist ${seqansimd_compile_blacklist} ${_travis_compile_blacklist})
        set(umesimd_compile_blacklist ${umesimd_compile_blacklist} ${_travis_compile_blacklist})
        message(STATUS "Don't compile ${seqansimd_compile_blacklist} on travis, because gcc<=8.x crashes occasionally")
    endif()

    # detect simd support by try-compile and try-run
    detect_simd_support()

    # exclude all simd extensions where a simple try-compile failed
    list_complement(_simd_intrinsics_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}" "${SEQAN_SIMD_COMPILER_SUPPORTS_INTRINSICS}")
    list_complement(_simd_seqansimd_blacklist "${SEQAN_SIMD_SUPPORTED_EXTENSIONS}" "${SEQAN_SIMD_COMPILER_SUPPORTS_SEQANSIMD}")
    set(seqansimd_compile_blacklist ${seqansimd_compile_blacklist} ${_simd_intrinsics_blacklist} ${_simd_seqansimd_blacklist})
    set(umesimd_compile_blacklist ${umesimd_compile_blacklist} ${_simd_intrinsics_blacklist})

    if (UMESIMD_FOUND AND umesimd_test_blacklist AND seqansimd_test_blacklist)
        message(STATUS "Disable test `${target}` for seqan-simd and ume-simd and the following simd extensions ${seqansimd_test_blacklist} and ${umesimd_test_blacklist}")
        message(STATUS "\tReason: ${reason_for_disabled_test}")
    elseif (UMESIMD_FOUND AND umesimd_test_blacklist)
        message(STATUS "Disable test `${target}` for ume-simd and the following simd extensions ${umesimd_test_blacklist}")
        message(STATUS "\tReason: ${reason_for_disabled_test}")
    elseif (seqansimd_test_blacklist)
        message(STATUS "Disable test `${target}` for seqan-simd and the following simd extensions ${seqansimd_test_blacklist}")
        message(STATUS "\tReason: ${reason_for_disabled_test}")
    endif()

    if (SEQAN_SIMD_UTILITY_VERBOSE)
        message(STATUS "seqansimd_compile_blacklist: ${seqansimd_compile_blacklist}")
        message(STATUS "umesimd_compile_blacklist: ${umesimd_compile_blacklist}")
        message(STATUS "seqansimd_test_blacklist: ${seqansimd_test_blacklist}")
        message(STATUS "umesimd_test_blacklist: ${umesimd_test_blacklist}")
    endif()

    add_simd_executables("${target}" "${seqansimd_compile_blacklist}")
    add_simd_tests("${target}" "${seqansimd_test_blacklist}")

    if (UMESIMD_FOUND)
        add_simd_executables("${target}_umesimd" "${umesimd_compile_blacklist}")
        add_simd_tests("${target}_umesimd" "${umesimd_test_blacklist}")
    endif()
endmacro(add_simd_platform_tests)
