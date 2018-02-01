// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// ==========================================================================
// generic SIMD interface for SSE3 / AVX2
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_H_

#include <seqan/basic.h>

// Check if more than simd back end is selected
#if SEQAN_SEQANSIMD_ENABLED && SEQAN_UMESIMD_ENABLED
    #error UME::SIMD and SEQAN::SIMD are both enabled, you can only use one SIMD back end.
#endif

// MSVC doesn't define SSE4 macros, even if instruction set is available (e.g.
// when AVX is defined)
#if defined(COMPILER_MSVC) && defined(__AVX__) && !defined(__SSE4_1__) && !defined(__SSE4_2__)
    #define __SSE4_1__ 1
    #define __SSE4_2__ 1
#endif

// Define global macro to check if simd instructions are enabled.
#if defined(__AVX512F__) || defined(__AVX2__) || (defined(__SSE4_1__) && defined(__SSE4_2__))
    #define SEQAN_SIMD_ENABLED
#else
    #undef SEQAN_SIMD_ENABLED
    #undef SEQAN_SEQANSIMD_ENABLED
    #undef SEQAN_UMESIMD_ENABLED
#endif

// SIMD operations make only sense on modern 64bit architectures.
#if SEQAN_IS_32_BIT
    // TODO(marehr): If we switch to jenkins, filter out these warnings
    #if !(defined(NDEBUG) || defined(SEQAN_ENABLE_TESTING))
        #pragma message("SIMD acceleration is only available on 64bit systems")
    #endif
    #undef SEQAN_SIMD_ENABLED
    #undef SEQAN_SEQANSIMD_ENABLED
    #undef SEQAN_UMESIMD_ENABLED
#endif // SEQAN_IS_32_BIT

// Fallback to seqan's simd implementation if nothing was specified.
#if defined(SEQAN_SIMD_ENABLED) && !defined(SEQAN_UMESIMD_ENABLED) && !defined(SEQAN_SEQANSIMD_ENABLED)
    #define SEQAN_SEQANSIMD_ENABLED
#endif

// Seqan's simd implementation supports only avx2 and sse4.
#if defined(SEQAN_SEQANSIMD_ENABLED) && !(defined(__AVX2__) || (defined(__SSE4_1__) && defined(__SSE4_2__)))
    #undef SEQAN_SIMD_ENABLED
    #undef SEQAN_SEQANSIMD_ENABLED
#endif

#if defined(SEQAN_SEQANSIMD_ENABLED) && (defined(COMPILER_MSVC) || defined(COMPILER_WINTEL))
    #error SEQAN::SIMD (vector extension) is not supported by msvc and windows intel compiler, try compiling with -DSEQAN_UMESIMD_ENABLED
#endif

// SIMD operations have severe performance issues on <= gcc4.9
#if defined(SEQAN_SEQANSIMD_ENABLED) && defined(COMPILER_GCC) && (__GNUC__ <= 4)
    // TODO(marehr): If we switch to jenkins, filter out these warnings
    #if !(defined(NDEBUG) || defined(SEQAN_ENABLE_TESTING))
        #pragma message("SIMD acceleration was disabled for <=gcc4.9, because of known performance issues " \
                        "https://github.com/seqan/seqan/issues/2017. Use a more recent gcc compiler.")
    #endif
    #undef SEQAN_SIMD_ENABLED
    #undef SEQAN_SEQANSIMD_ENABLED
    #undef SEQAN_UMESIMD_ENABLED
#endif // defined(COMPILER_GCC) && __GNUC__ <= 4

// Define maximal size of vector in byte.
#if defined(SEQAN_SEQANSIMD_ENABLED) && defined(__AVX512F__) && defined(COMPILER_GCC)
    // gcc compiler supports auto vectorization
    #define SEQAN_SIZEOF_MAX_VECTOR 64
#elif defined(SEQAN_SEQANSIMD_ENABLED) && defined(__AVX512F__)
    // TODO(marehr): If we switch to jenkins, filter out these warnings
    #if !(defined(NDEBUG) || defined(SEQAN_ENABLE_TESTING))
        #pragma message("SEQAN_SIMD doesn't support AVX512 (except gcc), thus falling back to AVX2 " \
                        "(we are using some back ported instruction for AVX2 which where introduced since AVX512)")
    #endif
    #define SEQAN_SIZEOF_MAX_VECTOR 32
#elif defined(__AVX512F__)
    #define SEQAN_SIZEOF_MAX_VECTOR 64
#elif defined(__AVX2__)
    #define SEQAN_SIZEOF_MAX_VECTOR 32
#elif defined(__SSE4_1__) && defined(__SSE4_2__)
    #define SEQAN_SIZEOF_MAX_VECTOR 16
#endif

#include "simd/simd_base.h"
#include "simd/simd_base_seqan_impl.h"

#if defined(SEQAN_SEQANSIMD_ENABLED)
    #if SEQAN_SIZEOF_MAX_VECTOR >= 16
    #include "simd/simd_base_seqan_impl_sse4.2.h"
    #endif // SEQAN_SIZEOF_MAX_VECTOR >= 16

    #if SEQAN_SIZEOF_MAX_VECTOR >= 32
    #include "simd/simd_base_seqan_impl_avx2.h"
    #endif // SEQAN_SIZEOF_MAX_VECTOR >= 32

    #if SEQAN_SIZEOF_MAX_VECTOR >= 64
    #include "simd/simd_base_seqan_impl_avx512.h"
    #endif // SEQAN_SIZEOF_MAX_VECTOR >= 64

    #include "simd/simd_base_seqan_interface.h"
#endif // defined(SEQAN_SEQANSIMD_ENABLED)

#if defined(SEQAN_UMESIMD_ENABLED)
    #include "simd/simd_base_umesimd_impl.h"
#endif

#endif // SEQAN_INCLUDE_SEQAN_SIMD_H_
