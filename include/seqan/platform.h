// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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

#ifndef SEQAN_PLATFORM_H
#define SEQAN_PLATFORM_H

#include <cinttypes>

#ifdef _MSC_VER
    #include "platform/platform_windows.h"
#elif __ICC
    #include "platform/platform_icc.h"
#else
    #include "platform/platform_gcc.h"
#endif

// NOTE(esiragusa): nvcc header must be included even if __CUDACC__ is not defined.
#include "platform/platform_nvcc.h"

/**
 * SEQAN_AUTO_PTR_NAME .... alias for the auto_ptr class template deprecated in C++11.
 *                          @deprecated use the RHS
 *
 * SEQAN_FORWARD_ARG and
 * SEQAN_FORWARD_CARG ..... macros to insert between argument type and name ...
 *                          @deprecated use the RHS
 *
 * SEQAN_FORWARD_RETURN ... or return type and function name to declare forwarding of variables
 *                          @deprecated use the RHS
 *
 * SEQAN_FORWARD .......... pass a variable as (of type T) as it was given to a function
 *                          @deprecated use the RHS
 *
 * SEQAN_MOVE ............. pass a variable to a function and never use it again
 *                          @deprecated use the RHS
 */

#define SEQAN_AUTO_PTR_NAME     unique_ptr
#define SEQAN_FORWARD_ARG       &&
#define SEQAN_FORWARD_CARG      &&
#define SEQAN_FORWARD_RETURN    &&
#define SEQAN_FORWARD(T, x)     std::forward<T>(x)
#define SEQAN_MOVE(x)           std::move(x)

// backwards compatibility
#define SEQAN_CXX11_STL 1
#define SEQAN_CXX11_STANDARD 1
#define SEQAN_CXX11_COMPLETE 1

// C++ restrict keyword, see e.g. platform_gcc.h
#ifndef SEQAN_RESTRICT
#define SEQAN_RESTRICT
#endif

// C++ branch hints
#ifndef SEQAN_LIKELY
#define SEQAN_LIKELY(x) (x)
#endif

#ifndef SEQAN_UNLIKELY
#define SEQAN_UNLIKELY(x) (x)
#endif

// A macro to eliminate warnings on GCC and Clang
#if defined(__GNUC__) || defined(__clang__)
#  define SEQAN_UNUSED __attribute__((unused))
#else
#  define SEQAN_UNUSED
#endif
// backwards compatibility
#define SEQAN_UNUSED_TYPEDEF SEQAN_UNUSED

// HAS_EXECINFO
// note that this is always set by FindSeqAn.cmake
// this is a fallback for non cmake environments
#ifndef SEQAN_HAS_EXECINFO
#ifdef PLATFORM_WINDOWS
    #define SEQAN_HAS_EXECINFO 0
#elif defined(__has_include)
    #if __has_include(<execinfo.h>)
        #define SEQAN_HAS_EXECINFO 1
    #else
        #define SEQAN_HAS_EXECINFO 0
    #endif
#else // assume that it is there
    #define SEQAN_HAS_EXECINFO 1
#endif
#endif

// ASYNCHRONOUS I/O
#ifndef SEQAN_ASYNC_IO
// FreeBSD only has proper support on 64Bit
#if defined(__FreeBSD__)
    #if defined(__x86_64__) || defined(__aarch64__) || defined(__ia64__) || defined(__ppc64__)
        #define SEQAN_ASYNC_IO 1
    #else
        #define SEQAN_ASYNC_IO 0
    #endif
// Clang (and future gcc) can detect it
#elif defined(__has_include) && defined(__unix__)
    #if __has_include(<aio.h>)
        #define SEQAN_ASYNC_IO 1
    #else
        #define SEQAN_ASYNC_IO 0
    #endif
// OpenBSD doesn't have it
#elif defined(__OpenBSD__)
    #define SEQAN_ASYNC_IO 0
// we assume the rest have it (Linux, OSX, Win)
#else
    #define SEQAN_ASYNC_IO 1
#endif
#endif //ndef SEQAN_ASYNC_IO

#endif
