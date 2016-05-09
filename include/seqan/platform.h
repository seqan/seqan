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

// ==========================================================================
// Define Used STD Library (e.g. the STL of the GNU compiler)
// ==========================================================================

#include <cstddef> // makes __GLIBCXX__ available
#include <ciso646> // makes _LIBCPP_VERSION available
#define STD_LIB_VS   defined(_MSC_VER)
#define STD_LIB_GNU  defined(__GLIBCXX__)
#define STD_LIB_LLVM defined(_LIBCPP_VERSION)

// ==========================================================================
// Define Compilers
// ==========================================================================

#define COMPILER_MSVC  defined(_MSC_VER) && !defined(__ICC) && !defined(__clang__)
#define COMPILER_GCC   defined(__GNUC__) && !defined(__ICC) && !defined(__clang__)
#define COMPILER_INTEL defined(__ICC)
#define COMPILER_CLANG defined(__clang__)

// ==========================================================================
// Define Integers
// ==========================================================================

/*!
 * @defgroup StandardIntegers Standard Integers
 * @brief Integers defined globally by the SeqAn library.
 *
 * For protability, SeqAn defines the integers in this group.
 *
 * @typedef StandardIntegers#__int64
 * @headerfile <seqan/platform.h>
 * @brief Signed 64-bit integer type.
 * @deprecated Use int64_t instead.
 *
 * @signature typedef (...) __int64;
 *
 * @typedef StandardIntegers#__uint64
 * @headerfile <seqan/platform.h>
 * @brief Unsigned 64-bit integer type.
 * @deprecated Use uint64_t instead.
 *
 * @signature typdef (...) __uint64;
 *
 * @typedef StandardIntegers#__int32
 * @headerfile <seqan/platform.h>
 * @brief Signed 32-bit integer type.
 * @deprecated Use int32_t instead.
 *
 * @signature typedef (...) __int32;
 *
 * @typedef StandardIntegers#__uint32
 * @headerfile <seqan/platform.h>
 * @brief Unsigned 32-bit integer type.
 * @deprecated Use uint32_t instead.
 *
 * @signature typdef (...) __uint32;
 *
 * @typedef StandardIntegers#__int16
 * @headerfile <seqan/platform.h>
 * @brief Signed 16-bit integer type.
 * @deprecated Use int16_t instead.
 *
 * @signature typedef (...) __int16;
 *
 * @typedef StandardIntegers#__uint16
 * @headerfile <seqan/platform.h>
 * @brief Unsigned 16-bit integer type.
 * @deprecated Use uint16_t instead.
 *
 * @signature typdef (...) __uint16;
 *
 * @typedef StandardIntegers#__int8
 * @headerfile <seqan/platform.h>
 * @brief Signed 8-bit integer type.
 * @deprecated Use int8_t instead.
 *
 * @signature typedef (...) __int8;
 *
 * @typedef StandardIntegers#__uint8
 * @headerfile <seqan/platform.h>
 * @brief Unsigned 8-bit integer type.
 * @deprecated Use uint8_t instead.
 *
 * @signature typdef (...) __uint8;
 */

typedef uint64_t __uint64; // nolint
typedef uint32_t __uint32; // nolint
typedef uint16_t __uint16; // nolint
typedef uint8_t __uint8;   // nolint

#if !(COMPILER_INTEL || COMPILER_MSVC)
typedef int64_t __int64;   // nolint
typedef int32_t __int32;   // nolint
typedef int16_t __int16;   // nolint
typedef int8_t __int8;     // nolint
#endif

#if !COMPILER_MSVC
#define finline __inline__
#else // !COMPILER_MSVC
#define finline __forceinline
#endif

// TODO(marehr): always define _FILE_OFFSET_BITS and _LARGEFILE_SOURCE
// if msvc doesn't supprt those flags, why not define them anyway.
#if !COMPILER_MSVC
    #ifndef _FILE_OFFSET_BITS
    #define _FILE_OFFSET_BITS 64
    #endif

    #ifndef _LARGEFILE_SOURCE
    #define _LARGEFILE_SOURCE
    #endif
#endif // !COMPILER_MSVC

/*!
 * @macro SEQAN_IS_64_BIT
 * @headerfile <seqan/platform.h>
 * @brief 1 if the architecture is 64 bit, 0 otherwise.
 *
 * @signature #define SEQAN_IS_64_BIT
 *
 * @macro SEQAN_IS_32_BIT
 * @headerfile <seqan/platform.h>
 * @brief 1 if the architecture is 32 bit, 0 otherwise.
 *
 * @signature #define SEQAN_IS_32_BIT
 */

// The symbols SEQAN_IS_64_BIT and SEQAN_IS_32_BIT can be used to check
// whether we are on a 32 bit or on a 64 bit machine.
#define SEQAN_IS_64_BIT defined(__amd64__) || defined(__x86_64__) || defined(__aarch64__) || defined(__ia64__) || defined(__ppc64__) || defined(_WIN64)
#define SEQAN_IS_32_BIT !SEQAN_IS_64_BIT

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

/**
 * The following macros need to be removed when cuda support will be dropped.
 */
#define SEQAN_FUNC inline
#define SEQAN_HOST_DEVICE
#define SEQAN_HOST
#define SEQAN_DEVICE
#define SEQAN_GLOBAL

// ==========================================================================
// C++ restrict keyword
// ==========================================================================
#if COMPILER_GCC || COMPILER_CLANG
#define SEQAN_RESTRICT  __restrict__
#else
#define SEQAN_RESTRICT
#endif

// ==========================================================================
// C++ branch hints
// ==========================================================================
#if COMPILER_GCC || COMPILER_CLANG || COMPILER_INTEL
#define SEQAN_LIKELY(expr) __builtin_expect(!!(expr), 1)
#else
#define SEQAN_LIKELY(x)    (x)
#endif

#if COMPILER_GCC || COMPILER_CLANG || COMPILER_INTEL
#define SEQAN_UNLIKELY(expr) __builtin_expect(!!(expr), 1)
#else
#define SEQAN_UNLIKELY(x)    (x)
#endif

// A macro to eliminate warnings on GCC and Clang
#if COMPILER_GCC || COMPILER_CLANG || COMPILER_INTEL
#define SEQAN_UNUSED __attribute__((unused))
#else
#define SEQAN_UNUSED
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
    #define SEQAN_HAS_EXECINFO !!__has_include(<execinfo.h>)
    #else // assume that it is there
    #define SEQAN_HAS_EXECINFO 1
    #endif
#endif

// ASYNCHRONOUS I/O
#ifndef SEQAN_ASYNC_IO
    // FreeBSD only has proper support on 64Bit
    #if defined(__FreeBSD__)
    #define SEQAN_ASYNC_IO SEQAN_IS_64_BIT
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
