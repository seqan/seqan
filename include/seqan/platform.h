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

// ==========================================================================
// Define Used STD Library (e.g. the STL of the GNU compiler)
// ==========================================================================

#include <cstddef> // makes __GLIBCXX__ available
#include <ciso646> // makes _LIBCPP_VERSION available

/*!
 * @macro STDLIB_VS
 * @headerfile <seqan/platform.h>
 * @brief The standard library implemented by Visual C++, if defined
 * @signature #define STDLIB_VS
 */
#ifdef _MSC_VER
#define STDLIB_VS
#endif

/*!
 * @macro STDLIB_GNU
 * @headerfile <seqan/platform.h>
 * @brief The standard library implemented by GNU/GCC, if defined
 * @signature #define STDLIB_GNU
 */
#ifdef __GLIBCXX__
#define STDLIB_GNU
#endif

/*!
 * @macro STDLIB_LLVM
 * @headerfile <seqan/platform.h>
 * @brief The standard library implemented by clang/llvm, if defined
 * @signature #define STDLIB_LLVM
 */
#ifdef _LIBCPP_VERSION
#define STDLIB_LLVM
#endif

// ==========================================================================
// Define Compilers
// ==========================================================================

/*!
 * @macro COMPILER_LINTEL
 * @headerfile <seqan/platform.h>
 * @brief The compiler is the intel compiler (icc), if defined
 * @signature #define COMPILER_LINTEL
 */
#if defined(__ICC)
#define COMPILER_LINTEL
#if __ICC < 1600
     #warning ICC versions older than 16 are not supported.
#endif
#endif

/*!
 * @macro COMPILER_WINTEL
 * @headerfile <seqan/platform.h>
 * @brief The compiler is the intel compiler for windows, if defined
 * @signature #define COMPILER_WINTEL
 */
#if defined(__ICL)
#define COMPILER_WINTEL
#if __ICL < 1600
     #warning Intel compiler (windows) versions older than 16 are not supported.
#endif
#endif

/*!
 * @macro COMPILER_CLANG
 * @headerfile <seqan/platform.h>
 * @brief The compiler is the llvm compiler (clang), if defined
 * @signature #define COMPILER_CLANG
 */
#if defined(__clang__)
#define COMPILER_CLANG
#define COMPILER_VERSION (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
#if COMPILER_VERSION < 30500
    #warning Clang versions older than 3.5.0 are not supported.
#endif
#undef COMPILER_VERSION
#endif

/*!
 * @macro COMPILER_MSVC
 * @headerfile <seqan/platform.h>
 * @brief The compiler is the microsoft visual studio compiler (msvc), if defined
 * @signature #define COMPILER_MSVC
 */
#if defined(_MSC_VER) && !defined(COMPILER_WINTEL) && !defined(COMPILER_CLANG)
#define COMPILER_MSVC
#if _MSC_VER < 1900
#error Visual Studio versions older than version 14 / "2015" are not supported.
#endif
#endif

/*!
 * @macro COMPILER_GCC
 * @headerfile <seqan/platform.h>
 * @brief The compiler is the gnu compiler (gcc), if defined
 * @signature #define COMPILER_GCC
 */
#if defined(__GNUC__) && !defined(COMPILER_LINTEL) && !defined(COMPILER_CLANG)
#define COMPILER_GCC
#define COMPILER_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if COMPILER_VERSION < 40901
    #warning GCC versions older than 4.9.1 are not supported.
#endif
#undef COMPILER_VERSION
#endif

// ==========================================================================
// C++ standard macros
// ==========================================================================

#ifndef STDLIB_VS // all Visual Studio >= 2015 compilers are c++14 by default
#if __cplusplus < 201300
    #error SeqAn requires C++14! You must compile your application with -std=c++14, -std=gnu++14 or -std=c++1y.
#endif
#endif

// ==========================================================================
// Platform Macros (Backwards Compatibility)
// ==========================================================================
/*!
 * @macro PLATFORM_GCC
 * @headerfile <seqan/platform.h>
 * @brief Defined if the compiler is GCC (or compatible).
 * @deprecated Use STDLIB_VS, STDLIB_GNU or STDLIB_LLVM to know which
 *     standard lib is currently used. Or use COMPILER_MSVC, COMPILER_GCC,
 *     COMPILER_LINTEL, COMPILER_WINTEL or COMPILER_CLANG to know which compiler
 *     is currently used.
 *
 * @signature #define PLATFORM_GCC
 */

#ifdef STDLIB_VS
#define PLATFORM_WINDOWS
#define PLATFORM_WINDOWS_VS
#else
#define PLATFORM_GCC
#endif

#if defined(PLATFORM_GCC) && defined(COMPILER_CLANG)
#define PLATFORM_CLANG
#endif

#if defined(PLATFORM_GCC) && defined(COMPILER_LINTEL)
#define PLATFORM_INTEL
#endif

#if defined(PLATFORM_GCC) && defined(COMPILER_GCC)
#define PLATFORM_GNU
#endif

// ==========================================================================
// Disable Warnings
// ==========================================================================

// Disable warning for identifer name truncation.  There is not much we can
// do about this.  Boost also has this problem and they chose to suppress
// it globally.  So did we.
//
// Documentation of C4503 from Microsoft:
//   https://msdn.microsoft.com/en-us/library/074af4b6%28v=vs.140%29.aspx
// Boost Warnings Guidelines:
//   https://svn.boost.org/trac/boost/wiki/Guidelines/WarningsGuidelines
// Intel compiler for windows also triggers this error:
//   seqan/pipe/pipe_base.h(263): warning #2586: 'bundle5' : decorated name
//   length exceeded, name was truncated
#if defined(COMPILER_MSVC) || defined(COMPILER_WINTEL)
#pragma warning( disable : 4503 )
#endif

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

#if !(defined(COMPILER_LINTEL) || defined(STDLIB_VS))
typedef int64_t __int64;   // nolint
typedef int32_t __int32;   // nolint
typedef int16_t __int16;   // nolint
typedef int8_t __int8;     // nolint
#endif

#if defined(COMPILER_MSVC) || defined(COMPILER_WINTEL)
#define finline __forceinline
#else
#define finline __inline__
#endif

// TODO(marehr): always define _FILE_OFFSET_BITS and _LARGEFILE_SOURCE
// if msvc doesn't supprt those flags, why not define them anyway.
#if !defined(COMPILER_MSVC)
    #ifndef _FILE_OFFSET_BITS
    #define _FILE_OFFSET_BITS 64
    #endif

    #ifndef _LARGEFILE_SOURCE
    #define _LARGEFILE_SOURCE
    #endif
#endif // !defined(COMPILER_MSVC)

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
#if defined(__amd64__) || defined(__x86_64__) || defined(__aarch64__) || defined(__ia64__) || defined(__ppc64__) || defined(_WIN64)
#define SEQAN_IS_64_BIT 1
#define SEQAN_IS_32_BIT 0
#else
#define SEQAN_IS_64_BIT 0
#define SEQAN_IS_32_BIT 1
#endif

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
#if defined(COMPILER_GCC) || defined(COMPILER_CLANG)
#define SEQAN_RESTRICT  __restrict__
#else
#define SEQAN_RESTRICT
#endif

// ==========================================================================
// C++ branch hints
// ==========================================================================
#if defined(COMPILER_GCC) || defined(COMPILER_CLANG) && !defined(STDLIB_VS) || defined(COMPILER_LINTEL)
#define SEQAN_LIKELY(expr) __builtin_expect(!!(expr), 1)
#define SEQAN_UNLIKELY(expr) __builtin_expect(!!(expr), 0)
#else
#define SEQAN_LIKELY(x)    (x)
#define SEQAN_UNLIKELY(x)    (x)
#endif

// A macro to eliminate warnings on GCC and Clang
#if defined(COMPILER_GCC) || defined(COMPILER_CLANG) || defined(COMPILER_LINTEL)
#define SEQAN_UNUSED __attribute__((unused))
#else
#define SEQAN_UNUSED
#endif
// backwards compatibility
#define SEQAN_UNUSED_TYPEDEF SEQAN_UNUSED

// eliminate fallthrough warnings
#define SEQAN_FALLTHROUGH
#if defined(__has_cpp_attribute)
    #if __has_cpp_attribute(fallthrough)
        #undef SEQAN_FALLTHROUGH
        #if __cplusplus < 201500 && defined(COMPILER_GCC)
            #define SEQAN_FALLTHROUGH [[gnu::fallthrough]];
        #elif __cplusplus < 201500 && defined(COMPILER_CLANG)
            #define SEQAN_FALLTHROUGH [[clang::fallthrough]];
        #else
            #define SEQAN_FALLTHROUGH [[fallthrough]];
        #endif
    #endif
#endif

// HAS_EXECINFO
// note that this is always set by seqan-config.cmake
// this is a fallback for non cmake environments
#ifndef SEQAN_HAS_EXECINFO
    #ifdef STDLIB_VS
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

// There is a bug in clang35 on FreeBSD that crashes the compiler.
// Since this is only triggered on some codepaths we define a macro
// here that we can later access.
#if defined(__FreeBSD__) && defined(COMPILER_CLANG)
#define COMPILER_VERSION (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
#if (COMPILER_VERSION >= 30500) && (COMPILER_VERSION < 30600)
    #define SEQAN_CLANG35_FREEBSD_BUG 1
#endif
#undef COMPILER_VERSION
#endif

#endif
