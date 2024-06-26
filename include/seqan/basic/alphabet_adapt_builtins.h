// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Adaptions of builting types such as bool, int, but also "builtin-level"
// user defined types such as wchar_t, int64_t, uint64_t to the alphabet
// concepts they are in.
// ==========================================================================

#ifndef SEQAN_INCLUDE_BASIC_ALPHABET_ADAPT_BUILTINS_H_
#define SEQAN_INCLUDE_BASIC_ALPHABET_ADAPT_BUILTINS_H_

#include <limits>

namespace seqan2 {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunctions MaxValue_, MinValue_
// ----------------------------------------------------------------------------

// We would want to have this here, however this is not possible with the
// current implementation.

// ----------------------------------------------------------------------------
// Metafunction BitsPerValue
// ----------------------------------------------------------------------------

template <>
struct BitsPerValue<bool>
{
    typedef int Type;
    enum { VALUE = 1 };
};

// ----------------------------------------------------------------------------
// Metafunction IsCharType
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should probably become a concept.

/*!
 * @mfn IsCharType
 * @headerfile <seqan/basic.h>
 *
 * @brief Return whether the argument is <tt>char</tt>, <tt>wchar_t</tt>, <tt>char const</tt>, or <tt>wchar_t
 *               const</tt>.
 *
 * @signature IsCharType<T>::Type;
 * @signature IsCharType<T>::VALUE;
 *
 * @tparam T Type to check type of.
 *
 * This metafunction is used to enable and disable templated adaptions of arrays to sequences for builtin character
 * types only.
 *
 * The return value is <tt>True</tt>/<tt>true</tt> for <tt>char</tt>, <tt>wchar_t</tt>, <tt>char const</tt>, and
 * <tt>wchar_t const</tt>.
 */

template <typename T>
struct IsCharType;

template <typename T>
struct IsCharType
{
    typedef False Type;
    enum { VALUE = 0 };
};

template <typename T>
struct IsCharType<T const>
    : IsCharType<T> {};

template <>
struct IsCharType<char>
{
    typedef True Type;
    enum { VALUE = 1 };
};

template <>
struct IsCharType<wchar_t>
{
    typedef True Type;
    enum { VALUE = 1 };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function gapValueImpl()                                               [char]
// ----------------------------------------------------------------------------

inline char const &
gapValueImpl(char *)
{
    static char const _gap = '-';
    return _gap;
}

inline char const &
gapValueImpl(char const *)
{
    static char const _gap = '-';
    return _gap;
}

// ----------------------------------------------------------------------------
// Function unknownValueImpl()                                           [char]
// ----------------------------------------------------------------------------

inline char const &
unknownValueImpl(char *)
{
    static char const _unknown = 'N';
    return _unknown;
}

inline char const &
unknownValueImpl(char const *)
{
    static char const _unknown = 'N';
    return _unknown;
}

// ----------------------------------------------------------------------------
// Function supremumValueImpl()
// ----------------------------------------------------------------------------

template <typename T>
[[deprecated("Use std::numeric_limits<T>::max() instead.")]]
inline T const &
supremumValueImpl(T *)
{
    static T const x = std::numeric_limits<T>::max();
    return x;
}

[[deprecated("Use std::numeric_limits<T>::max() instead.")]]
inline long double const &
supremumValueImpl(long double *)
{
    static long double const _value = std::numeric_limits<long double>::infinity( );
    return _value;
}

[[deprecated("Use std::numeric_limits<T>::max() instead.")]]
inline double const &
supremumValueImpl(double *)
{
    static double const _value = std::numeric_limits<double>::infinity( );
    return _value;
}

[[deprecated("Use std::numeric_limits<T>::max() instead.")]]
inline float const &
supremumValueImpl(float *)
{
    static float const _value = std::numeric_limits<float>::infinity( );
    return _value;
}

// ----------------------------------------------------------------------------
// Function infimumValueImpl()
// ----------------------------------------------------------------------------

template <typename T>
[[deprecated("Use std::numeric_limits<T>::min() instead.")]]
inline T const &
infimumValueImpl(T *)
{
    static T const x = std::numeric_limits<T>::min();
    return x;
}

[[deprecated("Use std::numeric_limits<T>::min() instead.")]]
inline float const &
infimumValueImpl(float *)
{
    static float const _value = -std::numeric_limits<float>::infinity( );
    return _value;
}

[[deprecated("Use std::numeric_limits<T>::min() instead.")]]
inline double const &
infimumValueImpl(double *)
{
    static double const _value = -std::numeric_limits<double>::infinity( );
    return _value;
}

[[deprecated("Use std::numeric_limits<T>::min() instead.")]]
inline long double const &
infimumValueImpl(long double *)
{
    static long double const _value = -std::numeric_limits<long double>::infinity( );
    return _value;
}

}  // namespace seqan2

#endif  // #ifndef SEQAN_INCLUDE_BASIC_ALPHABET_ADAPT_BUILTINS_H_
