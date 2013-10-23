// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Implementation of appendLexicalValue().
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_STREAM_STREAM_PUT_H_
#define CORE_INCLUDE_SEQAN_STREAM_STREAM_PUT_H_

namespace seqan {

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
// Metafunction IntegerFormatString_
// ----------------------------------------------------------------------------

// Return the format string for numbers.

template <typename TUnsigned, unsigned SIZE, typename T = void>
struct IntegerFormatString_;


template <typename T>
struct IntegerFormatString_<False, 1, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 1, T>::VALUE[] = "%hhi";


template <typename T>
struct IntegerFormatString_<True, 1, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 1, T>::VALUE[] = "%hhu";


template <typename T>
struct IntegerFormatString_<False, 2, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 2, T>::VALUE[] = "%hi";


template <typename T>
struct IntegerFormatString_<True, 2, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 2, T>::VALUE[] = "%hu";


template <typename T>
struct IntegerFormatString_<False, 4, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 4, T>::VALUE[] = "%i";


template <typename T>
struct IntegerFormatString_<True, 4, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 4, T>::VALUE[] = "%u";


template <typename T>
struct IntegerFormatString_<False, 8, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<False, 8, T>::VALUE[] = "%lli";


template <typename T>
struct IntegerFormatString_<True, 8, T>
{
    static const char VALUE[];
};
template <typename T>
const char IntegerFormatString_<True, 8, T>::VALUE[] = "%llu";


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function appendLexicalValue()
// ----------------------------------------------------------------------------

// Generic version for integers.

template <typename TTarget, typename TInteger>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TInteger> >, typename Size<TTarget>::Type)
appendLexicalValue(TTarget & target, TInteger i)
{
    // 1 byte has at most 3 decimal digits (plus 1 for the NULL character)
    char buffer[sizeof(TInteger) * 3 + 2];
    size_t len = snprintf(buffer, sizeof(buffer), IntegerFormatString_<typename Is<UnsignedIntegerConcept<TInteger> >::Type,
                          sizeof(TInteger)>::VALUE, i);
    write3(target, toRange(buffer + 0, buffer + len));
    return len;
}

// Special case for floats and doubles.

template <typename TTarget>
inline typename Size<TTarget>::Type
appendLexicalValue(TTarget & target, float source)
{
    char buffer[32];
    size_t len = snprintf(buffer, 32, "%g", source);
    write3(target, toRange(buffer + 0, buffer + len));
    return len;
}

template <typename TTarget>
inline typename Size<TTarget>::Type
appendLexicalValue(TTarget & target, double source)
{
    char buffer[32];
    size_t len = snprintf(buffer, 32, "%g", source);
    write3(target, toRange(buffer + 0, buffer + len));
    return len;
}

// NOTE(esiragusa): should be removed.
template <typename TTarget>
inline typename Size<TTarget>::Type
appendLexicalValue(TTarget & target, char c)
{
    writeValue(target, c);
    return 1;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_STREAM_STREAM_PUT_H_
