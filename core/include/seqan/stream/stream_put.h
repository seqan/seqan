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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Implementation of streamPut().
//
// We put this into its own header since the generic version requires some
// metafunctions and more involved code.
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
// Function streamPut()
// ----------------------------------------------------------------------------

// Forward for generic case.

template <typename TStream, typename TSource>
inline int
streamPut(TStream & stream, TSource const & source);

// Important special case of char.

template <typename TStreamSpec>
inline int
streamPut(Stream<TStreamSpec> & stream, char const c)
{
    return streamWriteChar(stream, c);
}

// Important special case of ::std::string.

template <typename TStreamSpec>
inline int
streamPut(Stream<TStreamSpec> & stream, ::std::string const & source)
{
    return (streamWriteBlock(stream, source.data(), source.length()) == source.length())  ?   0 : 1;
}

// Important special case of ::std::stringstream.

template <typename TStreamSpec>
inline int
streamPut(Stream<TStreamSpec> & stream, ::std::stringstream const & source)
{
    return streamPut(stream, source.str());
}

// Important special case of CharString.

template <typename TStreamSpec, typename TSpec>
inline int
streamPut(Stream<TStreamSpec> & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source)) == length(source))  ?   0 : 1;
}

// Special case for floats and doubles.

template <typename TStream>
int _streamPut(TStream & stream, float source, False const & /*tag*/)
{
    char buffer[32];
    size_t len = snprintf(buffer, 32, "%g", source);
    return (streamWriteBlock(stream, buffer, len) == len) ? 0 : 1;
}

template <typename TStream>
int _streamPut(TStream & stream, double source, False const & /*tag*/)
{
    char buffer[32];
    size_t len = snprintf(buffer, 32, "%g", source);
    return (streamWriteBlock(stream, buffer, len) == len) ? 0 : 1;
}

// Generic version for integers.

template <typename TStream, typename TInt>
SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TInt> >, int)
_streamPut(TStream & stream, TInt i, False const & /*tag*/)
{
    // 1 byte has at most 3 decimal digits (plus 1 for the NULL character)
    char buf[sizeof(TInt) * 3 + 2];
    size_t len = snprintf(buf, sizeof(buf), IntegerFormatString_<typename Is<UnsignedIntegerConcept<TInt> >::Type,
                          sizeof(TInt)>::VALUE, i);
    return (streamWriteBlock(stream, buf, len) == len) ? 0 : 1;
}

template <typename TStream>
int _streamPut(TStream & stream, unsigned long i, False const & /*tag*/)
{
    // 1 byte has at most 3 decimal digits (plus 1 for the NULL character)
    char buf[sizeof(unsigned long) * 3 + 2];
    size_t len = snprintf(buf, sizeof(buf), "%lu", i);
    return (streamWriteBlock(stream, buf, len) == len) ? 0 : 1;
}

template <typename TStream>
int _streamPut(TStream & stream, long i, False const & /*tag*/)
{
    // 1 byte has at most 3 decimal digits (plus 1 for the NULL character)
    char buf[sizeof(long) * 3 + 2];
    size_t len = snprintf(buf, sizeof(buf), "%li", i);
    return (streamWriteBlock(stream, buf, len) == len) ? 0 : 1;
}

template <typename TStream>
int _streamPut(TStream & stream, char c, False const & /*tag*/)
{
    return streamWriteChar(stream, c);
}

// Generic fallback version, based on stringstream.
//
// Parameters:
// (1) stream The stream to write to.
// (2) source The value to write to stream.
// (3) tag    True if the type of source is a sequence and False otherwise.

template <typename TStream, typename TSource>
SEQAN_FUNC_DISABLE_IF(Is<IntegerConcept<TSource> >, int)
_streamPut(TStream & stream, TSource const & source, False const & /*tag*/)
{
    ::std::stringstream s;

    s << source;
    if (s.fail())
        return s.fail();

    return streamPut(stream, s);
}

template <typename TStream, typename TSource>
inline int
_streamPut(TStream & target, TSource const & source, True const & /*tag*/)
{
    typename Iterator<TSource const, Standard>::Type it = begin(source, Standard());
    typename Iterator<TSource const, Standard>::Type itEnd = end(source, Standard());
    int res = 0;

    for (; it != itEnd && res == 0; ++it)
        res = streamPut(target, getValue(it));

    return res;
}

// Case: Character arrays.

template <typename TStream>
inline int
_streamPut(TStream & stream, char const * source, True const & /*tag*/)
{
    return (streamWriteBlock(stream, source, strlen(source)) == strlen(source)) ? 0 : 1;
}

// Case: Array.
// TODO(holtgrew): Requires atEnd(it) <==> *it == 0. Remove?

template <typename TStream, typename TSourceValue>
inline int
_streamPut(TStream & stream, TSourceValue const * source, True const & /*tag*/)
{
    int res = 0;
    for (; !atEnd(source) && res == 0; ++source)
        res = _streamWrite(stream, *source);
    return res;
}

// Function entry for generic version.

template <typename TStream, typename TSource>
inline int
streamPut(TStream & stream, TSource const & source)
{
    return _streamPut(stream, source, typename IsSequence<TSource const>::Type());
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_STREAM_STREAM_PUT_H_
