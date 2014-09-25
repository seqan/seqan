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
// Basic definitions for the stream module.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_BASIC_STREAM_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_BASIC_STREAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TDirection>
struct StreamIterator;

template <typename TValue, typename TTraits, typename TValue2>
inline void writeValue(std::basic_ostream<TValue, TTraits> &ostream, TValue2 val);

template <typename TValue, typename TTraits, typename TValue2>
inline void writeValue(std::ostreambuf_iterator<TValue, TTraits> &iter, TValue2 val);

template <typename TContainer, typename TValue>
inline void writeValue(Iter<TContainer, StreamIterator<Output> > &iter, TValue val);

template <typename TValue, typename TTraits>
inline bool atEnd(std::istreambuf_iterator<TValue, TTraits> const &it);

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
length(std::basic_string<TChar, TCharTraits, TAlloc> const & me);

/*!
 * @macro SEQAN_HAS_ZLIB
 * @headerfile <seqan/stream.h>
 * @brief Defined as 0 or 1, depending on zlib being available.
 *
 * @signature #define SEQAN_HAS_ZLIB 0  // or 1
 */

/*!
 * @macro SEQAN_HAS_BZIP2
 * @headerfile <seqan/stream.h>
 * @brief Defined as 0 or 1, depending on bzlib being available.
 *
 * @signature #define SEQAN_HAS_BZIP 0  // or 1
 */

/**
.Macro.SEQAN_HAS_ZLIB
..cat:Input/Output
..cat:From Outside
..signature:SEQAN_HAS_ZLIB
..summary:If set to 1 then zlib is available, i.e. including $<zlib.h>$ and linking against libz works.
..remarks:This flag is normally set from the outside by your build system using compiler flags.

.Macro.SEQAN_HAS_BZIP2
..cat:Input/Output
..cat:From Outside
..signature:SEQAN_HAS_BZLIB
..summary:If set to 1 then bzlib2 is available, i.e. including $<bzlib.h>$ and linking against libbzip2 works.
..remarks:This flag is normally set from the outside by your build system using compiler flags.
 */

// ============================================================================
// Tags
// ============================================================================

// ============================================================================
// Concepts
// ============================================================================

// --------------------------------------------------------------------------
// Concept StreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT(StreamConcept, (TStream))
{
};

// --------------------------------------------------------------------------
// Concept InputStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT_REFINE(InputStreamConcept, (TStream), (StreamConcept))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(InputStreamConcept)
    {}
};

// --------------------------------------------------------------------------
// Concept OutputStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT_REFINE(OutputStreamConcept, (TStream), (StreamConcept))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(OutputStreamConcept)
    {}
};

// --------------------------------------------------------------------------
// Concept BidirectionalStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT_REFINE(BidirectionalStreamConcept, (TStream), (InputStreamConcept)(OutputStreamConcept))
{};

// ============================================================================
// Exceptions
// ============================================================================

// ----------------------------------------------------------------------------
// Exception ParseError
// ----------------------------------------------------------------------------

struct ParseError : RuntimeError
{
    template <typename TString>
    ParseError(TString const &message):
        RuntimeError(message)
    {}
};

// ----------------------------------------------------------------------------
// Exception UnexpectedEnd
// ----------------------------------------------------------------------------

struct UnexpectedEnd : ParseError
{
    UnexpectedEnd():
        ParseError("Unexpected end of input.")
    {}
};

// ----------------------------------------------------------------------------
// Exception EmptyFieldError
// ----------------------------------------------------------------------------

struct EmptyFieldError : ParseError
{
    EmptyFieldError(std::string fieldName):
        ParseError(fieldName + " field was empty.")
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TObject, typename TDirection>
struct DirectionIterator :
    If<Is<StreamConcept<TObject> >,
       Iter<TObject, StreamIterator<TDirection> >,
       typename Iterator<TObject, Rooted>::Type>
{};

// --------------------------------------------------------------------------
// Metafunction BasicStream
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TTraits = std::char_traits<TValue> >
struct BasicStream :
    If<
        IsSameType<TDirection, Input>,
        std::basic_istream<TValue, TTraits>,
        typename If<
            IsSameType<TDirection, Output>,
            std::basic_ostream<TValue, TTraits>,
            std::basic_iostream<TValue, TTraits>
            >::Type
        >
{};

// --------------------------------------------------------------------------
// Metafunction IosOpenMode
// --------------------------------------------------------------------------

template <typename TDirection, typename TDummy = void>
struct IosOpenMode;


template <typename TDummy>
struct IosOpenMode<Input, TDummy>
{
    static const int VALUE;
};

template <typename TDummy>
struct IosOpenMode<Output, TDummy>
{
    static const int VALUE;
};

template <typename TDummy>
struct IosOpenMode<Bidirectional, TDummy>
{
    static const int VALUE;
};


template <typename TDummy>
const int IosOpenMode<Input, TDummy>::VALUE = std::ios::in;

template <typename TDummy>
const int IosOpenMode<Output, TDummy>::VALUE = std::ios::out;

template <typename TDummy>
const int IosOpenMode<Bidirectional, TDummy>::VALUE = std::ios::in | std::ios::out;


// --------------------------------------------------------------------------
// Metafunction MagicHeader
// --------------------------------------------------------------------------

template <typename TTag, typename T = void>
struct MagicHeader;

template <typename T>
struct MagicHeader<Nothing, T>
{
    static char const * VALUE;
};

template <typename T>
char const * MagicHeader<Nothing, T>::VALUE = NULL;

// --------------------------------------------------------------------------
// Metafunction FileFormatExtensions
// --------------------------------------------------------------------------

// TODO(weese:) rename FileFormatExtensions to FileTypeExtensions or FileExtensions
template <typename TFormat, typename T = void>
struct FileFormatExtensions;

template <typename T>
struct FileFormatExtensions<Nothing, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileFormatExtensions<Nothing, T>::VALUE[1] =
{
    ""
};              // default output extension

// ----------------------------------------------------------------------------
// Metafunction IntegerFormatString_
// ----------------------------------------------------------------------------
// Return the format string for numbers.

template <typename TUnsigned, unsigned SIZE, typename T = void>
struct IntegerFormatString_;


template <typename TUnsigned, typename T>
struct IntegerFormatString_<TUnsigned, 1, T> :
    IntegerFormatString_<TUnsigned, 2, T> {};


template <typename T>
struct IntegerFormatString_<False, 2, T>
{
    static const char VALUE[];
    typedef short Type;
};
template <typename T>
const char IntegerFormatString_<False, 2, T>::VALUE[] = "%hi%n";


template <typename T>
struct IntegerFormatString_<True, 2, T>
{
    static const char VALUE[];
    typedef unsigned short Type;
};
template <typename T>
const char IntegerFormatString_<True, 2, T>::VALUE[] = "%hu%n";


template <typename T>
struct IntegerFormatString_<False, 4, T>
{
    static const char VALUE[];
    typedef int Type;
};
template <typename T>
const char IntegerFormatString_<False, 4, T>::VALUE[] = "%i%n";


template <typename T>
struct IntegerFormatString_<True, 4, T>
{
    static const char VALUE[];
    typedef unsigned Type;
};
template <typename T>
const char IntegerFormatString_<True, 4, T>::VALUE[] = "%u%n";


template <typename T>
struct IntegerFormatString_<False, 8, T>
{
    static const char VALUE[];
    typedef __int64 Type;
};
template <typename T>
const char IntegerFormatString_<False, 8, T>::VALUE[] = "%lli%n";


template <typename T>
struct IntegerFormatString_<True, 8, T>
{
    static const char VALUE[];
    typedef __uint64 Type;
};
template <typename T>
const char IntegerFormatString_<True, 8, T>::VALUE[] = "%llu%n";

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeValue()
// ----------------------------------------------------------------------------

// resizable containers
template <typename TSequence, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TSequence> >, void)
writeValue(TSequence &cont, TValue val)
{
    appendValue(cont, val);
}

// Range
template <typename TIterator, typename TValue>
inline void
writeValue(Range<TIterator> &range, TValue val)
{
    assignValue(range.begin, val);
    ++range.begin;
}

// ----------------------------------------------------------------------------
// Function writeValue(Iter)
// ----------------------------------------------------------------------------

// resizable containers
template <typename TSequence, typename TSpec, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TSequence> >, void)
writeValue(Iter<TSequence, TSpec> &iter, TValue val)
{
    typedef Iter<TSequence, TSpec> TIter;

    TSequence &cont = container(iter);
    typename Position<TIter>::Type pos = position(iter);
    typename Size<TIter>::Type len = length(cont);

    if (pos < len)
    {
        assignValue(iter, val);
        ++iter;
    }
    else
    {
        if (pos > len)
            resize(cont, pos - 1);
        appendValue(cont, val);
        setPosition(iter, pos + 1);
    }
}

// non-resizable containers
template <typename TNoSequence, typename TSpec, typename TValue>
inline SEQAN_FUNC_DISABLE_IF(Is<ContainerConcept<TNoSequence> >, void)
writeValue(Iter<TNoSequence, TSpec> &iter, TValue val)
{
    SEQAN_ASSERT_LT(position(iter), length(container(iter)));

    assignValue(iter, val);
    ++iter;
}

template <typename TTargetValue, typename TValue>
inline void
writeValue(TTargetValue * iter, TValue val)
{
    *iter++ = val;
}

// ----------------------------------------------------------------------------
// Function _write(); Element-wise
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize, typename TIChunk, typename TOChunk>
inline void _write(TTarget &target, TFwdIterator &iter, TSize n, TIChunk, TOChunk)
{
    for (; n > (TSize)0; --n, ++iter)
        writeValue(target, getValue(iter));
}

// ----------------------------------------------------------------------------
// Function _write(); Chunked
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize, typename TIValue, typename TOValue>
inline void _write(TTarget &target, TFwdIterator &iter, TSize n, Range<TIValue*> *, Range<TOValue*> *)
{
    typedef Nothing* TNoChunking;
    typedef typename Size<TTarget>::Type TTargetSize;

    Range<TIValue*> ichunk;
    Range<TOValue*> ochunk;

    while (n != (TSize)0)
    {
        getChunk(ichunk, iter, Input());
        getChunk(ochunk, target, Output());

        TTargetSize minChunkSize = std::min((TTargetSize)length(ichunk), (TTargetSize)length(ochunk));

        if (SEQAN_UNLIKELY(minChunkSize == 0u))
        {
            reserveChunk(target, n, Output());
            reserveChunk(iter, n, Input());
            getChunk(ochunk, target, Output());
            getChunk(ichunk, iter, Input());
            minChunkSize = std::min((TTargetSize)length(ichunk), (TTargetSize)length(ochunk));
            if (SEQAN_UNLIKELY(minChunkSize == 0u))
            {
                _write(target, iter, n, TNoChunking(), TNoChunking());
                return;
            }
        }

        if (minChunkSize > (TTargetSize)n)
            minChunkSize = (TTargetSize)n;

        arrayCopyForward(ichunk.begin, ichunk.begin + minChunkSize, ochunk.begin);

        iter += minChunkSize;                      // advance input iterator
        advanceChunk(target, minChunkSize);
        n -= minChunkSize;
    }
}

// chunked, target is pointer (e.g. readRawPod)
template <typename TOValue, typename TFwdIterator, typename TSize>
inline SEQAN_FUNC_DISABLE_IF(IsSameType<typename Chunk<TFwdIterator>::Type, Nothing>, void)
write(TOValue *ptr, TFwdIterator &iter, TSize n)
{
    typedef Nothing* TNoChunking;
    typedef typename Size<TFwdIterator>::Type TSourceSize;
    typedef typename Chunk<TFwdIterator>::Type TIChunk;

    TIChunk ichunk;

    while (n != (TSize)0)
    {
        getChunk(ichunk, iter, Input());
        TSourceSize chunkSize = length(ichunk);

        if (SEQAN_UNLIKELY(chunkSize == 0u))
        {
            reserveChunk(iter, n, Input());
            getChunk(ichunk, iter, Input());
            TSourceSize chunkSize = length(ichunk);
            if (SEQAN_UNLIKELY(chunkSize == 0u))
            {
                _write(ptr, iter, n, TNoChunking(), TNoChunking());
                return;
            }
        }

        if (chunkSize > (TSourceSize)n)
            chunkSize = (TSourceSize)n;

        arrayCopyForward(ichunk.begin, ichunk.begin + chunkSize, ptr);

        iter += chunkSize;                          // advance input iterator
        ptr += chunkSize;
        n -= chunkSize;
    }
}

// chunked, source is pointer (e.g. readRawPod)
template <typename TTarget, typename TIValue, typename TSize>
inline SEQAN_FUNC_DISABLE_IF(IsSameType<typename Chunk<TTarget>::Type, Nothing>, void)
write(TTarget &target, TIValue *ptr, TSize n)
{
    typedef Nothing* TNoChunking;
    typedef typename Size<TTarget>::Type TTargetSize;
    typedef typename Chunk<TTarget>::Type TOChunk;

    TOChunk ochunk;

    while (n != (TSize)0)
    {
        getChunk(ochunk, target, Output());
        TTargetSize chunkSize = length(ochunk);

        if (SEQAN_UNLIKELY(chunkSize == 0u))
        {
            reserveChunk(target, n, Output());
            getChunk(ochunk, target, Output());
            chunkSize = length(ochunk);
            if (SEQAN_UNLIKELY(chunkSize == 0u))
            {
                _write(target, ptr, n, TNoChunking(), TNoChunking());
                return;
            }
        }

        if (chunkSize > (TTargetSize)n)
            chunkSize = (TTargetSize)n;

        arrayCopyForward(ptr, ptr + chunkSize, ochunk.begin);

        ptr += chunkSize;                      // advance input iterator
        advanceChunk(target, chunkSize);
        n -= chunkSize;
    }
}

// ----------------------------------------------------------------------------
// Function write(TValue *)
// ----------------------------------------------------------------------------
// NOTE(esiragusa): should it be defined for Streams and Containers?

//template <typename TTarget, typename TValue, typename TSize>
//inline SEQAN_FUNC_ENABLE_IF(Or<Is<OutputStreamConcept<TTarget> >, Is<ContainerConcept<TTarget> > >, void)
//write(TTarget &target, TValue *ptr, TSize n)
//{
//    typedef Range<TValue*>                          TRange;
//    typedef typename Iterator<TRange, Rooted>::Type TIterator;
//    typedef typename Chunk<TIterator>::Type*        TIChunk;
//    typedef typename Chunk<TTarget>::Type*          TOChunk;
//
//    TRange range(ptr, ptr + n);
//    TIterator iter = begin(range, Rooted());
//    _write(target, iter, n, TIChunk(), TOChunk());
//}

// ----------------------------------------------------------------------------
// Function write(Iterator<Input>)
// ----------------------------------------------------------------------------

//TODO(singer): Enable this!
template <typename TTarget, typename TFwdIterator, typename TSize>
//inline SEQAN_FUNC_ENABLE_IF(Or<Is<OutputStreamConcept<TTarget> >, Is<ContainerConcept<TTarget> > >, void)
inline SEQAN_FUNC_ENABLE_IF(And<
    Is<IntegerConcept<TSize> >,
    Is<Convertible<typename Value<TTarget>::Type,
                   typename Value<TFwdIterator>::Type> > >, void)
write(TTarget &target, TFwdIterator &iter, TSize n)
{
    typedef typename Chunk<TFwdIterator>::Type* TIChunk;
    typedef typename Chunk<TTarget>::Type*      TOChunk;

    _write(target, iter, n, TIChunk(), TOChunk());
}

// write for more complex values (defer to write of iterator value)
// used for Strings of Pairs
template <typename TTarget, typename TFwdIterator, typename TSize>
//inline SEQAN_FUNC_ENABLE_IF(Or<Is<OutputStreamConcept<TTarget> >, Is<ContainerConcept<TTarget> > >, void)
inline SEQAN_FUNC_ENABLE_IF(And<
    Is<IntegerConcept<TSize> >,
    Not< Is<Convertible<typename Value<TTarget>::Type,
                        typename Value<TFwdIterator>::Type> > > >, void)
write(TTarget &target, TFwdIterator &iter, TSize n)
{
    for (; n > (TSize)0; --n, ++iter)
    {
        write(target, *iter);
        writeValue(target, ' ');
    }
}

// ----------------------------------------------------------------------------
// Function write(TContainer)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TContainer> >, void)
write(TTarget &target, TContainer &cont)
{
    typename DirectionIterator<TContainer, Input>::Type iter = directionIterator(cont, Input());
    write(target, iter, length(cont));
}

template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TContainer const> >, void)
write(TTarget &target, TContainer const &cont)
{
    typename DirectionIterator<TContainer const, Input>::Type iter = directionIterator(cont, Input());
    write(target, iter, length(cont));
}

template <typename TTarget, typename TValue>
inline void
write(TTarget &target, TValue * ptr)
{
    write(target, ptr, length(ptr));
}

// ----------------------------------------------------------------------------
// Function appendNumber()
// ----------------------------------------------------------------------------
// Generic version for integers.

template <typename TTarget, typename TInteger>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TInteger> >, typename Size<TTarget>::Type)
appendNumber(TTarget & target, TInteger i)
{
    typedef IntegerFormatString_<typename Is<UnsignedIntegerConcept<TInteger> >::Type,
                                 xsizeof(TInteger)> TInt;

    // 1 byte has at most 3 decimal digits (plus 2 for '-' and the NULL character)
    char buffer[sizeof(TInteger) * 3 + 2];
    int offset;
    size_t len = snprintf(buffer, sizeof(buffer),
                          TInt::VALUE, static_cast<typename TInt::Type>(i), &offset);
    char *bufPtr = buffer;
    write(target, bufPtr, len);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendNumber(bool)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, bool source)
{
    writeValue(target, '0' + source);
    return 1;
}

// ----------------------------------------------------------------------------
// Function appendNumber(float)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, float source)
{
    char buffer[32];
    int offset;
    size_t len = snprintf(buffer, 32, "%g%n", source, &offset);
    char *bufPtr = buffer;
    write(target, bufPtr, len);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendNumber(double)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, double source)
{
    char buffer[32];
    int offset;
    size_t len = snprintf(buffer, 32, "%g%n", source, &offset);
    char *bufPtr = buffer;
    write(target, bufPtr, len);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendRawPod()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTarget>
inline typename Size<TTarget>::Type
appendRawPod(TTarget & target, TValue const & val)
{
    write(target, (unsigned char*)&val, sizeof(TValue));
    return sizeof(TValue);
}

// ----------------------------------------------------------------------------
// Function write(TNumber); write fundamental type
// ----------------------------------------------------------------------------

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue> >,
                                 Is<FundamentalConcept<TValue> > >, void)
write(TTarget &target, TValue &number)
{
    if (sizeof(TValue) == 1)
        writeValue(target, number);     // write chars as chars
    else
        appendNumber(target, number);
}

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue> >,
                                 Is<FundamentalConcept<TValue> > >, void)
write(TTarget &target, TValue const &number)
{
    if (sizeof(TValue) == 1)
        writeValue(target, number);     // write chars as chars
    else
        appendNumber(target, number);
}

// ----------------------------------------------------------------------------
// Function write(TNumber); write non-fundamental, convertible type
// ----------------------------------------------------------------------------

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue> >,
                                 Not<Is<FundamentalConcept<TValue> > > >, void)
write(TTarget &target, TValue &number)
{
    writeValue(target, number);
}

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue const> >,
                                 Not<Is<FundamentalConcept<TValue const> > > >, void)
write(TTarget &target, TValue const &number)
{
    writeValue(target, number);
}

// ----------------------------------------------------------------------------
// Function read(Iterator<Input>)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TSize> >, TSize)
read(TTarget &target, TFwdIterator &iter, TSize n)
{
    TSize i;
    for (i = 0; !atEnd(iter) && i < n; ++i, ++iter)
        writeValue(target, value(iter));
    return i;
}

// ----------------------------------------------------------------------------
// Function read(TContainer)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TContainer>
inline typename Size<TTarget>::Type
read(TTarget &target, TContainer &cont)
{
    typename DirectionIterator<TContainer, Input>::Type iter = directionIterator(cont, Input());
    return read(target, iter, length(cont));
}

// ----------------------------------------------------------------------------
// operator<<
// ----------------------------------------------------------------------------

template <typename TStream, typename TIterator>
inline TStream &
operator<<(TStream & target,
           Range<TIterator> const & source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}

}  // namespace seqean

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_BASIC_STREAM_H_
