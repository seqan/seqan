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
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Basic definitions for the stream module.
// ==========================================================================

#ifndef SEQAN_STREAM_STREAM_BASE_H_
#define SEQAN_STREAM_STREAM_BASE_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// --------------------------------------------------------------------------
// Compression Type Tags
// --------------------------------------------------------------------------

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

#if SEQAN_HAS_ZLIB
struct GZFile_;
typedef Tag<GZFile_> GZFile;
#endif

#if SEQAN_HAS_BZIP2
struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;
#endif

// --------------------------------------------------------------------------
// Direction Tags
// --------------------------------------------------------------------------

struct Input_;
typedef Tag<Input_> Input;

struct Output_;
typedef Tag<Output_> Output;

struct Bidirectional_;
typedef Tag<Bidirectional_> Bidirectional;

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction BasicStream
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection>
struct BasicStream:
    If<
        IsSameType<TDirection, Input>,
        std::basic_istream<TValue>,
        typename If<
            IsSameType<TDirection, Output>,
            std::basic_ostream<TValue>,
            std::basic_iostream<TValue>
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

// ============================================================================
// Concepts
// ============================================================================

// --------------------------------------------------------------------------
// Concept InputStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT(InputStreamConcept, (TStream))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(InputStreamConcept)
    {
    }
};

// --------------------------------------------------------------------------
// Concept OutputStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT(OutputStreamConcept, (TStream))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(OutputStreamConcept)
    {
    }
};

// --------------------------------------------------------------------------
// Concept BidirectionalStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT_REFINE(BidirectionalStreamConcept, (TStream), (InputStreamConcept)(OutputStreamConcept))
{
};

// ============================================================================
// Forwards
// ============================================================================
// TODO(esiragusa): remove this when chunking goes into basic.

template <typename TObject> struct Chunk;

template <typename TValue, typename TTraits, typename TValue2>
inline void writeValue(std::ostreambuf_iterator<TValue, TTraits> &iter, TValue2 val);

template <typename TValue, typename TTraits> inline bool
atEnd(std::istreambuf_iterator<TValue, TTraits> const &it);

// ============================================================================
// Functions
// ============================================================================
// TODO(esiragusa): not unique to streams - move them into basic
// TODO(esiragusa): tests

// ----------------------------------------------------------------------------
// Function writeValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TValue>
inline void writeValue(TContainer &cont, TValue val)
{
    appendValue(cont, val);
}

// ----------------------------------------------------------------------------
// Function writeValue(Iter)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TValue>
inline void writeValue(Iter<TContainer, TSpec> &iter, TValue val)
{
    typedef Iter<TContainer, TSpec> TIter;

    TContainer &cont = container(iter);
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

// ----------------------------------------------------------------------------
// Function _write(); Element-wise
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize, typename TIChunk, typename TOChunk>
inline void _write(TTarget &target, TFwdIterator &iter, TSize n, TIChunk, TOChunk)
{
    for (; n > (TSize)0; --n, ++iter)
        writeValue(target, value(iter));
}

// ----------------------------------------------------------------------------
// Function _write(); Chunked
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize, typename TIValue, typename TOValue>
inline void _write(TTarget &target, TFwdIterator &iter, TSize n, Range<TIValue*> *, Range<TOValue*> *)
{
    Range<TIValue*> ichunk;
    Range<TOValue*> ochunk;

    typename Size<TTarget>::Type minChunkSize;
    for (; n > (TSize)0; n -= minChunkSize)
    {
        ichunk = getChunk(iter, Input());
        minChunkSize = ichunk.end - ichunk.begin;
        SEQAN_ASSERT_GT(minChunkSize, 0u);

        reserveChunk(target, minChunkSize);
        ochunk = getChunk(end(target, Rooted()), Output());

        typename Size<TTarget>::Type olen = ochunk.end - ochunk.begin;

        if (minChunkSize > olen)
            minChunkSize = olen;

        if (minChunkSize > n)
            minChunkSize = n;

        ichunk.end = ichunk.begin + minChunkSize;

        for (; ichunk.begin != ichunk.end; ++ichunk.begin, ++ochunk.begin)
            assignValue(ochunk.begin, getValue(ichunk.begin));

        iter += minChunkSize;                      // advance input iterator
        advanceChunk(target, minChunkSize);
    }
}

// ----------------------------------------------------------------------------
// Function write(TValue *)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TValue, typename TSize>
inline void write(TTarget &target, TValue *ptr, TSize n)
{
    typedef Range<TValue*>                          TRange;
    typedef typename Iterator<TRange, Rooted>::Type TIterator;
    typedef typename Chunk<TIterator>::Type*        TIChunk;
    typedef typename Chunk<TTarget>::Type*          TOChunk;

    TRange range(ptr, ptr + n);
    TIterator iter = begin(range, Rooted());
    _write(target, iter, n, TIChunk(), TOChunk());
}

// ----------------------------------------------------------------------------
// Function write(Iterator<Input>)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize>
inline void write(TTarget &target, TFwdIterator &iter, TSize n)
{
    typedef typename Chunk<TFwdIterator>::Type* TIChunk;
    typedef typename Chunk<TTarget>::Type*      TOChunk;

    _write(target, iter, n, TIChunk(), TOChunk());
}

// ----------------------------------------------------------------------------
// Function read(Iterator<Input>)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize>
inline TSize read(TTarget &target, TFwdIterator &iter, TSize n)
{
    TSize i;
    for (i = 0; !atEnd(iter) && i < n; ++i, ++iter)
        writeValue(target, value(iter));
    return i;
}

// ----------------------------------------------------------------------------
// Function write(TContainer)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TContainer>
inline void write(TTarget &target, TContainer const &cont)
{
    typename Iterator<TContainer const, Rooted>::Type iter = begin(cont, Rooted());
    write(target, iter, length(cont));
}

// ----------------------------------------------------------------------------
// Function read(TContainer)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TContainer>
inline void read(TTarget &target, TContainer &cont)
{
    typename Iterator<TContainer, Rooted>::Type iter = begin(cont, Rooted());
    read(target, iter, length(cont));
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_H_
