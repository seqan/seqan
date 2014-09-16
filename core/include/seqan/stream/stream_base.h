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
// Tag Raw
// --------------------------------------------------------------------------

struct Raw_;
typedef Tag<Raw_> Raw;

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

// TODO(singer): remove this
//#if SEQAN_HAS_ZLIB
struct GZFile_;
typedef Tag<GZFile_> GZFile;
//#endif

//#if SEQAN_HAS_BZIP2
struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;
//#endif

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
// Compression Type Tags
// --------------------------------------------------------------------------

struct GZFile_;
typedef Tag<GZFile_> GZFile;

struct BgzfFile_;
typedef Tag<BgzfFile_> BgzfFile;

struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;

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


template <typename T>
struct MagicHeader<GZFile, T>
{
    static char const VALUE[3];
};

template <typename T>
char const MagicHeader<GZFile, T>::VALUE[3] = { 0x1f, '\x8b', 0x08 };  // gzip's magic number


template <typename T>
struct MagicHeader<BZ2File, T>
{
    static char const VALUE[3];
};

template <typename T>
char const MagicHeader<BZ2File, T>::VALUE[3] = { 0x42, 0x5a, 0x68 };  // bzip2's magic number


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


template <typename T>
struct FileFormatExtensions<GZFile, T>
{
    static char const * VALUE[3];
};

template <typename T>
char const * FileFormatExtensions<GZFile, T>::VALUE[3] =
{
    ".gz",      // default output extension
    ".Z",
    ".zip"
};


template <typename T>
struct FileFormatExtensions<BgzfFile, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileFormatExtensions<BgzfFile, T>::VALUE[2] =
{
    ".bgzf",      // default output extension
    ".bam"        // BAM files are bgzf compressed
};


template <typename T>
struct FileFormatExtensions<BZ2File, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileFormatExtensions<BZ2File, T>::VALUE[2] =
{
    ".bz2",      // default output extension
    ".bz"
};

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
// Forwards
// ============================================================================
// TODO(esiragusa): remove this when chunking goes into basic.

template <typename TDirection>
struct StreamIterator;

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

// streams
template <typename TContainer, typename TValue>
inline void writeValue(Iter<TContainer, StreamIterator<Output> > &iter, TValue val)
{
    setValue(iter, val);
    //goNext(iter);     // implicitly done by setValue above
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
            reserveChunk(target, n);
            reserveChunk(iter, n);
            getChunk(ichunk, iter, Input());
            getChunk(ochunk, target, Output());
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

// chunked, target is pointer (e.g. readRawByte)
template <typename TOValue, typename TFwdIterator, typename TSize>
inline void write(TOValue *ptr, TFwdIterator &iter, TSize n)
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
            reserveChunk(iter, n);
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

// chunked, source is pointer (e.g. readRawByte)
template <typename TTarget, typename TIValue, typename TSize>
inline void write(TTarget &target, TIValue *ptr, TSize n)
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
            reserveChunk(target, n);
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
    Is<IntegerConcept<TSize> >
    Is<Convertible<typename Value<TTarget>::Type,
                   typename Value<TFwdIterator>::Type> >, void)
write(TTarget &target, TFwdIterator &iter, TSize n)
{
    typedef typename Chunk<TFwdIterator>::Type* TIChunk;
    typedef typename Chunk<TTarget>::Type*      TOChunk;

    _write(target, iter, n, TIChunk(), TOChunk());
}

// ----------------------------------------------------------------------------
// Function write(TContainer)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TContainer> >, void)
write(TTarget &target, TContainer &cont)
{
    typename Iterator<TContainer, Rooted>::Type iter = begin(cont, Rooted());
    write(target, iter, length(cont));
}

template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TContainer const> >, void)
write(TTarget &target, TContainer const &cont)
{
    typename Iterator<TContainer const, Rooted>::Type iter = begin(cont, Rooted());
    write(target, iter, length(cont));
}

template <typename TTarget, typename TValue>
inline void
write(TTarget &target, TValue * ptr)
{
    write(target, ptr, length(ptr));
}

// ----------------------------------------------------------------------------
// Function write(StringSet)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSequence, typename TSpec>
inline void
write(TTarget &target, StringSet<TSequence, TSpec> &seqs)
{
    typedef typename Size<StringSet<TSequence, TSpec> >::Type TSize;
    for (TSize i = 0; i < length(seqs); ++i)
    {
        write(target, seqs[i]);
        writeValue(target, '\n');
    }
}

template <typename TTarget, typename TSequence, typename TSpec>
inline void
write(TTarget &target, StringSet<TSequence, TSpec> const &seqs)
{
    typedef typename Size<StringSet<TSequence, TSpec> const>::Type TSize;
    for (TSize i = 0; i < length(seqs); ++i)
    {
        write(target, seqs[i]);
        writeValue(target, '\n');
    }
}

template <typename TTarget, typename TSequence, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TSequence> >, void)
write(TTarget &target, String<TSequence, TSpec> &seqs)
{
    typedef typename Size<String<TSequence, TSpec> >::Type TSize;
    for (TSize i = 0; i < length(seqs); ++i)
    {
        write(target, seqs[i]);
        writeValue(target, '\n');
    }
}

template <typename TTarget, typename TSequence, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TSequence> >, void)
write(TTarget &target, String<TSequence, TSpec> const &seqs)
{
    typedef typename Size<String<TSequence, TSpec> const>::Type TSize;
    for (TSize i = 0; i < length(seqs); ++i)
    {
        write(target, seqs[i]);
        writeValue(target, '\n');
    }
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
inline void read(TTarget &target, TContainer &cont)
{
    typename Iterator<TContainer, Rooted>::Type iter = begin(cont, Rooted());
    read(target, iter, length(cont));
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_H_
