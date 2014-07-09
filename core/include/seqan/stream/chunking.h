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
// TODO(esiragusa): tests

#ifndef SEQAN_STREAM_CHUNKING_
#define SEQAN_STREAM_CHUNKING_

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Chunk
// --------------------------------------------------------------------------

// Chunking is not support for any object (default fallback).
template <typename TObject>
struct Chunk
{
    typedef Nothing Type;
};

// Chunk interface for rooted iterators
template <typename TContainer, typename TValue, typename TSpec>
struct Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >
{
    typedef Range<TValue*> Type;
};

template <typename TValue, typename TSpec>
struct Chunk<String<TValue, TSpec> >:
    Chunk<typename Iterator<String<TValue, TSpec>, Rooted>::Type> {};

template <typename TValue, typename TTraits>
struct Chunk<StreamBuffer<TValue, TTraits> >
{
    typedef Range<TValue*> Type;
};

template <typename TStream, typename TDirection>
struct Chunk<Iter<TStream, StreamIterator<Tag<TDirection> > > >:
    Chunk<typename Iter<TStream, StreamIterator<Tag<TDirection> > >::TStreamBuffer> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function reserveChunk()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TSize>
inline void reserveChunk(TIterator &, TSize)
{}

template <typename TValue, typename TSpec, typename TSize>
inline void reserveChunk(String<TValue, TSpec> &str, TSize size)
{
    reserve(str, length(str) + size);
}

template <typename TContainer, typename TSpec, typename TSize>
inline void reserveChunk(Iter<TContainer, TSpec> &iter, TSize size)
{
    typedef Iter<TContainer, TSpec> TIter;

    TContainer &cont = container(iter);
    typename Size<TIter>::Type newCap = length(cont) + size;

    if (newCap <= capacity(cont))
        return;

    typename Position<TIter>::Type pos = position(iter);
    reserve(cont, newCap);
    setPosition(iter, pos);
}

template <typename TStream, typename TDirection, typename TSize>
inline void reserveChunk(Iter<TStream, StreamIterator<TDirection> > &, TSize)
{}

// ----------------------------------------------------------------------------
// Function advanceChunk()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TSize>
inline void advanceChunk(TIterator &iter, TSize size)
{
    iter += size;
}

template <typename TContainer, typename TSpec, typename TSize>
inline void advanceChunk(Iter<TContainer, TSpec> &iter, TSize size)
{
    typedef Iter<TContainer, TSpec> TIter;

    TContainer &cont = container(iter);
    typename Position<TIter>::Type pos = position(iter);

    iter += size;
    if (pos > length(cont))
        _setLength(cont, pos);
}

template <typename TStream, typename TDirection, typename TSize>
inline void advanceChunk(Iter<TStream, StreamIterator<TDirection> > &iter, TSize size)
{
    iter.streamBuf->seekoff(size, std::ios_base::cur, (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
}

// extend target string size
template <typename TValue, typename TSpec, typename TSize>
inline void advanceChunk(String<TValue, TSpec> &str, TSize size)
{
    _setLength(str, length(str) + size);
}

// ----------------------------------------------------------------------------
// Function getChunk()
// ----------------------------------------------------------------------------

// StreamBuffer
template <typename TValue, typename TTraits>
inline typename Chunk<StreamBuffer<TValue, TTraits> >::Type
getChunk(StreamBuffer<TValue, TTraits> &buf, Input)
{
    return toRange(buf.gptr(), buf.egptr());
}

template <typename TValue, typename TTraits>
inline typename Chunk<StreamBuffer<TValue, TTraits> >::Type
getChunk(StreamBuffer<TValue, TTraits> &buf, Output)
{
    return toRange(buf.pptr(), buf.epptr());
}

// StreamIterator
template <typename TStream, typename TDirection>
inline typename Chunk<Iter<TStream, StreamIterator<Tag<TDirection> > > >::Type
getChunk(Iter<TStream, StreamIterator<Tag<TDirection> > > &iter, Tag<TDirection>)
{
    typedef typename Iter<TStream, StreamIterator<Input> >::TStreamBuffer TStreamBuffer;
    SEQAN_ASSERT(iter.streamBuf != NULL);
    return getChunk(*iter.streamBuf, Tag<TDirection>());
}

// AdaptorIterator
template <typename TContainer, typename TValue, typename TSpec>
inline typename Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >::Type
getChunk(Iter<TContainer, AdaptorIterator<TValue*, TSpec> > &rootedIter, Input)
{
    return toRange(hostIterator(rootedIter), end(container(rootedIter), Standard()));
}

template <typename TContainer, typename TValue, typename TSpec>
inline typename Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >::Type
getChunk(Iter<TContainer, AdaptorIterator<TValue*, TSpec> > &rootedIter, Output)
{
    TContainer &cont = container(rootedIter);
    return toRange(hostIterator(rootedIter), begin(cont, Standard()) + capacity(cont));
}

// SeqAn's strings
template <typename TValue, typename TSpec>
inline typename Chunk<String<TValue, TSpec> >::Type
getChunk(String<TValue, TSpec> &cont, Output)
{
    return toRange(end(cont, Standard()), begin(cont, Standard()) + capacity(cont));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_CHUNKING_
