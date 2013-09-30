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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Implementation of iterators working on streams.
// ==========================================================================

#ifndef SEQAN_STREAM_ITER_H_
#define SEQAN_STREAM_ITER_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

template <typename TSpec = void>
struct StreamIterator {};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class StreamIterator Iter
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
class Iter<TStream, StreamIterator<TSpec> >
{
public:
    TStream *stream;

    Iter():
        stream()
    {}

    Iter(TStream &stream):
        stream(&stream)
    {}

    Iter & operator=(Iter const &other)
    {
        stream = other.iter;
        return *this;
    }

    operator typename Value<TStream>::Type () const
    {
        return streamGet(*stream);
    }

    template <typename TValue2>
    TValue2 & operator=(TValue2 &val)
    {
        streamPut(*stream, val);
        return val;
    }

    template <typename TValue2>
    TValue2 const & operator=(TValue2 const &val)
    {
        streamPut(*stream, val);
        return val;
    }

    Iter& operator*() const
    {
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
struct Reference<Iter<TStream, StreamIterator<TSpec> > >
{
    typedef Iter<TStream, StreamIterator<TSpec> >  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
struct Value<Iter<TStream, StreamIterator<TSpec> > > : Value<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
struct Position<Iter<TStream, StreamIterator<TSpec> > > : Position<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
struct Difference<Iter<TStream, StreamIterator<TSpec> > > : Difference<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
struct Size<Iter<TStream, StreamIterator<TSpec> > > : Size<TStream> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
inline bool
atEnd(Iter<TStream, StreamIterator<TSpec> > const & iter)
{
    SEQAN_ASSERT_NEQ(iter.stream, NULL);
    return streamEof(*(iter.stream));
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
inline typename Position<Iter<TStream, StreamIterator<TSpec> > const>::Type
position(Iter<TStream, StreamIterator<TSpec> > const & iter)
{
    SEQAN_ASSERT_NEQ(iter.stream, NULL);
    return streamTell(*(iter.stream));
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
inline void
setPosition(Iter<TStream, StreamIterator<TSpec> > const & iter)
{
    SEQAN_ASSERT_NEQ(iter.stream, NULL);
    streamSeek(*(iter.stream));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ITER_H_
