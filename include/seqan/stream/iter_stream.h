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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Our own implementation of a streambuf_iterator. We could not use the STL's
// iterator as we need to have access to the underlying streambuf which is a
// private member of the STL iterator.
// ==========================================================================
// TODO(esiragusa): tests

#ifndef SEQAN_INCLUDE_SEQAN_STREAM_ITER_STREAM_H_
#define SEQAN_INCLUDE_SEQAN_STREAM_ITER_STREAM_H_

namespace seqan2 {

// ============================================================================
// Tags
// ============================================================================

template <typename TDirection>
struct StreamIterator {};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class StreamBuffer
// ----------------------------------------------------------------------------

/*!
 * @class StreamBuffer
 * @headerfile <seqan/stream.h>
 * @brief Reinterprets the std::basic_streambuf to grant access to protected member functions.
 *
 * @signature template <typename TValue[, typenam TTraits]>
 *            class StreamBuffer : public std::basic_streambuf<TValue, TTraits>;
 *
 * @tparam TValue  The value type of the stream buffer.
 * @tparam TTraits The traits to use, defaults to <tt>std::char_traits&lt;TValue&gt;</tt>.
 */
// TODO(holtgrew): Add documentation for member functions.

 // Unfortunately some of the most useful members of basic_streambuf are
 // protected, so we define a subclass to cast and access them
template <typename TValue, typename TTraits_ = std::char_traits<TValue>>
struct StreamBuffer : public std::basic_streambuf<TValue, TTraits_>
{
    using TTraits      = TTraits_;
    using TBasicStream = std::basic_streambuf<TValue, TTraits_>;

    using TBasicStream::eback;
    using TBasicStream::gptr;
    using TBasicStream::egptr;
    using TBasicStream::gbump;
    using TBasicStream::underflow;

    using TBasicStream::pbase;
    using TBasicStream::pptr;
    using TBasicStream::epptr;
    using TBasicStream::pbump;
    using TBasicStream::overflow;
};

// NOTE(rrahn): This is a wrapper for the StreamBuffer class.
// Since we usually work with std::basic_iostreams and their derivatives, we cannot simply cast a pointer to
// std::basic_filebuf to StreamBuffer to expose it's protected member functions.
// To do so, we only use the StreamBuffer to inherit from basic_streambuf (the base class of all buffer implementations.)
// and only expose the protected member functions as public functions. We then store the original basic_streambuf
// in this wrapper class and whenever access to the protected members is required we use a reinterpret_cast to convert
// basic_streambuf* into the public StreamBuffer*. The reinterpret_cast has zero overhead.
// This fixes an undetected error reported w/ the sanitizer option of gcc (https://github.com/seqan/seqan/issues/2104).
template <typename TValue, typename TTraits_ = std::char_traits<TValue>>
class StreamBufferWrapper
{
protected:
    // Having the throw in a separate function results in overflow()'s assembly using a jump instruction,
    // instead of inline construction and throw of the exception.
    [[noreturn]] void overflow_failure()
    {
        throw IOError{"StreamBuffer: Calling overflow() failed."};
    }

public:

    typedef std::basic_streambuf<TValue, TTraits_> TBasicStreamBuffer;
    typedef StreamBuffer<TValue, TTraits_>         TPubStreamBuffer_;
    typedef typename TPubStreamBuffer_::TTraits    TTraits;

    TBasicStreamBuffer * streamBuf{nullptr};

    StreamBufferWrapper() = default;

    explicit StreamBufferWrapper(TBasicStreamBuffer * _basicStreamBuf) : streamBuf(_basicStreamBuf)
    {}

    size_t chunkSize(Input)
    {
        return baseBuf()->egptr() - baseBuf()->gptr();
    }

    size_t chunkSize(Output)
    {
        return baseBuf()->epptr() - baseBuf()->pptr();
    }

    template <typename TOffset>
    void advanceChunk(TOffset ofs, Input)
    {
        baseBuf()->gbump(ofs);
    }

    template <typename TOffset>
    void advanceChunk(TOffset ofs, Output)
    {
        baseBuf()->pbump(ofs);
    }

    void reserveChunk(Input)
    {
        if (baseBuf()->gptr() == baseBuf()->egptr())
            baseBuf()->underflow();
    }

    void reserveChunk(Output)
    {
        if (baseBuf()->pptr() == baseBuf()->epptr())
            if (baseBuf()->overflow(EOF) == EOF)
                overflow_failure();
    }

    template <typename TOffset>
    typename std::streampos
    seekoff(TOffset ofs, std::ios_base::seekdir way, std::ios_base::openmode which)
    {
        return streamBuf->pubseekoff(ofs, way, which);
    }

    template <typename TOffset, typename TDirection>
    void goFurther(TOffset ofs, TDirection dir)
    {
        size_t left = chunkSize(dir);
        if (SEQAN_LIKELY((size_t)ofs <= left))
        {
            advanceChunk(ofs, dir);
            return;
        }

        while (true)
        {
            size_t adv = std::min((size_t)ofs, left);
            advanceChunk(adv, dir);
            ofs -= adv;
            if (ofs == 0)
                return;

            SEQAN_IF_CONSTEXPR (IsSameType<TDirection, Input>::VALUE)
                baseBuf()->underflow();
            else
                baseBuf()->overflow();
            left = chunkSize(dir);

            if (SEQAN_UNLIKELY(left == 0))
            {
                // if chunking isn't available try to seek
                typename TTraits::pos_type res = seekoff(ofs,
                                                         std::ios_base::cur,
                                                         (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);

                // if seek doesn't work manually skip characters (when reading)
                if (res == typename TTraits::pos_type(typename TTraits::off_type(-1)))
                {
                    SEQAN_IF_CONSTEXPR (IsSameType<TDirection, Input>::VALUE)
                    {
                        for (; ofs != 0; --ofs)
                            baseBuf()->sbumpc();
                    }
                    SEQAN_IF_CONSTEXPR (IsSameType<TDirection, Output>::VALUE)
                    {
                        for (; ofs != 0; --ofs)
                            baseBuf()->sputc('\0');
                    }
                }
                return;
            }
        }
    }

    TPubStreamBuffer_* baseBuf() const
    {
        return reinterpret_cast<TPubStreamBuffer_ *>(streamBuf);
    }
};

// ----------------------------------------------------------------------------
// Class StreamIterator
// ----------------------------------------------------------------------------

/*!
 * @class StreamIterator
 * @extends Iter
 * @brief Abstract base class for input and output stream iterators.
 *
 * @signature template <typename TStream, typename TDirection>
 *            class Iter<TStream, StreamIterator<TDirection> >;
 *
 * @tparam TStream    The @link StreamConcept @endlink to iterate over.
 * @tparam TDirection The iterator direction, one of the @link DirectionTags @endlink.
 */

// ----------------------------------------------------------------------------
// Class Input StreamIterator
// ----------------------------------------------------------------------------

/*!
 * @class InputStreamIterator Input StreamIterator
 * @extends StreamIterator
 * @brief @link Iter @endlink specialiazion for reading from @link StreamConcept streams @endlink.
 *
 * @signature template <typename TStream>
 *            class Iter<TStream, StreamIterator<Input> >;
 *
 * @tparam TStream    The @link StreamConcept @endlink to iterate over.
 */

template <typename TStream>
class Iter<TStream, StreamIterator<Input> >
{
public:
    typedef typename Value<TStream>::Type                    TValue;
    typedef std::basic_istream<TValue>                       TIStream;
    typedef std::basic_streambuf<TValue>                     TBasicBuffer;
    typedef StreamBufferWrapper<TValue>                      TStreamBufferWrapper;
    typedef typename TStreamBufferWrapper::TPubStreamBuffer_ TStreamBuffer;

    TStreamBufferWrapper streamBufWrapper{nullptr};

    /*!
     * @fn InputStreamIterator::Iter
     * @brief The constructors.
     *
     * @signature Iter::Iter();
     * @signature Iter::Iter(stream);
     * @signature Iter::Iter(streamBuffer);
     *
     * @param[in] stream    The <tt>TStream</tt> to read from.
     * @param[in] streamBuf A @link StreamBuffer @endlink to read from.
     *
     * Allows default construction, construction from stream, as well as from a @link StreamBuffer @endlink.
     */
    Iter() = default;

    Iter(TIStream & stream) : streamBufWrapper(stream.rdbuf())
    {
        // printf("streamBuf: %p\n", streamBuf);
        stream.exceptions(std::ios_base::badbit);
    }

    Iter(TBasicBuffer * buf) : streamBufWrapper(buf)
    {}
};

// ----------------------------------------------------------------------------
// Class StreamIterator
// ----------------------------------------------------------------------------

/*!
 * @class OutputStreamIterator Output StreamIterator
 * @extends StreamIterator
 * @brief @link Iter @endlink specialiazion for writing to @link StreamConcept streams @endlink.
 *
 * @signature template <typename TStream>
 *            class Iter<TStream, StreamIterator<Output> >;
 *
 * @tparam TStream    The @link StreamConcept @endlink to iterate over.
 */
template <typename TStream>
class Iter<TStream, StreamIterator<Output> >
{
public:
    typedef typename Value<TStream>::Type                    TValue;
    typedef std::basic_ostream<TValue>                       TOStream;
    typedef std::basic_streambuf<TValue>                     TBasicBuffer;
    typedef StreamBufferWrapper<TValue>                      TStreamBufferWrapper;
    typedef typename TStreamBufferWrapper::TPubStreamBuffer_ TStreamBuffer;

    TStreamBufferWrapper streamBufWrapper{nullptr};

    /*!
     * @fn Iter::Iter
     * @brief Constructor.
     *
     * @signature Iter::Iter()
     * @signature Iter::Iter(stream)
     * @signature Iter::Iter(streamBuf)
     *
     * @param[in] stream    The <tt>TStream</tt> to write to.
     * @param[in] streamBuf A @link StreamBuffer @endlink to write to.
     *
     * Allows default construction, construction from stream, as well as from a @link StreamBuffer @endlink.
     */
    Iter() = default;

    Iter(TOStream & stream) : streamBufWrapper(stream.rdbuf())
    {
        stream.exceptions(std::ios_base::badbit);
    }

    Iter(TBasicBuffer * buf) : streamBufWrapper(buf)
    {}

    template <typename TValue2>
    TValue2 & operator=(TValue2 &val)
    {
        setValue(*this, val);
        return val;
    }

    template <typename TValue2>
    TValue2 const & operator=(TValue2 const &val)
    {
        setValue(*this, val);
        return val;
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Chunk
// ----------------------------------------------------------------------------

/*!
 * @mfn StreamBuffer#Chunk
 * @brief Return chunk type for StreamBuffer
 *
 * @signature Chunk<TStreamBuffer>::Type;
 *
 * @tparam TStreamBuffer The StreamBuffer to query for its chunk type.
 * @return Type          The chunk type of the stream buffer.
 */

template <typename TValue, typename TTraits>
struct Chunk<StreamBuffer<TValue, TTraits> >
{
    typedef Range<TValue*> Type;
};

template <typename TStream, typename TDirection>
struct Chunk<Iter<TStream, StreamIterator<Tag<TDirection> > > >:
    Chunk<typename Iter<TStream, StreamIterator<Tag<TDirection> > >::TStreamBuffer> {};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

/*!
 * @mfn StreamBuffer#Reference
 * @brief Return reference for StreamBuffer.
 *
 * @signature Reference<TStreamBuffer>::Type;
 *
 * @tparam TStreamBuffer The StreamBuffer to query for its reference type.
 * @return Type          The reference type of the stream buffer.
 */

template <typename TStream>
struct Reference<Iter<TStream, StreamIterator<Input> > >:
    Value<Iter<TStream, StreamIterator<Input> > > {};

template <typename TStream>
struct Reference<Iter<TStream, StreamIterator<Input> > const>:
    Value<Iter<TStream, StreamIterator<Input> > > {};

template <typename TStream>
struct Reference<Iter<TStream, StreamIterator<Output> > >
{
    typedef Iter<TStream, StreamIterator<Output> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

/*!
 * @mfn StreamBuffer#GetValue
 * @brief Return get value for StreamBuffer.
 *
 * @signature GetValue<TStreamBuffer>::Type;
 *
 * @tparam TStreamBuffer The StreamBuffer to query for its get value type.
 * @return Type          The get value type of the stream buffer.
 */

template <typename TStream>
struct GetValue<Iter<TStream, StreamIterator<Input> > >:
    Reference<Iter<TStream, StreamIterator<Input> > const> {};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

/*!
 * @mfn StreamBuffer#Position
 * @brief Return position for StreamBuffer.
 *
 * @signature Position<TStreamBuffer>::Type;
 *
 * @tparam TStreamBuffer The StreamBuffer to query for its position type.
 * @return Type          The position type of the stream buffer.
 */

template <typename TStream, typename TDirection>
struct Position<Iter<TStream, StreamIterator<TDirection> > > : Position<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

/*!
 * @mfn StreamBuffer#Difference
 * @brief Return difference for StreamBuffer.
 *
 * @signature Difference<TStreamBuffer>::Type;
 *
 * @tparam TStreamBuffer The StreamBuffer to query for its difference type.
 * @return Type          The difference type of the stream buffer.
 */

template <typename TStream, typename TDirection>
struct Difference<Iter<TStream, StreamIterator<TDirection> > > : Difference<TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/*!
 * @mfn StreamBuffer#Size
 * @brief Return size for StreamBuffer.
 *
 * @signature Size<TStreamBuffer>::Type;
 *
 * @tparam TStreamBuffer The StreamBuffer to query for its size type.
 * @return Type          The size type of the stream buffer.
 */

template <typename TStream, typename TDirection>
struct Size<Iter<TStream, StreamIterator<TDirection> > > : Size<TStream> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function directionIterator()
// ----------------------------------------------------------------------------

/*!
 * @fn StreamConcept#directionIterator
 * @brief Returns direction iterator for Stream.
 *
 * @signature TDirIter directionIterator(stream, dirTag);
 *
 * @param[in] stream The @link StreamConcept @endlink object to compute iterator for.
 * @param[in] dirTag Direction tag, one of the @link DirectionTags @endlink.
 */

template <typename TStream, typename TDirection>
inline SEQAN_FUNC_ENABLE_IF(Is<StreamConcept<TStream> >, Iter<TStream, StreamIterator<TDirection> >)
directionIterator(TStream &stream, TDirection const &)
{
    return Iter<TStream, StreamIterator<TDirection> >(stream);
}

/*!
 * @fn ContainerConcept#directionIterator
 * @brief Returns direction iterator for a container.
 *
 * @signature TDirIter directionIterator(streamBuf, dirTag);
 *
 * @param[in] streamBuf The @link ContainerConcept container @endlink object to compute iterator for.
 * @param[in] dirTag    Direction tag, one of the @link DirectionTags @endlink.
 *
 * @return TDirIter The resulting @link ContainerConcept#DirectionIterator @endlink.
 */

template <typename TContainer, typename TDirection>
inline SEQAN_FUNC_DISABLE_IF(Is<StreamConcept<TContainer> >, typename Iterator<TContainer, Rooted>::Type)
directionIterator(TContainer &cont, TDirection const &)
{
    return begin(cont, Rooted());
}

// ----------------------------------------------------------------------------
// Function reserveChunk()
// ----------------------------------------------------------------------------

/*!
 * @fn StreamIterator#reserveChunk
 * @brief Reserve a chunk in the host of the StreamIterator
 *
 * @signature void reserveChunk(iter, len, dirTag);
 *
 * @param[in] iter   The @link StreamIterator @endlink object to reserve chunks for.
 * @param[in] len    The length of the chunk to reserve.
 * @param[in] dirTag Direction tag, one of @link DirectionTags#Input Input @endlink and @link
 *                   DirectionTags#Input Output @endlink .
 */

template <typename TStream, typename TDirection, typename TSize>
inline void reserveChunk(Iter<TStream, StreamIterator<TDirection> > &iter, TSize, Input dir)
{
    iter.streamBufWrapper.reserveChunk(dir);
}

template <typename TStream, typename TDirection, typename TSize>
inline void reserveChunk(Iter<TStream, StreamIterator<TDirection> > &iter, TSize, Output dir)
{
    iter.streamBufWrapper.reserveChunk(dir);
}

// ----------------------------------------------------------------------------
// Function advanceChunk()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Documentation missing below here.

template <typename TStream, typename TDirection, typename TSize>
inline void advanceChunk(Iter<TStream, StreamIterator<TDirection> > &iter, TSize size)
{
    iter.streamBufWrapper.advanceChunk(size, TDirection());
}

// ----------------------------------------------------------------------------
// Function getChunk()
// ----------------------------------------------------------------------------

// StreamBuffer
template <typename TChunk, typename TValue, typename TTraits>
inline void
getChunk(TChunk &result, StreamBuffer<TValue, TTraits> &buf, Input)
{
    return assignRange(result, buf.gptr(), buf.egptr());
}

template <typename TChunk, typename TValue, typename TTraits>
inline void
getChunk(TChunk &result, StreamBuffer<TValue, TTraits> &buf, Output)
{
    return assignRange(result, buf.pptr(), buf.epptr());
}

// StreamIterator
template <typename TChunk, typename TStream, typename TDirection>
inline void
getChunk(TChunk &result, Iter<TStream, StreamIterator<Tag<TDirection> > > &iter, Tag<TDirection>)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    getChunk(result, *iter.streamBufWrapper.baseBuf(), Tag<TDirection>());
}

// ----------------------------------------------------------------------------
// Function value() - Input
// ----------------------------------------------------------------------------

template <typename TStream>
inline typename Reference<Iter<TStream, StreamIterator<Input> > >::Type
value(Iter<TStream, StreamIterator<Input> > &iter)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    return iter.streamBufWrapper.baseBuf()->sgetc();
}
template <typename TStream>
inline typename Reference<Iter<TStream, StreamIterator<Input> > const>::Type
value(Iter<TStream, StreamIterator<Input> > const &iter)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    return iter.streamBufWrapper.baseBuf()->sgetc();
}

// ----------------------------------------------------------------------------
// Function value() - Ouput
// ----------------------------------------------------------------------------

template <typename TStream>
inline Iter<TStream, StreamIterator<Output> > &
value(Iter<TStream, StreamIterator<Output> > & iter)
{
    return iter;
}
template <typename TStream>
inline Iter<TStream, StreamIterator<Output> > const &
value(Iter<TStream, StreamIterator<Output> > const & iter)
{
    return iter;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue>
inline void
setValue(Iter<TStream, StreamIterator<Output> > & iter, TValue const &val)
{
    return setValue(const_cast<Iter<TStream, StreamIterator<Output> > const &>(iter), val);
}

template <typename TStream, typename TValue>
inline void
setValue(Iter<TStream, StreamIterator<Output> > const & iter, TValue const &val)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    iter.streamBufWrapper.baseBuf()->sputc((typename Value<Iter<TStream, StreamIterator<Output> > >::Type)val);
}

// ----------------------------------------------------------------------------
// Function writeValue()
// ----------------------------------------------------------------------------

// streams
template <typename TContainer, typename TValue>
inline void writeValue(Iter<TContainer, StreamIterator<Output> > &iter, TValue val)
{
    setValue(iter, val);
    //goNext(iter);     // implicitly done by setValue above
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
goNext(Iter<TStream, StreamIterator<Input> > & iter)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    iter.streamBufWrapper.baseBuf()->sbumpc();
}

template <typename TStream>
inline void
goNext(Iter<TStream, StreamIterator<Output> > &)
{
    // We do nothing here, as the stream is advanced by sputc whenever you assign
    // a value to the iterator with *iter= or setValue
}

// we intentionally don't return an iterator here, as the copied iterator wouldn't
// point to the position before the increment.
template <typename TContainer, typename TSpec>
inline void
operator++(Iter<TContainer, StreamIterator<Input> > & iter, int)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    iter.streamBufWrapper.baseBuf()->sbumpc();
}

// ----------------------------------------------------------------------------
// Function goFurther()
// ----------------------------------------------------------------------------

template <typename TStream, typename TOffset, typename TDirection>
inline void
goFurther(Iter<TStream, StreamIterator<TDirection> > &iter, TOffset ofs)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    iter.streamBufWrapper.goFurther(ofs, TDirection());
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection>
inline typename Position<Iter<TStream, StreamIterator<TDirection> > const>::Type
position(Iter<TStream, StreamIterator<TDirection> > const & iter)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    return iter.streamBufWrapper.baseBuf()->pubseekoff(0, std::ios_base::cur,
                                      (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TStream, typename TDirection, typename TPosition>
inline void
setPosition(Iter<TStream, StreamIterator<TDirection> > const & iter, TPosition pos)
{
    SEQAN_ASSERT(iter.streamBufWrapper.baseBuf() != nullptr);
    iter.streamBufWrapper.baseBuf()->pubseekpos(pos, (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TStream>
inline bool
atEnd(Iter<TStream, StreamIterator<Input> > const & iter)
{
    typedef typename Value<Iter<TStream, StreamIterator<Input> > >::Type TValue;
    typedef StreamBuffer<TValue> TStreamBuffer;

    if (SEQAN_UNLIKELY(iter.streamBufWrapper.baseBuf() == nullptr))
    {
        return true;
    }
    else
    {
        TStreamBuffer * const buf = iter.streamBufWrapper.baseBuf();
        if (SEQAN_LIKELY(buf->gptr() < buf->egptr()))
            return false;
        else
            return TStreamBuffer::TTraits::eq_int_type(buf->sgetc(), TStreamBuffer::TTraits::eof());
    }
}

template <typename TStream>
inline bool
atEnd(Iter<TStream, StreamIterator<Output> > const & iter)
{
    typedef typename Value<Iter<TStream, StreamIterator<Input> > >::Type TValue;
    typedef StreamBuffer<TValue> TStreamBuffer;

    if (SEQAN_UNLIKELY(iter.streamBufWrapper.baseBuf() == nullptr))
    {
        return true;
    }
    else
    {
        TStreamBuffer * const buf = iter.streamBufWrapper.baseBuf();
        if (SEQAN_LIKELY(buf->pptr() < buf->epptr()))
            return false;
        else
            return TStreamBuffer::TTraits::eq_int_type(buf->overflow(), TStreamBuffer::TTraits::eof());
    }
}

}  // namespace seqan2

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_STREAM_ITER_STREAM_H_
