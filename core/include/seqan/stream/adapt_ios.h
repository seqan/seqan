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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Adaptions for std::ios streams.
// ==========================================================================

#ifndef SEQAN_STREAM_ADAPT_IOS_H_
#define SEQAN_STREAM_ADAPT_IOS_H_

#define SEQAN_ASSERT_BADBIT(s) SEQAN_ASSERT_MSG(s.exceptions() | std::ios_base::badbit, \
        "The badbit exception is not set in the stream. Call either std::exceptions() or streamInit() on the stream.")

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================
// TODO(esiragusa): remove this when chunking goes into basic.

template <typename TObject> struct Chunk;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
struct Position<std::istreambuf_iterator<TValue, TTraits> >
{
    typedef typename std::istreambuf_iterator<TValue, TTraits>::difference_type Type;
};

template <typename TValue, typename TTraits>
struct Position<std::basic_istream<TValue, TTraits> >
{
    typedef typename std::basic_istream<TValue, TTraits>::pos_type Type;
};

template <typename TValue, typename TTraits>
struct Position<std::basic_ostream<TValue, TTraits> >
{
    typedef typename std::basic_ostream<TValue, TTraits>::pos_type Type;
};

template <typename TValue, typename TTraits>
struct Position<std::basic_iostream<TValue, TTraits> >
{
    typedef typename std::basic_iostream<TValue, TTraits>::pos_type Type;
};

template <typename TValue, typename TTraits>
struct Position<std::basic_fstream<TValue, TTraits> > :
    Position<std::basic_iostream<TValue, TTraits>  > {};

template <typename TValue, typename TTraits>
struct Position<std::basic_stringstream<TValue, TTraits> > :
    Position<std::basic_iostream<TValue, TTraits>  > {};


template <typename TValue, typename TTraits>
struct Position<std::basic_ifstream<TValue, TTraits> > :
    Position<std::basic_istream<TValue, TTraits>  > {};

template <typename TValue, typename TTraits>
struct Position<std::basic_istringstream<TValue, TTraits> > :
    Position<std::basic_istream<TValue, TTraits>  > {};


template <typename TValue, typename TTraits>
struct Position<std::basic_ofstream<TValue, TTraits> > :
    Position<std::basic_ostream<TValue, TTraits>  > {};

template <typename TValue, typename TTraits>
struct Position<std::basic_ostringstream<TValue, TTraits> > :
    Position<std::basic_ostream<TValue, TTraits>  > {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
struct Value<std::istreambuf_iterator<TValue, TTraits> >
{
    typedef typename std::istreambuf_iterator<TValue, TTraits>::char_type Type;
};

template <typename TValue, typename TTraits>
struct Value<std::basic_istream<TValue, TTraits> >
{
    typedef typename std::basic_istream<TValue, TTraits>::char_type Type;
};

template <typename TValue, typename TTraits>
struct Value<std::basic_ostream<TValue, TTraits> >
{
    typedef typename std::basic_ostream<TValue, TTraits>::char_type Type;
};

template <typename TValue, typename TTraits>
struct Value<std::basic_iostream<TValue, TTraits> >
{
    typedef typename std::basic_iostream<TValue, TTraits>::char_type Type;
};

template <typename TValue, typename TTraits>
struct Value<std::basic_fstream<TValue, TTraits> > :
    Value<std::basic_iostream<TValue, TTraits>  > {};

template <typename TValue, typename TTraits>
struct Value<std::basic_stringstream<TValue, TTraits> > :
    Value<std::basic_iostream<TValue, TTraits>  > {};


template <typename TValue, typename TTraits>
struct Value<std::basic_ifstream<TValue, TTraits> > :
    Value<std::basic_istream<TValue, TTraits>  > {};

template <typename TValue, typename TTraits>
struct Value<std::basic_istringstream<TValue, TTraits> > :
    Value<std::basic_istream<TValue, TTraits>  > {};


template <typename TValue, typename TTraits>
struct Value<std::basic_ofstream<TValue, TTraits> > :
    Value<std::basic_ostream<TValue, TTraits>  > {};

template <typename TValue, typename TTraits>
struct Value<std::basic_ostringstream<TValue, TTraits> > :
    Value<std::basic_ostream<TValue, TTraits>  > {};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
struct Reference<std::istreambuf_iterator<TValue, TTraits> >
{
    typedef typename std::istreambuf_iterator<TValue, TTraits>::char_type Type;
//    typedef typename std::istreambuf_iterator<TValue, TTraits>::reference Type;
};

// ----------------------------------------------------------------------------
// Metafunction Chunk
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
struct Chunk<std::basic_streambuf<TValue, TTraits> >
{
    typedef Range<TValue*> Type;
};

//template <typename TValue, typename TTraits>
//struct Chunk<std::istreambuf_iterator<TValue, TTraits> >:
//    Chunk<std::basic_streambuf<TValue, TTraits> > {};
//
//template <typename TValue, typename TTraits>
//struct Chunk<std::ostreambuf_iterator<TValue, TTraits> >:
//    Chunk<std::basic_streambuf<TValue, TTraits> > {};

// ----------------------------------------------------------------------------
// Metafunction DefaultOpenMode<std::>
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
struct DefaultOpenMode<std::basic_ifstream<TValue, TTraits> >
{
    enum { VALUE = OPEN_RDONLY };
};

template <typename TValue, typename TTraits>
struct DefaultOpenMode<std::basic_ofstream<TValue, TTraits> >
{
    enum { VALUE = OPEN_WRONLY | OPEN_CREATE };
};

// ============================================================================
// Concepts
// ============================================================================

template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_istream<TValue, TTraits>), (InputStreamConcept));
template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_ifstream<TValue, TTraits>), (InputStreamConcept));
template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_istringstream<TValue, TTraits>), (InputStreamConcept));

template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_ostream<TValue, TTraits>), (OutputStreamConcept));
template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_ofstream<TValue, TTraits>), (OutputStreamConcept));
template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_ostringstream<TValue, TTraits>), (OutputStreamConcept));

template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_iostream<TValue, TTraits>), (BidirectionalStreamConcept));
template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_fstream<TValue, TTraits>), (BidirectionalStreamConcept));
template <typename TValue, typename TTraits>
SEQAN_CONCEPT_IMPL((std::basic_stringstream<TValue, TTraits>), (BidirectionalStreamConcept));

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _getSTLStyleOrigin()
// ----------------------------------------------------------------------------

template <typename TPos>
inline std::ios_base::seekdir
_getSTLStyleOrigin(TPos origin)
{
    switch (origin)
    {
    case SEEK_CUR:
        return std::ios_base::cur;
    case SEEK_END:
        return std::ios_base::end;
    case SEEK_SET:
    default:
        return std::ios_base::beg;
    }
}

// ----------------------------------------------------------------------------
// Function _getSTLStyleOpenMode()
// ----------------------------------------------------------------------------

inline std::ios_base::openmode
_getSTLStyleOpenMode(int openMode)
{
    std::ios_base::openmode flags = std::ios_base::binary;

    if (openMode & OPEN_RDONLY)
        flags |= std::ios_base::in;
    if (openMode & OPEN_WRONLY)
    {
        flags |= std::ios_base::out;
        if (openMode & OPEN_APPEND)
            flags |= std::ios_base::app;
        else
            flags |= std::ios_base::trunc;
    }

    return flags;
}

// ----------------------------------------------------------------------------
// Function getChunk()
// ----------------------------------------------------------------------------

//template <typename TValue, typename TTraits, typename TDirection>
//inline typename Chunk<std::basic_streambuf<TValue, TTraits> >::Type
//getChunk(std::basic_streambuf<TValue, TTraits> const &buf, Tag<TDirection>)
//{
//    return getChunk(static_cast<StreamBuffer<TValue, TTraits> &>(buf), Tag<TDirection>());
//}

// ----------------------------------------------------------------------------
// Function writeValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits, typename TValue2>
inline void writeValue(std::ostreambuf_iterator<TValue, TTraits> &iter, TValue2 val)
{
    *iter = val;
    ++iter;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline bool atEnd(std::istreambuf_iterator<TValue, TTraits> const &it)
{
    return *it == TTraits::eof();
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------
// Deprecated?

template <typename TStream>
inline bool
streamEof(TStream & stream)
{
    // Peak sets eofbit if the next char is EOF.
    SEQAN_ASSERT_BADBIT(stream);
    stream.peek();
    SEQAN_ASSERT_NOT(stream.fail());
    return stream.eof();
}

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------
// Deprecated?

template <typename TStream>
inline typename Value<TStream>::Type
streamPeek(TStream & stream)
{
    typedef typename Value<TStream>::Type TValue;

    // Peak sets eofbit if the next char is EOF.
    SEQAN_ASSERT_BADBIT(stream);
    TValue val = (TValue)stream.peek();
    SEQAN_ASSERT_NOT(stream.fail());
    return val;
}

// ----------------------------------------------------------------------------
// Function streamGet()
// ----------------------------------------------------------------------------
// Deprecated?

template <typename TStream>
inline typename Value<TStream>::Type
streamGet(TStream & stream)
{
    typedef typename Value<TStream>::Type TValue;

    SEQAN_ASSERT_BADBIT(stream);
    SEQAN_ASSERT_NOT(streamEof(stream));
    TValue val = (TValue)stream.get();
    SEQAN_ASSERT_NOT(stream.fail());
    return val;
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------
// Deprecated?

template <typename TStream, typename TValue>
inline void
streamPut(TStream & stream, TValue const & val)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.put(val);
    SEQAN_ASSERT_NOT(stream.fail());
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------
// Deprecated?

template <typename TStream, typename TValue>
inline void
streamFlush(TStream & stream)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream << std::flush;
    SEQAN_ASSERT_NOT(stream.fail());
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------
// Deprecated?

// Input
template <typename TStream, typename TDelta, typename TPos>
inline SEQAN_FUNC_ENABLE_IF(
    And<    Is<InputStreamConcept<TStream> >,
        Not<Is<BidirectionalStreamConcept<TStream> > > >, void)
streamSeek(TStream & stream, TDelta delta, TPos origin)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.seekg(delta, _getSTLStyleOrigin(origin));
    SEQAN_ASSERT_NOT(stream.fail());
}

// Output
template <typename TStream, typename TDelta, typename TPos>
inline SEQAN_FUNC_ENABLE_IF(
    And<    Is<OutputStreamConcept<TStream> >,
        Not<Is<BidirectionalStreamConcept<TStream> > > >, void)
streamSeek(TStream & stream, TDelta delta, TPos origin)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.seekp(delta, _getSTLStyleOrigin(origin));
    SEQAN_ASSERT_NOT(stream.fail());
}

// Bidirectional
template <typename TStream, typename TDelta, typename TPos>
inline SEQAN_FUNC_ENABLE_IF(Is<BidirectionalStreamConcept<TStream> >, void)
streamSeek(TStream & stream, TDelta delta, TPos origin)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.seekp(delta, _getSTLStyleOrigin(origin));
    SEQAN_ASSERT_NOT(stream.fail());
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------
// Deprecated?

// Input
template <typename TStream>
inline SEQAN_FUNC_ENABLE_IF(
    And<    Is<InputStreamConcept<TStream> >,
        Not<Is<BidirectionalStreamConcept<TStream> > > >, typename Position<TStream>::Type)
streamTell(TStream & stream)
{
    typedef typename Position<TStream>::Type TPos;

    SEQAN_ASSERT_BADBIT(stream);
    TPos pos = stream.tellg();
    SEQAN_ASSERT_NOT(stream.fail());
    return pos;
}

// Output
template <typename TStream>
inline SEQAN_FUNC_ENABLE_IF(
    And<    Is<OutputStreamConcept<TStream> >,
        Not<Is<BidirectionalStreamConcept<TStream> > > >, typename Position<TStream>::Type)
streamTell(TStream & stream)
{
    typedef typename Position<TStream>::Type TPos;

    SEQAN_ASSERT_BADBIT(stream);
    TPos pos = stream.tellp();
    SEQAN_ASSERT_NOT(stream.fail());
    return pos;
}

// Bidirectional
template <typename TStream>
inline SEQAN_FUNC_ENABLE_IF(Is<BidirectionalStreamConcept<TStream> >, typename Position<TStream>::Type)
streamTell(TStream & stream)
{
    typedef typename Position<TStream>::Type  TPos;

    SEQAN_ASSERT_BADBIT(stream);
    SEQAN_ASSERT_EQ(stream.tellg(), stream.tellp());
    TPos pos = stream.tellg();
    SEQAN_ASSERT_NOT(stream.fail());
    return pos;
}

// ----------------------------------------------------------------------------
// Function streamInit()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline void streamInit(std::basic_ios<TValue, TTraits> & stream)
{
    stream.exceptions(std::ios_base::badbit);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline bool open(std::basic_fstream<TValue, TTraits> & stream, const char *fileName, int openMode)
{
    streamInit(stream);
    stream.open(fileName, _getSTLStyleOpenMode(openMode));
//    SEQAN_ASSERT_NOT(stream.fail());
    return stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool open(std::basic_ifstream<TValue, TTraits> & stream, const char *fileName, int openMode)
{
    streamInit(stream);
    stream.open(fileName, _getSTLStyleOpenMode(openMode));
//    SEQAN_ASSERT_NOT(stream.fail());
    return stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool open(std::basic_ofstream<TValue, TTraits> & stream, const char *fileName, int openMode)
{
    streamInit(stream);
    stream.open(fileName, _getSTLStyleOpenMode(openMode));
//    SEQAN_ASSERT_NOT(stream.fail());
    return stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool open(std::basic_fstream<TValue, TTraits> & stream, const char *fileName)
{
    return open(stream, fileName, DefaultOpenMode<std::basic_fstream<TValue, TTraits> >::VALUE);
}

template <typename TValue, typename TTraits>
inline bool open(std::basic_ifstream<TValue, TTraits> & stream, const char *fileName)
{
    return open(stream, fileName, DefaultOpenMode<std::basic_ifstream<TValue, TTraits> >::VALUE);
}

template <typename TValue, typename TTraits>
inline bool open(std::basic_ofstream<TValue, TTraits> & stream, const char *fileName)
{
    return open(stream, fileName, DefaultOpenMode<std::basic_ofstream<TValue, TTraits> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline bool close(std::basic_fstream<TValue, TTraits> & stream)
{
    stream.close();
    return !stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool close(std::basic_ifstream<TValue, TTraits> & stream)
{
    stream.close();
    return !stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool close(std::basic_ofstream<TValue, TTraits> & stream)
{
    stream.close();
    return !stream.is_open();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_IOS_H_
