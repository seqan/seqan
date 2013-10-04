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
        "The badbit exception is not set in the stream. Call either std::exceptions() or init() on the stream.")

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

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
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
struct HasStreamFeature<std::basic_ostream<TValue, TTraits>, IsInput> : False {};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<std::, IsOutput>
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
struct HasStreamFeature<std::basic_istream<TValue, TTraits>, IsOutput> : False {};

// ----------------------------------------------------------------------------
// Feature inheritance
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits, typename TSpec>
struct HasStreamFeature<std::basic_fstream<TValue, TTraits>, TSpec> :
    HasStreamFeature<std::basic_iostream<TValue, TTraits>, TSpec > {};

template <typename TValue, typename TTraits, typename TSpec>
struct HasStreamFeature<std::basic_stringstream<TValue, TTraits>, TSpec> :
    HasStreamFeature<std::basic_iostream<TValue, TTraits>, TSpec > {};


template <typename TValue, typename TTraits, typename TSpec>
struct HasStreamFeature<std::basic_ifstream<TValue, TTraits>, TSpec> :
    HasStreamFeature<std::basic_istream<TValue, TTraits>, TSpec > {};

template <typename TValue, typename TTraits, typename TSpec>
struct HasStreamFeature<std::basic_istringstream<TValue, TTraits>, TSpec> :
    HasStreamFeature<std::basic_istream<TValue, TTraits>, TSpec > {};


template <typename TValue, typename TTraits, typename TSpec>
struct HasStreamFeature<std::basic_ofstream<TValue, TTraits>, TSpec> :
    HasStreamFeature<std::basic_ostream<TValue, TTraits>, TSpec > {};

template <typename TValue, typename TTraits, typename TSpec>
struct HasStreamFeature<std::basic_ostringstream<TValue, TTraits>, TSpec> :
    HasStreamFeature<std::basic_ostream<TValue, TTraits>, TSpec > {};

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
    enum { VALUE = OPEN_WRONLY | OPEN_APPEND };
};

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
// Function streamEof()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline bool
streamEof(std::basic_istream<TValue, TTraits> & stream)
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

template <typename TValue, typename TTraits>
inline typename Value<std::basic_istream<TValue, TTraits> >::Type
streamPeek(std::basic_istream<TValue, TTraits> & stream)
{
    typedef typename Value<std::basic_istream<TValue, TTraits> >::Type TValue_;

    // Peak sets eofbit if the next char is EOF.
    SEQAN_ASSERT_BADBIT(stream);
    TValue_ val = stream.peek();
    SEQAN_ASSERT_NOT(stream.fail());
    return val;
}

// ----------------------------------------------------------------------------
// Function streamGet()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline typename Value<std::basic_istream<TValue, TTraits> >::Type
streamGet(std::basic_istream<TValue, TTraits> & stream)
{
    typedef typename Value<std::basic_istream<TValue, TTraits> >::Type TValue_;

    SEQAN_ASSERT_BADBIT(stream);
    SEQAN_ASSERT_NOT(streamEof(stream));
    TValue_ val = TTraits::to_char_type(stream.get());
    SEQAN_ASSERT_NOT(stream.fail());
    return val;
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits, typename TValue2>
inline void
streamPut(std::basic_ostream<TValue, TTraits> & stream, TValue2 val)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.put(val);
    SEQAN_ASSERT_NOT(stream.fail());
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline void
streamFlush(std::basic_ios<TValue, TTraits> & stream)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream << std::flush;
    SEQAN_ASSERT_NOT(stream.fail());
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits, typename TDelta, typename TPos>
inline void
streamSeek(std::basic_istream<TValue, TTraits> & stream, TDelta delta, TPos origin)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.seekg(delta, _getSTLStyleOrigin(origin));
    SEQAN_ASSERT_NOT(stream.fail());
}

template <typename TValue, typename TTraits, typename TDelta, typename TPos>
inline void
streamSeek(std::basic_ostream<TValue, TTraits> & stream, TDelta delta, TPos origin)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.seekp(delta, _getSTLStyleOrigin(origin));
    SEQAN_ASSERT_NOT(stream.fail());
}

template <typename TValue, typename TTraits, typename TDelta, typename TPos>
inline void
streamSeek(std::basic_iostream<TValue, TTraits> & stream, TDelta delta, TPos origin)
{
    SEQAN_ASSERT_BADBIT(stream);
    stream.seekp(delta, _getSTLStyleOrigin(origin));
    SEQAN_ASSERT_NOT(stream.fail());
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline typename Position<std::basic_istream<TValue, TTraits> >::Type
streamTell(std::basic_istream<TValue, TTraits> & stream)
{
    typedef typename Position<std::basic_istream<TValue, TTraits> >::Type   TPos;

    SEQAN_ASSERT_BADBIT(stream);
    TPos pos = stream.tellg();
    SEQAN_ASSERT_NOT(stream.fail());
    return pos;
}

template <typename TValue, typename TTraits>
inline typename Position<std::basic_ostream<TValue, TTraits> >::Type
streamTell(std::basic_ostream<TValue, TTraits> & stream)
{
    typedef typename Position<std::basic_ostream<TValue, TTraits> >::Type   TPos;

    SEQAN_ASSERT_BADBIT(stream);
    TPos pos = stream.tellp();
    SEQAN_ASSERT_NOT(stream.fail());
    return pos;
}

template <typename TValue, typename TTraits>
inline typename Position<std::basic_iostream<TValue, TTraits> >::Type
streamTell(std::basic_iostream<TValue, TTraits> & stream)
{
    typedef typename Position<std::basic_iostream<TValue, TTraits> >::Type  TPos;

    SEQAN_ASSERT_BADBIT(stream);
    SEQAN_ASSERT_EQ(stream.tellg(), stream.tellp());
    TPos pos = stream.tellg();
    SEQAN_ASSERT_NOT(stream.fail());
    return pos;
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline void
init(std::basic_ios<TValue, TTraits> & stream)
{
    stream.exceptions(std::ios_base::badbit);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline bool
open(std::basic_fstream<TValue, TTraits> & stream, const char *fileName, int openMode)
{
    init(stream);
    stream.open(fileName, _getSTLStyleOpenMode(openMode));
//    SEQAN_ASSERT_NOT(stream.fail());
    return stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool
open(std::basic_ifstream<TValue, TTraits> & stream, const char *fileName, int openMode)
{
    init(stream);
    stream.open(fileName, _getSTLStyleOpenMode(openMode));
//    SEQAN_ASSERT_NOT(stream.fail());
    return stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool
open(std::basic_ofstream<TValue, TTraits> & stream, const char *fileName, int openMode)
{
    init(stream);
    stream.open(fileName, _getSTLStyleOpenMode(openMode));
//    SEQAN_ASSERT_NOT(stream.fail());
    return stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool
open(std::basic_fstream<TValue, TTraits> & stream, const char *fileName)
{
    return open(stream, fileName, DefaultOpenMode<std::basic_fstream<TValue, TTraits> >::VALUE);
}

template <typename TValue, typename TTraits>
inline bool
open(std::basic_ifstream<TValue, TTraits> & stream, const char *fileName)
{
    return open(stream, fileName, DefaultOpenMode<std::basic_ifstream<TValue, TTraits> >::VALUE);
}

template <typename TValue, typename TTraits>
inline bool
open(std::basic_ofstream<TValue, TTraits> & stream, const char *fileName)
{
    return open(stream, fileName, DefaultOpenMode<std::basic_ofstream<TValue, TTraits> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline bool
close(std::basic_fstream<TValue, TTraits> & stream)
{
    stream.close();
    return !stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool
close(std::basic_ifstream<TValue, TTraits> & stream)
{
    stream.close();
    return !stream.is_open();
}

template <typename TValue, typename TTraits>
inline bool
close(std::basic_ofstream<TValue, TTraits> & stream)
{
    stream.close();
    return !stream.is_open();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_IOS_H_
