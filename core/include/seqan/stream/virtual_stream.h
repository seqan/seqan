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
// Virtual stream class to automatically compress/decompress files/streams.
// It adapts the zipstream and bzip2stream classes (Jonathan de Halleux, 2003)
// http://www.codeproject.com/Articles/4457/zipstream-bzip2stream-iostream-wrappers-for-the-zl
// ==========================================================================

#ifndef SEQAN_STREAM_VIRTUAL_STREAM_
#define SEQAN_STREAM_VIRTUAL_STREAM_

#if SEQAN_HAS_ZLIB
#include "zipstream/zipstream.hpp"
#endif

#if SEQAN_HAS_BZIP2
#include "zipstream/bzip2stream.hpp"
#endif

namespace seqan {

// ============================================================================
// Tags, Enums
// ============================================================================

// --------------------------------------------------------------------------
// TagList CompressedFileTypes
// --------------------------------------------------------------------------

typedef
#if SEQAN_HAS_ZLIB
    TagList<GZFile,
#endif
#if SEQAN_HAS_BZIP2
    TagList<BZ2File,
#endif
    TagList<Nothing>
#if SEQAN_HAS_BZIP2
    >
#endif
#if SEQAN_HAS_ZLIB
    >
#endif
    CompressedFileTypes;  // if TagSelector is set to -1, the file format is auto-detected

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction MagicHeader
// --------------------------------------------------------------------------

template <typename TTag, typename T = void>
struct MagicHeader;


template <typename T>
struct MagicHeader<GZFile, T>
{
    static unsigned char const VALUE[3];
};

template <typename T>
unsigned char const MagicHeader<GZFile, T>::VALUE[3] = { 0x1f, 0x8b, 0x08 };  // gzip's magic number


template <typename T>
struct MagicHeader<BZ2File, T>
{
    static unsigned char const VALUE[3];
};

template <typename T>
unsigned char const MagicHeader<BZ2File, T>::VALUE[3] = { 0x42, 0x5a, 0x68 };  // bzip2's magic number


template <typename T>
struct MagicHeader<Nothing, T>
{
    static unsigned char const *VALUE;
};

template <typename T>
unsigned char const *MagicHeader<Nothing, T>::VALUE = NULL;


// TODO(weese:) The following defines makes the old guessFormat functions in file_format_mmap.h obsolete. Disable them!
template <typename T>
struct MagicHeader<Fasta, T>
{
    static unsigned char const VALUE[1];
};

template <typename T>
unsigned char const MagicHeader<Fasta, T>::VALUE[1] = { '>' }; // Fasta's first character


template <typename T>
struct MagicHeader<Fastq, T>
{
    static unsigned char const VALUE[1];
};

template <typename T>
unsigned char const MagicHeader<Fastq, T>::VALUE[1] = { '@' };  // Fastq's first character


// --------------------------------------------------------------------------
// Metafunction FileFormatExtensions
// --------------------------------------------------------------------------

// TODO(weese:) rename FileFormatExtensions to FileTypeExtensions or FileExtensions

template <typename T>
struct FileFormatExtensions<GZFile, T>
{
    static char const * VALUE[3];
};

template <typename T>
char const * FileFormatExtensions<GZFile, T>::VALUE[3] = {
    ".gz",      // default output extension
    ".Z",
    ".zip" };


template <typename TTag, typename T>
struct FileFormatExtensions;

template <typename T>
struct FileFormatExtensions<BZ2File, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileFormatExtensions<BZ2File, T>::VALUE[2] = {
    ".bz2",      // default output extension
    ".bz" };


template <typename T>
struct FileFormatExtensions<Nothing, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileFormatExtensions<Nothing, T>::VALUE[1] = {
    "" };       // default output extension

// --------------------------------------------------------------------------
// Metafunction VirtualStreamSwitch_
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TFormatTag>
struct VirtualStreamSwitch_
{
    typedef Nothing Type;
};

template <typename TValue>
struct VirtualStreamSwitch_<TValue, Input, GZFile>
{
    typedef zlib_stream::basic_zip_istream<TValue> Type;
};

template <typename TValue>
struct VirtualStreamSwitch_<TValue, Output, GZFile>
{
    typedef zlib_stream::basic_zip_ostream<TValue> Type;
};

template <typename TValue>
struct VirtualStreamSwitch_<TValue, Input, BZ2File>
{
    typedef bzip2_stream::basic_bzip2_istream<TValue> Type;
};

template <typename TValue>
struct VirtualStreamSwitch_<TValue, Output, BZ2File>
{
    typedef bzip2_stream::basic_bzip2_ostream<TValue> Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class VirtualStreamFactoryContext_
// --------------------------------------------------------------------------

template <typename TVirtualStream>
struct VirtualStreamFactoryContext_;

// --------------------------------------------------------------------------
// Class VirtualStreamContext_
// --------------------------------------------------------------------------

// a compressed stream lives in the VirtualStreamContext_ and provides a basic_streambuf

// generic subclass with virtual destructor
template <typename TValue, typename TDirection, typename TFormatTag = void>
struct VirtualStreamContext_:
    VirtualStreamContext_<TValue, TDirection>
{
    typename VirtualStreamSwitch_<TValue, TDirection, TFormatTag>::Type stream;

    template <typename TObject>
    VirtualStreamContext_(TObject &object):
        stream(object)
    {
        this->streamBuf = stream.rdbuf();
    }

    ~VirtualStreamContext_();
};

// special case: no compression, we simply forward the file stream
template <typename TValue, typename TDirection>
struct VirtualStreamContext_<TValue, TDirection, Nothing>:
    VirtualStreamContext_<TValue, TDirection>
{
    template <typename TObject>
    VirtualStreamContext_(TObject &object)
    {
        this->streamBuf = object.rdbuf();
    }
};

// base class
template <typename TValue, typename TDirection>
struct VirtualStreamContext_<TValue, TDirection>
{
    std::basic_streambuf<TValue> *streamBuf;
    VirtualStreamContext_(): streamBuf() {}
    virtual ~VirtualStreamContext_() {}
};

// --------------------------------------------------------------------------
// Class VirtualStream
// --------------------------------------------------------------------------

// The VirtualStream class handles a file or input stream and auto-detects data
// compression from file name or stream.
// We inherit from std::basic_Xstream to provide the convenient stream interface.
template <typename TValue, typename TDirection>
class VirtualStream: public BasicStream<TValue, TDirection>::Type
{
public:
    typedef std::fstream                                    TFile;
    typedef typename BasicStream<TValue, TDirection>::Type  TBasicStream;
    typedef std::basic_streambuf<TValue>                    TStreamBuffer;
    typedef typename BasicStream<TValue, TDirection>::Type  TStream;
    typedef VirtualStreamContext_<TValue, TDirection>       TVirtualStreamContext;

    TFile                   file;
    TStreamBuffer           *streamBuf;
    TVirtualStreamContext   *context;

    VirtualStream():
        TBasicStream(NULL),
        streamBuf(),
        context()
    {}

    VirtualStream(TStreamBuffer &streamBuf):
        TBasicStream(NULL),
        streamBuf(streamBuf),
        context()
    {}

    VirtualStream(TStream &stream):
        TBasicStream(NULL),
        streamBuf(),
        context()
    {
        open(*this, stream);
    }

    VirtualStream(const char *fileName, int openMode):
        TBasicStream(NULL),
        streamBuf(),
        context()
    {
        open(*this, fileName, openMode);
    }

    ~VirtualStream()
    {
        close(*this);
    }


    operator TStreamBuffer*() const
    {
        return streamBuf;
    }

    operator bool() const
    {
        return streamBuf != NULL;
    }

    void _init()
    {
        this->init(streamBuf);
    }

    TStreamBuffer* rdbuf() const
    {
        return streamBuf;
    }
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection>
struct Value<VirtualStream<TValue, TDirection> >
{
    typedef TValue Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection>
struct Position<VirtualStream<TValue, TDirection> >:
    Position<typename VirtualStream<TValue, TDirection>::TFile> {};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection>
struct Iterator<VirtualStream<TValue, TDirection>, Standard>
{
    typedef Iter<VirtualStream<TValue, TDirection>, StreamIterator<TDirection> >    Type;
};

// ----------------------------------------------------------------------------
// Concepts
// ----------------------------------------------------------------------------

template <typename TValue>
SEQAN_CONCEPT_IMPL((VirtualStream<TValue, Input>), (InputStreamConcept));

template <typename TValue>
SEQAN_CONCEPT_IMPL((VirtualStream<TValue, Output>), (OutputStreamConcept));

template <typename TValue>
SEQAN_CONCEPT_IMPL((VirtualStream<TValue, Bidirectional>), (BidirectionalStreamConcept));

// --------------------------------------------------------------------------
// Class VirtualStreamFactoryContext_
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection>
struct VirtualStreamFactoryContext_<VirtualStream<TValue, TDirection> >
{
    typedef VirtualStream<TValue, TDirection>       TVirtualStream;
    typedef typename TVirtualStream::TStream        TStream;

    TStream &stream;
    VirtualStreamFactoryContext_(TStream &stream):
        stream(stream) {}
};

template <typename TVirtualStream>
struct Value<VirtualStreamFactoryContext_<TVirtualStream> >
{
    typedef typename TVirtualStream::TVirtualStreamContext *Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function tagApply()
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TFormat>
inline VirtualStreamContext_<TValue, TDirection>*
tagApply(VirtualStreamFactoryContext_<VirtualStream<TValue, TDirection> > &ctx, Tag<TFormat>)
{
    return new VirtualStreamContext_<TValue, TDirection, Tag<TFormat> >(ctx.stream);
}


template <typename TContext>
inline typename Value<TContext>::Type
tagApply(TContext &, TagSelector<>)
{
    return typename Value<TContext>::Type();
}

template <typename TContext, typename TTagList>
inline typename Value<TContext>::Type
tagApply(TContext &ctx, TagSelector<TTagList> &format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
        return tagApply(ctx, TFormatTag());

    return tagApply(ctx, static_cast<typename TagSelector<TTagList>::Base &>(format));
}

// --------------------------------------------------------------------------
// Function flush()
// --------------------------------------------------------------------------

template<
	typename Elem, 
	typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT >
inline void
flush(zlib_stream::basic_zip_istream<Elem,Tr,ElemA,ByteT,ByteAT> &)
{}

template<
	typename Elem, 
	typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT >
inline void
flush(zlib_stream::basic_zip_ostream<Elem,Tr,ElemA,ByteT,ByteAT> &stream)
{
    stream.zflush();
}

template<
	typename Elem, 
	typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT >
inline void
flush(bzip2_stream::basic_bzip2_istream<Elem,Tr,ElemA,ByteT,ByteAT> &)
{}

template<
	typename Elem, 
	typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT >
inline void
flush(bzip2_stream::basic_bzip2_ostream<Elem,Tr,ElemA,ByteT,ByteAT> &stream)
{
    stream.zflush();
}

// --------------------------------------------------------------------------
// Function guessFormat()
// --------------------------------------------------------------------------

// read first bytes of a file/stream and compare with file format's magic header
template <typename TStream, typename TFormat_>
inline bool
guessFormat(TStream &istream, Tag<TFormat_>)
{
    typedef Tag<TFormat_> TFormat;

    if (MagicHeader<TFormat>::VALUE == NULL)
        return true;

    bool match = true;

    // check magic header
    unsigned i;
    for (i = 0; i != sizeof(MagicHeader<TFormat>::VALUE) / sizeof(char); ++i)
    {
        int c = (int)istream.get();
        if (c != MagicHeader<TFormat>::VALUE[i])
        {
            match = false;
            if (c != EOF)
                ++i;
            break;
        }
    }

    // unget all read characters
    for (; i > 0; --i)
        istream.unget();

    return match;
}

// --------------------------------------------------------------------------
// Function open()
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TStream>
inline bool
open(VirtualStream<TValue, TDirection> &stream, TStream &fileStream)
{
    typedef VirtualStream<TValue, TDirection> TVirtualStream;
    typedef typename TVirtualStream::TFile TFile;
    typedef typename TVirtualStream::TStreamBuffer TStreamBuffer;

    // detect compression type from file extension
    TagSelector<CompressedFileTypes> fileType;
    guessFormat(fileStream, fileType);

    VirtualStreamFactoryContext_<TVirtualStream> ctx(fileStream);

    // create a new (un)zipper buffer
    stream.context = tagApply(ctx, fileType);
    if (stream.context == NULL)
        return false;
    stream.streamBuf = stream.context->streamBuf;

    // reset our outer stream interface
    stream._init();
    return true;
}

template <typename TValue, typename TDirection>
inline bool
open(VirtualStream<TValue, TDirection> &stream, const char *fileName, int openMode)
{
    typedef VirtualStream<TValue, TDirection> TVirtualStream;
    typedef typename TVirtualStream::TFile TFile;
    typedef typename TVirtualStream::TStreamBuffer TStreamBuffer;
    
    if (!open(stream.file, fileName, openMode))
        return false;

    // detect compression type from file extension
    TagSelector<CompressedFileTypes> fileType;
    guessFormatFromFilename(fileName, fileType);

    VirtualStreamFactoryContext_<TVirtualStream> ctx(stream.file);

    // create a new (un)zipper buffer
    stream.context = tagApply(ctx, fileType);
    if (stream.context == NULL)
    {
        close(stream.file);
        return false;
    }
    stream.streamBuf = stream.context->streamBuf;

    // reset our outer stream interface
    stream._init();
    return true;
}

template <typename TValue, typename TDirection, typename TFormatTag>
VirtualStreamContext_<TValue, TDirection, TFormatTag>::~VirtualStreamContext_()
{
    flush(this->stream);
}

// --------------------------------------------------------------------------
// Function close()
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection>
inline bool
close(VirtualStream<TValue, TDirection> &stream)
{
    delete stream.context;
    stream.context = NULL;
    stream.streamBuf = NULL;
    return close(stream.file);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_VIRTUAL_STREAM_
