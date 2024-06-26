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
// ==========================================================================
// Virtual stream class to automatically compress/decompress files/streams.
// It adapts the zipstream and bzip2stream classes (Jonathan de Halleux, 2003)
// https://www.codeproject.com/Articles/4457/zipstream-bzip2stream-iostream-wrappers-for-the-zl
// and our own implementation of a parallel bgzfstream.
// ==========================================================================

#ifndef SEQAN_STREAM_VIRTUAL_STREAM_
#define SEQAN_STREAM_VIRTUAL_STREAM_

namespace seqan2 {

// ============================================================================
// Forwards
// ============================================================================

template <typename TStream>
struct StreamFormat;

// ============================================================================
// Tags, Enums
// ============================================================================

// --------------------------------------------------------------------------
// TagList CompressedFileTypes
// --------------------------------------------------------------------------

typedef
#if SEQAN_HAS_ZLIB
    TagList<BgzfFile,
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
    >
#endif
    CompressedFileTypes;  // if TagSelector is set to -1, the file format is auto-detected

// ============================================================================
// Metafunctions
// ============================================================================

#if SEQAN_HAS_ZLIB

template <typename Elem, typename Tr, typename ElemA, typename ByteT, typename ByteAT>
struct Value<basic_bgzf_istream<Elem, Tr, ElemA, ByteT, ByteAT> > :
    Value<std::basic_istream<Elem, Tr> > {};

template <typename Elem, typename Tr, typename ElemA, typename ByteT, typename ByteAT>
struct Position<basic_bgzf_istream<Elem, Tr, ElemA, ByteT, ByteAT> > :
    Position<std::basic_istream<Elem, Tr> > {};


template <typename Elem, typename Tr, typename ElemA, typename ByteT, typename ByteAT>
struct Value<basic_bgzf_ostream<Elem, Tr, ElemA, ByteT, ByteAT> > :
    Value<std::basic_ostream<Elem, Tr> > {};

template <typename Elem, typename Tr, typename ElemA, typename ByteT, typename ByteAT>
struct Position<basic_bgzf_ostream<Elem, Tr, ElemA, ByteT, ByteAT> > :
    Position<std::basic_ostream<Elem, Tr> > {};


template <typename Elem, typename Tr, typename ElemA, typename ByteT, typename ByteAT>
SEQAN_CONCEPT_IMPL((basic_bgzf_istream<Elem, Tr, ElemA, ByteT, ByteAT>), (InputStreamConcept));

template <typename Elem, typename Tr, typename ElemA, typename ByteT, typename ByteAT>
SEQAN_CONCEPT_IMPL((basic_bgzf_ostream<Elem, Tr, ElemA, ByteT, ByteAT>), (OutputStreamConcept));

#endif

// --------------------------------------------------------------------------
// Metafunction VirtualStreamSwitch_
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TFormatTag>
struct VirtualStreamSwitch_
{
    typedef Nothing Type;
};

#if SEQAN_HAS_ZLIB
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
struct VirtualStreamSwitch_<TValue, Input, BgzfFile>
{
    typedef basic_bgzf_istream<TValue> Type;
};

template <typename TValue>
struct VirtualStreamSwitch_<TValue, Output, BgzfFile>
{
    typedef basic_bgzf_ostream<TValue> Type;
};

#endif

#if SEQAN_HAS_BZIP2

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
#endif

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

// base class
template <typename TValue, typename TTraits>
struct VirtualStreamContextBase_
{
    std::basic_streambuf<TValue, TTraits> *streamBuf;
    VirtualStreamContextBase_(): streamBuf() {}
    virtual ~VirtualStreamContextBase_() {}
};

// generic subclass with virtual destructor
template <typename TValue, typename TDirection, typename TTraits, typename TFormatTag = void>
struct VirtualStreamContext_:
    VirtualStreamContextBase_<TValue, TTraits>
{
    typename VirtualStreamSwitch_<TValue, TDirection, TFormatTag>::Type stream;

    template <typename TObject>
    VirtualStreamContext_(TObject &object):
        stream(object)
    {
        this->streamBuf = stream.rdbuf();
    }
};

// special case: no compression, we simply forward the file stream
template <typename TValue, typename TDirection, typename TTraits>
struct VirtualStreamContext_<TValue, TDirection, TTraits, Nothing>:
    VirtualStreamContextBase_<TValue, TTraits>
{
    template <typename TObject>
    VirtualStreamContext_(TObject &object)
    {
        this->streamBuf = object.rdbuf();
    }
};

// --------------------------------------------------------------------------
// Class VirtualStream
// --------------------------------------------------------------------------
// The VirtualStream class auto-detects data compression from file name or stream.
// It inherits from std::basic_Xstream to provide the convenient stream interface.
// It accepts a file or stream.

/*!
 * @class VirtualStream
 * @implements StreamConcept
 * @headerfile <seqan/stream.h>
 * @brief Provides seamless (de)compression for another @link StreamConcept stream @endlink.
 *
 * @signature template <typename TValue, typename TDirection, typename TTraits>
 *            class VirtualStream;
 *
 * @tparam TValue       The stream value type, e.g. <tt>char</tt>.
 * @tparam TDirection   The stream direction, one of @link DirectionTags @endlink.
 * @tparam TTraits      The stream traits, defaults to <tt>std::char_traits&lt;TValue&gt;</tt>.
 */

template <typename TValue, typename TDirection, typename TTraits = std::char_traits<TValue> >
class VirtualStream: public BasicStream<TValue, TDirection, TTraits>::Type
{
public:
    typedef typename BasicStream<TValue, TDirection, TTraits>::Type TStream;                // the stream base class we expose
    typedef std::basic_fstream<TValue, TTraits>                     TFile;                  // if a real file should be opened
//    typedef FileStream<TValue, TDirection>                          TFile;                  // if a real file should be opened
    typedef BufferedStream<TStream, TDirection>                     TBufferedStream;        // if input stream is not buffered
    typedef std::basic_streambuf<TValue, TTraits>                   TStreamBuffer;          // the streambuf to use
    typedef VirtualStreamContextBase_<TValue, TTraits>              TVirtualStreamContext;  // the owner of the streambuf
    typedef typename StreamFormat<VirtualStream>::Type              TFormat;                // detected stream format

    TFile                   file;
    TBufferedStream         bufferedStream;
    TStreamBuffer           *streamBuf;
    TVirtualStreamContext   *context;
    TFormat                 format;

    /*!
     * @fn VirtualStream::VirtualStream
     * @brief Default constructor and construction from stream, stream buffer, or filename.
     *
     * @signature VirtualStream::VirtualStream();
     * @signature VirtualStream::VirtualStream(stream);
     * @signature VirtualStream::VirtualStream(streamBuf);
     * @signature VirtualStream::VirtualStream(fileName, openMode);
     *
     * @param[in] stream    The @link StreamConcept stream @endlink to attach to.
     * @param[in] streamBuf The @link StreamBuffer stream buffer @endlink to attach to.
     * @param[in] fileName  Path to the file to open. Type: <tt>char const *</tt>.
     * @param[in] openMode  The open mode. Type: <tt>int</tt>.
     */
    VirtualStream():
        TStream(NULL),
        streamBuf(),
        context()
    {}

    VirtualStream(TStreamBuffer &streamBuf):
        TStream(NULL),
        streamBuf(streamBuf),
        context()
    {}

    VirtualStream(TStream &stream):
        TStream(NULL),
        streamBuf(),
        context()
    {
        open(*this, stream);
    }

    VirtualStream(const char *fileName,
                  int openMode = DefaultOpenMode<VirtualStream>::VALUE):
        TStream(NULL),
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

    /*!
     * @fn VirtualStream::getFileExtensions
     * @brief Static function that returns a list of allowed file format extension.
     *
     * @signature TExtensionVector getFileExtensions()
     *
     * @return TExtensionVector A <tt>std::vector&lt;std::string&gt;</tt> with the allowed file extensions.
     */
    static std::vector<std::string>
    getFileExtensions()
    {
        std::vector<std::string> extensions;
        _getFileExtensions(extensions, TFormat());
        return extensions;
    }
};

// ----------------------------------------------------------------------------
// Metafunction StreamFormat
// ----------------------------------------------------------------------------

/*!
 * @defgroup StreamFormats Stream Formats
 * @brief Tags for identifying stream formats.
 */

/*!
 * @mfn VirtualStream#StreamFormat
 * @brief Metafunction for retrieving the format type of a @link VirtualStream @endlink.
 *
 * @signature StreamFormat<TStream>::Type;
 *
 * @tparam TStream  The stream file type to query for its file format type.
 * @return Type     The resulting @link VirtualStream @endlink file formats type.
 */

template <typename TValue, typename TDirection, typename TTraits>
struct StreamFormat<VirtualStream<TValue, TDirection, TTraits> >
{
    typedef TagSelector<CompressedFileTypes> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TTraits>
struct Value<VirtualStream<TValue, TDirection, TTraits> >
{
    typedef TValue Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TTraits>
struct Position<VirtualStream<TValue, TDirection, TTraits> >:
    Position<typename VirtualStream<TValue, TDirection, TTraits>::TFile> {};

// ----------------------------------------------------------------------------
// Metafunction Iterator<Standard>
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TTraits>
struct Iterator<VirtualStream<TValue, TDirection, TTraits>, TDirection>
{
    typedef Iter<VirtualStream<TValue, TDirection, TTraits>, StreamIterator<TDirection> > Type;
};

// --------------------------------------------------------------------------
// Metafunction DefaultOpenMode<Input>
// --------------------------------------------------------------------------

template <typename TValue, typename TDummy, typename TTraits>
struct DefaultOpenMode<VirtualStream<TValue, Input, TTraits>, TDummy>
{
    enum { VALUE = OPEN_RDONLY };
};

// --------------------------------------------------------------------------
// Metafunction DefaultOpenMode<Output>
// --------------------------------------------------------------------------

template <typename TValue, typename TDummy, typename TTraits>
struct DefaultOpenMode<VirtualStream<TValue, Output, TTraits>, TDummy>
{
    enum { VALUE = OPEN_WRONLY | OPEN_CREATE };
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

template <typename TValue, typename TDirection, typename TTraits>
struct VirtualStreamFactoryContext_<VirtualStream<TValue, TDirection, TTraits> >
{
    typedef VirtualStream<TValue, TDirection, TTraits>  TVirtualStream;
    typedef typename TVirtualStream::TStream            TStream;

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
// NOTE(esiragusa): private / impl.

template <typename TValue, typename TDirection, typename TTraits, typename TFormat>
inline VirtualStreamContextBase_<TValue, TTraits> *
tagApply(VirtualStreamFactoryContext_<VirtualStream<TValue, TDirection, TTraits> > &ctx, Tag<TFormat>)
{
    return new VirtualStreamContext_<TValue, TDirection, TTraits, Tag<TFormat> >(ctx.stream);
}

// ----------------------------------------------------------------------------
// _guessFormat wrapper
// ----------------------------------------------------------------------------

template <typename TValue, typename TStream, typename TCompressionType>
inline bool _guessFormat(VirtualStream<TValue, Input> &, TStream &fileStream, TCompressionType &compressionType)
{
    // peek the first character to initialize the underlying streambuf
    fileStream.rdbuf()->sgetc();
    return guessFormatFromStream(fileStream, compressionType);
}

template <typename TValue, typename TStream, typename TCompressionType>
inline bool _guessFormat(VirtualStream<TValue, Output> &, TStream &, TCompressionType &)
{
    return true;
}

// ----------------------------------------------------------------------------
// Function _getUncompressedBasename()
// ----------------------------------------------------------------------------

// single format
template <typename TFilename, typename TFormat>
inline typename Prefix<TFilename const>::Type
_getUncompressedBasename(TFilename const & fileName, TFormat const & format)
{
    return getBasename(fileName, format);
}

// make sure not to only cut the ".bgzf" extension and not ".bam"
template <typename TFilename>
inline typename Prefix<TFilename const>::Type
_getUncompressedBasename(TFilename const & fileName, BgzfFile const &)
{
    typedef typename Value<TFilename>::Type                                     TValue;
    typedef ModifiedString<TFilename const, ModView<FunctorLowcase<TValue> > >    TLowcase;

    TLowcase lowcaseFileName(fileName);

    if (endsWith(lowcaseFileName, ".bgzf"))
        return prefix(fileName, length(fileName) - 5);

    if (endsWith(lowcaseFileName, ".gz"))
        return prefix(fileName, length(fileName) - 3);

    return prefix(fileName, length(fileName));
}

// TagSelector
template <typename TFilename>
inline typename Prefix<TFilename const>::Type
_getUncompressedBasename(TFilename const & fname, TagSelector<> const & format)
{
    return getBasename(fname, format);
}

template <typename TFilename, typename TTagList>
inline typename Prefix<TFilename const>::Type
_getUncompressedBasename(TFilename const & fname, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        return _getUncompressedBasename(fname, TFormat());
    else
        return _getUncompressedBasename(fname, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// --------------------------------------------------------------------------
// Function open()
// --------------------------------------------------------------------------

/*!
 * @fn VirtualStream#open
 * @brief Open a VirtualStream.
 *
 * @signature bool open(stream, fileStream);
 * @signature bool open(stream, fileName, openMode);
 *
 * @param[in,out] stream        The VirtualStream to open.
 * @param[in,out  fileStream    File stream to attach to. Type: <tt>std::fstream</tt>.
 * @param[in]     fileName      Path to the file to open. Type: <tt>char const *</tt>.
 * @param[in]     openMode      The open mode. Type: <tt>int</tt>.
 * @return bool <tt>true</tt> in the case of success, <tt>false</tt> otherwise.
 */

template <typename TValue, typename TDirection, typename TTraits, typename TStream, typename TCompressionType>
inline SEQAN_FUNC_DISABLE_IF(IsPointer<TStream>, bool)
open(VirtualStream<TValue, TDirection, TTraits> &stream, TStream &fileStream, TCompressionType & compressionType)
{
    SEQAN_ASSERT_MSG(stream.context == NULL, "VirtualStream: close() must be called before re-opening.");

    typedef VirtualStream<TValue, TDirection, TTraits> TVirtualStream;
    typedef typename TVirtualStream::TBufferedStream TBufferedStream;

    // peek the first character to initialize the underlying streambuf (for in_avail)
    SEQAN_IF_CONSTEXPR (IsSameType<TDirection, Input>::VALUE)  // Only getc if input stream.
        fileStream.rdbuf()->sgetc();

    SEQAN_IF_CONSTEXPR (IsSameType<TDirection, Input>::VALUE && !IsSameType<TStream, TBufferedStream>::VALUE)
    {
        if (fileStream.rdbuf()->in_avail() < 2)
        {
            stream.bufferedStream.setStream(fileStream);
            return open(stream, stream.bufferedStream, compressionType);
        }
    }

    VirtualStreamFactoryContext_<TVirtualStream> ctx(fileStream);

    // try to detect/verify format
    if (!_guessFormat(stream, fileStream, compressionType))
        return false;

    // create a new (un)zipper buffer
    stream.context = tagApply(ctx, compressionType);
    if (stream.context == NULL)
        return false;

    SEQAN_ASSERT(stream.context->streamBuf != NULL);
    stream.streamBuf = stream.context->streamBuf;

    // reset our outer stream interface
    stream._init();
    return true;
}

template <typename TValue, typename TDirection, typename TTraits, typename TStream, typename TCompressionType>
inline SEQAN_FUNC_DISABLE_IF(IsPointer<TStream>, bool)
open(VirtualStream<TValue, TDirection, TTraits> &stream, TStream &fileStream, TCompressionType const & compressionType)
{
    assign(stream.format, compressionType);
    return open(stream, fileStream, stream.format);
}

template <typename TValue, typename TStream>
inline SEQAN_FUNC_DISABLE_IF(IsPointer<TStream>, bool)
open(VirtualStream<TValue, Input> &stream, TStream &fileStream)
{
    // detect compression type from file extension
    assign(stream.format, typename StreamFormat<VirtualStream<TValue, Input> >::Type());
    return open(stream, fileStream, stream.format);
}

template <typename TValue, typename TDirection, typename TTraits>
inline bool
open(VirtualStream<TValue, TDirection, TTraits> &stream,
     const char *fileName,
     int openMode = DefaultOpenMode<VirtualStream<TValue, TDirection, TTraits> >::VALUE)
{
    SEQAN_ASSERT_MSG(stream.context == NULL, "VirtualStream: close() must be called before re-opening.");

    typedef VirtualStream<TValue, TDirection, TTraits> TVirtualStream;

    if (!open(stream.file, fileName, openMode))
        return false;

    // detect compression type from file extension
    assign(stream.format, typename StreamFormat<TVirtualStream>::Type());

    if (IsSameType<TDirection, Input>::VALUE && _isPipe(fileName))
        open(stream, stream.file, stream.format);               // read from a pipe (without file extension)
    else
        guessFormatFromFilename(fileName, stream.format);       // read/write from/to a file (with extension)

    VirtualStreamFactoryContext_<TVirtualStream> ctx(stream.file);

    // create a new (un)zipper buffer
    stream.context = tagApply(ctx, stream.format);
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

// --------------------------------------------------------------------------
// Function close()
// --------------------------------------------------------------------------

/*!
 * @fn VirtualStream#close
 * @brief Close a VirtualStream.
 *
 * @signature bool close(stream);
 *
 * @param[in,out] stream The VirtualStream to close.
 * @return bool <tt>true</tt> in the case of success, <tt>false</tt> otherwise.
 */

template <typename TValue, typename TDirection, typename TTraits>
inline bool
close(VirtualStream<TValue, TDirection, TTraits> &stream)
{
    delete stream.context;
    stream.context = NULL;
    stream.streamBuf = NULL;
    assign(stream.format, typename StreamFormat<VirtualStream<TValue, TDirection, TTraits> >::Type());
    return !stream.file.is_open() || close(stream.file);
}

// ----------------------------------------------------------------------------
// Function format()
// ----------------------------------------------------------------------------

/*!
 * @fn VirtualStream#format
 * @brief Return the format of a VirtualStream.
 *
 * @signature TFormat format(stream);
 *
 * @param[in] stream The VirtualStream to check.
 * @return TFormat   The type as returned from @link VirtualStream#StreamFormat @endlink.
 */

template <typename TValue, typename TDirection, typename TTraits>
inline typename StreamFormat<VirtualStream<TValue, TDirection, TTraits> >::Type &
format(VirtualStream<TValue, TDirection, TTraits> &stream)
{
    return stream.format;
}

}  // namespace seqan2

#endif  // #ifndef SEQAN_STREAM_VIRTUAL_STREAM_
