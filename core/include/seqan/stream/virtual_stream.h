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
// Virtual stream class to subsume different implementations, e.g. for raw
// and compressed streams.
// ==========================================================================

#ifndef SEQAN_STREAM_VIRTUAL_STREAM_
#define SEQAN_STREAM_VIRTUAL_STREAM_

#include <functional>

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================


// ============================================================================
// Classes
// ============================================================================


template <typename TTag, typename T>
struct FileExtensions<TTag, T>;

template <typename T>
struct FileExtensions<BZ2File, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileExtensions<BZ2File, T>::VALUE[2] = {
    ".bz2",      // default output extension
    ".bzip2" };


template <typename T>
struct FileExtensions<GZFile, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileExtensions<GZFile, T>::VALUE[2] = {
    ".gz",      // default output extension
    ".gzip" };



typedef
    TagList<GZFile,
    TagList<BZ2File,
    TagList<Raw> > > SeqFormats;  // if TagSelector is set to -1, the file format is auto-detected

typedef TagSelector<SeqFormats> CompressedFileTypes;



template <typename TStream, typename TStreamBuffer>
struct StreamFactoryContext
{
    typedef TStreamBuffer TResultType;

    TStream &stream;
    TStream(TStream &stream):
        stream(stream) {}
};


template <typename TValue, typename TDirection>
class VirtualStream
{
public:
    typedef std::basic_streambuf<TValue> TStreamBuffer;
    typedef If<IsSameType<TDirection, Input>, std::basic_istream<TValue>, std::basic_ostream<TValue> > TStream;
    typedef TStreamFactoryFunctor<TStreamBuffer>

    std::fstream    file;
    TStreamBuffer   *streamBuf;
    bool            ownerOfStreamBuf;

    VirtualStream():
        streamBuf()
    {}

    VirtualStream(TStreamBuffer &streamBuf):
        streamBuf(streamBuf),
        ownerOfStreamBuf(false)
    {}

    VirtualStream(TStream &stream):
        streamBuf()
    {
        open(*this, stream);
    }

    VirtualStream(const char *fileName, int openMode):
        streamBuf(),
        ownerOfStreamBuf(false)
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
};

template <typename TValue, typename TDirection>
inline bool
close(VirtualStream<TValue, TDirection> &stream)
{
    if (ownerOfStreamBuf)
        delete streamBuf;

    stream.streamBuf = NULL;
    return true;
}

template <typename TValue, typename TDirection>
inline bool
open(VirtualStream<TValue, TDirection> &stream, const char *fileName, int openMode)
{
    if (!open(stream.file, fileName, openMode))
        return false;

    TagSelector<CompressedFileTypes> fileType;
    guessFormatFromFilename(fileName, fileType);

    if (value(fileType) == Find<CompressedFileTypes, Raw>::VALUE)
    {
        streamBuf = file.rdbuf();
        ownerOfStreamBuf = false;
    }
    else
    {
        StreamFactoryContext<std::fstream, TStreamBuffer> ctx(file);
        streamBuf = tagApply(ctx, fileType);    // create a new unzipper buffer
        if (streamBuf == NULL)
        {
            close(stream.file);
            return false;
        }
        ownerOfStreamBuf = true;
    }
}


template <typename TStream, typename TStreamBuffer>
inline TStreamBuffer*
tagApply(StreamFactoryContext<TStream, TStreamBuffer> &ctx, GZFile)
{
    return new basic_unzip_streambuf<TValue>(ctx.stream, 15, 4096, 4096);
}

template <typename TStream, typename TStreamBuffer>
inline TStreamBuffer*
tagApply(StreamFactoryContext<TStream, TStreamBuffer> &ctx, BZ2File)
{
    return new basic_unbzip2_streambuf<TValue>(ctx.stream, 0, false, 4096, 4096);
}






template <typename TContext>
inline typename TContext::TResultType
tagApply(TContext const &, TagSelector<>)
{
    return TContext::TResultType();
}

template <typename TContext, typename TTagList>
inline typename TContext::TResultType
tagApply(TContext &ctx, TagSelector<TTagList> &format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
    {
        return tagApply(ctx, TFormatTag()))
    }
    return tagApply(ctx, static_cast<typename TagSelector<TTagList>::Base &>(format));
}



}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_VIRTUAL_STREAM_
