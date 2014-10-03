// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Base class for high level access to formatted files. It supports different
// file formats, compression, auto detection from file extensions and content
// as well as reading/writing from/to stdin/stdout.
// ==========================================================================

#ifndef SEQAN_STREAM_SMART_FILE_H_
#define SEQAN_STREAM_SMART_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSmartFile>
struct FileFormat;

template <typename TSmartFile, typename TStorageSpec>
struct SmartFileContext;

// ============================================================================
// Classes
// ============================================================================

template <typename TObject, typename TStorageSpec>
struct StorageSwitch
{
    typedef TObject Type;
};

template <typename TObject, typename TSpec>
struct StorageSwitch<TObject, Dependent<TSpec> >
{
    typedef TObject * Type;
};

// ----------------------------------------------------------------------------
// Class SmartFile
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec = void>
struct SmartFile
{
    typedef VirtualStream<char, TDirection>                             TStream;
    typedef typename Iterator<TStream, TDirection>::Type                TIter;
    typedef typename FileFormat<SmartFile>::Type                        TFileFormat;
    typedef typename SmartFileContext<SmartFile, Owner<> >::Type        TOwnerContext;
    typedef typename SmartFileContext<SmartFile, Dependent<> >::Type    TDependentContext;

    TStream             stream;
    TIter               iter;
    TFileFormat         format;
    TOwnerContext       data_context;
    TDependentContext   context;

    SmartFile() :
        context(data_context)
    {}

    // filename based c'tors
    explicit
    SmartFile(const char *fileName, int openMode = DefaultOpenMode<SmartFile>::VALUE) :
        context(data_context)
    {
        _open(*this, fileName, openMode, True());
    }

    // stream based c'tors
    template <typename TValue>
    explicit
    SmartFile(std::basic_istream<TValue> &istream,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Input>, IsSameType<TValue, char> >)) :
        context(data_context)
    {
        _open(*this, istream, _mapFileFormatToCompressionFormat(format), True());
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TValue, typename TFormat>
    SmartFile(std::basic_ostream<TValue> &ostream,
              TFormat const &format,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Output>, IsSameType<TValue, char> >)) :
        context(data_context)
    {
        bool success = open(*this, ostream, format);
        ignoreUnusedVariableWarning(dummy);
        ignoreUnusedVariableWarning(success);
        SEQAN_ASSERT(success);
    }

    // now everything given another context
    explicit
    SmartFile(TDependentContext &otherCtx) :
        context(otherCtx)
    {}

    explicit
    SmartFile(SmartFile &other) :
        context(other.context)
    {}

    SmartFile(TDependentContext &otherCtx, const char *fileName, int openMode = DefaultOpenMode<SmartFile>::VALUE) :
        context(otherCtx)
    {
        _open(*this, fileName, openMode, True());
    }

    template <typename TValue>
    explicit
    SmartFile(TDependentContext &otherCtx,
              std::basic_istream<TValue> &istream,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Input>, IsSameType<TValue, char> >)) :
        context(otherCtx)
    {
        _open(*this, istream, _mapFileFormatToCompressionFormat(format), True());
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TValue, typename TFormat>
    SmartFile(TDependentContext &otherCtx,
              std::basic_ostream<TValue> &ostream,
              TFormat const &format,
              SEQAN_CTOR_ENABLE_IF(And<IsSameType<TDirection, Output>, IsSameType<TValue, char> >)) :
        context(otherCtx)
    {
        bool success = open(*this, ostream, format);
        ignoreUnusedVariableWarning(dummy);
        ignoreUnusedVariableWarning(success);
        SEQAN_ASSERT(success);
    }

    ~SmartFile()
    {
        close(*this);
    }

    operator TDependentContext & ()
    {
        return context;
    }

    static std::vector<std::string>
    getFileFormatExtensions()
    {
        std::vector<std::string> extensions;

        _getCompressionExtensions(extensions,
                                  TFileFormat(),
                                  typename FileFormat<TStream>::Type(),
//                                  true);
                                  IsSameType<TDirection, Output>::VALUE);
        return extensions;
    }
};

// ----------------------------------------------------------------------------
// Concepts
// ----------------------------------------------------------------------------

template <typename TFileType, typename TSpec>
SEQAN_CONCEPT_IMPL((SmartFile<TFileType, Input, TSpec>), (InputStreamConcept));

template <typename TFileType, typename TSpec>
SEQAN_CONCEPT_IMPL((SmartFile<TFileType, Output, TSpec>), (OutputStreamConcept));

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DirectionIterator
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
struct DirectionIterator<SmartFile<TFileType, TDirection, TSpec>, TDirection>
{
    typedef typename SmartFile<TFileType, TDirection, TSpec>::TIter Type;
};

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TSmartFile, typename TStorageSpec>
struct SmartFileContext
{
    typedef Nothing Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
struct Value<SmartFile<TFileType, TDirection, TSpec> > :
    Value<typename SmartFile<TFileType, TDirection, TSpec>::TStream> {};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
struct Position<SmartFile<TFileType, TDirection, TSpec> > :
    Position<typename SmartFile<TFileType, TDirection, TSpec>::TStream> {};

// ----------------------------------------------------------------------------
// Metafunction DefaultOpenMode
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec, typename TDummy>
struct DefaultOpenMode<SmartFile<TFileType, TDirection, TSpec>, TDummy> :
    DefaultOpenMode<typename SmartFile<TFileType, TDirection, TSpec>::TStream, TDummy> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function directionIterator()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline typename SmartFile<TFileType, TDirection, TSpec>::TIter &
directionIterator(SmartFile<TFileType, TDirection, TSpec> & file, TDirection const &)
{
    return file.iter;
}

// ----------------------------------------------------------------------------
// Function format()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline typename FileFormat<SmartFile<TFileType, TDirection, TSpec> >::Type &
format(SmartFile<TFileType, TDirection, TSpec> & file)
{
    return file.format;
}

// ----------------------------------------------------------------------------
// Function setFormat()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec, typename TFormat>
inline void
setFormat(SmartFile<TFileType, TDirection, TSpec> & file, TFormat format)
{
    assign(file.format, format);
}

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TSpec>
inline bool guessFormat(SmartFile<TFileType, Input, TSpec> & file)
{
    return guessFormatFromStream(file.stream, file.format);
}

template <typename TFileType, typename TSpec>
inline bool guessFormat(SmartFile<TFileType, Output, TSpec> &)
{
    return true;
}

// --------------------------------------------------------------------------
// _mapFileFormatToCompressionFormat
// --------------------------------------------------------------------------

template <typename TFormat>
inline Nothing
_mapFileFormatToCompressionFormat(TFormat)
{
    return Nothing();
}

template <typename TCompressedFileTypes>
inline void
_mapFileFormatToCompressionFormat(TagSelector<TCompressedFileTypes> &,
                                  TagSelector<void> const &)
{}

template <typename TCompressedFileTypes, typename TTagList>
inline void
_mapFileFormatToCompressionFormat(TagSelector<TCompressedFileTypes> & result,
                                  TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;
    typedef typename TagSelector<TTagList>::Base TBase;

    if (isEqual(format, TFormat()))
        assign(result, _mapFileFormatToCompressionFormat(TFormat()));
    else
        _mapFileFormatToCompressionFormat(result, static_cast<TBase const &>(format));
}

template <typename TFileFormatList>
inline TagSelector<CompressedFileTypes>
_mapFileFormatToCompressionFormat(TagSelector<TFileFormatList> format)
{
    TagSelector<CompressedFileTypes> compressionType;
    _mapFileFormatToCompressionFormat(compressionType, format);
    return compressionType;
}


template <typename TSmartFile, typename TFormat>
inline void
_checkThatStreamOutputFormatIsSet(TSmartFile const &, TFormat const &)
{}

template <typename TFileType, typename TSpec, typename TFileFormatList>
inline void
_checkThatStreamOutputFormatIsSet(SmartFile<TFileType, Output, TSpec> const &, TagSelector<TFileFormatList> const &format)
{
    if (value(format) < 0)
        SEQAN_FAIL("SmartFile: File format not specified, use setFormat().");
}

// --------------------------------------------------------------------------
// Function open(stream)
// --------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec,
          typename TStream, typename TCompressionFormat, typename TThrowExceptions>
inline bool _open(SmartFile<TFileType, TDirection, TSpec> & file,
                  TStream &stream,
                  TCompressionFormat const &compressionFormat,
                  TThrowExceptions = True())
{
    if (!open(file.stream, stream, compressionFormat))
    {
        if (TThrowExceptions::VALUE)
            SEQAN_THROW(UnknownFileFormat());
        return false;
    }

    if (!guessFormat(file))
    {
        if (TThrowExceptions::VALUE)
            SEQAN_THROW(UnknownFileFormat());
        return false;
    }

    file.iter = directionIterator(file.stream, TDirection());

    return true;
}

template <typename TFileType, typename TDirection, typename TSpec, typename TStream>
inline bool open(SmartFile<TFileType, TDirection, TSpec> & file,
                 TStream &stream)
{
    return _open(file, stream, _mapFileFormatToCompressionFormat(file.format), False());
}

template <typename TFileType, typename TSpec, typename TStream, typename TFormat_>
inline bool open(SmartFile<TFileType, Output, TSpec> & file,
                 TStream &stream,
                 Tag<TFormat_> format)
{
    setFormat(file, format);
    return _open(file, stream, _mapFileFormatToCompressionFormat(file.format), False());
}

template <typename TFileType, typename TSpec, typename TStream, typename TFormats>
inline bool open(SmartFile<TFileType, Output, TSpec> & file,
                 TStream &stream,
                 TagSelector<TFormats> format)
{
    setFormat(file, format);
    return _open(file, stream, _mapFileFormatToCompressionFormat(file.format), False());
}

// ----------------------------------------------------------------------------
// Function open(fileName)
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec, typename TThrowExceptions>
inline bool _open(SmartFile<TFileType, TDirection, TSpec> & file,
                  const char *fileName,
                  int openMode = DefaultOpenMode<SmartFile<TFileType, TDirection, TSpec> >::VALUE,
                  TThrowExceptions = True())
{
    if (!open(file.stream, fileName, openMode))
    {
        if (TThrowExceptions::VALUE)
            SEQAN_THROW(FileOpenError(fileName));
        return false;
    }

    if (IsSameType<TDirection, Input>::VALUE && _isPipe(fileName))
    {
        if (!guessFormat(file))
        {
            if (TThrowExceptions::VALUE)                        // read from a pipe (without file extension)
                SEQAN_THROW(UnknownFileFormat());
            return false;
        }
    }
    else
    {
        Prefix<const char *>::Type basename = _getUncompressedBasename(fileName, format(file.stream));
        if (!guessFormatFromFilename(basename, file.format))    // read/write from/to a file (with extension)
        {
            close(file.stream);
            if (TThrowExceptions::VALUE)
                SEQAN_THROW(UnknownExtensionError(fileName));
            return false;
        }
    }

    file.iter = directionIterator(file.stream, TDirection());
    return true;
}

template <typename TFileType, typename TDirection, typename TSpec>
inline bool open(SmartFile<TFileType, TDirection, TSpec> & file,
                 const char *fileName,
                 int openMode = DefaultOpenMode<SmartFile<TFileType, TDirection, TSpec> >::VALUE)
{
    return _open(file, fileName, openMode, False());
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline bool close(SmartFile<TFileType, TDirection, TSpec> & file)
{
    file.iter = typename DirectionIterator<SmartFile<TFileType, TDirection, TSpec>, TDirection>::Type();
    return close(file.stream);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename SmartFile<TFileType, TDirection, TSpec>::TStream> >, bool)
atEnd(SmartFile<TFileType, TDirection, TSpec> const & file)
{
    return atEnd(file.iter);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TSpec>
inline typename Position<SmartFile<TFileType, Output, TSpec> >::Type
position(SmartFile<TFileType, Output, TSpec> & file)
{
    return file.stream.tellp();
}

template <typename TFileType, typename TSpec>
inline typename Position<SmartFile<TFileType, Input, TSpec> >::Type
position(SmartFile<TFileType, Input, TSpec> & file)
{
    return file.stream.tellg();
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TSpec, typename TPosition>
inline bool
setPosition(SmartFile<TFileType, Output, TSpec> & file, TPosition pos)
{
    return (TPosition)file.stream.rdbuf()->pubseekpos(pos, std::ios_base::out) == pos;
}

template <typename TFileType, typename TSpec, typename TPosition>
inline bool
setPosition(SmartFile<TFileType, Input, TSpec> & file, TPosition pos)
{
    return (TPosition)file.stream.rdbuf()->pubseekpos(pos, std::ios_base::in) == pos;
}

// ----------------------------------------------------------------------------
// Function context()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline typename SmartFileContext<SmartFile<TFileType, TDirection, TSpec>, Dependent<> >::Type &
context(SmartFile<TFileType, TDirection, TSpec> & file)
{
    return file.context;
}

template <typename TFileType, typename TDirection, typename TSpec>
inline typename SmartFileContext<SmartFile<TFileType, TDirection, TSpec>, Dependent<> >::Type const &
context(SmartFile<TFileType, TDirection, TSpec> const & file)
{
    return file.context;
}

// ----------------------------------------------------------------------------
// Function getFileFormatExtensions()
// ----------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Function _getFileFormatExtensions()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TFormat_, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    Tag<TFormat_> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly,
    Nothing)
{
    typedef Tag<TFormat_> TFormat;

    std::vector<std::string> compressionExtensions;
    _getFileFormatExtensions(compressionExtensions, compress, primaryExtensionOnly);

    unsigned len = (primaryExtensionOnly)? 1 : sizeof(FileFormatExtensions<TFormat>::VALUE) / sizeof(char*);
    for (unsigned i = 0; i < len; ++i)
        for (unsigned j = 0; j < compressionExtensions.size(); ++j)
        {
            size_t jj = (j == 0)? compressionExtensions.size() - 1 : j - 1;    // swap first and last compression extension
            appendValue(stringSet, (std::string)FileFormatExtensions<TFormat>::VALUE[i] + compressionExtensions[jj]);
        }
}

template <typename TStringSet, typename TFormat_, typename TCompressionFormats, typename TCompression_>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    Tag<TFormat_> const & formatTag,
    TCompressionFormats const &,
    bool primaryExtensionOnly,
    Tag<TCompression_>)
{
    _getFileFormatExtensions(stringSet, formatTag, primaryExtensionOnly);
}

template <typename TStringSet, typename TFormat_, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    Tag<TFormat_> const & formatTag,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly)
{
    _getCompressionExtensions(stringSet, formatTag, compress, primaryExtensionOnly, _mapFileFormatToCompressionFormat(formatTag));
}

template <typename TStringSet, typename TTag, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    TagList<TTag, void> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly = false)
{
    _getCompressionExtensions(stringSet, TTag(), compress, primaryExtensionOnly);
}

template <typename TStringSet, typename TTag, typename TSubList, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    TagList<TTag, TSubList> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly = false)
{
    _getCompressionExtensions(stringSet, TTag(), compress, primaryExtensionOnly);
    _getCompressionExtensions(stringSet, TSubList(), compress, primaryExtensionOnly);
}

template <typename TStringSet, typename TTagList, typename TCompressionFormats>
inline void
_getCompressionExtensions(
    TStringSet &stringSet,
    TagSelector<TTagList> const & /*formatTag*/,
    TCompressionFormats const & compress,
    bool primaryExtensionOnly = false)
{
    _getCompressionExtensions(stringSet, TTagList(), compress, primaryExtensionOnly);
}



template <typename TFileType, typename TDirection, typename TSpec>
static std::vector<std::string>
getFileFormatExtensions(SmartFile<TFileType, TDirection, TSpec> const & file)
{
    return file.getFileFormatExtensions();
}


}  // namespace seqan

#endif // SEQAN_STREAM_SMART_FILE_H_
