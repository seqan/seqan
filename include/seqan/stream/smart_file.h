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

// TODO(holtgrew): Document or mark as internal by adding underscore. Add an explaining sentence anyway.

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

// TODO(holtgrew): Dave/Enrico, please double-check, especially the context stuff.
// TODO(holtgrew): We definitely need extensive explanation and tutorials for this, especially the context stuff.

/*!
 * @class SmartFile
 * @implements InputStreamConcept
 * @implements OutputStreamConcept
 * @headerfile <seqan/stream.h>
 * @brief Base class for high-level I/O functions.
 *
 * @signature template <typename TFileType, typename TDirection[, typename TSpec]>
 *            struct SmartFile;
 *
 * @tparam TFileType  A type specifying the file format.
 * @tparam TDirection The direction of the file, one of @link DirectionTags#Input Input
 *                    @endlink and @link DirectionTags#Output @endlink.
 * @tparam TSpec      A tag for the specialization, defauls to <tt>void</tt>.
 *
 * It depends on <tt>TDirection</tt> whether it implements @link InputStreamConcept @endlink or
 * @link OutputStreamConcept @endlink.
 *
 * SmartFile bundles the necessary information for easily building (aka metaprogramming) easy-to use high-level I/O
 * functionality:
 *
 * <ul>
 * <li>An underlying @link VirtualStream @endlink that is used for the I/O and allows attaching to
 *     already open streams as well as opening and owning the stream.  Further, transparent
 *     access to compressed files is possible.</li>
 * <li>A @link StreamConcept#DirectionIterator direction iterator @endlink for reading/writing.</li>
 * <li>A <b>context</b> used for the I/O.  The context object can be retrieved using @link SmartFile#context context
 *     @endlink.  This context can be used for storing contig and sample names, in the case of VCF I/O, for example,
 *     @link NameStoreCache @endlink objects, and buffers.  Each SmartFile has exactly one context.</li>
 * <li>The context type can be constructed both as dependent and independent using @link SmartFile#SmartFileContext
 *     SmartFileContext @endlink.  Each SmartFile stores an independent context and a dependent context.  When
 *     constructing with another SmartFile or another SmartFile's context, the dependent context of this SmartFile is
 *     set to depend on the dependent context of the other Smart File.  Otherwise, it is set to depend on this
 *     SmartFile's independent context.  This allows chaining of contexts, e.g. for updating the contig names of an
 *     output smart file from an input smart file.</li>
 * </ul>
 */

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

    // TODO(holtgrew): Only documented first few constructors, must be extended after describing the contexts concept.
    /*!
     * @fn SmartFile::SmartFile
     * @brief Provides default construction.
     *
     * @signature SmartFile::SmartFile();
     * @signature SmartFile::SmartFile(fileName[, openMode]);
     * @signature SmartFile::SmartFile(stream);
     * @signature SmartFile::SmartFile(other);
     * @signature SmartFile::SmartFile(otherContext);
     * @signature SmartFile::SmartFile(otherContext, fileName[, openMode]);
     * @signature SmartFile::SmartFile(otherContext, stream);
     *
     * @param[in] fileName     Path to file to open, <tt>char const *</tt>.
     * @param[in] openMode     Optionally, the file open mode, default obtained from <tt>TDirection</tt>.
     * @param[in] stream       A <tt>std::basic_istream&lt;&gt;</tt> to read from or <tt>std::basic_ostream&lt;&gt;</tt>
     *                         to write to, depending on <tt>TDirection</tt>.
     * @param[in] other        A second SmartFile, this SmartFile's dependent context will depend on <tt>other</tt>'s
     *                         dependent context.
     * @param[in] otherContext The dependent context of another SmartFile, this SmartFile's dependent context will depend on <tt>otherContext</tt>.
     *
     * @throw IOError The variants that accept the <tt>fileName</tt> parameter throw an exception of this
     *                type in case opening the file fails.
     */
    SmartFile() : context(data_context)
    {}

    // filename based c'tors
    explicit SmartFile(const char * fileName, int openMode = DefaultOpenMode<SmartFile>::VALUE) :
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

    /*!
     * @fn SmartFile::operatorTDependentContext SmartFile::operator TDependentContext
     * @brief Allows conversion to a dependent context for the SmartFile
     * @signature TDependentContext & SmartFile::operator TDependentContext();
     *
     * @return TDependentContext The dependent context of this SmartFile.
     */

    operator TDependentContext & ()
    {
        return context;
    }

    /*!
     * @fn SmartFile::getFileExtensions
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

/*!
 * @mfn SmartFile#SmartFileContext
 * @brief Returns the context type for a SmartFile.
 *
 * @signature SmartFileContext<TSmartFile, TStorageSpec>::Type
 *
 * @tparam TSmartFile   The SmartFile to query.
 * @tparam TStorageSpec The storage specification, passed as specialization to any @link StringSet @endlink
 *                      contained in the context.
 * @tparam Type         The resulting smart file context type.
 */

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
// Helper Function _throwIf().
// ----------------------------------------------------------------------------

// Helper functions that reduce number of "uncaught exception" false positives in static analysis tools.
template <typename TException> void _throwIf(TException const & e, True const & /*tag*/) { SEQAN_THROW(e); }
template <typename TException> void _throwIf(TException const & /*e*/, False const & /*tag*/) { /*nop*/ }

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

// TODO(holtgrew): Document me, including demo snippet showing isEqual()...

template <typename TFileType, typename TDirection, typename TSpec>
inline typename FileFormat<SmartFile<TFileType, TDirection, TSpec> >::Type &
format(SmartFile<TFileType, TDirection, TSpec> & file)
{
    return file.format;
}

// ----------------------------------------------------------------------------
// Function setFormat()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document me, including demo snippet.

template <typename TFileType, typename TDirection, typename TSpec, typename TFormat>
inline void
setFormat(SmartFile<TFileType, TDirection, TSpec> & file, TFormat format)
{
    assign(file.format, format);
}

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document me, including demo snippet.

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
// Function open()
// --------------------------------------------------------------------------

// TODO(holtgrew): Complete documentation, the combination of formats etc are a bit hard to grasp quickly.

/*!
 * @fn SmartFile#open
 * @brief Open a SmartFile
 */

template <typename TFileType, typename TDirection, typename TSpec,
          typename TStream, typename TCompressionFormat, typename TThrowExceptions>
inline bool _open(SmartFile<TFileType, TDirection, TSpec> & file,
                  TStream &stream,
                  TCompressionFormat const &compressionFormat,
                  TThrowExceptions = True())
{
    if (!open(file.stream, stream, compressionFormat))
    {
        _throwIf(UnknownFileFormat(), TThrowExceptions());
        return false;
    }

    if (!guessFormat(file))
    {
        _throwIf(UnknownFileFormat(), TThrowExceptions());
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
        _throwIf(FileOpenError(fileName), TThrowExceptions());
        return false;
    }

    if (IsSameType<TDirection, Input>::VALUE && _isPipe(fileName))
    {
        if (!guessFormat(file))
        {
            // read from a pipe (without file extension)
            _throwIf(UnknownFileFormat(), TThrowExceptions());
            return false;
        }
    }
    else
    {
        Prefix<const char *>::Type basename = _getUncompressedBasename(fileName, format(file.stream));
        if (!guessFormatFromFilename(basename, file.format))    // read/write from/to a file (with extension)
        {
            close(file.stream);
            _throwIf(UnknownExtensionError(fileName), TThrowExceptions());
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

/*!
 * @fn SmartFile#close
 * @brief Close a SmartFile.
 *
 * @signature bool close(file);
 *
 * @param[in,out] file The SmartFile to close.
 * @return bool <tt>true</tt> in the case of success, <tt>false</tt> otherwise.
 */

template <typename TFileType, typename TDirection, typename TSpec>
inline bool close(SmartFile<TFileType, TDirection, TSpec> & file)
{
    setFormat(file, typename FileFormat<SmartFile<TFileType, TDirection, TSpec> >::Type());
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

/*!
 * @fn SmartFile#context
 * @brief Return the SmartFile's dependent context object.
 *
 * @signature TDependentContext & context(smartFile);
 *
 * @param[in,out] smartFile The SmartFile to query for its context.
 *
 * @return TDependentContext The dependent context, type as returned from @link SmartFile#SmartFileContext @endlink.
 */

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

// --------------------------------------------------------------------------
// Function _getCompressionFormatExtensions()
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
    _getFileExtensions(compressionExtensions, compress, primaryExtensionOnly);

    unsigned len = (primaryExtensionOnly)? 1 : sizeof(FileExtensions<TFormat>::VALUE) / sizeof(char*);
    for (unsigned i = 0; i < len; ++i)
        for (unsigned j = 0; j < compressionExtensions.size(); ++j)
        {
            size_t jj = (j == 0)? compressionExtensions.size() - 1 : j - 1;    // swap first and last compression extension
            appendValue(stringSet, (std::string)FileExtensions<TFormat>::VALUE[i] + compressionExtensions[jj]);
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
    _getFileExtensions(stringSet, formatTag, primaryExtensionOnly);
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

// ----------------------------------------------------------------------------
// Function getFileExtensions()
// ----------------------------------------------------------------------------

/*!
 * @fn SmartFile#getFileExtensions
 * @brief Static function that returns a list of allowed file format extension.
 *
 * @signature TExtensionVector getFileExtensions(smartFile)
 *
 * @param[in] smartFile The SmartFile to query.
 * @return TExtensionVector A <tt>std::vector&lt;std::string&gt;</tt> with the allowed file extensions.
 *
 * This is a shortcut to @link SmartFile#getFileExtensions @endlink.
 */

template <typename TFileType, typename TDirection, typename TSpec>
static std::vector<std::string>
getFileExtensions(SmartFile<TFileType, TDirection, TSpec> const & file)
{
    return file.getFileExtensions();
}

}  // namespace seqan

#endif // SEQAN_STREAM_SMART_FILE_H_
