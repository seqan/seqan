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

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SmartFile
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec = void>
struct SmartFile;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TSmartFile>
struct FileFormats;

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
// Function setFormat()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec, typename TFormat>
inline void setFormat(SmartFile<TFileType, TDirection, TSpec> & file, TFormat format)
{
    assign(file.format, format);
}

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TSpec>
inline bool guessFormat(SmartFile<TFileType, Input, TSpec> & file)
{
    return guessFormat(file.stream, file.format);
}

template <typename TFileType, typename TSpec>
inline bool guessFormat(SmartFile<TFileType, Output, TSpec> &)
{
    return true;
}

// --------------------------------------------------------------------------
// _mapFileFormatToCompressionFormat
// --------------------------------------------------------------------------

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


// --------------------------------------------------------------------------
// Function open(stream)
// --------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec, typename TStream>
inline bool open(SmartFile<TFileType, TDirection, TSpec> & file,
                 TStream &stream)
{
    typedef typename SmartFile<TFileType, TDirection, TSpec>::TIter TIter;

    if (!open(file.stream, stream, _mapFileFormatToCompressionFormat(file.format)))
        return false;

    if (!guessFormat(file))
        return false;

    file.iter = TIter(file.stream);

    return true;
}

template <typename TFileType, typename TSpec, typename TStream, typename TFormat>
inline bool open(SmartFile<TFileType, Output, TSpec> & file,
                 TStream &stream,
                 TFormat format)
{
    setFormat(file, format);
    return open(file, stream);
}

// ----------------------------------------------------------------------------
// Function open(fileName)
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline bool open(SmartFile<TFileType, TDirection, TSpec> & file,
                 const char *fileName,
                 int openMode = DefaultOpenMode<SmartFile<TFileType, TDirection, TSpec> >::VALUE)
{
    typedef typename SmartFile<TFileType, TDirection, TSpec>::TIter TIter;

    if (!guessFormatFromFilename(fileName, file.format))
        return false;

    if (!open(file.stream, fileName, openMode))
        return false;

    file.iter = TIter(file.stream);

    return true;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline bool close(SmartFile<TFileType, TDirection, TSpec> & file)
{
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

}  // namespace seqan

#endif // SEQAN_STREAM_SMART_FILE_H_
