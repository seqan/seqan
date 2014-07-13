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

#ifndef SEQAN_BAM_IO_BAM_FILE_H_
#define SEQAN_BAM_IO_BAM_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TDirection, typename TSpec>
struct BamFile;

// ============================================================================
// Typedefs
// ============================================================================

template <typename TBamFile>
struct BamFileFormats
{
#if SEQAN_HAS_ZLIB
    typedef TagSelector<
                TagList<Bam,
                TagList<Sam
                > >
            > Type;
#else
    typedef Sam Type;
#endif
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class BamFile
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec = void>
struct BamFile
{
    typedef VirtualStream<char, TDirection>                 TStream;
    typedef typename Iterator<TStream, TDirection>::Type    TIter;
    typedef typename BamFileFormats<BamFile>::Type          TBamFileFormats;

    typedef StringSet<CharString>                           TNameStore;
    typedef NameStoreCache<TNameStore>                      TNameStoreCache;
    typedef BamIOContext<TNameStore, TNameStoreCache>       TBamIOContext;

    TNameStore      nameStore;
    TNameStoreCache nameStoreCache;
    TBamIOContext   bamIOCtx;
    TBamIOContext   *ctx;

    TStream         stream;
    TIter           iter;
    TBamFileFormats format;

    BamFile() :
        nameStoreCache(nameStore),
        bamIOCtx(nameStore, nameStoreCache),
        ctx(&bamIOCtx),
        iter(stream)
    {}

    BamFile(TBamIOContext &_ctx) :
        nameStoreCache(nameStore),
        ctx(&_ctx),
        iter(stream)
    {}

    BamFile(const char *fileName, int openMode = DefaultOpenMode<BamFile>::VALUE) :
        nameStoreCache(nameStore),
        bamIOCtx(nameStore, nameStoreCache),
        ctx(&bamIOCtx),
        iter(stream)
    {
        if (!open(*this, fileName, openMode))
            throw IOError(std::string("Could not open file ") + fileName);
    }

    BamFile(TBamIOContext &_ctx, const char *fileName, int openMode = DefaultOpenMode<BamFile>::VALUE) :
        nameStoreCache(nameStore),
        ctx(&_ctx),
        iter(stream)
    {
        if (!open(*this, fileName, openMode))
            throw IOError(std::string("Could not open file ") + fileName);
    }

    ~BamFile()
    {
        close(*this);
    }

    operator TBamIOContext & ()
    {
        return *ctx;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultOpenMode
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TDummy>
struct DefaultOpenMode<BamFile<TDirection, TSpec>, TDummy> :
    DefaultOpenMode<typename BamFile<TDirection, TSpec>::TStream, TDummy> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setFormat()
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TFormat>
inline void setFormat(BamFile<TDirection, TSpec> & file, TFormat)
{
    assign(file.format, TFormat());
}

// ----------------------------------------------------------------------------
// Function open(fileName)
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
inline bool open(BamFile<TDirection, TSpec> & file,
                 const char *fileName,
                 int openMode = DefaultOpenMode<BamFile<TDirection, TSpec> >::VALUE)
{
    typedef typename BamFile<TDirection, TSpec>::TIter TIter;

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

template <typename TDirection, typename TSpec>
inline bool close(BamFile<TDirection, TSpec> & file)
{
    return close(file.stream);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename BamFile<TDirection, TSpec>::TStream> >, bool)
atEnd(BamFile<TDirection, TSpec> const & file)
{
    return atEnd(file.iter);
}

// ----------------------------------------------------------------------------
// Function read(); BamHeader
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
inline void
readRecord(BamHeader & /* header */,
           BamIOContext<TNameStore, TNameStoreCache> & /* context */,
           TForwardIter & /* iter */,
           TagSelector<> const & /* format */)
{}

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TTagList>
inline void
readRecord(BamHeader & header,
           BamIOContext<TNameStore, TNameStoreCache> & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
        readRecord(header, context, iter, TFormatTag());
    else
        readRecord(header, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient BamFile variant
template <typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename BamFile<TDirection, TSpec>::TStream> >, void)
readRecord(BamHeader & header, BamFile<TDirection, TSpec> & file)
{
    readRecord(header, *file.ctx, file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function read(); BamAlignmentRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
inline void
readRecord(BamAlignmentRecord & /* record */,
           BamIOContext<TNameStore, TNameStoreCache> & /* context */,
           TForwardIter & /* iter */,
           TagSelector<> const & /* format */)
{}

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TTagList>
inline void
readRecord(BamAlignmentRecord & record,
           BamIOContext<TNameStore, TNameStoreCache> & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
        readRecord(record, context, iter, TFormatTag());
    else
        readRecord(record, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient BamFile variant
template <typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename BamFile<TDirection, TSpec>::TStream> >, void)
readRecord(BamAlignmentRecord & record, BamFile<TDirection, TSpec> & file)
{
    readRecord(record, *file.ctx, file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); BamHeader
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TTarget, typename TNameStore, typename TNameStoreCache>
inline void
write(TTarget & /* target */,
      BamHeader const & /* header */,
      BamIOContext<TNameStore, TNameStoreCache> & /* context */,
      TagSelector<> const & /* format */)
{}

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TTagList>
inline void
write(TTarget & target,
      BamHeader const & header,
      BamIOContext<TNameStore, TNameStoreCache> & context,
      TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
        write(target, header, context, TFormatTag());
    else
        write(target, header, context, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient BamFile variant
template <typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename BamFile<TDirection, TSpec>::TStream> >, void)
write(BamFile<TDirection, TSpec> & file, BamHeader & header)
{
    write(file.iter, header, *file.ctx, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); BamAlignmentRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TTarget, typename TNameStore, typename TNameStoreCache>
inline void
write(TTarget & /* target */,
      BamAlignmentRecord & /* record */,
      BamIOContext<TNameStore, TNameStoreCache> & /* context */,
      TagSelector<> const & /* format */)
{}

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TTagList>
inline void
write(TTarget & target,
      BamAlignmentRecord & record,
      BamIOContext<TNameStore, TNameStoreCache> & context,
      TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == LENGTH<TTagList>::VALUE - 1)
        write(target, record, context, TFormatTag());
    else
        write(target, record, context, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename BamFile<TDirection, TSpec>::TStream> >, void)
write(BamFile<TDirection, TSpec> & file, BamAlignmentRecord & record)
{
    write(file.iter, record, *file.ctx, file.format);
}

}  // namespace seqan

#endif // SEQAN_BAM_IO_BAM_FILE_H_
