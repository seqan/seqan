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
// Smart file for reading/writing files in SAM or BAM format.
// ==========================================================================

#ifndef SEQAN_BAM_IO_BAM_FILE_H_
#define SEQAN_BAM_IO_BAM_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

typedef SmartFile<Bam, Input>   BamFileIn;
typedef SmartFile<Bam, Output>  BamFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TOwnerSpec>
struct SmartFileContext<SmartFile<Bam, TDirection, TSpec>, TOwnerSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    typedef BamIOContext<TNameStore, TNameStoreCache, TOwnerSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormats<SmartFile<Bam, TDirection, TSpec> >
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

// --------------------------------------------------------------------------
// _mapBamFormatToCompressionFormat
// --------------------------------------------------------------------------

inline BgzfFile
_mapFileFormatToCompressionFormat(Bam)
{
    return BgzfFile();
}

inline Nothing
_mapFileFormatToCompressionFormat(Sam)
{
    return Nothing();
}

// ----------------------------------------------------------------------------
// Function read(); BamHeader
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(BamHeader & /* header */,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
           TForwardIter & /* iter */,
           TagSelector<> const & /* format */)
{
    SEQAN_FAIL("BamFileIn: File format not specified.");
}

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
readRecord(BamHeader & header,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(header, context, iter, TFormat());
    else
        readRecord(header, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient BamFile variant
template <typename TSpec>
inline void
readRecord(BamHeader & header, SmartFile<Bam, Input, TSpec> & file)
{
    readRecord(header, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function read(); BamAlignmentRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(BamAlignmentRecord & /* record */,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
           TForwardIter & /* iter */,
           TagSelector<> const & /* format */)
{
    SEQAN_FAIL("BamFileIn: File format not specified.");
}

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
readRecord(BamAlignmentRecord & record,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(record, context, iter, TFormat());
    else
        readRecord(record, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient BamFile variant
template <typename TSpec>
inline void
readRecord(BamAlignmentRecord & record, SmartFile<Bam, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); BamHeader
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
write(TTarget & /* target */,
      BamHeader const & /* header */,
      BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
      TagSelector<> const & /* format */)
{
    SEQAN_FAIL("BamFileOut: File format not specified.");
}

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
write(TTarget & target,
      BamHeader const & header,
      BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
      TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        write(target, header, context, TFormat());
    else
        write(target, header, context, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient BamFile variant
template <typename TSpec>
inline void
writeRecord(SmartFile<Bam, Output, TSpec> & file, BamHeader & header)
{
    write(file.iter, header, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function write(); BamAlignmentRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
write(TTarget & /* target */,
      BamAlignmentRecord & /* record */,
      BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
      TagSelector<> const & /* format */)
{
    SEQAN_FAIL("BamFileOut: File format not specified.");
}

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
write(TTarget & target,
      BamAlignmentRecord & record,
      BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
      TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        write(target, record, context, TFormat());
    else
        write(target, record, context, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
writeRecord(SmartFile<Bam, Output, TSpec> & file, BamAlignmentRecord & record)
{
    write(file.iter, record, context(file), file.format);
}

}  // namespace seqan

#endif // SEQAN_BAM_IO_BAM_FILE_H_
