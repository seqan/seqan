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

/*!
 * @class BamFileIn
 * @headerfile <seqan/bam_io.h>
 * @brief Class that provides an easy to use interface for reading SAM and BAM files.
 *
 * @signature class BamFileIn;
 *
 * @section Example
 *
 * Read SAM or BAM files.
 *
 * @include demos/bam_io/bam_file_in.cpp
 *
 * The output is as follows:
 *
 * @include demos/bam_io/bam_file_in.cpp.stdout
 */

/*!
 * @fn BamFileIn::BamFileIn
 * @brief Constructor
 *
 * @signature BamFileIn::BamFileIn([fileName[, openMode]]);
 *
 * @param[in] fileName The path to the SAM or BAM file to load, <tt>char const *</tt>.
 * @param[in] openMode The open mode. Type: <tt>int</tt>.
 */

typedef SmartFile<Bam, Input> BamFileIn;

/*!
 * @class BamFileOut
 * @headerfile <seqan/bam_io.h>
 * @brief Class that provides an easy to use interface for writing SAM and BAM files.
 *
 * @signature class BamFileOut;
 */

/*!
 * @fn BamFileIn::BamFileOut
 * @brief Constructor
 *
 * @signature BamFileOut::BamFileOut([fileName[, openMode]]);
 *
 * @param[in] fileName The path to the SAM or BAM file to write, <tt>char const *</tt>.
 * @param[in] openMode The open mode. Type: <tt>int</tt>.
 */

typedef SmartFile<Bam, Output> BamFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct SmartFileContext<SmartFile<Bam, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    typedef BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<SmartFile<Bam, TDirection, TSpec> >
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
// Function _mapBamFormatToCompressionFormat()
// --------------------------------------------------------------------------

inline BgzfFile
_mapFileFormatToCompressionFormat(Bam)
{
    return BgzfFile();
}

// ----------------------------------------------------------------------------
// Function readRecord(); BamHeader
// ----------------------------------------------------------------------------

/*!
 * @fn BamFileIn#readRecord
 * @brief Read one @link BamAlignmentHeader @endlink or @link BamAlignmentRecord @endlink from a @link BamFileIn @endlink object.
 *
 * @signature int readRecord(header, bamFileIn);
 * @signature int readRecord(record, bamFileIn);
 *
 * @param[out]   header     The @link BamAlignmentHeader @endlink to read the header information into. Of type
 *                          @link BamAlignmentHeader @endlink.
 * @param[out]   record     The @link BamAlignmentRecord @endlink to read the next alignment record into. Of type
 *                          @link BamAlignmentRecord @endlink.
 * @param[in,out] bamFileIn The @link BamFileIn @endlink object to read from.
 */

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
// Function readRecord(); BamAlignmentRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TBuffer, typename TForwardIter>
inline void
_readBamRecord(TBuffer & /* rawRecord */, TForwardIter & /* iter */, TagSelector<> const & /* format */)
{
    SEQAN_FAIL("BamFileIn: File format not specified.");
}

template <typename TBuffer, typename TForwardIter, typename TTagList>
inline void
_readBamRecord(TBuffer & rawRecord, TForwardIter & iter, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        _readBamRecord(rawRecord, iter, TFormat());
    else
        _readBamRecord(rawRecord, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord(); BamAlignmentRecord
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

template <typename TRecords, typename TSpec, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And<IsSameType<typename Value<TRecords>::Type, BamAlignmentRecord>,
                                IsInteger<TSize> >, TSize)
readBatch(TRecords & records, SmartFile<Bam, Input, TSpec> & file, TSize maxRecords)
{
    String<CharString> & buffers = context(file).buffers;
    if ((TSize)length(buffers) < maxRecords)
    {
        resize(buffers, maxRecords, Exact());
        resize(records, maxRecords, Exact());
    }

    TSize numRecords = 0;
    for (; numRecords < maxRecords && !atEnd(file.iter); ++numRecords)
        _readBamRecord(buffers[numRecords], file.iter, file.format);

//    SEQAN_OMP_PRAGMA(parallel for)
    for (int i = 0; i < (int)numRecords; ++i)
    {
        CharIterator bufIter = begin(buffers[i]);
        readRecord(records[i], context(file), bufIter, file.format);
    }
    return numRecords;
}

// ----------------------------------------------------------------------------
// Function writeRecord(); BamHeader
// ----------------------------------------------------------------------------

/*!
 * @fn BamFileOut#writeRecord
 * @brief Write one @link BamAlignmentHeader @endlink or @link BamAlignmentRecord @endlink to a @link BamFileOut @endlink object.
 *
 * @signature int writeRecord(bamFileOut, header);
 * @signature int writeRecord(bamFileOut, record);
 *
 * @param[in,out] bamFileOut    The @link BamFileOut @endlink object to write to.
 * @param[in]     header        The @link BamAlignmentHeader @endlink to write out.
 * @param[in]     record        The @link BamAlignmentRecord @endlink to write out.
*/

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
writeRecord(SmartFile<Bam, Output, TSpec> & file, BamHeader const & header)
{
    write(file.iter, header, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); BamAlignmentRecord
// ----------------------------------------------------------------------------

// support for dynamically chosen file formats
template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
write(TTarget & /* target */,
      BamAlignmentRecord const & /* record */,
      BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
      TagSelector<> const & /* format */)
{
    SEQAN_FAIL("BamFileOut: File format not specified.");
}

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
write(TTarget & target,
      BamAlignmentRecord const & record,
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
writeRecord(SmartFile<Bam, Output, TSpec> & file, BamAlignmentRecord const & record)
{
    write(file.iter, record, context(file), file.format);
}

template <typename TSpec, typename TRecords>
inline SEQAN_FUNC_ENABLE_IF(IsSameType<typename Value<TRecords>::Type, BamAlignmentRecord>, void)
writeRecords(SmartFile<Bam, Output, TSpec> & file, TRecords const & records)
{
    String<CharString> & buffers = context(file).buffers;
    if (length(buffers) < length(records))
        resize(buffers, length(records));

    SEQAN_OMP_PRAGMA(parallel for)
    for (int i = 0; i < (int)length(records); ++i)
    {
        clear(buffers[i]);
        write(buffers[i], records[i], context(file), file.format);
    }
    for (int i = 0; i < (int)length(records); ++i)
        write(file.iter, buffers[i]);
}

}  // namespace seqan

#endif // SEQAN_BAM_IO_BAM_FILE_H_
