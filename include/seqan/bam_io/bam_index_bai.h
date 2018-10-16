// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// The index building algorithm is based on the samtools implementation which
// is under the MIT License:

/* The MIT License

   Copyright (c) 2008-2018 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BamIndex
// ----------------------------------------------------------------------------

/*!
 * @class BamIndex
 * @headerfile <seqan/bam_io.h>
 *
 * @brief Access to BAM indices.
 *
 * @signature template <typename TSpec>
 *            class BamIndex;
 *
 * This is an abstract class; don't use it itself but its specializations.
 *
 * @see BamFileIn
 */

template <typename TSpec>
class BamIndex;

// ----------------------------------------------------------------------------
// Tag Bai
// ----------------------------------------------------------------------------

struct Bai_;
typedef Tag<Bai_> Bai;

// ----------------------------------------------------------------------------
// Helper Class BaiBamIndexBinData_
// ----------------------------------------------------------------------------

// Store the information of a bin.

struct BaiBamIndexBinData_
{
    String<Pair<uint64_t, uint64_t> > chunkBegEnds;
};

// ----------------------------------------------------------------------------
// Spec BAI BamIndex
// ----------------------------------------------------------------------------

/*!
 * @class BaiBamIndex
 * @headerfile <seqan/bam_io.h>
 * @extends BamIndex
 * @brief Access to BAI (samtools-style).
 *
 * @signature template <>
 *            class BamIndex<Bai>;
 */

/*!
 * @fn BaiBamIndex::BamIndex
 * @brief Constructor.
 *
 * @signature BamIndex::BamIndex();
 *
 * @section Remarks
 *
 * Only the default constructor is provided.
 */

template <>
class BamIndex<Bai>
{
public:
    typedef std::map<uint32_t, BaiBamIndexBinData_> TBinIndex_;
    typedef String<uint64_t> TLinearIndex_;

    uint64_t _unalignedCount;

    // 1<<14 is the size of the minimum bin.
    static const int32_t BAM_LIDX_SHIFT = 14;

    String<TBinIndex_> _binIndices;
    String<TLinearIndex_> _linearIndices;

    BamIndex() : _unalignedCount(std::numeric_limits<uint64_t>::max())
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function viewRecords()
// ----------------------------------------------------------------------------

/*!
 * @fn BamFileIn#viewRecords
 * @brief Extract all reads overlapping a region from a BamFileIn mimicking samtools view.
 *
 * You provide a 1-based region <tt>[regionStart, regionEnd]</tt> on the reference
 * <tt>refID</tt> and the functions extracts all reads that overlap this region
 * with at least one base. The overlapping records are **appended** in the
 * <tt>resultContainer</tt>, which must be a container over @link BamAlignmentRecord @endlink.
 * This function extracts the same records as samtools view would do when the
 * same region is specified.
 * The result container remains empty if an error occurred or no reads were found.
 *
 * @signature bool viewRecords(resultContainer, bamFileIn, bamIndex, refID, regionStart, regionEnd, index);
 *
 * @param[in,out] resultContainer Container of @link BamAlignmentRecord @endlink 's to store the alignments.
 * @param[in,out] bamFileIn       The @link BamFileIn @endlink to extract reads from.
 * @param[in]     bamIndex        The @link BamIndex @endlink over the `bamFileIn`.
 * @param[in]     refID           The reference id to extract reads from (<tt>int32_t</tt>).
 * @param[in]     regionStart     The begin of the region to extract reads from (<tt>int32_t</tt>).
 * @param[in]     regionEnd       The end of the region to extract reads from (<tt>int32_t</tt>).
 *
 * @throw ParseError when reading the bamFile or index fails.
 * @throw std::logical_error when the input arguments are invalid.
 *
 * @section Remarks
 * This function fails if <tt>refID</tt>/<tt>regionStart</tt> are invalid.
 */

//! See SAM specifications https://samtools.github.io/hts-specs/SAMv1.pdf page 16.
//! This function expects 0-based coordinates and includes the positions beg and end.
static inline void
_baiReg2bins(std::vector<uint16_t> & list, uint32_t const beg, uint32_t const end)
{
    SEQAN_ASSERT_LEQ(beg, end);
    SEQAN_ASSERT_LT(end, 536870912u); // BAM format only allows reference of size < 2^29

    list.clear();

    unsigned k;
    list.push_back(0);
    for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list.push_back(k);
    for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list.push_back(k);
    for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list.push_back(k);
    for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list.push_back(k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push_back(k);
}

inline bool
_getMinFileOffset(uint64_t & minFileOffset,
                  BamIndex<Bai> const & bamIndex,
                  int32_t const refId,
                  uint32_t const position)
{
    uint32_t windowIdx = position >> 14; // Linear index consists of 16kb windows.

    if (windowIdx >= length(bamIndex._linearIndices[refId]))
        return false;

    minFileOffset = bamIndex._linearIndices[refId][windowIdx];

    return true;
}

template <typename TContainer, typename TSpec>
inline
SEQAN_FUNC_ENABLE_IF(And<IsSameType<typename Value<TContainer>::Type, BamAlignmentRecord>,
                     Not<IsSameType<typename Value<TContainer>::Type, TContainer>>>, void)
viewRecords(TContainer & resultContainer,
            FormattedFile<Bam, Input, TSpec> & bamFile,
            BamIndex<Bai> const & bamIndex,
            int32_t const refId,
            uint32_t regionStart,
            uint32_t regionEnd)
{
    typedef typename std::map<uint32_t, BaiBamIndexBinData_>::const_iterator TMapIter;
    typedef String<Pair<uint64_t, uint64_t>>                                 TChunk;
    typedef Iterator<TChunk const, Rooted>::Type                             TChunkIter;

    // Sanity checks
    // -------------------------------------------------------------------------
    if (!isEqual(format(bamFile), Bam()))
        SEQAN_THROW(std::logic_error("You attempt to use a BAM format specific functionality on a non-BAM format (function viewRecord)."));

    if (refId < 0 ||
        static_cast<decltype(length(contigNames(context(bamFile))))>(refId) > length(contigNames(context(bamFile))))
        SEQAN_THROW(std::logic_error("Invalid reference identifier 'refId' passed to function viewRecords."));

    if (static_cast<decltype(length(bamIndex._binIndices))>(refId) >= length(bamIndex._binIndices))
        SEQAN_THROW(std::logic_error("BAM index bin directory does not match the given reference identifier 'refId'. "
                   "Maybe your BAM index is corrupted or refers to a different BAM file."));

    if (regionEnd > static_cast<uint32_t>(contigLengths(context(bamFile))[refId]))
        SEQAN_THROW(std::logic_error("The 'regionEnd' parameter point to a position greater than the length of the reference."));

    if (regionStart > regionEnd)
        SEQAN_THROW(std::logic_error("Invalid region specified. Parameter regionStart was greater than regionEnd."));

    if (regionStart == 0 || regionEnd == 0)
        SEQAN_THROW(std::logic_error("Invalid region specified. Region is expected to be 1-based but 0 was passed as an argument."));

    // Adjust region
    // -------------------------------------------------------------------------
    // We want to mimic samtools view. Since SAM/BAM files are 1-based we expect
    // the given region to be also 1-based, to be consistent with samtools.
    // Since SeqAn works 0-based we need to subtract 1.
    // (Underflow won't happen because of prior sanity checks)
    --regionStart;
    --regionEnd;

    // Retrieve the candidate bin identifiers for [regionStart, regionEnd).
    // -------------------------------------------------------------------------
    std::vector<uint16_t> candidateBins;
    _baiReg2bins(candidateBins, regionStart, regionEnd); // 0-based, closed interval

    // Retrieve the smallest required offset from the linear index.
    // -------------------------------------------------------------------------
    uint64_t linearMinOffset;
    if (!_getMinFileOffset(linearMinOffset, bamIndex, refId, regionStart))
        linearMinOffset = 0; // if linear index query was unsuccessful, value 0 never limits chunk seeking

    // Iterate over chunks in candidate bins and extract records
    // -------------------------------------------------------------------------
    BamAlignmentRecord record;
    for (uint16_t bin : candidateBins)
    {
        TMapIter mIt = bamIndex._binIndices[refId].find(bin);

        if (mIt == bamIndex._binIndices[refId].end()) // empty bins are not stored
            continue;

        for (TChunkIter chunkIt = begin(mIt->second.chunkBegEnds, Rooted()); !atEnd(chunkIt); goNext(chunkIt))
        {
            if (chunkIt->i2 < linearMinOffset) // i2 is the end-file-offset of the chunk
                continue; // use linearMinOffset to discard unnecessary chunks

            setPosition(bamFile, chunkIt->i1); // jump to candidate chunk

            while(static_cast<uint64_t>(position(bamFile)) != chunkIt->i2) // this includes end of file
            {
                readRecord(record, bamFile);

                if (record.beginPos == -1 /*invalid pos*/ || static_cast<uint32_t>(record.beginPos) > regionEnd)
                    break; // chunk of consecutive reads is already out of bounds

                // max(1, length) in if clause because if read is unmapped, the endPos calculation is wrong.
                // (Could be also done by checking if (unmapped) but might be less efficient than max.
                if (static_cast<uint32_t>(record.beginPos) <= regionEnd &&
                    (record.beginPos + std::max(1u, getAlignmentLengthInRef(record)) - 1) >= regionStart)
                    appendValue(resultContainer, record, Generous());
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function jumpToRegion()
// ----------------------------------------------------------------------------

/*!
 * @fn BamFileIn#jumpToRegion
 * @brief Seek in BamFileIn using an index.
 *
 * You provide a 1-based region <tt>[regionStart, regionEnd]</tt> on the reference
 * <tt>refID</tt> that you want to jump to and the function jumps to the first
 * alignment (bam record) whose **start position** lies inside this region, if any.
 *
 * Note: The current implementation only considers the read start positions.
 *       Therefore, we do guarantee that reads overlapping the region are included.
 *       To account for this limitation you may want to choose the start of your
 *       your region with an appropriate offset (e.g. start - length_of_read) or
 *       use viewRecords().
 *
 * @signature bool jumpToRegion(bamFileIn, hasAlignments, refID, regionStart, regionEnd, index);
 *
 * @param[in,out] bamFileIn     The @link BamFileIn @endlink to jump with.
 * @param[out]    hasAlignments A <tt>bool</tt> that is set true if the region <tt>[regionStart, regionEnd]</tt> has any
 *                              at least one overlapping alignment (bam record).
 * @param[in]     refID         The reference id to jump to (<tt>int32_t</tt>).
 * @param[in]     regionStart   The begin of the region to jump to (<tt>int32_t</tt>).
 * @param[in]     regionEnd     The end of the region to jump to (<tt>int32_t</tt>).
 * @param[in]     index         The @link BamIndex @endlink to use for the jumping.
 *
 * @return bool `false` if an error occurred while seeking, otherwise `true`.
 *
 * @throw ParseError when reading the bamFile or index fails.
 * @throw std::logical_error when the input arguments are invalid.
 *
 * @section Remarks
 *
 * This function fails if <tt>refID</tt>/<tt>regionStart</tt> are invalid.
 */

template <typename TSpec>
inline bool
jumpToRegion(FormattedFile<Bam, Input, TSpec> & bamFile,
             bool & hasAlignments,
             int32_t const refId,
             uint32_t regionStart,
             uint32_t regionEnd,
             BamIndex<Bai> const & index)
{
    hasAlignments = false;

    // Sanity checks
    // -------------------------------------------------------------------------
    // Sanity checks
    // -------------------------------------------------------------------------
    if (!isEqual(format(bamFile), Bam()))
        SEQAN_THROW(std::logic_error("You attempt to use a BAM format specific functionality on a non-BAM format (function viewRecord)."));

    if (refId < 0 ||
        static_cast<decltype(length(contigNames(context(bamFile))))>(refId) > length(contigNames(context(bamFile))))
        SEQAN_THROW(std::logic_error("Invalid reference identifier 'refId' passed to function viewRecords."));

    if (static_cast<decltype(length(index._binIndices))>(refId) >= length(index._binIndices))
        SEQAN_THROW(std::logic_error("BAM index bin directory does not match the given reference identifier 'refId'. "
                   "Maybe your BAM index is corrupted or refers to a different BAM file."));

    if (regionEnd > static_cast<uint32_t>(contigLengths(context(bamFile))[refId]))
        SEQAN_THROW(std::logic_error("The 'regionEnd' parameter point to a position greater than the length of the reference."));
        // correct end positions that exceed the reference length for release build.

    if (regionStart > regionEnd)
        SEQAN_THROW(std::logic_error("Invalid region specified. Parameter regionStart was greater than regionEnd."));

    if (regionStart == 0 || regionEnd == 0)
        SEQAN_THROW(std::logic_error("Invalid region specified. Region is expected to be 1-based but 0 was passed as an argument."));

    // Adjust region
    // -------------------------------------------------------------------------
    // We want to be consistent with viewRecords and thus expect 1-based coordinates.
    // Since SeqAn works 0-based we need to subtract 1.
    --regionStart;
    --regionEnd;

    // Retrieve the smallest required offset from the linear index.
    // -------------------------------------------------------------------------
    uint64_t linearMinOffset;
    if (!_getMinFileOffset(linearMinOffset, index, refId, regionStart)) // 0-based coordinates.
        return false;

    // Now only the minimum offset was found but it is not ensured that the
    // alignment at the start of the offset lies inside the region.
    // Therefore we must advance further or until the region ends.
    setPosition(bamFile, linearMinOffset);
    BamAlignmentRecord record;

    while (!atEnd(bamFile))
    {
        uint64_t seek_pos = position(bamFile); // store current pos to reset bamFile if necessary

        readRecord(record, bamFile);

        if (record.beginPos == -1 /*invalid pos*/ || static_cast<uint32_t>(record.beginPos) > regionEnd)
            break; // no records can be found anymore in coordinate sorted file

        if (static_cast<uint32_t>(record.beginPos) >= regionStart) // we do not support overlap but the start position of a record
        {
            hasAlignments = true;
            setPosition(bamFile, seek_pos);
            break;
        }
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function jumpToOrphans()
// ----------------------------------------------------------------------------

/*!
 * @fn BamFileIn#jumpToOrphans
 * @brief Seek to orphans block in BamFileIn using an index.
 *
 * @signature bool jumpToOrphans(bamFileIn, hasAlignments, index);
 *
 * @param[in,out] bamFileIn      The @link BamFileIn @endlink object to jump with.
 * @param[out]    hasAlignments  A <tt>bool</tt> that is set to true if there are any orphans.
 * @param[in]     index          The @link BamIndex @endlink to use for jumping.
 */

template <typename TSpec>
bool jumpToOrphans(FormattedFile<Bam, Input, TSpec> & bamFile,
                   bool & hasAlignments,
                   BamIndex<Bai> const & index)
{
    if (!isEqual(format(bamFile), Bam()))
        return false;

    hasAlignments = false;

    // Search linear indices for the largest entry of all references.
    uint64_t aliOffset = std::numeric_limits<uint64_t>::max();
    for (int i = length(index._linearIndices) - 1; i >= 0; --i)
        if (!empty(index._linearIndices[i]))
        {
            aliOffset = back(index._linearIndices[i]);
            break;
        }
    if (aliOffset == std::numeric_limits<uint64_t>::max())
        return false;  // No offset found.

    // Get index of the first orphan alignment by seeking from linear index bucket.
    BamAlignmentRecord record;
    uint64_t offset = std::numeric_limits<uint64_t>::max();
    uint64_t result = 0;
    if (!setPosition(bamFile, aliOffset))
        return false;  // Error while seeking.
    while (!atEnd(bamFile))
    {
        result = position(bamFile);
        readRecord(record, bamFile);
        if (record.rID == -1)
        {
            // Found alignment.
            hasAlignments = true;
            offset = result;
            break;
        }
    }

    // Jump back to the first alignment.
    if (offset != std::numeric_limits<uint64_t>::max())
    {
        if (!setPosition(bamFile, offset))
            return false;  // Error while seeking.
    }

    // Finding no orphan alignment is not an error, hasAilgnments is false then.
    return true;
}

// ----------------------------------------------------------------------------
// Function getUnalignedCount()
// ----------------------------------------------------------------------------

/*!
 * @fn BamIndex#getUnalignedCount
 * @brief Query index for number of unaligned reads.
 *
 * @signature uint64_t getUnalignedCount(index);
 *
 * @param[in] index     Index to query.
 * @return    uint64_t  The number of unaligned reads.
 */

inline uint64_t
getUnalignedCount(BamIndex<Bai> const & index)
{
    return index._unalignedCount;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn BamIndex#open
 * @brief Load a BAM index from a given file name.
 * @signature bool open(index, filename);

 * @param[in,out] index    Target data structure.
 * @param[in]     filename Path to file to load. Types: char const *
 *
 * @return        bool     Returns <tt>true</tt> on success, false otherwise.
 */

inline bool
open(BamIndex<Bai> & index, char const * filename)
{
    std::fstream fin(filename, std::ios::binary | std::ios::in);
    if (!fin.good())
        return false;  // Could not open file.

    // Read magic number.
    CharString buffer;
    resize(buffer, 4);
    fin.read(&buffer[0], 4);
    if (!fin.good())
        return false;
    if (buffer != "BAI\1")
        return false;  // Magic number is wrong.

    int32_t nRef = 0;
    fin.read(reinterpret_cast<char *>(&nRef), 4);
    if (!fin.good())
        return false;
    enforceLittleEndian(nRef);

    resize(index._linearIndices, nRef);
    resize(index._binIndices, nRef);

    for (int i = 0; i < nRef; ++i)  // For each reference.
    {
        // Read bin index.
        int32_t nBin = 0;
        fin.read(reinterpret_cast<char *>(&nBin), 4);
        if (!fin.good())
            return false;
        enforceLittleEndian(nBin);
        index._binIndices[i].clear();
        BaiBamIndexBinData_ data;
        for (int j = 0; j < nBin; ++j)  // For each bin.
        {
            clear(data.chunkBegEnds);

            uint32_t bin = 0;
            fin.read(reinterpret_cast<char *>(&bin), 4);
            if (!fin.good())
                return false;
            enforceLittleEndian(bin);

            int32_t nChunk = 0;
            fin.read(reinterpret_cast<char *>(&nChunk), 4);
            if (!fin.good())
                return false;
            enforceLittleEndian(nChunk);
            reserve(data.chunkBegEnds, nChunk);
            for (int k = 0; k < nChunk; ++k)  // For each chunk;
            {
                uint64_t chunkBeg = 0;
                uint64_t chunkEnd = 0;
                fin.read(reinterpret_cast<char *>(&chunkBeg), 8);
                fin.read(reinterpret_cast<char *>(&chunkEnd), 8);
                if (!fin.good())
                    return false;
                enforceLittleEndian(chunkBeg);
                enforceLittleEndian(chunkEnd);
                appendValue(data.chunkBegEnds, Pair<uint64_t>(chunkBeg, chunkEnd));
            }

            // Copy bin data into index.
            index._binIndices[i][bin] = data;
        }

        // Read linear index.
        int32_t nIntv = 0;
        fin.read(reinterpret_cast<char *>(&nIntv), 4);
        if (!fin.good())
            return false;
        enforceLittleEndian(nIntv);
        clear(index._linearIndices[i]);
        reserve(index._linearIndices[i], nIntv);
        for (int j = 0; j < nIntv; ++j)
        {
            uint64_t ioffset = 0;
            fin.read(reinterpret_cast<char *>(&ioffset), 8);
            if (!fin.good())
                return false;
            enforceLittleEndian(ioffset);
            appendValue(index._linearIndices[i], ioffset);
        }
    }

    if (!fin.good())
        return false;

    // Read (optional) number of alignments without coordinate.
    uint64_t nNoCoord = 0;
    fin.read(reinterpret_cast<char *>(&nNoCoord), 8);
    if (!fin.good())
    {
        fin.clear();
        nNoCoord = 0;
    }
    enforceLittleEndian(nNoCoord);
    index._unalignedCount = nNoCoord;

    return true;
}

// TODO(holtgrew): This is only here because of the read() function with TSequence in old file.h.

inline bool
open(BamIndex<Bai> & index, char * filename)
{
    return open(index, static_cast<char const *>(filename));
}

// ---------------------------------------------------------------------------
// Function save()
// ---------------------------------------------------------------------------

/*!
 * @fn BamIndex#save
 * @brief Save a BamIndex object.
 *
 * @signature bool save(baiIndex, baiFileName);
 *
 * @param[in] baiIndex    The BamIndex to write out.
 * @param[in] baiFileName The name of the BAI file to write to.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> otherwise.
 */

inline bool save(BamIndex<Bai> const & index, char const * baiFilename)
{
    // Open output file.
    std::ofstream out(baiFilename, std::ios::binary | std::ios::out);

    SEQAN_ASSERT_EQ(length(index._binIndices), length(index._linearIndices));

    // Write header.
    out.write("BAI\1", 4);
    int32_t numRefSeqs = length(index._binIndices);
    enforceLittleEndian(numRefSeqs);
    out.write(reinterpret_cast<char *>(&numRefSeqs), 4);

    // Write out indices.
    typedef BamIndex<Bai> const                TBamIndex;
    typedef TBamIndex::TBinIndex_ const        TBinIndex;
    typedef TBinIndex::const_iterator          TBinIndexIter;
    typedef TBamIndex::TLinearIndex_           TLinearIndex;
    for (unsigned i = 0; i < length(index._binIndices); ++i)
    {
        TBinIndex const & binIndex = index._binIndices[i];
        TLinearIndex const & linearIndex = index._linearIndices[i];

        // Write out binning index.
        int32_t numBins = binIndex.size();
        enforceLittleEndian(numBins);
        out.write(reinterpret_cast<char *>(&numBins), 4);
        for (TBinIndexIter itB = binIndex.begin(), itBEnd = binIndex.end(); itB != itBEnd; ++itB)
        {
            // Write out bin id.
            uint32_t tmp = itB->first;
            enforceLittleEndian(tmp);
            out.write(reinterpret_cast<char const *>(&tmp), 4);
            // Write out number of chunks.
            uint32_t numChunks = length(itB->second.chunkBegEnds);
            enforceLittleEndian(numChunks);
            out.write(reinterpret_cast<char *>(&numChunks), 4);
            // Write out all chunks.
            typedef Iterator<String<Pair<uint64_t> > const, Rooted>::Type TChunkIter;
            for (TChunkIter itC = begin(itB->second.chunkBegEnds); !atEnd(itC); goNext(itC))
            {
                uint64_t tmp1 = itC->i1;
                enforceLittleEndian(tmp1);
                uint64_t tmp2 = itC->i2;
                enforceLittleEndian(tmp2);
                out.write(reinterpret_cast<char const *>(&tmp1), 8);
                out.write(reinterpret_cast<char const *>(&tmp2), 8);
            }
        }

        // Write out linear index.
        int32_t numIntervals = length(linearIndex);
        enforceLittleEndian(numIntervals);
        out.write(reinterpret_cast<char *>(&numIntervals), 4);
        typedef Iterator<String<uint64_t> const, Rooted>::Type TLinearIndexIter;
        for (TLinearIndexIter it = begin(linearIndex, Rooted()); !atEnd(it); goNext(it))
        {
            uint64_t tmp = *it;
            enforceLittleEndian(tmp);
            out.write(reinterpret_cast<char const *>(&tmp), 8);
        }
    }

    // Write the number of unaligned reads if set.
    //std::cerr << "UNALIGNED\t" << index._unalignedCount << std::endl;
    if (index._unalignedCount != std::numeric_limits<uint64_t>::max())
    {
        uint64_t tmp = index._unalignedCount;
        enforceLittleEndian(tmp);
        out.write(reinterpret_cast<char const *>(&tmp), 8);
    }

    return out.good();  // false on error, true on success.
}


inline void _baiAddAlignmentChunkToBin(BamIndex<Bai> & index,
                                       uint32_t currBin,
                                       uint32_t currOffset,
                                       uint64_t prevOffset)
{
    // If this is not the first reference sequence then add previous reference data.
    Pair<uint64_t> newChunk(currOffset, prevOffset);

    // If no interest exists yet for this bin, create one and store alignment chunk.
    BamIndex<Bai>::TBinIndex_::iterator binIter = back(index._binIndices).find(currBin);
    if (binIter == back(index._binIndices).end())
    {
        BaiBamIndexBinData_ binData;
        appendValue(binData.chunkBegEnds, newChunk);
        back(index._binIndices).insert(std::make_pair(currBin, binData));
    }
    else
    {
        // Otherwise, just append alignment chunk.
        appendValue(binIter->second.chunkBegEnds, newChunk);
    }
}

// ---------------------------------------------------------------------------
// Function build()
// ---------------------------------------------------------------------------
// TODO(dadi): uncomment when BamIndex.build index is fixed. DOX commented out
/*
 * @fn BamIndex#build
 * @brief Create a BamIndex from BAM file.
 *
 * @signature bool build(baiIndex, bamFileName);
 *
 * @param[out] baiIndex    The BamIndex to build into.
 * @param[in]  bamFileName Path to the BAM file to build an index for.  Type: <tt>char const *</tt>.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> otherwise.
 */
inline bool build(BamIndex<Bai> & index, char const * bamFilename)
{
    // SEQAN_FAIL("This does not work yet!");

    index._unalignedCount = 0;
    clear(index._binIndices);
    clear(index._linearIndices);

    // Open BAM file for reading.
    BamFileIn bamFile;
    if (!open(bamFile, bamFilename))
        return false;  // Could not open BAM file.

    // Read BAM header.
    BamHeader header;
    readHeader(header, bamFile);

    uint32_t numRefSeqs = length(contigNames(context(bamFile)));

    // Scan over BAM file and create index.
    BamAlignmentRecord record;
    uint32_t currBin    = std::numeric_limits<uint32_t>::max();
    uint32_t prevBin    = std::numeric_limits<uint32_t>::max();
    int32_t currRefId   = BamAlignmentRecord::INVALID_REFID;
    int32_t prevRefId   = BamAlignmentRecord::INVALID_REFID;
    uint64_t currOffset = position(bamFile);
    uint64_t prevOffset = currOffset;
    int32_t prevPos     = std::numeric_limits<int32_t>::min();

    while (!atEnd(bamFile))
    {
        // Load next record.
        readRecord(record, bamFile);

        // Check ordering.
        if (prevRefId == record.rID && prevPos > record.beginPos)
            return false;

        // The reference sequence changed, close bins for previous reference.
        if (prevRefId != record.rID)
        {
            if (prevRefId != BamAlignmentRecord::INVALID_REFID)
            {
                _baiAddAlignmentChunkToBin(index, currBin, currOffset, prevOffset);

                // Add an index for all empty references between prevRefId (excluded) and record.rID (included).
                for (int i = prevRefId + 1; i < record.rID; ++i)
                {
                    BamIndex<Bai>::TBinIndex_ binIndex;
                    appendValue(index._binIndices, binIndex);
                    BamIndex<Bai>::TLinearIndex_ linearIndex;
                    appendValue(index._linearIndices, linearIndex);
                }

                // Update bin book keeping.
                currOffset = prevOffset;
                currBin    = record.bin;
                prevBin    = record.bin;
                currRefId  = record.rID;
            }
            else
            {
                // Otherwise, this is the first pass.  Create an index for all empty references up to and including
                // current refId.
                for (int i = 0; i < record.rID; ++i)
                {
                    BamIndex<Bai>::TBinIndex_ binIndex;
                    appendValue(index._binIndices, binIndex);
                    BamIndex<Bai>::TLinearIndex_ linearIndex;
                    appendValue(index._linearIndices, linearIndex);
                }
            }

            // Update reference book keeping.
            prevRefId = record.rID;
            prevBin = std::numeric_limits<int32_t>::min();
        }

        // If the alignment's reference id is valid and its bin is not a leaf.
        if (record.rID >= 0 && record.bin < 4681)
        {
            int32_t beginOffset = record.beginPos >> BamIndex<Bai>::BAM_LIDX_SHIFT;
            int32_t endPos      = getAlignmentLengthInRef(record);
            int32_t endOffset   = (endPos - 1) >> BamIndex<Bai>::BAM_LIDX_SHIFT;

            // Resize linear index if necessary.
            unsigned oldSize = length(index._linearIndices);
            unsigned newSize = endOffset + 1;
            if (oldSize < newSize)
                resize(back(index._linearIndices), newSize, 0);

            // Store offset.
            for (int i = beginOffset + 1; i <= endOffset; ++i)
                if (back(index._linearIndices)[i] == 0u)
                    back(index._linearIndices)[i] = prevOffset;
        }

        // Handle the case if we changed to a new BAI bin.
        if (record.bin != prevBin)
        {
            // If not first bin of reference, save previous bin data.
            if (currBin != std::numeric_limits<uint32_t>::max())
                _baiAddAlignmentChunkToBin(index, currBin, currOffset, prevOffset);

            // Update markers.
            currOffset = prevOffset;
            currBin    = record.bin;
            prevBin    = record.bin;
            currRefId  = record.rID;

            // If the reference id is invalid then break out.
            if (currRefId < 0)
                break;
        }

        // Make sure that the current file pointer is beyond prevOffset.
        if (position(bamFile) <= static_cast<int64_t>(prevOffset))
            return false;  // Calculating offsets failed.

        // Update prevOffset and prevPos.
        prevOffset = position(bamFile);
        prevPos    = record.beginPos;
    }

    // Count remaining unaligned records.
    while (!atEnd(bamFile))
    {
        SEQAN_ASSERT_GT(index._unalignedCount, 0u);

        readRecord(record, bamFile);
        if (record.rID >= 0)
            return false;  // Could not read record.

        index._unalignedCount += 1;
    }

    // After loading all alignments, if any data was read, perform checks.
    if (currRefId >= 0)
    {
        // Store last alignment chunk to its bin and then write last reference entry with data.
        _baiAddAlignmentChunkToBin(index, currBin, currOffset, prevOffset);

        // Finally, write any empty references remaining at end of file.
        SEQAN_ASSERT_GEQ(numRefSeqs, length(index._binIndices));
        BamIndex<Bai>::TBinIndex_ binIndex;
        resize(index._binIndices, numRefSeqs, binIndex);  // TODO(holtgrew): binIndex is unnecessary if resize used T() as default value as vector.resize() does.
    }

    // Merge small bins if possible.
    // SEQAN_FAIL("TODO: Merge bins!");
    return true;
}


}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
