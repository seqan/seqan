// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_TABIX_IO_TABIX_INDEX_TBI_H_
#define INCLUDE_SEQAN_TABIX_IO_TABIX_INDEX_TBI_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Class TabixIndexBinData_
// ----------------------------------------------------------------------------

// Store the information of a bin.

struct TabixIndexBinData_
{
    String<Pair<__uint64, __uint64> > chunkBegEnds;
};

// ----------------------------------------------------------------------------
// Class TabixIndex
// ----------------------------------------------------------------------------

/*!
 * @class TabixIndex
 * @headerfile <seqan/tabix_io.h>
 * @brief Access to Tabix indexed files.
 *
 * @signature class TabixIndex;
 */

/*!
 * @fn TabixIndex::TabixIndex
 * @brief Constructor.
 *
 * @signature TabixIndex::TabixIndex();
 *
 * @section Remarks
 *
 * Only the default constructor is provided.
 */

class TabixIndex
{
public:
    typedef std::map<__uint32, TabixIndexBinData_> TBinIndex_;
    typedef String<__uint64> TLinearIndex_;
    typedef StringSet<CharString, Owner<ConcatDirect<> > > TNameStore;

    __int32 format;             // Format (0: generic; 1: SAM; 2: VCF)
    __int32 colSeq;             // Column for the sequence name
    __int32 colBeg;             // Column for the start of a region
    __int32 colEnd;             // Column for the end of a region
    __int32 meta;               // Leading character for comment lines
    __int32 skip;               // # lines to skip at the beginning
    __uint64 unalignedCount;    // # unmapped reads without coordinates set

    // 1<<14 is the size of the minimum bin.
    static const __int32 BAM_LIDX_SHIFT = 14;

    String<TBinIndex_>          _binIndices;
    String<TLinearIndex_>       _linearIndices;
    TNameStore                  _nameStore;
    NameStoreCache<TNameStore>  _nameStoreCache;

    TabixIndex() :
        format(0),
        colSeq(1),
        colBeg(2),
        colEnd(3),
        meta('#'),
        skip(0),
        unalignedCount(maxValue<__uint64>())
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _tbiReg2bins()
// ----------------------------------------------------------------------------

static inline void
_tbiReg2bins(String<__uint16> & list, __uint32 beg, __uint32 end)
{
    unsigned k;
    if (beg >= end) return;
    if (end >= 1u<<29) end = 1u<<29;
    --end;
    appendValue(list, 0);
    for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) appendValue(list, k);
    for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) appendValue(list, k);
    for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) appendValue(list, k);
    for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) appendValue(list, k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) appendValue(list, k);
}

// ----------------------------------------------------------------------------
// Function jumpToRegion()
// ----------------------------------------------------------------------------

/*!
 * @fn TabixIndex#jumpToRegion
 * @brief Seek in VcfFileIn or BedFileIn using an index.
 *
 * You provide a region <tt>[pos, posEnd)</tt> on the reference <tt>refID</tt> that you want to jump to and the function
 * jumps to the first alignment in this region, if any.
 *
 * @signature bool jumpToRegion(bamFileIn, hasAlignments, refID, pos, posEnd, index);
 *
 * @param[in,out] bamFileIn     The @link BamFileIn @endlink to jump with.
 * @param[out]    hasAlignments A <tt>bool</tt> that is set true if the region <tt>[pos, posEnd)</tt> has any
 *                              alignments.
 * @param[in]     refID         The reference id to jump to (<tt>__int32</tt>).
 * @param[in]     pos           The begin of the region to jump to (<tt>__int32</tt>).
 * @param[in]     posEnd        The end of the region to jump to (<tt>__int32</tt>).
 * @param[in]     index         The @link TabixIndex @endlink to use for the jumping.
 *
 * @return bool true if seeking was successful, false if not.
 *
 * @section Remarks
 *
 * This function fails if <tt>refID</tt>/<tt>pos</tt> are invalid.
 */

template <typename TFormat, typename TSpec, typename TName>
inline bool
jumpToRegion(FormattedFile<TFormat, Input, TSpec> & file,
             TName seqName,
             __int32 pos,
             __int32 posEnd,
             TabixIndex const & index)
{
    if (!isEqual(format(bamFile), Bam()))
        return false;

    hasAlignments = false;
    if (refId < 0)
        return false;  // Cannot seek to invalid reference.
    if (static_cast<unsigned>(refId) >= length(index._binIndices))
        return false;  // Cannot seek to invalid reference.

    // ------------------------------------------------------------------------
    // Compute offset in BGZF file.
    // ------------------------------------------------------------------------
    __uint64 offset = MaxValue<__uint64>::VALUE;

    // Retrieve the candidate bin identifiers for [pos, posEnd).
    String<__uint16> candidateBins;
    _tbiReg2bins(candidateBins, pos, posEnd);

    // Retrieve the smallest required offset from the linear index.
    unsigned windowIdx = pos >> 14;  // Linear index consists of 16kb windows.
    __uint64 linearMinOffset = 0;
    if (windowIdx >= length(index._linearIndices[refId]))
    {
        // TODO(holtgrew): Can we simply always take case 1?

        // This is the case were we want to jump in a non-existing window.
        //
        // If there are no linear indices for this reference then we use the linear min offset of the next
        // reference that has an linear index.
        if (empty(index._linearIndices[refId]))
        {
            for (unsigned i = refId; i < length(index._linearIndices); ++i)
            {
                if (!empty(index._linearIndices[i]))
                {
                    linearMinOffset = front(index._linearIndices[i]);
                    if (linearMinOffset != 0u)
                        break;
                    for (unsigned j = 1; j < length(index._linearIndices[i]); ++j)
                    {
                        if (index._linearIndices[i][j] > linearMinOffset)
                        {
                            linearMinOffset = index._linearIndices[i][j];
                            break;
                        }
                    }
                    if (linearMinOffset != 0u)
                        break;
                }
            }
        }
        else
        {
            linearMinOffset = back(index._linearIndices[refId]);
        }
    }
    else
    {
        linearMinOffset = index._linearIndices[refId][windowIdx];
    }

    // Combine candidate bins and smallest required offset from linear index into candidate offset.
    typedef std::set<__uint64> TOffsetCandidates;
    TOffsetCandidates offsetCandidates;
    typedef typename Iterator<String<__uint16>, Rooted>::Type TCandidateIter;
    for (TCandidateIter it = begin(candidateBins, Rooted()); !atEnd(it); goNext(it))
    {
        typedef typename std::map<__uint32, TabixIndexBinData_>::const_iterator TMapIter;
        TMapIter mIt = index._binIndices[refId].find(*it);
        if (mIt == index._binIndices[refId].end())
            continue;  // Candidate is not in index!

        typedef typename Iterator<String<Pair<__uint64, __uint64> > const, Rooted>::Type TBegEndIter;
        for (TBegEndIter it2 = begin(mIt->second.chunkBegEnds, Rooted()); !atEnd(it2); goNext(it2))
            if (it2->i2 >= linearMinOffset)
                offsetCandidates.insert(it2->i1);
    }

    // Search through candidate offsets, find rightmost possible.
    //
    // Note that it is not necessarily the first.
    //
    // TODO(holtgrew): Can this be optimized similar to how bamtools does it?
    typedef typename TOffsetCandidates::const_iterator TOffsetCandidateIter;
    BamAlignmentRecord record;
    for (TOffsetCandidateIter candIt = offsetCandidates.begin(); candIt != offsetCandidates.end(); ++candIt)
    {
        setPosition(bamFile, *candIt);

        readRecord(record, bamFile);

        // std::cerr << "record.beginPos == " << record.beginPos << "\n";
        // __int32 endPos = record.beginPos + getAlignmentLengthInRef(record);
        if (record.rID != refId)
            continue;  // Wrong contig.
        if (!hasAlignments || record.beginPos <= pos)
        {
            // Found a valid alignment.
            hasAlignments = true;
            offset = *candIt;
        }

        if (record.beginPos >= posEnd)
            break;  // Cannot find overlapping any more.
    }

    if (offset != MaxValue<__uint64>::VALUE)
        setPosition(bamFile, offset);

    // Finding no overlapping alignment is not an error, hasAlignments is false.
    return true;
}

// ----------------------------------------------------------------------------
// Function getUnalignedCount()
// ----------------------------------------------------------------------------

/*!
 * @fn TabixIndex#getUnalignedCount
 * @brief Query index for number of unaligned reads.
 *
 * @signature __uint64 getUnalignedCount(index);
 *
 * @param[in] index     Index to query.
 * @return    __uint64  The number of unaligned reads.
 */

inline __uint64
getUnalignedCount(TabixIndex const & index)
{
    return index.unalignedCount;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn TabixIndex#open
 * @brief Load a BAM index from a given file name.
 * @signature bool open(index, filename);

 * @param[in,out] index    Target data structure.
 * @param[in]     filename Path to file to load. Types: char const *
 *
 * @return        bool     Returns <tt>true</tt> on success, false otherwise.
 */

inline bool
open(TabixIndex & index, char const * filename)
{
    typedef VirtualStream<char, Input> TInStream;
    
    TInStream tbi;
    if (!open(tbi, filename))
        return false;  // Could not open file.

    DirectionIterator<TInStream, Input>::Type iter = directionIterator(tbi, Input());

    // Read magic header.
    String<char, Array<4> > magic;
    read(magic, iter, 4);
    if (magic != "TBI\1")
        SEQAN_THROW(ParseError("Not in TBI format."));
    
    // Read parameters.
    __int32 nRef = 0;
    readRawPod(nRef, iter);
    readRawPod(index.format, iter);
    readRawPod(index.colSeq, iter);
    readRawPod(index.colBeg, iter);
    readRawPod(index.colEnd, iter);
    readRawPod(index.meta, iter);
    readRawPod(index.skip, iter);

    // Read concatenated names.
    __int32 lNm = 0;
    CharString tmp;
    readRawPod(lNm, iter);
    read(tmp, iter, lNm);

    // Split concatenated names at \0's.
    clear(index._nameStore);
    strSplit(index._nameStore, tmp, EqualsChar<'\0'>(), true, nRef - 1);

    clear(index._linearIndices);
    clear(index._binIndices);
    resize(index._linearIndices, nRef);
    resize(index._binIndices, nRef);

    TabixIndexBinData_ data;
    for (int i = 0; i < nRef; ++i)  // For each reference.
    {
        // Read bin index.
        __int32 nBin = 0;
        readRawPod(nBin, iter);

        for (int j = 0; j < nBin; ++j)  // For each bin.
        {
            __uint32 bin = 0;
            __int32 nChunk = 0;
            readRawPod(bin, iter);
            readRawPod(nChunk, iter);

            resize(data, nChunk);
            for (int k = 0; k < nChunk; ++k)  // For each chunk;
            {
                readRawPod(data.chunkBegEnds[k].i1, iter);
                readRawPod(data.chunkBegEnds[k].i2, iter);
            }

            // Copy bin data into index.
            index._binIndices[i][bin] = data;
        }

        // Read linear index.
        __int32 nIntv = 0;
        readRawPod(nIntv, iter);

        resize(index._linearIndices[i], nIntv);
        for (int j = 0; j < nIntv; ++j)
            readRawPod(index._linearIndices[i], iter);
    }


    // Read (optional) number of alignments without coordinate.
    if (!atEnd(iter))
        readRawPod(index._unalignedCount, iter);
    else
        index._unalignedCount = maxValue<__uint64>();

    return true;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_TABIX_IO_TABIX_INDEX_TBI_H_
