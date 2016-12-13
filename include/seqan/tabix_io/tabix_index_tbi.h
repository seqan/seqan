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
// (Read-only) Tabix index support.
//
// A Tabix index (Heng Li) allows one to randomly seek in a tab-seperated genome
// related file, e.g. VCF, GFF, SAM, BED, etc. The corresponding file only
// needs to be sorted by chromosomal position in advance and optionally
// compressed with 'bgzip'. The resulting file must be indexed with 'tabix'.
//
// TODOs:
//  - clean jumpToRegion(), I simply adapted the one from bam_index.h
//  - implement Tabix index creation (and BAM index creation as well)
// ==========================================================================

#ifndef INCLUDE_SEQAN_TABIX_IO_TABIX_INDEX_TBI_H_
#define INCLUDE_SEQAN_TABIX_IO_TABIX_INDEX_TBI_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

class TabixIndex;
bool open(TabixIndex & index, char const * filename);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Class TabixIndexBinData_
// ----------------------------------------------------------------------------

// Store the information of a bin.

struct TabixIndexBinData_
{
    String<Pair<uint64_t, uint64_t> > chunkBegEnds;
};

// ----------------------------------------------------------------------------
// Helper Class TabixRecord_
// ----------------------------------------------------------------------------

struct TabixRecord_
{
    CharString refName;
    int32_t posBeg;
    int32_t posEnd;
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
    typedef std::map<uint32_t, TabixIndexBinData_> TBinIndex_;
    typedef String<uint64_t> TLinearIndex_;
    typedef StringSet<CharString, Owner<ConcatDirect<> > > TNameStore;

    int32_t format;             // Format (0: generic; 1: SAM; 2: VCF)
    int32_t colSeq;             // Column for the sequence name
    int32_t colBeg;             // Column for the start of a region
    int32_t colEnd;             // Column for the end of a region
    int32_t meta;               // Leading character for comment lines
    int32_t skip;               // # lines to skip at the beginning
    uint64_t unalignedCount;    // # unmapped reads without coordinates set

    // 1<<14 is the size of the minimum bin.
    static const int32_t BAM_LIDX_SHIFT = 14;

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
        unalignedCount(maxValue<uint64_t>()),
        _nameStoreCache(_nameStore)
    {}

    TabixIndex(char const * fileName) :
        format(0),
        colSeq(1),
        colBeg(2),
        colEnd(3),
        meta('#'),
        skip(0),
        unalignedCount(maxValue<uint64_t>()),
        _nameStoreCache(_nameStore)
    {
        if (!open(*this, fileName))
            SEQAN_THROW(FileOpenError(fileName));
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _tbiReg2bins()
// ----------------------------------------------------------------------------

static inline void
_tbiReg2bins(String<uint16_t> & list, uint32_t beg, uint32_t end)
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
// Function _readTabixRecord()
// ----------------------------------------------------------------------------

template <typename TIter>
bool _readTabixRecord(TabixRecord_ & record, CharString & buffer, TIter & iter, TabixIndex const & index)
{
    clear(record.refName);
    record.posBeg = 0;
    record.posEnd = 0;

    // Skip comment lines.
    while (!atEnd(iter) && *iter == (char)index.meta)
        skipLine(iter);

    if (atEnd(iter))
        return false;
    
    // Extract columns.
    int32_t maxCol = std::max(index.colSeq, std::max(index.colBeg, index.colEnd));
    for (int32_t col = 1; col <= maxCol; ++col)
    {
        // Read column.
        clear(buffer);
        if (col < maxCol)
        {
            readUntil(buffer, iter, IsTab());
            skipOne(iter, IsTab());
        }
        else
        {
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
        }
        
        if (col == index.colSeq)
            record.refName = buffer;
        else if (col == index.colBeg)
            record.posBeg = lexicalCast<int32_t>(buffer);
        else if (col == index.colEnd)
            record.posEnd = lexicalCast<int32_t>(buffer);
    }

    if (index.colEnd == 0 || index.colEnd == index.colBeg)
        record.posEnd = record.posBeg + 1;
    
    // all text-based file formats are 1-based (we use 0-based positions internally)
    --record.posBeg;
    --record.posEnd;
    
    // Go to next line.
    skipLine(iter);
    return true;
}

// ----------------------------------------------------------------------------
// Function jumpToRegion()
// ----------------------------------------------------------------------------

/*!
 * @fn TabixIndex#jumpToRegion
 * @brief Seek in a tab-separated genome related file using a Tabix index.
 *
 * You provide a region <tt>[posBeg, posEnd)</tt> on the contig <tt>refName</tt> that you want to jump to and the function
 * jumps to the first entry in this region, if any.
 *
 * @signature bool jumpToRegion(fileIn, hasEntries, refName, posBeg, posEnd, index[, firstMatch]);
 *
 * @param[in,out] fileIn        The @link VcfFileIn @endlink, @link GffFileIn @endlink, or @link BedFileIn @endlink to jump with.
 * @param[out]    hasEntries    A <tt>bool</tt> that is set true if the region <tt>[posBeg, posEnd)</tt> has any
 *                              entries.
 * @param[in]     refName       The reference name to jump to.
 * @param[in]     posBeg        The begin of the region to jump to (<tt>int32_t</tt>).
 * @param[in]     posEnd        The end of the region to jump to (<tt>int32_t</tt>).
 * @param[in]     index         The @link TabixIndex @endlink to use for the jumping.
 * @param[in]     firstMatch    A <tt>bool</tt>, if <tt>true</tt> (default) this function seeks to the first
 *                              overlapping record. Otherwise, the function potentially stops before the first
 *                              overlapping record.
 *
 * @return bool true if seeking was successful, false if not.
 *
 * @section Remarks
 *
 * This function fails if <tt>refName</tt>/<tt>pos</tt> wasn't found.
 */

template <typename TFileFormat, typename TSpec, typename TName>
inline bool
jumpToRegion(FormattedFile<TFileFormat, Input, TSpec> & fileIn,
             bool & hasEntries,
             TName const & refName,
             int32_t posBeg,
             int32_t posEnd,
             TabixIndex const & index,
             bool firstMatch = true)
{
    hasEntries = false;

    // Get id of given contig name
    unsigned refId = 0;
    if (!getIdByName(refId, index._nameStoreCache, refName))
        return false;

    // ------------------------------------------------------------------------
    // Compute offset in BGZF file.
    // ------------------------------------------------------------------------
    uint64_t offset = MaxValue<uint64_t>::VALUE;

    // Retrieve the candidate bin identifiers for [posBeg, posEnd).
    String<uint16_t> candidateBins;
    _tbiReg2bins(candidateBins, posBeg, posEnd);

    // Retrieve the smallest required offset from the linear index.
    unsigned windowIdx = posBeg >> 14;  // Linear index consists of 16kb windows.
    uint64_t linearMinOffset = 0;
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
    typedef std::set<uint64_t> TOffsetCandidates;
    TOffsetCandidates offsetCandidates;
    typedef typename Iterator<String<uint16_t>, Rooted>::Type TCandidateIter;
    for (TCandidateIter it = begin(candidateBins, Rooted()); !atEnd(it); goNext(it))
    {
        typedef typename std::map<uint32_t, TabixIndexBinData_>::const_iterator TMapIter;
        TMapIter mIt = index._binIndices[refId].find(*it);
        if (mIt == index._binIndices[refId].end())
            continue;  // Candidate is not in index!

        typedef typename Iterator<String<Pair<uint64_t, uint64_t> > const, Rooted>::Type TBegEndIter;
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

    CharString buffer;
    TabixRecord_ record;
    for (TOffsetCandidateIter candIt = offsetCandidates.begin(); candIt != offsetCandidates.end(); ++candIt)
    {
        setPosition(fileIn, *candIt);

        if (!_readTabixRecord(record, buffer, fileIn.iter, index))
            break;

        if (record.refName != refName)
            continue;  // Wrong contig.
        
        if (!hasEntries || record.posBeg <= posBeg)
        {
            // Found a valid record.
            hasEntries = true;
            offset = *candIt;
        }

        if (record.posBeg >= posEnd)
            break;  // Cannot find overlapping any more.
    }

    if (offset != MaxValue<uint64_t>::VALUE)
    {
        setPosition(fileIn, offset);
        
        if (firstMatch)
        {
            // skip to the first overlapping record
            while (!atEnd(fileIn) && !(posBeg < record.posEnd && record.posBeg < posEnd))
            {
                offset = position(fileIn);
                _readTabixRecord(record, buffer, fileIn.iter, index);
                if (record.posBeg >= posEnd)
                {
                    hasEntries = false;
                    return true;
                }
            }
            setPosition(fileIn, offset);
        }
    }

    // Finding no overlapping records is not an error, hasEntries is false.
    return true;
}

// ----------------------------------------------------------------------------
// Function getUnalignedCount()
// ----------------------------------------------------------------------------

/*!
 * @fn TabixIndex#getUnalignedCount
 * @brief Query index for number of unaligned reads.
 *
 * @signature uint64_t getUnalignedCount(index);
 *
 * @param[in] index     Index to query.
 * @return    uint64_t  The number of unaligned reads.
 */

inline uint64_t
getUnalignedCount(TabixIndex const & index)
{
    return index.unalignedCount;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn TabixIndex#open
 * @brief Load a Tabix index from file.
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
    int32_t nRef = 0;
    readRawPod(nRef, iter);
    readRawPod(index.format, iter);
    readRawPod(index.colSeq, iter);
    readRawPod(index.colBeg, iter);
    readRawPod(index.colEnd, iter);
    readRawPod(index.meta, iter);
    readRawPod(index.skip, iter);

    // Read concatenated names.
    int32_t lNm = 0;
    CharString tmp;
    readRawPod(lNm, iter);
    read(tmp, iter, lNm);
    // Trim last terminating '\0's as they confuse strSplit() below
    while (!empty(tmp) && back(tmp) == '\0')
        resize(tmp, length(tmp) - 1);

    // Split concatenated names at \0's.
    clear(index._nameStore);
    strSplit(index._nameStore, tmp, EqualsChar<'\0'>(), true, nRef - 1);
    refresh(index._nameStoreCache);

    clear(index._linearIndices);
    clear(index._binIndices);
    resize(index._linearIndices, nRef);
    resize(index._binIndices, nRef);

    TabixIndexBinData_ data;
    for (int i = 0; i < nRef; ++i)  // For each reference.
    {
        // Read bin index.
        int32_t nBin = 0;
        readRawPod(nBin, iter);

        for (int j = 0; j < nBin; ++j)  // For each bin.
        {
            uint32_t bin = 0;
            int32_t nChunk = 0;
            readRawPod(bin, iter);
            readRawPod(nChunk, iter);

            resize(data.chunkBegEnds, nChunk);
            for (int k = 0; k < nChunk; ++k)  // For each chunk;
            {
                readRawPod(data.chunkBegEnds[k].i1, iter);
                readRawPod(data.chunkBegEnds[k].i2, iter);
            }

            // Copy bin data into index.
            index._binIndices[i][bin] = data;
        }

        // Read linear index.
        int32_t nIntv = 0;
        readRawPod(nIntv, iter);

        resize(index._linearIndices[i], nIntv);
        for (int j = 0; j < nIntv; ++j)
            readRawPod(index._linearIndices[i][j], iter);
    }


    // Read (optional) number of alignments without coordinate.
    if (!atEnd(iter))
        readRawPod(index.unalignedCount, iter);
    else
        index.unalignedCount = maxValue<uint64_t>();

    return true;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_TABIX_IO_TABIX_INDEX_TBI_H_
