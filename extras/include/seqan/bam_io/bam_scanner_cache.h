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

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_SCANNER_CACHE_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_SCANNER_CACHE_H_

#ifdef SEQAN_CXX11_STANDARD
#include <functional>
#include <unordered_map>
#else
#include <tr1/functional>
#include <tr1/unordered_map>
#endif

namespace seqan {

#ifdef SEQAN_CXX11_STANDARD
using namespace std;
#else
using namespace std::tr1;
#endif

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BamScannerCacheKey_
{
    __int32  rID;
    __int32  beginPos;

    bool operator== (BamScannerCacheKey_ const &other) const
    {
        return rID == other.rID && beginPos == other.beginPos;
    }
};

struct BamScannerCacheSearchKey_    
{
    BamScannerCacheKey_ cacheKey;
    __uint16 flags;
    __uint16 flagsMask;
};

struct BamScannerCacheHash_ :
    std::unary_function<BamScannerCacheKey_, size_t>
{    
    size_t operator()(BamScannerCacheKey_ const &v) const
    {
        return hash<__int32>()(v.rID) ^ hash<__int32>()(v.beginPos);
    }
};

class BamScannerCache
{
public:
    // The Key is a pair of (genomic pos, name) where genomic pos is a pair of (rId, pos).
    typedef String<BamAlignmentRecord> TRecords;
    typedef Size<TRecords>::Type TRecordId;
    typedef BamScannerCacheKey_ TKey;

    // A mapping from the key type to the BamAlignmentRecord at this position.
    typedef unordered_multimap<TKey, TRecordId, BamScannerCacheHash_> TMap;
    typedef TMap::const_iterator TMapIter;

    TRecords            records;
    String<TRecordId>   unusedIds;
    TMap                map;
    BamAlignmentRecord  tmpRecord;

    static const TRecordId INVALID_ID = (TRecordId)-1;
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

inline void
insertRecord(BamScannerCache &cache, seqan::BamAlignmentRecord const &record)
{
    int _id;
    if (empty(cache.unusedIds))
    {
        _id = length(cache.records);
        appendValue(cache.records, record);
    }
    else
    {
        _id = back(cache.unusedIds);
        eraseBack(cache.unusedIds);
        cache.records[_id] = record;
    }

    BamScannerCache::TKey key = { record.rID, record.beginPos };
    cache.map.insert(std::make_pair(key, _id));
}   

template <typename TQName1, typename TQName2>
inline bool
_qNamesEqual(TQName1 const &name1, TQName2 const &name2)
{
    unsigned len1 = length(name1);
    unsigned len2 = length(name2);

    if (len1 != len2)
        return false;

    if (len1 > 2 && len2 > 2)
    {
        if (name1[len1 - 2] == ':' || name1[len1 - 2] == '/')
            return prefix(name1, len1 - 1) == prefix(name2, len2 - 1);
    }
    return name1 == name2;
}


// search a certain segment that might refer to lastRecord
// if recursively search for a intermediate segments that might refer to lastRecord
inline bool
_recursivelyFindSegmentGraph(
    String<BamAlignmentRecord> &records,
    BamScannerCacheSearchKey_ &searchKey,
    unsigned segmentNo,
    BamScannerCache &cache)
{
    typedef BamScannerCache::TMapIter TMapIter;
    
    // search for segment using contigId, position
    std::pair<TMapIter, TMapIter> range = cache.map.equal_range(searchKey.cacheKey);
    for (TMapIter iter = range.first; iter != range.second;)
    {
        TMapIter it = iter;
        ++iter;
        BamAlignmentRecord &record = cache.records[it->second];

        // ... additionally check flags and qName
        if ((record.flag & searchKey.flagsMask) == searchKey.flags &&
            _qNamesEqual(record.qName, records[0].qName))
        {
            if (record.rNextId == records[0].rID && record.pNext == records[0].beginPos)
            {
                resize(records, segmentNo + 2);
                records[segmentNo + 1] = records[0];
                records[segmentNo] = record;
                appendValue(cache.unusedIds, it->second);
                cache.map.erase(it);
                return true;
            }
        }
    }

    for (; range.first != range.second; ++range.first)
    {
        BamAlignmentRecord &record = cache.records[range.first->second];

        // ... additionally check flags and qName
        if ((record.flag & searchKey.flagsMask) == searchKey.flags &&
            _qNamesEqual(record.qName, records[0].qName))
        {
            BamScannerCacheSearchKey_ searchKey = {
                { record.rNextId, record.pNext },
                (record.flag & BAM_FLAG_MULTIPLE) | ((record.flag & BAM_FLAG_NEXT_RC) >> 1),
                BAM_FLAG_MULTIPLE | BAM_FLAG_RC
            };
            if (_recursivelyFindSegmentGraph(records, searchKey, segmentNo + 1, cache))
            {
                records[segmentNo] = record;
                appendValue(cache.unusedIds, range.first->second);
                cache.map.erase(range.first);
                return true;
            }
        }
    }
    
    return false;
}


inline int
readMultiRecords(String<BamAlignmentRecord> &records, BamStream &bamStream, BamScannerCache &cache)
{
    typedef BamScannerCache::TMapIter TMapIter;

    if (empty(records))
        resize(records, 1);

    while (!atEnd(bamStream))
    {
        // read next record
        BamAlignmentRecord &record = records[0];
        int result = readRecord(record, bamStream);
        if (result != 0)
            return result;

        // is this a single-end read or single alignment?
        if (!hasFlagMultiple(record) ||
            record.rID == BamAlignmentRecord::INVALID_REFID ||
            record.beginPos == BamAlignmentRecord::INVALID_POS ||
            record.rNextId == BamAlignmentRecord::INVALID_REFID ||
            record.pNext == BamAlignmentRecord::INVALID_POS)
        {
            resize(records, 1);
            return 0;
        }

        if (record.rID < record.rNextId || (record.rID == record.rNextId && record.beginPos < record.pNext))
        {
            // store record to retrieve it later
            insertRecord(cache, record);
            continue;
        }

        // search mates in case of multiple templates
        BamScannerCacheSearchKey_ searchKey = {
            { record.rNextId, record.pNext },
            (record.flag & BAM_FLAG_MULTIPLE) | ((record.flag & BAM_FLAG_NEXT_RC) >> 1),
            BAM_FLAG_MULTIPLE | BAM_FLAG_RC
        };
        if (_recursivelyFindSegmentGraph(records, searchKey, 0, cache))
        {
            // order paired-end reads by first and second segment in the template
//            if (length(records) == 2)
//            {
//                if (hasFlagLast(records[0]))
//                    std::swap(records[0], records[1]);
//            }
            return 0;
        }
        else
        {
            // we could get here if both mates align to same position and we are the first
            // hence, insert our record to be retrieved by the second.
            if (records[0].beginPos != records[0].pNext)
            {
                std::cerr << "WARNING: Mate could not be found for:\n";
                write2(std::cerr, records[0], bamStream.bamIOContext, seqan::Sam());
            }
            insertRecord(cache, records[0]);
        }
    }
    clear(records);
    return 0;
}

}  // namespace seqan;

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_SCANNER_CACHE_H_
