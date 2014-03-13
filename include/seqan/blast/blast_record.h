// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Module for handling NCBI Blast I/O and E-Value computation
// ==========================================================================


#ifndef SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_
#define SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_

namespace seqan {



//TODO: doc

template <typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned>
struct BlastMatch
{
    TQId            qId;
    TSId            sId;

    TPos            qStart;
    TPos            qEnd;
    TPos            sStart;
    TPos            sEnd;

    TPos            sLength;

    TAlign          align;

    short           qFrameShift;
    short           sFrameShift;

    long            score;
    TPos            aliLength;
    TPos            identities;
    TPos            positives;
    TPos            mismatches;
    TPos            gaps;
    TPos            gapOpenings;

    double          eVal;
    double          bitScore;

    BlastMatch() :
        qId(TQId()), sId(TSId()), qStart(0), qEnd(0), sStart(0), sEnd(0),
        sLength(0), qFrameShift(0), sFrameShift(0), score(0), aliLength(0),
        identities(0), positives(0), mismatches(0), gaps(0), gapOpenings(0),
        eVal(0), bitScore(0)
    {}

    BlastMatch(TQId  & _qId, TSId _sId) :
        qId(_qId), sId(_sId), qStart(0), qEnd(0), sStart(0), sEnd(0),
        sLength(0), qFrameShift(0), sFrameShift(0), score(0), aliLength(0),
        identities(0), positives(0), mismatches(0), gaps(0), gapOpenings(0),
        eVal(0), bitScore(0)
    {}

//     BlastMatch(TQId && _qId, TSId && _sId) :
//         qId(std::move(_qId)), sId(std::move(_sId)), qStart(0), qEnd(0), sStart(0), sEnd(0),
//         sLength(0), qFrameShift(0), sFrameShift(0), score(0), aliLength(0),
//         identities(0), positives(0), mismatches(0), gaps(0), gapOpenings(0),
//         eVal(0), bitScore(0)
//     {}

    inline bool operator==(BlastMatch const & bm2) const
    {
        #ifndef SEQAN_CXX11_STANDARD //C++98
        if (qId != bm2.qId)
            return false;
        if (sId != bm2.sId)
            return false;
        if (qStart != bm2.qStart)
            return false;
        if (qEnd != bm2.qEnd)
            return false;
        if (sStart != bm2.sStart)
            return false;
        if (sEnd != bm2.sEnd)
            return false;
        if (align != bm2.align)
            return false;
        if (score != bm2.score)
            return false;
        if (aliLength != bm2.aliLength)
            return false;
        if (identities != bm2.identities)
            return false;
        if (positives != bm2.positives)
            return false;
        if (mismatches != bm2.mismatches)
            return false;
        if (gaps != bm2.gaps)
            return false;
        if (gapOpenings != bm2.gapOpenings)
            return false;
        if (eVal != bm2.eVal)
            return false;
        if (bitScore != bm2.bitScore)
            return false;
        return true;
        #else
        return std::tie(qId,
                        sId,
                        qStart,
                        qEnd,
                        sStart,
                        sEnd,
                        align,
                        score,
                        aliLength,
                        identities,
                        positives,
                        mismatches,
                        gaps,
                        gapOpenings,
                        eVal,
                        bitScore)
            == std::tie(bm2.qId,
                        bm2.sId,
                        bm2.qStart,
                        bm2.qEnd,
                        bm2.sStart,
                        bm2.sEnd,
                        bm2.align,
                        bm2.score,
                        bm2.aliLength,
                        bm2.identities,
                        bm2.positives,
                        bm2.mismatches,
                        bm2.gaps,
                        bm2.gapOpenings,
                        bm2.eVal,
                        bm2.bitScore);
        #endif
    }

    inline void operator=(BlastMatch const & bm2) const
    {
        #ifndef SEQAN_CXX11_STANDARD //C++98
        //TODO
        #else
        std::make_tuple(qId,
                        sId,
                        qStart,
                        qEnd,
                        sStart,
                        sEnd,
                        align,
                        score,
                        aliLength,
                        identities,
                        positives,
                        mismatches,
                        gaps,
                        gapOpenings,
                        eVal,
                        bitScore)
            = std::tie(bm2.qId,
                        bm2.sId,
                        bm2.qStart,
                        bm2.qEnd,
                        bm2.sStart,
                        bm2.sEnd,
                        bm2.align,
                        bm2.score,
                        bm2.aliLength,
                        bm2.identities,
                        bm2.positives,
                        bm2.mismatches,
                        bm2.gaps,
                        bm2.gapOpenings,
                        bm2.eVal,
                        bm2.bitScore);
        #endif
    }

    inline bool operator< (BlastMatch const & bm2) const
    {
        if (qId >= bm2.qId)
            return false;
//         if (sId >= bm2.sId)
//             return false;
        if (bitScore >= bm2.bitScore)
            return false;
        return true;
    }
};

template <typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned>
struct BlastRecord
{
    typedef         BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;
    TDbName         dbName;
    unsigned long   dbTotalLength;
    unsigned int    dbNumberOfSeqs;

    TQId            qId;
    unsigned long   qLength;
    std::list<TBlastMatch>  matches;
//     String<TBlastMatch> matches;

    BlastRecord() :
        dbName(TDbName()), dbTotalLength(0), dbNumberOfSeqs(0), qId(TQId()),
        qLength(0), matches()
    {}

    BlastRecord(TDbName const & _dbName) :
        dbName(_dbName), dbTotalLength(0), dbNumberOfSeqs(0), qId(TQId()),
        qLength(0), matches()
    {}
    BlastRecord(TDbName const & _dbName, TQId & _qId) :
        dbName(_dbName), dbTotalLength(0), dbNumberOfSeqs(0), qId(_qId),
        qLength(0), matches()
    {}

//     BlastRecord(TDbName const & _dbName, TQId && _qId) :
//         dbName(_dbName), dbTotalLength(0), dbNumberOfSeqs(0), qId(std::move(_qId)),
//         qLength(0)
//     {}

//     BlastRecord(TDbName && _dbName, TQId && _qId) :
//         dbName(_dbName), dbTotalLength(0), dbNumberOfSeqs(0), qId(_qId),
//         qLength(0)
//     {}


};

}

#endif // SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_
