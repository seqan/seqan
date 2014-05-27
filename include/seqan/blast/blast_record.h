// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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

/*!
 * @class BlastMatch
 * @headerfile <seqan/blast.h>
 * @signature struct BlastMatch<TQId, TSId, TPos, TAlign> { ... };
 * @brief An object to hold data members of a blast-match.
 *
 * @tparam TQId  Type of qId, defaults to @link CharString @endlink
 * @tparam TSId  Type of sId, defaults to @link CharString @endlink
 * @tparam TPos  Position type of the sequences, defaults to <tt>uint32_t</tt>
 * @tparam TAlign Type of the @link Align @endlink member, defaults to
 * <tt>Align<CharString, ArrayGaps></tt>
 *
 * @var TQId BlastMatch::qId;
 * @brief verbose Id of the query
 *
 * @var TSId BlastMatch::sId;
 * @brief verbose Id of the subject
 *
 * @var long BlastMatch::score;
 * @brief score of the alignment
 *
 * @var TPos BlastMatch::qStart;
 * @brief start of the query sequence
 *
 * @var TPos BlastMatch::qEnd;
 * @brief end of the query sequence
 *
 * @var TPos BlastMatch::sStart;
 * @brief start of the subject sequence
 *
 * @var TPos BlastMatch::sEnd;
 * @brief start of the subject sequence
 *
 * @var TPos BlastMatch::sLength;
 * @brief length of the subject sequence
 *
 * @var TPos BlastMatch::aliLength;
 * @brief length of the alignment
 *
 * @var TPos BlastMatch::identities;
 * @brief no of identical positions in alignment
 *
 * @var TPos BlastMatch::positives;
 * @brief no of positive scoring positions in alignment
 *
 * @var TPos BlastMatch::mismatches;
 * @brief no of non-identical positions in alignment
 *
 * @var TPos BlastMatch::gaps;
 * @brief no of gap-characters in alignment
 *
 * @var TPos BlastMatch::gapOpenings;
 * @brief number of contiguous gap stretches in alignment
 *
 * @var double BlastMatch::eVal;
 * @brief e-value of the alignment
 *
 * @var double BlastMatch::bitScore;
 * @brief bit-score of the alignment
 *
 * @var char BlastMatch::qFrameShift;
 * @brief one out of { -3, -2, -1, +1, +2, +3 } where the absolute value -1 is
 * shift of the translation frame and a negative sign indicates the reverse
 * complement strand [query sequence]
 *
 * @var char BlastMatch::sFrameShift;
 * @brief one out of { -3, -2, -1, +1, +2, +3 } where the absolute value -1 is
 * shift of the translation frame and a negative sign indicates the reverse
 * complement strand [subject sequence]
 *
 * @var TAlign BlastMatch::align;
 * @brief @link Align @endlink object of the alignment
 *
 * @fn BlastMatch::operator<
 * @signature bool operator< (BlastMatch const & bm2) const
 * @brief only qId and bit-score are compared, so that matches are sorted
 * by query sequence and then per bit-score
 *
 * TODO sees
 *
 */

template <typename TQId = CharString,
          typename TSId = CharString,
          typename TPos = uint32_t,
          typename TAlign = Align<CharString, ArrayGaps>>
struct BlastMatch
{
    TQId            qId;
    TSId            sId;

    long            score;

    TPos            qStart;
    TPos            qEnd;
    TPos            sStart;
    TPos            sEnd;

    TPos            sLength;

    TPos            aliLength;
    TPos            identities;
    TPos            positives;
    TPos            mismatches;
    TPos            gaps;
    TPos            gapOpenings;

    double          eVal;
    double          bitScore;

    signed char     qFrameShift;
    signed char     sFrameShift;

    TAlign          align;

    BlastMatch() :
        qId(TQId()), sId(TSId()), score(0), qStart(0), qEnd(0), sStart(0),
        sEnd(0), sLength(0),  aliLength(0), identities(0), positives(0),
        mismatches(0), gaps(0), gapOpenings(0), eVal(0), bitScore(0),
        qFrameShift(1), sFrameShift(1)
    {}

    BlastMatch(TQId const & _qId, TSId const & _sId) :
        qId(_qId), sId(_sId), score(0), qStart(0), qEnd(0), sStart(0),
        sEnd(0), sLength(0),  aliLength(0), identities(0), positives(0),
        mismatches(0), gaps(0), gapOpenings(0), eVal(0), bitScore(0),
        qFrameShift(1), sFrameShift(1)
    {}

    BlastMatch(TQId && _qId, TSId && _sId) :
        qId(std::move(_qId)), sId(std::move(_sId)), score(0), qStart(0),
        qEnd(0), sStart(0), sEnd(0), sLength(0),  aliLength(0), identities(0),
        positives(0), mismatches(0), gaps(0), gapOpenings(0), eVal(0),
        bitScore(0), qFrameShift(1), sFrameShift(1)
    {}

    inline bool operator==(BlastMatch const & bm2) const
    {
        return std::tie(qId,
                        sId,
                        qStart,
                        qEnd,
                        sStart,
                        sEnd,
//                         align,
                        qFrameShift,
                        sFrameShift,
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
                        bm2.qFrameShift,
                        bm2.sFrameShift,
//                         bm2.align,
                        bm2.score,
                        bm2.aliLength,
                        bm2.identities,
                        bm2.positives,
                        bm2.mismatches,
                        bm2.gaps,
                        bm2.gapOpenings,
                        bm2.eVal,
                        bm2.bitScore);
    }

    inline bool operator< (BlastMatch const & bm2) const
    {
        if (qId >= bm2.qId)
            return false;
        if (bitScore >= bm2.bitScore)
            return false;
        return true;
    }
};

/*!
 * @class BlastRecord
 * @headerfile <seqan/blast.h>
 * @signature struct BlastRecord<TDbName, TQId, TSId, TPos, TAlign> { ... };
 * @brief A record of blast-matches (belonging to one query).
 *
 * @tparam TDbName  Type of dbName, defaults to @link CharString @endlink
 * @tparam TQId  Type of qId, defaults to @link CharString @endlink
 * @tparam TSId  Type of sId, defaults to @link CharString @endlink
 * @tparam TPos  Position type of the sequences, defaults to <tt>uint32_t</tt>
 * @tparam TAlign Type of the @link Align @endlink member, defaults to
 * <tt>Align<CharString, ArrayGaps></tt>
 *
 * @typedef BlastRecord::TBLASTMatch
 * @signature typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;
 * @brief type of the contained matches
 *
 * @var TDbName BlastRecord::dbName;
 * @brief verbose name of the database
 *
 * @var uint64_t BlastRecord::dbTotalLength;
 * @brief summed sequence length of the database
 *
 * @var uint32_t BlastRecord::dbNumberOfSeqs;
 * @brief number of sequences in the database
 *
 * @var TQId BlastRecord::qId;
 * @brief verbose Id of the query
 *
 * @var TPos BlastRecord::qLength;
 * @brief length of the query sequence
 *
 * @var std::list<TBlastMatch> BlastRecord::matches;
 * @brief length of the query sequence
 */

template <typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TPos = uint32_t,
          typename TAlign = Align<CharString, ArrayGaps>>
struct BlastRecord
{
    typedef         BlastMatch<TQId, TSId, TPos, TAlign> TBlastMatch;
    TDbName         dbName;
    uint64_t        dbTotalLength;
    uint32_t        dbNumberOfSeqs;

    TQId            qId;
    TPos            qLength;
    std::list<TBlastMatch>  matches;

    BlastRecord() :
        dbName(TDbName()), dbTotalLength(0), dbNumberOfSeqs(0), qId(TQId()),
        qLength(0), matches()
    {}

    BlastRecord(TDbName const & _dbName) :
        dbName(_dbName), dbTotalLength(0), dbNumberOfSeqs(0), qId(TQId()),
        qLength(0), matches()
    {}

    BlastRecord(TDbName const & _dbName, TQId const &_qId) :
        dbName(_dbName), dbTotalLength(0), dbNumberOfSeqs(0), qId(_qId),
        qLength(0), matches()
    {}

    BlastRecord(TDbName && _dbName, TQId && _qId) :
        dbName(std::move(_dbName)), dbTotalLength(0), dbNumberOfSeqs(0),
        qId(std::move(_qId)), qLength(0)
    {}
};

}

#endif // SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_
