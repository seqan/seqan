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
 * You should set the following members manually: @link BlastMatch::qId @endlink, @link BlastMatch::sId @endlink,
 * @link BlastMatch::qLength @endlink, @link BlastMatch::sLength @endlink,
 * @link BlastMatch::qFrameShift @endlink and @link BlastMatch::sFrameShift @endlink .
 *
 * If you then also set a valid @link BlastMatch::align @endlink you can
 * let the other members be computed by
 * @link BlastMatch#calcStatsAndScore @endlink and
 * @link BlastMatch#calcBitScoreAndEValue @endlink.
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
 * @var TPos BlastMatch::qLength;
 * @brief length of the query sequences
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
 * @var double BlastMatch::eValue;
 * @brief e-value of the alignment
 *
 * @var double BlastMatch::bitScore;
 * @brief bit-score of the alignment
 *
 * @var char BlastMatch::qFrameShift;
 * @brief one out of { -3, -2, -1, +1, +2, +3 } where the absolute value -1 is
 * shift of the translation frame and a negative sign indicates the reverse
 * complement strand [query sequence, only applies for BlastFormatProgram ==
 * TBLASTN | TBLASTX]
 *
 * @var char BlastMatch::sFrameShift;
 * @brief one out of { -3, -2, -1, +1, +2, +3 } where the absolute value -1 is
 * shift of the translation frame and a negative sign indicates the reverse
 * complement strand [subject sequence, only applies for BlastFormatProgram ==
 * BLASTX | TBLASTX]
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

    long            score         = 0;

    TPos            qStart        = 0;
    TPos            qEnd          = 0;
    TPos            sStart        = 0;
    TPos            sEnd          = 0;

    TPos            qLength       = 0;
    TPos            sLength       = 0;

    TPos            aliLength     = 0;
    TPos            identities    = 0;
    TPos            positives     = 0;
    TPos            mismatches    = 0;
    TPos            gaps          = 0;
    TPos            gapOpenings   = 0;

    double          eValue        = 0;
    double          bitScore      = 0;

    signed char     qFrameShift   = 1;
    signed char     sFrameShift   = 1;

    TAlign          align;

    BlastMatch() :
        qId(TQId()), sId(TSId())
    {}

    BlastMatch(TQId const & _qId, TSId const & _sId) :
        qId(_qId), sId(_sId)
    {}

    BlastMatch(TQId && _qId, TSId && _sId) :
        qId(std::move(_qId)), sId(std::move(_sId))
    {}

    inline bool operator==(BlastMatch const & bm2) const
    {
        return std::tie(qId,
                        sId,
                        qStart,
                        qEnd,
                        sStart,
                        sEnd,
                        qLength,
                        sLength,
// align too big to compare
//                         align,
                        qFrameShift,
                        sFrameShift,
                        score,
                        aliLength,
                        identities,
                        positives,
                        mismatches,
                        gaps,
                        gapOpenings
// scores have rounding errors
//                         eValue,
//                         bitScore
                       )
            == std::tie(bm2.qId,
                        bm2.sId,
                        bm2.qStart,
                        bm2.qEnd,
                        bm2.sStart,
                        bm2.sEnd,
                        bm2.qLength,
                        bm2.sLength,
                        bm2.qFrameShift,
                        bm2.sFrameShift,
//                         bm2.align,
                        bm2.score,
                        bm2.aliLength,
                        bm2.identities,
                        bm2.positives,
                        bm2.mismatches,
                        bm2.gaps,
                        bm2.gapOpenings
// scores have rounding errors
//                         bm2.eValue,
//                         bm2.bitScore
                       );
    }

    inline bool operator< (BlastMatch const & bm2) const
    {
        if (bitScore >= bm2.bitScore)
            return true;
        //TODO check this; comparison should be with numeric id, not strings
//         if (qId <= bm2.qId)
//             return true;
        return false;
    }
};

/*!
 * @class BlastRecord
 * @headerfile <seqan/blast.h>
 * @signature struct BlastRecord<TDbName, TQId, TSId, TPos, TAlign> { ... };
 * @brief A record of blast-matches (belonging to one query).

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
 * @var TQId BlastRecord::qId;
 * @brief verbose Id of the query
 *
 * @var TPos BlastRecord::qLength;
 * @brief length of the query sequence
 *
 * @var std::list<TBlastMatch> BlastRecord::matches;
 * @brief length of the query sequence
 */

template <typename TQId = CharString,
          typename TSId = CharString,
          typename TPos = uint32_t,
          typename TAlign = Align<CharString, ArrayGaps>>
struct BlastRecord
{
    typedef         BlastMatch<TQId, TSId, TPos, TAlign> TBlastMatch;

    TQId            qId;
    TPos            qLength;
    std::vector<TBlastMatch>  matches;

    BlastRecord() :
        qId(TQId()), qLength(0), matches()
    {}

    BlastRecord(TQId const &_qId) :
        qId(_qId), qLength(0), matches()
    {}

    BlastRecord(TQId && _qId) :
        qId(std::move(_qId)), qLength(0)
    {}
};


/*!
 * @class BlastDbSpecs
 * @headerfile <seqan/blast.h>
 * @signature struct BlastRecord<TDbName> { ... };
 * @brief A record of blast-matches (belonging to one query).
 *
 * @tparam TDbName  Type of dbName, defaults to @link CharString @endlink
 *
 * @var TDbName BlastDbSpecs::dbName;
 * @brief verbose name of the database
 *
 * @var uint64_t BlastDbSpecs::dbTotalLength;
 * @brief summed sequence length of the database
 *
 * @var uint32_t BlastDbSpecs::dbNumberOfSeqs;
 * @brief number of sequences in the database
 */

template <typename TDbName = CharString>
struct BlastDbSpecs
{
    TDbName         dbName;
    uint64_t        dbTotalLength;
    uint32_t        dbNumberOfSeqs;

    BlastDbSpecs() :
        dbName(), dbTotalLength(0), dbNumberOfSeqs(0)
    {}

    BlastDbSpecs(TDbName const & _dbName) :
        dbName(_dbName), dbTotalLength(0), dbNumberOfSeqs(0)
    {}

    BlastDbSpecs(TDbName && _dbName) :
        dbName(std::move(_dbName)), dbTotalLength(0), dbNumberOfSeqs(0)
    {}
};

}

#endif // SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_
