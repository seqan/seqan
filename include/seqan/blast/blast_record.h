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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Module for handling NCBI Blast I/O and E-Value computation
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_
#define SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_

namespace seqan
{

/*!
 * @class BlastMatch
 * @implements AssignableConcept
 * @implements CopyConstructibleConcept
 * @implements DefaultConstructibleConcept
 * @implements EqualityComparableConcept
 * @implements LessThanComparableConcept
 * @headerfile <seqan/blast.h>
 * @signature struct BlastMatch<TAlign, TPos, TQId, TSId> { ... };
 * @brief An data structure to hold a blast match, also known as high-scoring segment pair (HSP)
 *
 * You should set the following members manually: @link BlastMatch::qId @endlink, @link BlastMatch::sId @endlink,
 * @link BlastMatch::qLength @endlink, @link BlastMatch::sLength @endlink,
 * @link BlastMatch::qFrameShift @endlink and @link BlastMatch::sFrameShift @endlink.
 *
 * If you then also set a valid @link BlastMatch::align @endlink you can
 * let the other members be computed by
 * @link BlastMatch#computeAlignmentStats @endlink and
 * @link BlastMatch#computeBitScore @endlink, @link BlastMatch#computeEValue @endlink.
 *
 * @tparam TAlign Type of the @link Align @endlink member, defaults to
 * <tt>Align<CharString, ArrayGaps></tt>
 * @tparam TPos  Position type of the sequences, defaults to <tt>__uint32</tt>
 * @tparam TQId  Type of qId, defaults to std::string
 * @tparam TSId  Type of sId, defaults to std::string
 */

template <typename TAlign_ = Align<CharString, ArrayGaps>,
          typename TPos_ = __uint32,
          typename TQId_ = std::string,
          typename TSId_ = std::string>
struct BlastMatch
{
    typedef TAlign_ TAlign;
    typedef TPos_ TPos;
    typedef TQId_ TQId;
    typedef TSId_ TSId;

    /*!
     * @var TQId BlastMatch::qId;
     * @brief The verbose Id of the query.
     *
     * @var TSId BlastMatch::sId;
     * @brief The verbose Id of the subject.
     */
    TQId            qId;
    TSId            sId;

    /*!
     * @var TPos BlastMatch::qStart;
     * @brief The start of the alignment on the (possibly translated) query sequence.
     *
     * @var TPos BlastMatch::qEnd;
     * @brief The end of the alignment on the (possibly translated) query sequence.
     *
     * @var TPos BlastMatch::sStart;
     * @brief The start of the alignment on the (possibly translated) subject sequence.
     *
     * @var TPos BlastMatch::sEnd;
     * @brief The end of the alignment on the (possibly translated) subject sequence.
     */
    TPos            qStart        = 0;
    TPos            qEnd          = 0;
    TPos            sStart        = 0;
    TPos            sEnd          = 0;

    /*!
     * @var TPos BlastMatch::qLength;
     * @brief The length of the original query sequence (possibly before translation).
     *
     * @var TPos BlastMatch::sLength;
     * @brief The length of the original subject sequence (possibly before translation).
     */
    TPos            qLength       = 0;
    TPos            sLength       = 0;
    /*!
     * @var char BlastMatch::qFrameShift;
     * @brief An indicator for query frame and query strand.
     *
     * one out of { -3, -2, -1, +1, +2, +3 } where the <tt>absolute value - 1</tt> is
     * the shift of the translation frame and a negative sign indicates the reverse
     * complement strand [query sequence, only applies for BlastFormatProgram ==
     * BLASTX | TBLASTX]; -1 implies reverse complement for BLASTN.
     *
     * @var char BlastMatch::sFrameShift;
     * @brief An indicator for subject frame and subject strand.
     *
     * one out of { -3, -2, -1, +1, +2, +3 } where the <tt>absolute value - 1</tt> is
     * the shift of the translation frame and a negative sign indicates the reverse
     * complement strand [subject sequence, only applies for BlastFormatProgram ==
     * TBLASTN | TBLASTX].
     */
    __int8          qFrameShift   = 1;
    __int8          sFrameShift   = 1;

    /*!
     * @var double BlastMatch::eValue;
     * @brief The e-value of the alignment.
     *
     * @var double BlastMatch::bitScore;
     * @brief The bit-score of the alignment.
     */
    double          eValue        = 0;
    double          bitScore      = 0;

    /*!
     * @var AlignmentStats BlastMatch::alignStats
     * @brief An @link AlignmentStats @endlink object holding further stats of the alignment.
     */
    AlignmentStats  alignStats;

    /*!
     * @var TAlign BlastMatch::align;
     * @brief @link Align @endlink object of the alignment.
     */
    TAlign          align;

    /*!
     * @fn BlastMatch::BlastMatch()
     * @brief Constructor, can be called with arguments for qId and sId.
     * @signature BlastMatch::BlastMatch()
     * BlastMatch::BlastMatch(qId, sId)
     */
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
                        sFrameShift
//                         alignStats
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
                        bm2.sFrameShift
//                         bm2.align,
//                         bm2.alignStats
// scores have rounding errors
//                         bm2.eValue,
//                         bm2.bitScore
                       );
    }

    // copy, move and assign implicitly

    /*!
     * @fn BlastMatch::operator<
     * @brief The comparison operator (for sorting by bit-score).
     * @signature bool BlastMatch::operator< (BlastMatch const & bm2) const
     *
     * To facilitate fast sorting of matches in a @link BlastRecord @endlink, only the bit-score is compared. Also
     * large bit-score are sorted to front (i.e. operator< on BlastMatch checks operator>= on the bitScores).
     */

    inline bool operator< (BlastMatch const & bm2) const
    {
        return (bitScore >= bm2.bitScore);
    }

    inline void _clear()
    {
        clear(qId);
        clear(sId);

        qStart        = 0;
        qEnd          = 0;
        sStart        = 0;
        sEnd          = 0;
        qLength       = 0;
        sLength       = 0;

        eValue        = 0;
        bitScore      = 0;
        qFrameShift   = 1;
        sFrameShift   = 1;
        clear(alignStats);
        clear(align.data_rows);
    }

    inline void _maxInitialize()
    {
        qId           = "not init";
        sId           = "not init";

        qStart        = std::numeric_limits<TPos>::max();
        qEnd          = std::numeric_limits<TPos>::max();
        sStart        = std::numeric_limits<TPos>::max();
        sEnd          = std::numeric_limits<TPos>::max();

        qLength       = std::numeric_limits<TPos>::max();
        sLength       = std::numeric_limits<TPos>::max();

        eValue        = std::numeric_limits<double>::max();
        bitScore      = std::numeric_limits<double>::max();

        qFrameShift   = std::numeric_limits<int8_t>::max();
        sFrameShift   = std::numeric_limits<int8_t>::max();

        alignStats.numGaps              = std::numeric_limits<unsigned>::max();
        alignStats.numGapOpens          = std::numeric_limits<unsigned>::max();
        alignStats.numGapExtensions     = std::numeric_limits<unsigned>::max();
        alignStats.numInsertions        = std::numeric_limits<unsigned>::max();
        alignStats.numDeletions         = std::numeric_limits<unsigned>::max();
        alignStats.numMatches           = std::numeric_limits<unsigned>::max();
        alignStats.numMismatches        = std::numeric_limits<unsigned>::max();
        alignStats.numPositiveScores    = std::numeric_limits<unsigned>::max();
        alignStats.numNegativeScores    = std::numeric_limits<unsigned>::max();
        alignStats.alignmentLength      = std::numeric_limits<unsigned>::max();
        alignStats.alignmentSimilarity  = std::numeric_limits<float>::max();
        alignStats.alignmentIdentity    = std::numeric_limits<float>::max();
        alignStats.alignmentScore       = std::numeric_limits<unsigned>::max();

        clear(align.data_rows);
    }
};

inline bool
_memberIsSet(CharString const & in)
{
    return (in != "not init");
}

template <typename TNumber>
inline SEQAN_FUNC_ENABLE_IF(Is<NumberConcept<TNumber> >, bool)
_memberIsSet(TNumber const & in)
{
    return (in != std::numeric_limits<TNumber>::max());
}

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign>
inline void
clear(BlastMatch<TAlign, TPos, TQId, TSId> & match)
{
    match._clear();
}

/*!
 * @class BlastRecord
 * @implements FormattedFileRecordConcept
 * @headerfile <seqan/blast.h>
 * @signature struct BlastRecord<TMatch> { ... };
 * @brief A record of blast-matches (belonging to one query).
 *
 * @tparam TMatch Specialization of @link BlastMatch @endlink
 */

template <typename TBlastMatch_ = BlastMatch<>>
struct BlastRecord
{
    /*!
     * @typedef BlastRecord::TBlastMatch
     * @signature typedef TBlastMatch_ TBlastMatch;
     * @brief type of the contained matches
     */
    typedef TBlastMatch_                TBlastMatch;
    typedef typename TBlastMatch::TQId  TQId;
    typedef typename TBlastMatch::TPos  TPos;

    /*!
     * @var TQId BlastRecord::qId;
     * @brief verbose Id of the query
     */
    TQId            qId;

    /*!
     * @var TPos BlastRecord::qLength;
     * @brief length of the query sequence
     */
    TPos            qLength;

    /*!
     * @var std::list<TBlastMatch> BlastRecord::matches;
     * @brief list of the contained matches
     */
    std::list<TBlastMatch>  matches;

    /*!
     * @fn BlastRecord::BlastRecord()
     * @brief constructor, can be passed the qId
     * @signature BlastRecord::BlastRecord()
     * BlastRecord::BlastRecord(qid)
     */
    BlastRecord() :
        qId(TQId()), qLength(0), matches()
    {}

    BlastRecord(TQId const & _qId) :
        qId(_qId), qLength(0), matches()
    {}

    BlastRecord(TQId && _qId) :
        qId(std::move(_qId)), qLength(0), matches()
    {}

    // copy, move and assign implicitly
};

template <typename TMatch>
inline void
clear(BlastRecord<TMatch> & blastRecord)
{
    clear(blastRecord.qId);
    blastRecord.qLength = 0;
    clear(blastRecord.matches);
}

}

#endif // SEQAN_EXTRAS_BLAST_BLAST_RECORD_H_
