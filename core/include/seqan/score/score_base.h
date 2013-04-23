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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

// TODO(holtgrew): Should the public interface for the class Score not be defined here?

#ifndef SEQAN_SCORE_SCORE_BASE_H_
#define SEQAN_SCORE_SCORE_BASE_H_

namespace SEQAN_NAMESPACE_MAIN {

/**
.Class.Score:
..cat:Miscellaneous
..summary:A scoring scheme.
..signature:Score<TValue, TSpec>
..param.TValue:The value type.
...default:int
..param.TSpec:The specializing type.
...default:$Tag.Simple$
..include:seqan/score.h
*/
template <typename TValue = int, typename TSpec = Simple>
class Score;


///.Metafunction.Value.param.T.type:Class.Score
template <typename TValue, typename TSpec>
struct Value<Score<TValue, TSpec> > {
    typedef TValue Type;
};

// --------------------------------------------------------------------------
// Metafunction ConsensusScoreSequenceEntry
// --------------------------------------------------------------------------

/**
.Metafunction.Score#SequenceEntryForScore
..cat:Alignments
..class:Class.Score
..signature:SequenceEntryForScore<TScore, TSequence>::Type
..summary:Returns representation type for a character of a position in a sequence.
..description:This is used for unified interfaces for position dependent and independent scores.
..param.TScore:The score type to use.
...type:Class.Score
..param.TSequence:The underlying sequence of the alignments or gaps.
...type:Concept.SequenceConcept
..return:The type to use for the representation of sequence entries.
..see:Metafunction.Score#SequenceEntryForScore
..see:Function.Score#sequenceEntryForScore
..include:seqan/score.h
*/

template <typename TScore, typename TSequence>
struct SequenceEntryForScore
{
    typedef typename Value<TSequence>::Type Type;
};

// --------------------------------------------------------------------------
// Function sequenceEntryForScore
// --------------------------------------------------------------------------

/**
.Function.Score#sequenceEntryForScore
..summary:Helper function for element access, depending on score type.
..cat:Alignments
..signature:sequenceEntryForScore(scoringScheme, seq, pos)
..param.scoringScheme:The scoring scheme to get the representation for.
...type:Class.Score
..param.seq:The sequence to get the representation for.
...type:Concept.SequenceConcept
..param.pos:The position of the character.
..return:
Representation of the character $seq[pos]$ to be used for the given scoring scheme.
The resulting type is $SequenceEntryForScore<TScore, TSequence>::Type$.
*/

// TODO(rmaerker): Check if using iterator instead would be more efficient than subscript operator.
template <typename TScore, typename TSequence, typename TPosition>
inline typename Value<TSequence>::Type
sequenceEntryForScore(TScore const & /*scoringScheme*/, TSequence const & seq, TPosition pos)
{
    return seq[pos];
}

/**
.Function.scoreGapOpenHorizontal
..class:Class.Score
..cat:Scoring
..signature:scoreGapOpenHorizontal(score, entry1, entry2)
..param.score:The scoring scheme to use.
...type:Class.Score
..param.entry1:Entry in sequence one.
...type:Metafunction.Score#SequenceEntryForScore
..param.entry2:Entry in sequence two.
...type:Metafunction.Score#SequenceEntryForScore
..summary:Returns the score for opening a horizontal gap after $entry1$.
..see:Function.Score#sequenceEntryForScore
..see:Class.ConsensusScoreSequenceEntry
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapOpen(me);
}


/**
.Function.scoreGapOpenVertical
..class:Class.Score
..cat:Scoring
..signature:scoreGapOpenVertical(score, entry1, entry2)
..param.score:The scoring scheme to use.
...type:Class.Score
..param.entry1:Entry in sequence one.
...type:Metafunction.Score#SequenceEntryForScore
..param.entry2:Entry in sequence two.
...type:Metafunction.Score#SequenceEntryForScore
..summary:Returns the score for opening a vertical gap after $entry2$.
..see:Function.Score#sequenceEntryForScore
..see:Class.ConsensusScoreSequenceEntry
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapOpenVertical(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapOpen(me);
}


/**
.Function.scoreGapExtendHorizontal
..cat:Scoring
..signature:scoreGapExtendHorizontal(score, entry1, entry2)
..param.score:The scoring scheme to use.
...type:Class.Score
..param.entry1:Entry in sequence one.
...type:Metafunction.Score#SequenceEntryForScore
..param.entry2:Entry in sequence two.
...type:Metafunction.Score#SequenceEntryForScore
..summary:Returns the score for extending a horizontal gap after $entry1$.
..see:Function.Score#sequenceEntryForScore
..see:Class.ConsensusScoreSequenceEntry
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}


/**
.Function.scoreGapExtendVertical
..cat:Scoring
..signature:scoreGapExtendVertical(score, entry1, entry2)
..param.score:The scoring scheme to use.
...type:Class.Score
..param.entry1:Entry in sequence one.
...type:Metafunction.Score#SequenceEntryForScore
..param.entry2:Entry in sequence two.
...type:Metafunction.Score#SequenceEntryForScore
..summary:Returns the score for extending a vertical gap after $entry2$.
..see:Function.Score#sequenceEntryForScore
..see:Class.ConsensusScoreSequenceEntry
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapExtendVertical(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}


/**
.Function.scoreGapHorizontal
..cat:Scoring
..signature:scoreGapHorizontal(score, entry1, entry2)
..param.score:The scoring scheme to use.
...type:Class.Score
..param.entry1:Entry in sequence one.
...type:Metafunction.Score#SequenceEntryForScore
..param.entry2:Entry in sequence two.
...type:Metafunction.Score#SequenceEntryForScore
..summary:Returns the score for a horizontal gap after $entry1$.
..see:Function.Score#sequenceEntryForScore
..see:Class.ConsensusScoreSequenceEntry
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapHorizontal(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGap(me);
}


/**
.Function.scoreGapVertical
..cat:Scoring
..signature:scoreGapVertical(score, entry1, entry2)
..param.score:The scoring scheme to use.
...type:Class.Score
..param.entry1:Entry in sequence one.
...type:Metafunction.Score#SequenceEntryForScore
..param.entry2:Entry in sequence two.
...type:Metafunction.Score#SequenceEntryForScore
..summary:Returns the score for a vertical gap after $entry2$.
..see:Function.Score#sequenceEntryForScore
..see:Class.ConsensusScoreSequenceEntry
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapVertical(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGap(me);
}


/**
.Function.score
..cat:Scoring
..signature:score(score, entry1, entry2)
..param.score:The scoring scheme to use.
...type:Class.Score
..param.entry1:Entry in sequence one.
...type:Metafunction.Score#SequenceEntryForScore
..param.entry2:Entry in sequence two.
...type:Metafunction.Score#SequenceEntryForScore
..summary:Returns the score for aligning the characters $valH$ and $valV$.
..see:Function.Score#sequenceEntryForScore
..see:Class.ConsensusScoreSequenceEntry
..include:seqan/score.h
 */
template <typename TValue, typename TSpec, typename TSeqHVal, typename TSeqVVal>
inline TValue
score(Score<TValue, TSpec> const & me, TSeqHVal valH, TSeqVVal valV) {
    SEQAN_CHECKPOINT;
    if (valH == valV)
        return scoreMatch(me);
    else
        return scoreMismatch(me);
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SCORE_SCORE_BASE_H_
