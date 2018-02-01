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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_SCORE_SCORE_SIMD_WRAPPER_H_
#define INCLUDE_SEQAN_SCORE_SCORE_SIMD_WRAPPER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TScoreSpec>
struct ScoreSimdWrapper;

// ----------------------------------------------------------------------------
// Class Score<SimdWrapper>
// ----------------------------------------------------------------------------

/*!
 * @class ScoreSimdWrapper
 * @headerfile <seqan/score.h>
 * @extends SimpleScore
 * @brief A wrapper class to extend score calsses for SIMD usage.
 *
 * @signature template<typename TScoreVec, typename TScoreSpec>
 *            class Score<TScoreVec, ScoreSimdWrapper<TScoreSpec> >;
 *
 * @tparam TScoreVec    The simd vector type.
 * @tparam TScoreSpec   The score type to be wrapped.
 *
 * Extends a @link SimpleScore @endlink or a @link MatrixScore @endlink for the use with simd vectors.
 */

template <typename TScoreVec, typename TScore>
class Score<TScoreVec, ScoreSimdWrapper<TScore> >
{
public:
    using TVecValue = typename Value<TScoreVec>::Type;
    using TBaseScoreSpec = typename Spec<TScore>::Type;
    using TBaseScore = Score<typename IfC<sizeof(TVecValue) <= 2,
                                          int32_t,
                                          typename IfC<sizeof(TVecValue) == 8, int64_t, TVecValue>::Type
                                         >::Type,
                             TBaseScoreSpec>;

    // We can be either a SimpleScore or a ScoreMatrix.
    TScoreVec data_match        = createVector<TScoreVec>(0);
    TScoreVec data_mismatch     = createVector<TScoreVec>(-1);
    TScoreVec data_gap_extend   = createVector<TScoreVec>(-1);
    TScoreVec data_gap_open     = createVector<TScoreVec>(-1);

    TBaseScore _baseScore;   // Only needed for the ScoreMatrix data table.

    // Default Constructor.
    Score()
    {}

    // Constructing from base score.
    template <typename TScoreVal2, typename TScoreSpec2>
    Score(Score<TScoreVal2, TScoreSpec2> const & pScore,
          SEQAN_CTOR_ENABLE_IF(And<IsSameType<TScoreSpec2, Simple>, IsSameType<TScoreSpec2, TBaseScoreSpec> >)) :
            data_match(createVector<TScoreVec>(scoreMatch(pScore))),
            data_mismatch(createVector<TScoreVec>(scoreMismatch(pScore))),
            data_gap_extend(createVector<TScoreVec>(scoreGapExtend(pScore))),
            data_gap_open(createVector<TScoreVec>(scoreGapOpen(pScore)))
    {
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TScoreVal2, typename TScoreSpec2>
    Score(Score<TScoreVal2, TScoreSpec2> const & pScore,
          SEQAN_CTOR_ENABLE_IF(And<IsScoreMatrix_<TScoreSpec2>, IsSameType<TScoreSpec2, TBaseScoreSpec> >)) :
            data_gap_extend(createVector<TScoreVec>(scoreGapExtend(pScore))),
            data_gap_open(createVector<TScoreVec>(scoreGapOpen(pScore))),
            _baseScore(pScore)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    //TODO(rrahn): implement the assignment operator.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function baseScore();
// ----------------------------------------------------------------------------

template <typename TValue, typename TScore>
inline auto &
underlying(Score<TValue, ScoreSimdWrapper<TScore> > & me)
{
    return me._baseScore;
}

template <typename TValue, typename TScore>
inline auto const &
underlying(Score<TValue, ScoreSimdWrapper<TScore> > const & me)
{
    return me._baseScore;
}

// ----------------------------------------------------------------------------
// Function score(); SimpleScore Wrapper
// ----------------------------------------------------------------------------

template <unsigned LENGTH>
struct VectorLength_
{};

template <typename TValue, typename TScore, typename TSeqHVal, typename TSeqVVal>
inline SEQAN_FUNC_DISABLE_IF(IsScoreMatrix_<TScore>, TValue)
score(Score<TValue, ScoreSimdWrapper<TScore> > const & me, TSeqHVal const & valH, TSeqVVal const & valV)
{
    return blend(scoreMismatch(me), scoreMatch(me), cmpEq(valH, valV));
}

// ----------------------------------------------------------------------------
// Function score(); ScoreMatrix Wrapper
// ----------------------------------------------------------------------------

template <typename TValue, typename TScore, typename TVal1, typename TVal2>
inline SEQAN_FUNC_ENABLE_IF(IsScoreMatrix_<TScore>, TValue)
score(Score<TValue, ScoreSimdWrapper<TScore> > const & sc, TVal1 const & val1, TVal2 const & val2)
{
    return gather(&sc._baseScore.data_tab[0], val1 + val2);
}

}

#endif  // #ifndef INCLUDE_SEQAN_SCORE_SCORE_SIMD_WRAPPER_H_
