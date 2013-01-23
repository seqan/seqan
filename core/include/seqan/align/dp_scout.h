// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// The dp scout is a structure that stores the current maximal score and its
// host position in the underlying dp-matrix.
// This class can be overloaded to implement different behaviors of tracking
// the maximal score, e.g., for the split breakpoint computation.
// ==========================================================================

#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================



// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPScoutState_
// ----------------------------------------------------------------------------

template <typename TSpec>
class DPScoutState_;

template <>
class DPScoutState_<Default> : public Nothing  // empty member optimization
{};

// ----------------------------------------------------------------------------
// Class DPScout_
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSpec>
class DPScout_;

/**
 * The default implementation of the dp scout simply stores one maximum
 * and its corresponding position.
 *
 * The state must be a Nothing and is left untouched and unused.
 */
template <typename TScoreValue>
class DPScout_<TScoreValue, Default>
{
public:
    TScoreValue _maxScore;          // The maximal score.
    unsigned int _maxHostPosition;  // The corrsponding host position within the underlying dp-matrix.

    DPScout_() : _maxScore(MinValue<TScoreValue>::VALUE), _maxHostPosition(0) {}

    DPScout_(DPScoutState_<Default> const & /*state*/) :
        _maxScore(MinValue<TScoreValue>::VALUE), _maxHostPosition(0) {}

    DPScout_(DPScout_ const & other) :
        _maxScore(other._maxScore), _maxHostPosition(other._maxHostPosition) {}

    DPScout_ & operator=(DPScout_ const & other)
    {
        if (this != &other)
        {
            _maxScore = other._maxScore;
            _maxHostPosition = other._maxHostPosition;
        }
        return *this;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForAlignmentAlgorithm_
// ----------------------------------------------------------------------------

// Given an alignment algorithm tag such as GlobalAlignment_ or LocalAlignment_, returns the specialization tag for the
// corresponding DPScout_ specialization.

template <typename TAlignmentAlgorithm>
struct ScoutSpecForAlignmentAlgorithm_
{
    typedef Default Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _scoutBestScore()
// ----------------------------------------------------------------------------
/**
 * Tracks the new score, if it is the new maximum.
 */
template <typename TScoreValue, typename TSpec, typename TTraceMatrixNavigator>
inline void
_scoutBestScore(DPScout_<TScoreValue, TSpec> & dpScout, TScoreValue const & score,
                TTraceMatrixNavigator const & navigator, bool isLastColumn = false,
                bool isLastRow = false)
{
    (void)isLastColumn;
    (void)isLastRow;

    if (score > dpScout._maxScore)
    {
        dpScout._maxScore = score;
        dpScout._maxHostPosition = position(navigator);
    }
}

// ----------------------------------------------------------------------------
// Function maxScore()
// ----------------------------------------------------------------------------
/**
 * Returns the current maximal score.
 */
template <typename TScoreValue, typename TScoutSpec>
inline TScoreValue const
maxScore(DPScout_<TScoreValue, TScoutSpec> const & dpScout)
{
    return dpScout._maxScore;
}

// ----------------------------------------------------------------------------
// Function maxHostPosition()
// ----------------------------------------------------------------------------

/**
 * Returns the host position that holds the current maximum score.
 */

template <typename TScoreValue, typename TScoutSpec>
inline unsigned int
maxHostPosition(DPScout_<TScoreValue, TScoutSpec> const & dpScout)
{
    return dpScout._maxHostPosition;
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_
