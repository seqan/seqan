// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// The terminator specialization of dp scout offers the possibility to have
// dp generation stop, if specified criteria are met.
// To do this, define HasTerminationCriterium_<> for your algorithm and
// implement a DPScoutState for your terminator specialization. In your
// overloaded _scoutBestScore() or _computeCell() you can call
// terminateScout() on your Scout to have DP-generation stop.
// see dp_scout_xdrop.h for an example.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Terminator_
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Terminator_;

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

template <typename TDPCell, typename TSpec>
class DPScout_;

// The default implementation of the dp scout simply stores one maximum
// and its corresponding position.
//
// The state must be a Nothing and is left untouched and unused.
template <typename TDPCell, typename TSpec>
class DPScout_
{
public:
    using TScoreValue = typename Value<TDPCell>::Type;

    TDPCell _maxScore         = TDPCell();
    uint32_t _maxHostPosition = DPCellDefaultInfinity<TScoreValue>::VALUE; // The corresponding host position within the underlying dp-matrix.

    DPScout_() = default;

    DPScout_(DPScoutState_<Default> const & /*state*/) {}
};

// Terminator_ Specialization
template <typename TDPCell, typename TSpec>
class DPScout_<TDPCell, Terminator_<TSpec> >
    : public DPScout_<TDPCell, Default>
{
public:
    typedef DPScout_<TDPCell, Default>  TParent;

    DPScoutState_<Terminator_<TSpec> > * state = nullptr;
    bool terminationCriteriumMet               = false;

    DPScout_() = default;

    DPScout_(DPScoutState_<Terminator_<TSpec> > & pState) :
        DPScout_<TDPCell, Default>(),
        state(&pState)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForAlignmentAlgorithm_
// ----------------------------------------------------------------------------

// Given an alignment algorithm tag such as GlobalAlignment_ or LocalAlignment_, returns the specialization tag for the
// corresponding DPScout_ specialization.

template <typename TAlignmentAlgorithm, typename TScoutState>
struct ScoutSpecForAlignmentAlgorithm_
{
    typedef If<HasTerminationCriterium_<TAlignmentAlgorithm>,
               Terminator_<>,
               Default> Type;
};

// ----------------------------------------------------------------------------
// Metafunction ScoutStateSpecForScout_
// ----------------------------------------------------------------------------

// Given an dp scout this meta-function returns the appropriate specialization for the scout state.

template <typename TScout>
struct ScoutStateSpecForScout_
{
    typedef Default Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _scoutBestScore()
// ----------------------------------------------------------------------------

// Tracks the new score, if it is the new maximum.
template <typename TDPCell, typename TSpec, typename TTraceMatrixNavigator,
          typename TIsLastColumn, typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, TSpec> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    if (_scoreOfCell(activeCell) > _scoreOfCell(dpScout._maxScore))
    {
        dpScout._maxScore = activeCell;
        dpScout._maxHostPosition = position(navigator);
    }
}

// TODO(rmaerker): Why is this needed?
template <typename TDPCell, typename TSpec, typename TTraceMatrixNavigator, typename TIsLastColumn>
inline void
_scoutBestScore(DPScout_<TDPCell, TSpec> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/)
{
    _scoutBestScore(dpScout, activeCell, navigator, TIsLastColumn(), False());
}

// TODO(rmaerker): Why is this needed?
template <typename TDPCell, typename TSpec, typename TTraceMatrixNavigator>
inline void
_scoutBestScore(DPScout_<TDPCell, TSpec> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator)
{
    _scoutBestScore(dpScout, activeCell, navigator, False(), False());
}

// ----------------------------------------------------------------------------
// Function maxScore()
// ----------------------------------------------------------------------------

// Returns the current maximal score.
template <typename TDPCell, typename TScoutSpec>
inline typename Value<TDPCell>::Type const
maxScore(DPScout_<TDPCell, TScoutSpec> const & dpScout)
{
    return _scoreOfCell(dpScout._maxScore);
}

// ----------------------------------------------------------------------------
// Function maxHostPosition()
// ----------------------------------------------------------------------------

// Returns the host position that holds the current maximum score.
template <typename TDPCell, typename TScoutSpec>
inline unsigned int
maxHostPosition(DPScout_<TDPCell, TScoutSpec> const & dpScout)
{
    return dpScout._maxHostPosition;
}

// ----------------------------------------------------------------------------
// Function _terminationCriteriumIsMet()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline bool
_terminationCriteriumIsMet(DPScout_<TDPCell, Terminator_<TSpec> > const & scout)
{
    return scout.terminationCriteriumMet;
}

// ----------------------------------------------------------------------------
// Function terminateScout()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
terminateScout(DPScout_<TDPCell, Terminator_<TSpec> > & scout)
{
    scout.terminationCriteriumMet = true;
}

// ----------------------------------------------------------------------------
// Function _preInitScoutHorizontal()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
_preInitScoutHorizontal(DPScout_<TDPCell, TSpec> const & /*scout*/)
{
    // no-op.
}

// ----------------------------------------------------------------------------
// Function _reachedHorizontalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec, typename TIter>
constexpr inline bool
_reachedHorizontalEndPoint(DPScout_<TDPCell, TSpec> const & /*scout*/,
                           TIter const & /*hIt*/)
{
    return false;
}

// ----------------------------------------------------------------------------
// Function _preInitScoutVertical()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
_preInitScoutVertical(DPScout_<TDPCell, TSpec> const & /*scout*/)
{
    // no-op.
}

// ----------------------------------------------------------------------------
// Function _reachedVerticalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec, typename TIter>
constexpr inline bool
_reachedVerticalEndPoint(DPScout_<TDPCell, TSpec> const & /*scout*/,
                         TIter const & /*iter*/)
{
    return false;
}

// ----------------------------------------------------------------------------
// Function _nextHorizontalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
_nextHorizontalEndPos(DPScout_<TDPCell, TSpec> const & /*scout*/)
{
    // no-op.
}

// ----------------------------------------------------------------------------
// Function _nextVerticalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
_nextVerticalEndPos(DPScout_<TDPCell, TSpec> const & /*scout*/)
{
    // no-op.
}

// ----------------------------------------------------------------------------
// Function _incHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
_incHorizontalPos(DPScout_<TDPCell, TSpec> const & /*scout*/)
{
    // no-op.
}

// ----------------------------------------------------------------------------
// Function _incVerticalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
_incVerticalPos(DPScout_<TDPCell, TSpec> const & /*scout*/)
{
    // no-op.
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_
