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

#ifndef INCLUDE_SEQAN_DP_PARALLEL_DP_PARALLEL_SCOUT_SIMD_H_
#define INCLUDE_SEQAN_DP_PARALLEL_DP_PARALLEL_SCOUT_SIMD_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPScoutState_; DPTiled
// ----------------------------------------------------------------------------

// The overloaded DPScoutState which simply stores the pointers to the corresponding buffer.
template <typename TBuffer, typename TThreadContext, typename TSimdSpec>
class DPScoutState_<DPTiled<TBuffer, TThreadContext, TSimdSpec> > :
    public DPScoutState_<DPTiled<TBuffer, TThreadContext, void> >,
    public DPScoutState_<TSimdSpec>
{
public:

    DPScoutState_() = default;

    DPScoutState_(TBuffer & horBuffer, TBuffer & verBuffer) :
        DPScoutState_<DPTiled<TBuffer, TThreadContext, void> >(horBuffer, verBuffer),
        DPScoutState_<TSimdSpec>()
    {}

    DPScoutState_(TBuffer & horBuffer, TBuffer & verBuffer, TThreadContext && pThreadContext) :
        DPScoutState_<DPTiled<TBuffer, TThreadContext, void> >(horBuffer, verBuffer, std::move(pThreadContext)),
        DPScoutState_<TSimdSpec>()
    {}
};

// ----------------------------------------------------------------------------
// Class DPScout_; DPTiled
// ----------------------------------------------------------------------------

// Overloaded DPScout to store the corresponding buffer for the current dp tile.
template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSimdSpec>
class DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > > :
    public DPScout_<TDPCell, SimdAlignmentScout<TSimdSpec>>
{
public:
    using TBase = DPScout_<TDPCell, SimdAlignmentScout<TSimdSpec> >;

    DPScoutState_<DPTiled<TBuffer, TThreadContext, TSimdSpec> > state;
    size_t   horizontalPos;
    size_t   verticalPos;
    bool  forceTracking;

    DPScout_(DPScoutState_<DPTiled<TBuffer, TThreadContext, TSimdSpec> > & state,
             bool const pForceTracking) :
        TBase(static_cast<DPScoutState_<TSimdSpec>&>(state)),
        state(state),
        forceTracking(pForceTracking)
    {}

    DPScout_(DPScoutState_<DPTiled<TBuffer, TThreadContext, TSimdSpec> > & state) : DPScout_(state, false)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForSimdAlignment_
// ----------------------------------------------------------------------------

template<typename TAlignmentAlgorithm, typename TBuffer, typename TThreadContext>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm,
                                       DPScoutState_<DPTiled<TBuffer, TThreadContext, SimdAlignEqualLength> > >
{
    using Type = DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<SimdAlignEqualLength> >;
};

template<typename TAlignmentAlgorithm, typename TBuffer, typename TThreadContext, typename TTraits>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm,
                                       DPScoutState_<DPTiled<TBuffer,
                                                             TThreadContext,
                                                             SimdAlignVariableLength<TTraits> > > >
{
    using Type = DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > >;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isTrackingEnabled()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSimdSpec>
inline bool
isTrackingEnabled(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > > const & dpScout,
                  True const & /*unused*/,
                  True const & /*unused*/)
{
    // TODO(rrahn): Implement me!
    return (dpScout.forceTracking);
}

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSimdSpec>
inline bool
isTrackingEnabled(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > > const & dpScout,
                  True const & /*unused*/,
                  False const & /*unused*/)
{
    // TODO(rrahn): Implement me!
    return (dpScout.forceTracking);
}

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSimdSpec>
inline bool
isTrackingEnabled(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > > const & dpScout,
                  False const & /*unused*/,
                  True const & /*unused*/)
{
    // TODO(rrahn): Implement me!
    return (dpScout.forceTracking);
}

// ----------------------------------------------------------------------------
// Function _scoutBestScore()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSimdSpec,
          typename TTraceMatrixNavigator,
          typename TIsLastColumn,
          typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > > & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & isLastColumn,
                TIsLastRow const & isLastRow)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec>>>::TBase;
    _scoutBestScore(static_cast<TScoutBase&>(dpScout), activeCell, navigator, isLastColumn, isLastRow);
}

// ----------------------------------------------------------------------------
// Function maxHostCoordinate()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSimdSpec,
typename TDimension>
inline auto
maxHostCoordinate(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > > const & dpScout,
                  TDimension const dimension)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > >::TBase;
    return maxHostCoordinate(static_cast<TScoutBase const &>(dpScout), dimension);
}

// ----------------------------------------------------------------------------
// Function _setSimdLane()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSimdSpec,
typename TPosition>
inline void
_setSimdLane(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > > & dpScout,
             TPosition const pos)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<TSimdSpec> > >::TBase;
    _setSimdLane(static_cast<TScoutBase&>(dpScout), pos);
}

// ----------------------------------------------------------------------------
// Function _preInitScoutHorizontal()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits>
inline void
_preInitScoutHorizontal(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits>>>>::TBase;
    _preInitScoutHorizontal(static_cast<TScoutBase&>(scout));
    scout.horizontalPos = 0;
}

// ----------------------------------------------------------------------------
// Function _preInitScoutVertical()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits>
inline void
_preInitScoutVertical(DPScout_<TDPCell,
                               DPTiled<TBuffer,
                                       TThreadContext,
                                       SimdAlignmentScout<SimdAlignVariableLength<TTraits>>>> & scout)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits>>>>::TBase;
    _preInitScoutVertical(static_cast<TScoutBase&>(scout));
    scout.verticalPos = 0;
}

// ----------------------------------------------------------------------------
// Function _reachedHorizontalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits, typename TIter>
inline bool
_reachedHorizontalEndPoint(DPScout_<TDPCell,
                                    DPTiled<TBuffer,
                                            TThreadContext,
                                            SimdAlignmentScout<SimdAlignVariableLength<TTraits>>>> & scout,
                           TIter const & hIt)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits>>>>::TBase;
    return _reachedHorizontalEndPoint(static_cast<TScoutBase&>(scout), hIt);
}

// ----------------------------------------------------------------------------
// Function _reachedVerticalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits, typename TIter>
inline bool
_reachedVerticalEndPoint(DPScout_<TDPCell,
                                  DPTiled<TBuffer,
                                          TThreadContext,
                                          SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout,
                         TIter const & vIt)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    return _reachedVerticalEndPoint(static_cast<TScoutBase&>(scout), vIt);
}

// ----------------------------------------------------------------------------
// Function _nextHorizontalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits>
inline void
_nextHorizontalEndPos(DPScout_<TDPCell,
                               DPTiled<TBuffer,
                                       TThreadContext,
                                       SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _nextHorizontalEndPos(static_cast<TScoutBase&>(scout));
}

// ----------------------------------------------------------------------------
// Function _nextVerticalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits>
inline void
_nextVerticalEndPos(DPScout_<TDPCell,
                             DPTiled<TBuffer,
                                     TThreadContext,
                                     SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _nextVerticalEndPos(static_cast<TScoutBase&>(scout));
}

// ----------------------------------------------------------------------------
// Function _incHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits>
inline void
_incHorizontalPos(DPScout_<TDPCell,
                           DPTiled<TBuffer,
                                   TThreadContext,
                                   SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _incHorizontalPos(static_cast<TScoutBase&>(scout));
    ++scout.horizontalPos;
}

// ----------------------------------------------------------------------------
// Function _incVerticalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraits>
inline void
_incVerticalPos(DPScout_<TDPCell,
                         DPTiled<TBuffer,
                                 TThreadContext,
                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell,
                                         DPTiled<TBuffer,
                                                 TThreadContext,
                                                 SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _incVerticalPos(static_cast<TScoutBase&>(scout));
    ++scout.verticalPos;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_DP_PARALLEL_DP_PARALLEL_SCOUT_SIMD_H_
