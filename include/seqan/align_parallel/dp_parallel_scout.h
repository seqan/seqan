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

#ifndef INCLUDE_SEQAN_DP_PARALLEL_DP_PARALLEL_SCOUT_H_
#define INCLUDE_SEQAN_DP_PARALLEL_DP_PARALLEL_SCOUT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPTileBuffer
// ----------------------------------------------------------------------------

// The structure owning the horizontal/vertical buffer.
template <typename TDPCellBuff, typename TBuffer = String<TDPCellBuff> >
struct DPTileBuffer
{
    TBuffer horizontalBuffer;
    TBuffer verticalBuffer;
};

// ----------------------------------------------------------------------------
// Tag DPTiled<TBuffer>
// ----------------------------------------------------------------------------

// Tag used to subclass DPScoutState and DPScout.
// T represents the buffer type.
template <typename TBuffer, typename TThreadContext = Default, typename TSimdSpec = void>
struct DPTiled;

// ----------------------------------------------------------------------------
// Class DPScoutState_; DPTiled
// ----------------------------------------------------------------------------

// The overloaded DPScoutState which simply stores the pointers to the corresponding buffer.
template <typename TBuffer, typename TThreadContext>
class DPScoutState_<DPTiled<TBuffer, TThreadContext, void> >
{
public:

    using TDPCell = typename Value<typename Value<TBuffer>::Type, 1>::Type;

    TBuffer* ptrHorBuffer = nullptr;
    TBuffer* ptrVerBuffer = nullptr;
    TThreadContext threadContext{};

    DPScoutState_() = default;

    DPScoutState_(TBuffer & horBuffer, TBuffer & verBuffer) :
        ptrHorBuffer(&horBuffer),
        ptrVerBuffer(&verBuffer)
    {}

    DPScoutState_(TBuffer & horBuffer, TBuffer & verBuffer, TThreadContext pThreadContext) :
        ptrHorBuffer(&horBuffer),
        ptrVerBuffer(&verBuffer),
        threadContext(std::move(pThreadContext))
    {}
};

// ----------------------------------------------------------------------------
// Class DPScout_; DPTiled
// ----------------------------------------------------------------------------

// Overloaded DPScout to store the corresponding buffer for the current dp tile.
template <typename TDPCell, typename TBuffer, typename TThreadContext>
class DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, void> > :
    public DPScout_<TDPCell, Default>
{
public:
    using TBase = DPScout_<TDPCell, Default>;

    DPScoutState_<DPTiled<TBuffer, TThreadContext, void> > state;

    size_t   horizontalPos;
    size_t   verticalPos;
    bool     forceTracking;

    DPScout_(DPScoutState_<DPTiled<TBuffer, TThreadContext, void> > state,
             bool pForceTracking = false) :
        TBase(),
        state(state),
        forceTracking(pForceTracking)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForSimdAlignment_
// ----------------------------------------------------------------------------

template<typename TAlignmentAlgorithm, typename TThreadContext, typename TBuffer>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<DPTiled<TBuffer, TThreadContext, void> > >
{
    using Type = DPTiled<TBuffer, TThreadContext, void>;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSpec,
          typename TIsLastColumn,
          typename TIsLastRow>
inline bool
isTrackingEnabled(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, TSpec> > const & /*dpScout*/,
                  TIsLastColumn const & /*unused*/,
                  TIsLastRow const & /*unused*/)
{
    return false;
}

template <typename TDPCell, typename TBuffer, typename TThreadContext>
inline bool
isTrackingEnabled(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, void> > const & dpScout,
                  True const & /*unused*/,
                  True const & /*unused*/)
{
    return (dpScout.forceTracking || (dpScout.state.threadContext.task._lastHBlock &&
                                       dpScout.state.threadContext.task._lastVBlock));
}

template <typename TDPCell, typename TBuffer, typename TThreadContext>
inline bool
isTrackingEnabled(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, void> > const & dpScout,
                  True const & /*unused*/,
                  False const & /*unused*/)
{
    return (dpScout.forceTracking || dpScout.state.threadContext.task._lastHBlock);
}

template <typename TDPCell, typename TBuffer, typename TThreadContext>
inline bool
isTrackingEnabled(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, void> > const & dpScout,
                  False const & /*unused*/,
                  True const & /*unused*/)
{
    return (dpScout.forceTracking || dpScout.state.threadContext.task._lastVBlock);
}

// ----------------------------------------------------------------------------
// Function _scoutBestScore()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer,
          typename TTraceMatrixNavigator,
          typename TIsLastColumn,
          typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, DPTiled<TBuffer, Default, void> > & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & isLastColumn,
                TIsLastRow const & isLastRow)
{
    using TBaseScout = typename DPScout_<TDPCell, DPTiled<TBuffer, Default, void> >::TBase;
    _scoutBestScore(static_cast<TBaseScout&>(dpScout), activeCell, navigator, isLastColumn, isLastRow);
}

// Tracks the new score, if it is the new maximum.
template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TTraceMatrixNavigator,
          typename TIsLastColumn, typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, void> > & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & isLastColumn,
                TIsLastRow const & isLastRow)
{
    using TBaseScout = typename DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, void> >::TBase;
    if (isTrackingEnabled(dpScout, isLastColumn, isLastRow))
        _scoutBestScore(static_cast<TBaseScout&>(dpScout), activeCell, navigator, isLastColumn, isLastRow);
}

// ----------------------------------------------------------------------------
// Function _preInitScoutHorizontal()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSpec>
inline void
_preInitScoutHorizontal(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, TSpec> > & scout)
{
    scout.horizontalPos = 0;
}

// ----------------------------------------------------------------------------
// Function _preInitScoutVertical()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSpec>
inline void
_preInitScoutVertical(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, TSpec> > & scout)
{
    scout.verticalPos = 0;
}

// ----------------------------------------------------------------------------
// Function _incHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSpec>
inline void
_incHorizontalPos(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, TSpec> > & scout)
{
    ++scout.horizontalPos;
}

// ----------------------------------------------------------------------------
// Function _incVerticalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TThreadContext, typename TSpec>
inline void
_incVerticalPos(DPScout_<TDPCell, DPTiled<TBuffer, TThreadContext, TSpec> > & scout)
{
    ++scout.verticalPos;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_DP_PARALLEL_DP_PARALLEL_SCOUT_H_
