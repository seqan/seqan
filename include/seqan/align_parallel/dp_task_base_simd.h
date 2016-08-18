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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_BASE_SIMD_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_BASE_SIMD_H_

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
template <typename TBuffer, typename TSimdSpec>
class DPScoutState_<DPTiled<TBuffer, TSimdSpec> > :
    public DPScoutState_<DPTiled<TBuffer, void> >,
    public DPScoutState_<TSimdSpec>
{
public:

    DPScoutState_() = default;

    DPScoutState_(TBuffer & horBuffer, TBuffer & verBuffer) :
        DPScoutState_<DPTiled<TBuffer, void> >(horBuffer, verBuffer),
        DPScoutState_<TSimdSpec>()
    {}
};

// ----------------------------------------------------------------------------
// Class DPScout_; DPTiled
// ----------------------------------------------------------------------------

// Overloaded DPScout to store the corresponding buffer for the current dp tile.
template <typename TDPCell, typename TBuffer, typename TSimdSpec>
class DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<TSimdSpec> > > :
    public DPScout_<TDPCell, SimdAlignmentScout<TSimdSpec> >
{
public:
    using TBase = DPScout_<TDPCell, SimdAlignmentScout<TSimdSpec> >;

    DPScoutState_<DPTiled<TBuffer, TSimdSpec> > state = {};

    DPScout_() = default;

    DPScout_(DPScoutState_<DPTiled<TBuffer, TSimdSpec> > & state) :
        TBase(static_cast<DPScoutState_<TSimdSpec>&>(state)),
        state(state)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForSimdAlignment_
// ----------------------------------------------------------------------------

template<typename TAlignmentAlgorithm, typename TBuffer>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<DPTiled<TBuffer, SimdAlignEqualLength> > >
{
    using Type = DPTiled<TBuffer, SimdAlignmentScout<SimdAlignEqualLength> >;
};

template<typename TAlignmentAlgorithm, typename TBuffer, typename TTraits>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<DPTiled<TBuffer, SimdAlignVariableLength<TTraits> > > >
{
    using Type = DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > >;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _scoutBestScore()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TSimdSpec,
          typename TTraceMatrixNavigator,
          typename TIsLastColumn,
          typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<TSimdSpec> > > & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<TSimdSpec> > >::TBase;
    _scoutBestScore(static_cast<TScoutBase&>(dpScout), activeCell, navigator, TIsLastColumn(), TIsLastRow());
}

// ----------------------------------------------------------------------------
// Function _preInitScoutHorizontal()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits>
inline void
_preInitScoutHorizontal(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _preInitScoutHorizontal(static_cast<TScoutBase&>(scout));
}

// ----------------------------------------------------------------------------
// Function _preInitScoutVertical()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits>
inline void
_preInitScoutVertical(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _preInitScoutVertical(static_cast<TScoutBase&>(scout));
}

// ----------------------------------------------------------------------------
// Function _reachedHorizontalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits, typename TIter>
inline bool
_reachedHorizontalEndPoint(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout,
                           TIter const & hIt)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    return _reachedHorizontalEndPoint(static_cast<TScoutBase&>(scout), hIt);
}

// ----------------------------------------------------------------------------
// Function _reachedVerticalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits, typename TIter>
inline bool
_reachedVerticalEndPoint(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout,
                         TIter const & vIt)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    return _reachedVerticalEndPoint(static_cast<TScoutBase&>(scout), vIt);
}

// ----------------------------------------------------------------------------
// Function _nextHorizontalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits>
inline void
_nextHorizontalEndPos(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _nextHorizontalEndPos(static_cast<TScoutBase&>(scout));
}

// ----------------------------------------------------------------------------
// Function _nextVerticalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits>
inline void
_nextVerticalEndPos(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _nextVerticalEndPos(static_cast<TScoutBase&>(scout));
}

// ----------------------------------------------------------------------------
// Function _incHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits>
inline void
_incHorizontalPos(DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _incHorizontalPos(static_cast<TScoutBase&>(scout));
}

// ----------------------------------------------------------------------------
// Function _incVerticalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TBuffer, typename TTraits>
inline void
_incVerticalPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout)
{
    using TScoutBase = typename DPScout_<TDPCell, DPTiled<TBuffer, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > >::TBase;
    _incVerticalPos(static_cast<TScoutBase&>(scout));
}


namespace impl
{
template <typename TDPCell, typename TTrace,
          typename TTasks,
          typename TBuffer,
          typename TPos,
          typename TFunc>
inline void
loadIntoSimd(Pair<TDPCell, TTrace> & target,
             TTasks const & tasks,
             TBuffer const & buffer,
             TPos const pos,
             TFunc && f,
             AffineGaps const & /*unsused*/)
{
    using TSimdVec = typename Value<TDPCell>::Type;
    using TVecVal = typename Value<TSimdVec>::Type;
    using TPair = typename std::decay<decltype(buffer[0][0])>::type;

    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreHorVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVerVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> traceVec;

    auto zipCont = makeZipView(tasks, scoreVec, scoreHorVec, scoreVerVec, traceVec);

//    std::for_each(begin(zipCont, Standard()), end(zipCont, Standard()),
    std::for_each(begin(zipCont), end(zipCont),
            [&](auto tuple)
            {
                auto blockId = f(std::get<0>(tuple));
                TPair val = (length(buffer[blockId]) > pos) ? buffer[blockId][pos] : TPair{};
                // We might access values out of bounds here.
                std::get<1>(tuple) = val.i1._score;
                std::get<2>(tuple) = val.i1._horizontalScore;
                std::get<3>(tuple) = val.i1._verticalScore;
                std::get<4>(tuple) = val.i2;
            });

    target.i1._score = load<TSimdVec>(&scoreVec[0]);
    target.i1._horizontalScore = load<TSimdVec>(&scoreHorVec[0]);
    target.i1._verticalScore = load<TSimdVec>(&scoreVerVec[0]);
    target.i2 = load<TSimdVec>(&traceVec[0]);
}

template <typename TBuffer,
          typename TTasks,
          typename TDPCell, typename TTrace,
          typename TPos,
          typename TFunc>
inline void
storeIntoBuffer(TBuffer & buffer,
                TTasks const & tasks,
                Pair<TDPCell, TTrace> const & source,
                TPos const pos,
                TFunc && f,
             	AffineGaps const & /*unsused*/)
{
    using TSimdVec = typename Value<TDPCell>::Type;
    using TVecVal = typename Value<TSimdVec>::Type;

    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreHorVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVerVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> traceVec;

    storeu(&scoreVec[0], source.i1._score);
    storeu(&scoreHorVec[0], source.i1._horizontalScore);
    storeu(&scoreVerVec[0], source.i1._verticalScore);
    storeu(&traceVec[0], source.i2);

    auto zipCont = makeZipView(tasks, scoreVec, scoreHorVec, scoreVerVec, traceVec);

    //    std::for_each(begin(zipCont, Standard()), end(zipCont, Standard()),
    std::for_each(begin(zipCont), end(zipCont),
            [&](auto tuple)
            {
                auto blockId = f(std::get<0>(tuple));
                if (length(buffer[blockId]) > pos)
                {
                    auto& pair = buffer[blockId][pos];
                    pair.i1._score = std::get<1>(tuple);
                    pair.i1._horizontalScore = std::get<2>(tuple);
                    pair.i1._verticalScore = std::get<3>(tuple);
                    pair.i2 = std::get<4>(tuple);
                }
            });
}

template <typename TTasks, typename TBufferSet, typename TFunc>
inline auto
gatherSimdBuffer(bool & allSameLength,
                 TTasks const & tasks,
                 TBufferSet const & buffer,
                 TFunc && getBlockId)
{
    using TBuffer     = typename Value<TBufferSet>::Type;
    using TBuffValue  = typename Value<TBuffer>::Type;
    using TDPCell     = typename Value<TBuffValue, 1>::Type;
    using TDPCellSpec = typename Spec<TDPCell>::Type;

    // TODO(rrahn): Pass simd type from outside.
    using TSimdVec = typename SimdVector<short>::Type;
    using TDPCellSimd = DPCell_<TSimdVec, TDPCellSpec>;
    using TTraceValue = typename TraceBitMap_<TSimdVec>::Type;
    using TBufferValue = Pair<TDPCellSimd, TTraceValue>;

    // Check for valid simd length.
    SEQAN_ASSERT_EQ(LENGTH<TSimdVec>::VALUE, length(tasks));

    String<TBufferValue, Alloc<OverAligned> > simdSet;

    auto maxLength = length(buffer[getBlockId(tasks[0])]);
    std::for_each(begin(tasks) + 1, end(tasks),
                  [&](auto& task)
                  {
                      auto len = length(buffer[getBlockId(task)]);
                      if (len != maxLength)
                          allSameLength = false;
                      maxLength = (len > maxLength) ? len : maxLength;
                  });

    resize(simdSet, maxLength, Exact());
    for (unsigned i = 0; i < length(simdSet); ++i)
    {
        loadIntoSimd(simdSet[i], tasks, buffer, i, std::forward<TFunc>(getBlockId), TDPCellSpec());
    }
    return simdSet;
}

template <typename TBufferSet,
          typename TTasks,
          typename TBufferValue, typename TSpec,
          typename TFunc>
inline void
scatterSimdBuffer(TBufferSet & buffer,
                  TTasks const & tasks,
                  String<TBufferValue, TSpec> const & simdSet,
                  TFunc && getBlockId)
{
    using TBuffer     = typename Value<TBufferSet>::Type;
    using TBuffValue  = typename Value<TBuffer>::Type;
    using TDPCell     = typename Value<TBuffValue, 1>::Type;  // Buffer value is Pair<DPCell, TTraceValue>
    using TDPCellSpec = typename Spec<TDPCell>::Type;

    // Check for valid simd length.
    for (unsigned i = 0; i < length(simdSet); ++i)
    {
        storeIntoBuffer(buffer, tasks, simdSet[i], i, std::forward<TFunc>(getBlockId), TDPCellSpec());
    }
}

template <typename TTasks,
          typename TDPCell, typename TTraceValue, typename TScoreMat, typename TTraceMat,
          typename TBuffer,
          typename TSetSeqH,
          typename TSetSeqV,
          typename TScore,
          typename TBand,
          typename TDPConfig>
inline void computeSimdBatch(TTasks const & tasks,
                             DPContext<TDPCell, TTraceValue, TScoreMat, TTraceMat> & dpContext,
                             DPScoutState_<DPTiled<TBuffer, void> > & state,
                             TSetSeqH const & seqH,
                             TSetSeqV const & seqV,
                             TScore const & scoringScheme,
                             TBand const & band,
                             bool const allSameLength,
                             TDPConfig const & config)
{
    using TSeqH = typename Value<TSetSeqH>::Type;
    using TSeqV = typename Value<TSetSeqV>::Type;
    using TSimdVec = typename Value<TDPCell>::Type;

    // Prepare sequence set.
    StringSet<TSeqH, Dependent<> > depSetH;
    StringSet<TSeqV, Dependent<> > depSetV;
    for (auto& task : tasks)
    {
        appendValue(depSetH, seqH[task->_col]);
        appendValue(depSetV, seqV[task->_row]);
    }

    // Dummy trace set.
    StringSet<String<Nothing> > trace;  // We need to instantiate it, but it will not be used.

    // Create a SIMD scoring scheme.
    Score<TSimdVec, ScoreSimdWrapper<TScore> > simdScoringScheme(scoringScheme);

    // Preapare and run alingment.
    String<TSimdVec, Alloc<OverAligned> > stringSimdH;
    String<TSimdVec, Alloc<OverAligned> > stringSimdV;

    if (allSameLength)
    {
        DPScoutState_<DPTiled<TBuffer, SimdAlignEqualLength> > simdScoutState(*state.ptrHorBuffer, *state.ptrVerBuffer);
        _prepareSimdAlignment(stringSimdH, stringSimdV, depSetH, depSetV, simdScoutState);
        _computeAlignment(dpContext, trace, simdScoutState, stringSimdH, stringSimdV, simdScoringScheme,
                          band, TDPConfig(), false, false);
    }
    else
    {
        using TSimdScoutTrait = SimdAlignVariableLengthTraits<TSimdVec, decltype(depSetH), decltype(depSetV)>;
        DPScoutState_<DPTiled<TBuffer, SimdAlignVariableLength<TSimdScoutTrait> > > simdScoutState(*state.ptrHorBuffer, *state.ptrVerBuffer);
        _prepareSimdAlignment(stringSimdH, stringSimdV, depSetH, depSetV, simdScoutState);

        simdScoutState.dimV = length(stringSimdV);
        simdScoutState.isLocalAlignment = IsLocalAlignment_<TDPConfig>::VALUE;
        simdScoutState.right = IsFreeEndGap_<TDPConfig, DPLastColumn>::VALUE;
        simdScoutState.bottom = IsFreeEndGap_<TDPConfig, DPLastRow>::VALUE;

        _computeAlignment(dpContext, trace, simdScoutState, stringSimdH, stringSimdV, simdScoringScheme,
                          band, TDPConfig(), false, false);
    }
}

}  // namespace impl

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_BASE_SIMD_H_
