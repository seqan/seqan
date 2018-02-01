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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_UTIL_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_UTIL_H_

namespace seqan
{
namespace impl
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// Helper meta-function to extract the correct DP Property.
template <typename TAlgotrithm>
struct AlgorithmProperty
{
    template <typename TTask>
    inline static bool
    isTrackingEnabled(TTask const & tile)
    {
        return isLastColumn(tile) && isLastRow(tile);
    }
};

template <typename TFreeEndGaps>
struct AlgorithmProperty<GlobalAlignment_<TFreeEndGaps>>
{
    template <typename TTask>
    inline static bool
    isTrackingEnabled(TTask const & tile)
    {
        return (IsFreeEndGap_<TFreeEndGaps, DPLastColumn>::VALUE && inLastColumn(tile)) ||
               (IsFreeEndGap_<TFreeEndGaps, DPLastRow>::VALUE && inLastRow(tile)) ||
               (inLastColumn(tile) && inLastRow(tile));
    }
};

template <typename TSpec>
struct AlgorithmProperty<LocalAlignment_<TSpec>>
{
    template <typename TTask>
    inline static bool
    isTrackingEnabled(TTask const & /*tile*/)
    {
        return true;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function computeTile()
// ----------------------------------------------------------------------------

// Wrapper function to call alignment core for the specific block.
template <typename TScoreValue, typename TTraceValue, typename TScoreMatHost, typename TTraceMatHost,
          typename TDPScout,
          typename TSequenceH,
          typename TSequenceV,
          typename TScoringScheme,
          typename TDPSettings>
inline void
computeTile(DPContext<TScoreValue, TTraceValue, TScoreMatHost, TTraceMatHost> & dpContext,
            TDPScout & scout,
            TSequenceH const & seqH,
            TSequenceV const & seqV,
            TScoringScheme const & scoringScheme,
            TDPSettings const & /*settings*/)
{
    using TDPTraits = typename TDPSettings::TTraits;

    using TScoreMatrixSpec = typename DefaultScoreMatrixSpec_<typename TDPTraits::TAlgorithmType>::Type;

    using TDPScoreMatrix = DPMatrix_<TScoreValue, TScoreMatrixSpec, TScoreMatHost>;
    using TDPTraceMatrix = DPMatrix_<TTraceValue, FullDPMatrix, TTraceMatHost>;

    using TDPScoreMatrixNavigator = DPMatrixNavigator_<TDPScoreMatrix, DPScoreMatrix, NavigateColumnWise>;
    using TDPTraceMatrixNavigator = DPMatrixNavigator_<TDPTraceMatrix, DPTraceMatrix<typename TDPTraits::TTracebackType>, NavigateColumnWise>;

    using TDPProfile = DPProfile_<typename TDPTraits::TAlgorithmType,
                                  typename TDPTraits::TGapType,
                                  typename TDPTraits::TTracebackType,
                                  Parallel>;

    // Setup the score and trace matrix.
    TDPScoreMatrix dpScoreMatrix;
    TDPTraceMatrix dpTraceMatrix;

    setLength(dpScoreMatrix, +DPMatrixDimension_::HORIZONTAL, length(seqH) + 1);
    setLength(dpScoreMatrix, +DPMatrixDimension_::VERTICAL, length(seqV) + 1);

    setLength(dpTraceMatrix, +DPMatrixDimension_::HORIZONTAL, length(seqH) + 1);
    setLength(dpTraceMatrix, +DPMatrixDimension_::VERTICAL, length(seqV) + 1);

    // Resue the buffer from the cache.
    setHost(dpScoreMatrix, getDpScoreMatrix(dpContext));
    setHost(dpTraceMatrix, getDpTraceMatrix(dpContext));

    resize(dpScoreMatrix);
    // We do not need to allocate the memory for the trace matrix if the traceback is disabled.
    if /*constexpr*/(IsTracebackEnabled_<typename TDPTraits::TTracebackType>::VALUE)
    {
        static_assert(std::is_same<typename TDPTraits::TTracebackType, TracebackOff>::value, "Traceback not implemented!");
        resize(dpTraceMatrix);
    }

    // Initialize the navigators.
    TDPScoreMatrixNavigator dpScoreMatrixNavigator{dpScoreMatrix, DPBandConfig<BandOff>{}};
    TDPTraceMatrixNavigator dpTraceMatrixNavigator{dpTraceMatrix, DPBandConfig<BandOff>{}};

    // Execute the alignment.
    _computeAlignmentImpl(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator, seqH, seqV,
                          scoringScheme, DPBandConfig<BandOff>{}, TDPProfile(), NavigateColumnWise{});
}

#ifdef SEQAN_SIMD_ENABLED
// Some utility functions.
template <typename TTasks,
          typename TScoreValueScalar,
          typename TScoreValueSimd>
inline auto
doComputeOffset(TTasks const &tasks,
                TScoreValueScalar const & /*scalarScore*/,
                TScoreValueSimd const & /*simdScore*/)
{
    String<TScoreValueScalar> offset;
    resize(offset, length(tasks), std::numeric_limits<TScoreValueScalar>::min(), Exact());

    size_t pos = 0;

    for (auto task : tasks)
    {
            offset[pos] = front(context(*task).tileBuffer.horizontalBuffer[column(*task)]).i1._score;
        ++pos;
    }

    return offset;
}

template <typename TTasks,
          typename TScoreValue>
inline auto
doComputeOffset(TTasks const &tasks,
                TScoreValue const & /*scalarScore*/,
                TScoreValue const & /*simdScore*/)
{
    String<TScoreValue> offset;
    resize(offset, length(tasks), 0, Exact());
    return offset;
}

template <typename TTasks,
          typename TTaskTraits>
inline auto
computeOffset(TTasks const &tasks, TTaskTraits const & /*traits*/)
{
    using TDPSettings       = typename TTaskTraits::TDPSettings;
    using TScoreValueScalar = typename Value<typename TDPSettings::TScoringScheme>::Type;
    using TScoreValueSimd   = typename Value<typename TDPSettings::TSimdScoringScheme>::Type;
    using TDPSimdValue      = typename Value<TScoreValueSimd>::Type;

    return doComputeOffset(tasks, TScoreValueScalar{}, TDPSimdValue{});
}

template <typename TDPCell, typename TTrace,
          typename TTasks,
          typename TPos,
          typename TFunc,
          typename TOffset>
inline void
loadIntoSimd(Pair<TDPCell, TTrace> & target,
             TTasks const & tasks,
             TPos const pos,
             TFunc && getBuffer,
             TOffset const & offset,
             LinearGaps const & /*unsused*/)
{
    using TSimdVec = typename Value<TDPCell>::Type;
    using TVecVal = typename Value<TSimdVec>::Type;

    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> traceVec;

    auto zipCont = makeZipView(tasks, scoreVec, traceVec, offset);

    std::for_each(begin(zipCont), end(zipCont),
                  [&, getBuffer = std::move(getBuffer)](auto tuple)
                  {
                      auto & buffer = *getBuffer(*std::get<0>(tuple));
                      auto val = (length(buffer) > pos) ? buffer[pos] : typename std::decay<decltype(buffer[0])>::type{};

                      // We might access values out of bounds here.
                      std::get<1>(tuple) = static_cast<TVecVal>(val.i1._score - std::get<3>(tuple));
                      std::get<2>(tuple) = val.i2;
                  });

    target.i1._score = load<TSimdVec>(&scoreVec[0]);
    target.i2 = load<TSimdVec>(&traceVec[0]);
}

template <typename TDPCell, typename TTrace,
          typename TTasks,
          typename TPos,
          typename TFunc,
          typename TOffset>
inline void
loadIntoSimd(Pair<TDPCell, TTrace> & target,
             TTasks const & tasks,
             TPos const pos,
             TFunc && getBuffer,
             TOffset const & offset,
             AffineGaps const & /*unsused*/)
{
    using TSimdVec = typename Value<TDPCell>::Type;
    using TVecVal = typename Value<TSimdVec>::Type;

    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreHorVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVerVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> traceVec;

    auto zipCont = makeZipView(tasks, scoreVec, scoreHorVec, scoreVerVec, traceVec, offset);

    std::for_each(begin(zipCont), end(zipCont),
                  [&, getBuffer = std::move(getBuffer)](auto tuple)
                  {
                      auto & buffer = *getBuffer(*std::get<0>(tuple));
                      auto val = (length(buffer) > pos) ? buffer[pos] : typename std::decay<decltype(buffer[0])>::type{};
                      using TDPCellVar = decltype(val.i1);
                      using TDPCell16 = DPCell_<TVecVal, AffineGaps>;

                      // We might access values out of bounds here.
                      std::get<1>(tuple) = static_cast<TVecVal>(val.i1._score - std::get<5>(tuple));

                      std::get<2>(tuple) =
                        (val.i1._horizontalScore <= DPCellDefaultInfinity<TDPCellVar>::VALUE) ?
                            DPCellDefaultInfinity<TDPCell16>::VALUE :
                            static_cast<TVecVal>(val.i1._horizontalScore - std::get<5>(tuple));
                      std::get<3>(tuple) =
                        (val.i1._verticalScore <= DPCellDefaultInfinity<TDPCellVar>::VALUE) ?
                        DPCellDefaultInfinity<TDPCell16>::VALUE :
                        static_cast<TVecVal>(val.i1._verticalScore - std::get<5>(tuple));
                      std::get<4>(tuple) = val.i2;
                  });

    target.i1._score = load<TSimdVec>(&scoreVec[0]);
    target.i1._horizontalScore = load<TSimdVec>(&scoreHorVec[0]);
    target.i1._verticalScore = load<TSimdVec>(&scoreVerVec[0]);
    target.i2 = load<TSimdVec>(&traceVec[0]);
}

template <typename TTasks,
          typename TDPCell, typename TTrace,
          typename TPos,
          typename TFunc,
          typename TOffset>
inline void
storeIntoBuffer(TTasks & tasks,
                Pair<TDPCell, TTrace> const & source,
                TPos const pos,
                TFunc && getBuffer,
                TOffset const & offset,
                LinearGaps const & /*unsused*/)
{
    using TSimdVec = typename Value<TDPCell>::Type;
    using TVecVal = typename Value<TSimdVec>::Type;

    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> scoreVec;
    alignas(sizeof(TSimdVec)) std::array<TVecVal, LENGTH<TSimdVec>::VALUE> traceVec;

    storeu(&scoreVec[0], source.i1._score);
    storeu(&traceVec[0], source.i2);

    auto zipCont = makeZipView(tasks, scoreVec, traceVec, offset);

    std::for_each(begin(zipCont), end(zipCont),
                  [&, getBuffer = std::move(getBuffer)] (auto tuple)
                  {
                      auto & buffer = *getBuffer(*std::get<0>(tuple));
                      if (length(buffer) > pos)
                      {
                          auto & pair = buffer[pos];
                          pair.i1._score = std::get<1>(tuple) + std::get<3>(tuple);
                          pair.i2 = std::get<2>(tuple);
                      }
                  });
}

template <typename TTasks,
          typename TDPCell, typename TTrace,
          typename TPos,
          typename TFunc,
          typename TOffset>
inline void
storeIntoBuffer(TTasks & tasks,
                Pair<TDPCell, TTrace> const & source,
                TPos const pos,
                TFunc && getBuffer,
                TOffset const & offset,
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

    auto zipCont = makeZipView(tasks, scoreVec, scoreHorVec, scoreVerVec, traceVec, offset);

    std::for_each(begin(zipCont), end(zipCont),
                  [&, getBuffer = std::move(getBuffer)](auto tuple)
                  {
                      auto & buffer = *getBuffer(*std::get<0>(tuple));
                      if (length(buffer) > pos)
                      {
                          auto & pair = buffer[pos];
                          pair.i1._score = std::get<1>(tuple) + std::get<5>(tuple);
                          pair.i1._horizontalScore = std::get<2>(tuple) + std::get<5>(tuple);
                          pair.i1._verticalScore = std::get<3>(tuple) + std::get<5>(tuple);
                          pair.i2 = std::get<4>(tuple);
                      }
                  });
}

template <typename TTasks,
          typename TFunc,
          typename TOffset,
          typename TExecTraits>
inline auto
gatherSimdBuffer(TTasks const & tasks,
                 TFunc && getBuffer,
                 TOffset const & offset,
                 TExecTraits const & /*traits*/)
{
    // Check for valid simd length.
    SEQAN_ASSERT_EQ(LENGTH<typename TExecTraits::TScoreValue>::VALUE, length(tasks));

    String<typename TExecTraits::TBufferValue, Alloc<OverAligned> > simdSet;

    auto maxLength = length(*getBuffer(*tasks[0]));
    std::for_each(begin(tasks, Standard()) + 1, end(tasks, Standard()),
                  [&](auto & task)
                  {
                      auto len = length(*getBuffer(*task));
                      maxLength = (len > maxLength) ? len : maxLength;
                  });

    resize(simdSet, maxLength, Exact());
    for (unsigned i = 0; i < length(simdSet); ++i)
    {
        loadIntoSimd(simdSet[i], tasks, i, std::forward<TFunc>(getBuffer), offset, typename TExecTraits::TGapType());
    }
    return simdSet;
}

template <typename TTasks,
          typename TBufferValue, typename TSpec,
          typename TFunc,
          typename TOffset,
          typename TExecTraits>
inline void
scatterSimdBuffer(TTasks & tasks,
                  String<TBufferValue, TSpec> const & simdSet,
                  TFunc && getBuffer,
                  TOffset const & offset,
                  TExecTraits const & /*traits*/)
{
    for (unsigned i = 0; i < length(simdSet); ++i)
    {
        storeIntoBuffer(tasks, simdSet[i], i, std::forward<TFunc>(getBuffer), offset, typename TExecTraits::TGapType());
    }
}

// Compute tasks as simd alignment.
template <typename TDPCell, typename TTraceValue, typename TScoreMat, typename TTraceMat,
          typename TTasks,
          typename TSimdBufferH,
          typename TSimdBufferV,
          typename TDPLocal,
          typename TOffset,
          typename TExecTraits>
inline void
computeSimdBatch(DPContext<TDPCell, TTraceValue, TScoreMat, TTraceMat> & cache,
                 TSimdBufferH                                          & bufferH,
                 TSimdBufferV                                          & bufferV,
                 TTasks                                                & tasks,
                 TDPLocal                                              & dpLocal,
                 TOffset                                               & offset,
                 TExecTraits                                     const & /*traits*/)
{
    // Now what?
    using TSeqH    = typename TExecTraits::TSeqH;
    using TSeqV    = typename TExecTraits::TSeqV;
    using TSimdVec = typename TExecTraits::TScoreValue;

    // Prepare sequence set.
    StringSet<TSeqH, Dependent<> > depSetH;
    StringSet<TSeqV, Dependent<> > depSetV;
    bool allSameLength = true;
    auto ptrTask = tasks[0];
    auto lenH = length(context(*ptrTask).seqHBlocks[column(*ptrTask)]);
    auto lenV = length(context(*ptrTask).seqVBlocks[row(*ptrTask)]);

    for (auto ptrTask : tasks)
    {
        appendValue(depSetH, context(*ptrTask).seqHBlocks[column(*ptrTask)]);
        appendValue(depSetV, context(*ptrTask).seqVBlocks[row(*ptrTask)]);
        if (lenH != length(context(*ptrTask).seqHBlocks[column(*ptrTask)]) ||
            lenV != length(context(*ptrTask).seqVBlocks[row(*ptrTask)]))
        {
            allSameLength = false;
        }
    }

    // Dummy trace set.
    StringSet<String<Nothing> > trace;  // We need to instantiate it, but it will not be used.

    // We can compute with one simd score, but might collect them here.
    auto const & scoringScheme = context(*tasks[0]).dpSettings.simdScoringScheme;

    // Preapare and run alingment.
    String<TSimdVec, Alloc<OverAligned> > stringSimdH;
    String<TSimdVec, Alloc<OverAligned> > stringSimdV;

    if (allSameLength)
    {
        using TScoutState = DPScoutState_<DPTiled<TSimdBufferH, Default, SimdAlignEqualLength>>;
        TScoutState scoutState(bufferH, bufferV);
        _prepareSimdAlignment(stringSimdH, stringSimdV, depSetH, depSetV, scoutState);

        using TScoutSpec = typename ScoutSpecForAlignmentAlgorithm_<typename TExecTraits::TAlgorithmType, TScoutState>::Type;
        using TDPScout = DPScout_<TDPCell, TScoutSpec>;

        TDPScout dpScout(scoutState);
        // We rather want to set
        computeTile(cache, dpScout, stringSimdH, stringSimdV, scoringScheme, context(*tasks[0]).dpSettings);

        // Now we need to run the scout check for all tasks.

        // We want to get the state here from the scout.
        for (size_t pos = 0; pos < length(tasks); ++pos)
        {
            auto & task = *tasks[pos];
            if (AlgorithmProperty<typename TExecTraits::TAlgorithmType>::isTrackingEnabled(task))
            {
                // TODO(rrahn): Implement the interface.
                // TODO(rrahn): Make it a member function of a policy so that we don't have to implement the specifics here
                _setSimdLane(dpScout, pos);
                auto & taskContext = context(task);
                updateMax(intermediate(dpLocal, taskContext.alignmentId),
                          {maxScoreAt(dpScout) + offset[pos], 0u},
                          column(task),
                          row(task));
            }
        }
    }
    else
    {
        using TDPSettings = std::decay_t<decltype(context(*tasks[0]).dpSettings)>;
        using TDPTraits = typename TDPSettings::TTraits;

        using TDPProfile = DPProfile_<typename TDPTraits::TAlgorithmType,
                                      typename TDPTraits::TGapType,
                                      typename TDPTraits::TTracebackType,
                                      Parallel>;

        using TSimdScoutTrait = SimdAlignVariableLengthTraits<TSimdVec,
                                                              decltype(depSetH),
                                                              decltype(depSetV),
                                                              TDPProfile>;
        using TScoutState = DPScoutState_<DPTiled<TSimdBufferH, Default, SimdAlignVariableLength<TSimdScoutTrait>>>;

        String<size_t> lengthsH;
        String<size_t> lengthsV;

        TScoutState scoutState(bufferH, bufferV);
        _prepareSimdAlignment(stringSimdH, stringSimdV, depSetH, depSetV, lengthsH, lengthsV, scoutState);

        using TScoutSpec = typename ScoutSpecForAlignmentAlgorithm_<typename TExecTraits::TAlgorithmType, TScoutState>::Type;
        using TDPScout = DPScout_<TDPCell, TScoutSpec>;

        TDPScout dpScout(scoutState);
        computeTile(cache, dpScout, stringSimdH, stringSimdV, scoringScheme, context(*tasks[0]).dpSettings);
        // We want to get the state here from the scout.
        for (size_t pos = 0; pos < length(tasks); ++pos)
        {
            auto & task = *tasks[pos];
            if (AlgorithmProperty<typename TExecTraits::TAlgorithmType>::isTrackingEnabled(task))
            {
                // TODO(rrahn): Implement the interface.
                // TODO(rrahn): Make it a member function of a policy so that we don't have to implement the specifics here
                _setSimdLane(dpScout, pos);
                auto & taskContext = context(task);
                updateMax(intermediate(dpLocal, taskContext.alignmentId),
                          {maxScoreAt(dpScout) + offset[pos], 0u},
                          column(task),
                          row(task));
            }
        }
    }
}
#endif // SEQAN_SIMD_ENABLED
}  // namespace impl
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_UTIL_H_
