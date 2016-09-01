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

template <typename TSimdVec, typename TTasks, typename TBufferSet, typename TFunc>
inline auto
gatherSimdBuffer(TTasks const & tasks,
                 TBufferSet const & buffer,
                 TFunc && getBlockId)
{
    using TBuffer     = typename Value<TBufferSet>::Type;
    using TBuffValue  = typename Value<TBuffer>::Type;
    using TDPCell     = typename Value<TBuffValue, 1>::Type;
    using TDPCellSpec = typename Spec<TDPCell>::Type;

    // TODO(rrahn): Pass simd type from outside.
    using TDPCellSimd = DPCell_<TSimdVec, TDPCellSpec>;
    using TTraceValue = typename TraceBitMap_<TSimdVec>::Type;
    using TBufferValue = Pair<TDPCellSimd, TTraceValue>;

    // Check for valid simd length.
    SEQAN_ASSERT_EQ(LENGTH<TSimdVec>::VALUE, length(tasks));

    String<TBufferValue, Alloc<OverAligned> > simdSet;

    auto maxLength = length(buffer[getBlockId(tasks[0])]);
//    std::for_each(begin(tasks) + 1, end(tasks),
//                  [&](auto& task)
//                  {
//                      auto len = length(buffer[getBlockId(task)]);
//                      if (len != maxLength)
//                          allSameLength = false;
//                      maxLength = (len > maxLength) ? len : maxLength;
//                  });

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

template <typename TDPCell, typename TTraceValue, typename TScoreMat, typename TTraceMat,
          typename TStateThreadContext,
          typename TSimdBufferH,
          typename TSimdBufferV,
          typename TSetSeqH,
          typename TSetSeqV,
          typename TScore,
          typename TBand,
          typename TDPConfig>
inline void
computeSimdBatch(DPContext<TDPCell, TTraceValue, TScoreMat, TTraceMat> & dpContext,
                 TStateThreadContext stateThreadContext,
                 TSimdBufferH & bufferH,
                 TSimdBufferV & bufferV,
                 TSetSeqH const & seqH,
                 TSetSeqV const & seqV,
                 TScore const & scoringScheme,
                 TBand const & band,
                 TDPConfig const & config)
{
    using TSeqH = typename Value<TSetSeqH>::Type;
    using TSeqV = typename Value<TSetSeqV>::Type;
    using TSimdVec = typename Value<TDPCell>::Type;

    // Prepare sequence set.
    StringSet<TSeqH, Dependent<> > depSetH;
    StringSet<TSeqV, Dependent<> > depSetV;
    bool allSameLength = true;
    auto lenH = length(seqH[stateThreadContext.mTask[0]->_col]);
    auto lenV = length(seqV[stateThreadContext.mTask[0]->_row]);
    for (auto& task : stateThreadContext.mTask)
    {
//        if (stateThreadContext.mTask[0]->_col == 88 && stateThreadContext.mTask[0]->_row == 8)
//        {
//            std::cout << "SeqH " << seqH[task->_col] << "\n";
//            std::cout << "SeqV " << seqV[task->_row] << "\n";
//        }
        appendValue(depSetH, seqH[task->_col]);
        appendValue(depSetV, seqV[task->_row]);
        if (lenH != length(seqH[task->_col]) || lenV != length(seqV[task->_row]))
            allSameLength = false;
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
        using TScoutState = DPScoutState_<DPTiled<TSimdBufferH, TStateThreadContext, SimdAlignEqualLength> >;

        TScoutState scoutState(bufferH, bufferV, std::move(stateThreadContext));
        _prepareSimdAlignment(stringSimdH, stringSimdV, depSetH, depSetV, scoutState);
        _computeAlignment(dpContext, trace, scoutState, stringSimdH, stringSimdV, simdScoringScheme,
                          band, TDPConfig(), false, false);
    }
    else
    {
        using TSimdScoutTrait = SimdAlignVariableLengthTraits<TSimdVec, decltype(depSetH), decltype(depSetV)>;
        using TScoutState = DPScoutState_<DPTiled<TSimdBufferH, TStateThreadContext, SimdAlignVariableLength<TSimdScoutTrait> > >;

        TScoutState scoutState(bufferH, bufferV, std::move(stateThreadContext));
        _prepareSimdAlignment(stringSimdH, stringSimdV, depSetH, depSetV, scoutState);

        scoutState.dimV = length(stringSimdV);
        scoutState.isLocalAlignment = IsLocalAlignment_<TDPConfig>::VALUE;
        scoutState.right = IsFreeEndGap_<TDPConfig, DPLastColumn>::VALUE;
        scoutState.bottom = IsFreeEndGap_<TDPConfig, DPLastRow>::VALUE;

        _computeAlignment(dpContext, trace, scoutState, stringSimdH, stringSimdV, simdScoringScheme,
                          band, TDPConfig(), false, false);
    }
//
//    if (stateThreadContext.mTask[0]->_col == 88 && stateThreadContext.mTask[0]->_row == 8)
//    {
//        for (unsigned j = 0; j < length(stateThreadContext.mTask); ++j)
//        {
//            std::cout << "\n\nSeqH ";
//            for (unsigned i = 0; i < length(stringSimdH); ++i)
//            {
//                std::cout << static_cast<Dna>(stringSimdH[i][j]);
//            }
//            std::cout << "\nSeqV ";
//            for (unsigned i = 0; i < length(stringSimdV); ++i)
//            {
//                std::cout << static_cast<Dna>(stringSimdV[i][j]);
//            }
//            std::cout << "\n";
//        }
//    }
}

}  // namespace impl

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_BASE_SIMD_H_
