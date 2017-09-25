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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_H_

namespace seqan
{
// ============================================================================
// Forwards
// ============================================================================

std::atomic_flag _debug_cout_flag = ATOMIC_FLAG_INIT;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSeqHBlocks,
          typename TSeqVBlocks,
          typename TTileBuffer,
          typename TDPSettings,
          typename TEvent = WavefrontTaskEvent>
struct WavefrontAlignmentContext
{
    size_t              mAlignmentId{0};
    TSeqHBlocks const & mSeqHBlocks;
    TSeqVBlocks const & mSeqVBlocks;
    TTileBuffer       & mTileBuffer;
    TDPSettings const & mDPSettings;
    TEvent            * mEventPtr{nullptr};
};

template <typename TAlignmentContext>
class WavefrontTask
{
public:

    using TContext  = TAlignmentContext;

    TContext &              mContext;

    std::array<WavefrontTask*, 2>  mSuccessor{nullptr, nullptr};
    size_t                  mCol{0};
    size_t                  mRow{0};
    std::atomic<size_t>     mRefCount{0};

    bool                    mLastTileH{false};
    bool                    mLastTileV{false};


    //-------------------------------------------------------------------------
    // Constructor
    WavefrontTask() = delete;

    WavefrontTask(TContext & context, std::array<WavefrontTask*, 2> successor,
                  size_t const col,
                  size_t const row,
                  size_t const refCount,
                  bool const lastTileH,
                  bool const lastTileV) :
                mContext(context),
                mSuccessor(std::move(successor)),
                mCol(col), mRow(row),
            	mRefCount(refCount),
                mLastTileH(lastTileH), mLastTileV(lastTileV)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TContext>
struct TaskExecutionTraits;

template <typename ...TArgs>
struct TaskExecutionTraits<WavefrontAlignmentContext<TArgs...>>
{
    using TaskContext_      = WavefrontAlignmentContext<TArgs...>;

    using TSeqHBlocks       = typename std::decay<decltype(std::declval<TaskContext_>().mSeqHBlocks)>::type;
    using TSeqVBlocks       = typename std::decay<decltype(std::declval<TaskContext_>().mSeqVBlocks)>::type;
    using TWavefrontBuffer  = typename std::decay<decltype(std::declval<TaskContext_>().mTileBuffer)>::type;
    using TDPSettings       = typename std::decay<decltype(std::declval<TaskContext_>().mDPSettings)>::type;

    using TTileBuffer       = typename std::decay<decltype(std::declval<TWavefrontBuffer>().horizontalBuffer[0])>::type;
    using TDPScoutState     = DPScoutState_<DPTiled<TTileBuffer>>;

    // Sequence types.

    using TSeqH = typename Value<TSeqHBlocks>::Type;
    using TSeqV = typename Value<TSeqVBlocks>::Type;

    // DPTrait type forwarding.
    using TDPTraits         = typename TDPSettings::TTraits;
    using TScoreValue       = typename Value<typename TDPSettings::TScoringScheme>::Type;
    using TAlgorithmType    = typename TDPTraits::TAlgorithmType;
    using TTracebackType    = typename TDPTraits::TTracebackType;
    using TGapType          = typename TDPTraits::TGapType;

    // Wavefront Alignment Context.
    using TDPCell           = DPCell_<TScoreValue, TGapType>;

    using TScoutSpec        = typename ScoutSpecForAlignmentAlgorithm_<TAlgorithmType, TDPScoutState>::Type;
    using TDPScout          = DPScout_<TDPCell, TScoutSpec>;
};


template <typename TWavefrontAlignmentContextConcept>
struct SimdTaskExecutionTraits : public TaskExecutionTraits<TWavefrontAlignmentContextConcept>
{
    using TBase = TaskExecutionTraits<TWavefrontAlignmentContextConcept>;

    using TScoreValue = typename TBase::TDPSettings::TScoreValueSimd;
    using TDPCell     = DPCell_<TScoreValue, typename TBase::TGapType>;
    using TTraceValue = typename TraceBitMap_<TScoreValue>::Type;
    using TBufferValue = Pair<TDPCell, TTraceValue>;
};

// ============================================================================
// Functions
// ============================================================================

template <typename ...TArgs>
inline void
setRefCount(WavefrontTask<TArgs...> & me, size_t const count)
{
    me.mRefCount.store(count, std::memory_order_relaxed);
}

template <typename ...TArgs>
inline unsigned
decrementRefCount(WavefrontTask<TArgs...> & me)
{
    return --me.mRefCount;
}

template <typename ...TArgs>
inline unsigned
incrementRefCount(WavefrontTask<TArgs...> & me)
{
    return ++me.mRefCount;
}

template <typename TTask>
inline auto
column(TTask const & task)
{
    return task.mCol;
}

template <typename TTask>
inline auto
row(TTask const & task)
{
    return task.mRow;
}

template <typename TTask>
inline bool
inLastColumn(TTask const & task)
{
    return task.mLastTileH;
}

template <typename TTask>
inline auto
inLastRow(TTask const & task)
{
    return task.mLastTileV;
}

template <typename TTask>
inline bool
isLastTask(TTask const & task)
{
    return inLastColumn(task) && inLastRow(task);
}

template <typename TTask>
inline auto &
successor(TTask & task)
{
    return task.mSuccessor;
}

template <typename TTask>
inline auto const &
successor(TTask const & task)
{
    return task.mSuccessor;
}

template <typename TTask>
inline auto &
context(TTask & task)
{
    return task.mContext;
}

template <typename TTask>
inline auto const &
context(TTask const & task)
{
    return task.mContext;
}

template <typename TAlgorithm, typename TTask>
inline bool
isTrackTile(TTask const & task)
{
    return isLastColumn(task) && isLastRow(task);
}

template <typename TTask>
inline bool
isTrackTile(TTask const & task)
{
    return isLastColumn(task) && isLastRow(task);
}

template <typename TTask, typename TDPLocalData>
inline void
executeScalar(TTask & task, TDPLocalData & dpLocal)
{
    using TExecTraits = TaskExecutionTraits<typename TTask::TContext>;

    auto& taskContext = context(task);
    // Load the cache from the local data.
    auto & dpCache = cache(dpLocal, taskContext.mAlignmentId);
    auto & buffer = taskContext.mTileBuffer;

    // Capture the buffer.
    typename TExecTraits::TDPScoutState scoutState(buffer.horizontalBuffer[column(task)],
                                                   buffer.verticalBuffer[row(task)]);  // Task local

    typename TExecTraits::TDPScout scout(scoutState);

    // DEBUG: Remove!
//        auto bufHBegin = _taskContext.getTileBuffer().horizontalBuffer[_col];
//        auto bufVBegin = _taskContext.getTileBuffer().verticalBuffer[_row];

    impl::computeTile(dpCache, scout,
                      taskContext.mSeqHBlocks[column(task)],
                      taskContext.mSeqVBlocks[row(task)],
                      taskContext.mDPSettings.mScoringScheme,
                      taskContext.mDPSettings);

    // We want to get the state here from the scout.
    if(impl::AlgorithmProperty<typename TExecTraits::TAlgorithmType>::isTrackingEnabled(task))
    {
        // TODO(rrahn): Implement the interface.
        // TODO(rrahn): Make it a member function of a policy so that we don't have to implement the specifics here
        updateMax(intermediate(dpLocal, taskContext.mAlignmentId),
                  {maxScore(scout), maxHostPosition(scout)},
                  column(task),
                  row(task));
    }


    // DEBUG: Remove!
//        _taskContext.getDebugBuffer().matrix[_col][_row] = {bufHBegin, bufVBegin,
//                                                            _taskContext.getTileBuffer().horizontalBuffer[_col],
//                                                            _taskContext.getTileBuffer().verticalBuffer[_row],
//                                                            _col, _row, false};
    // TODO(rrahn): Add traceback later.
//    if (IsTracebackEnabled_<typename TTaskConfig::TDPConfig>::VALUE)
//    {
//        swap(getDpTraceMatrix(dpContext), pTls.mLocalTraceStore.localScalarTraceMatrix());
//        _taskContext.getTraceProxy().insert(std::make_pair(_col, _row),
//                                            TTraceProxyValue{&pTls.mLocalTraceStore,
//                                                             {0, static_cast<uint16_t>(length(pTls.mLocalTraceStore.mScalarTraceVec) - 1)},
//                                                              0});
//    }
//    #ifdef DP_ALIGN_STATS
//        ++serialCounter;
//    #endif
}

template <typename TBuffer>
inline void
printSimdBuffer(TBuffer const & buffer, size_t const l)
{
    for (auto simdHolder : buffer)
    {
        std::cout << "<";
        unsigned i = 0;
        for (; i < l - 1; ++i)
        {
            std::cout << simdHolder.i1._score[i] << ", ";
        }
        std::cout << simdHolder.i1._score[i] << ">\n";
    }
}

#ifdef SEQAN_SIMD_ENABLED
template <typename TTasks, typename TDPLocalData>
inline void
executeSimd(TTasks & tasks, TDPLocalData & dpLocal)
{
//    for (auto task : tasks)
//    {
//        executeScalar(*task, dpLocal);
//    }
    using TTask = typename std::remove_pointer<typename Value<TTasks>::Type>::type;
    using TExecTraits = SimdTaskExecutionTraits<typename TTask::TContext>;

//    auto& taskContext = context(task);
//    // Load the cache from the local data.
//    auto & dpCache = cache(dpLocal, taskContext.mAlignmentId);
//    auto & buffer = taskContext.mTileBuffer;
//
//    // Capture the buffer.
//    typename TExecTraits::TDPScoutState scoutState(buffer.horizontalBuffer[column(task)],
//                                                   buffer.verticalBuffer[row(task)]);  // Task local
//
//    typename TExecTraits::TDPScout scout(scoutState);
    //offset version:

    auto offset = impl::computeOffset(tasks, TExecTraits{});

    // Has to be adapted to take the correct buffer from the corresponding task.
    auto simdBufferH = impl::gatherSimdBuffer(tasks,
                                              [](auto& task)
                                              {
                                                  return &context(task).mTileBuffer.horizontalBuffer[column(task)];
                                              },
                                              offset,
                                              TExecTraits{});
    auto simdBufferV = impl::gatherSimdBuffer(tasks,
                                              [](auto& task)
                                              {
                                                  return &context(task).mTileBuffer.verticalBuffer[row(task)];
                                              },
                                              offset,
                                              TExecTraits{});

    // Does not really make sense.
    auto & cache = simdCache(dpLocal, 0);
    // Run alignment.
    impl::computeSimdBatch(cache, simdBufferH, simdBufferV, tasks, dpLocal, offset, TExecTraits{});

    // Write back into buffer.
    impl::scatterSimdBuffer(tasks,
                            simdBufferH,
                            [](auto & task)
                            {
                                return &context(task).mTileBuffer.horizontalBuffer[column(task)];
                            },
                            offset,
                            TExecTraits{});
    impl::scatterSimdBuffer(tasks,
                            simdBufferV,
                            [](auto & task)
                            {
                                return &context(task).mTileBuffer.verticalBuffer[row(task)];
                            },
                            offset,
                            TExecTraits{});
}
#endif  // namespace seqan

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_H_
