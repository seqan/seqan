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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_H_

namespace seqan
{
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Context used per task. Access information like the infixes of the sequences for this block and other.
template <typename TSeqHBlocks,
          typename TSeqVBlocks,
          typename TTileBuffer,
          typename TDPSettings,
          typename TEvent = WavefrontTaskEvent>
struct WavefrontAlignmentContext
{
    size_t              alignmentId{0};
    TSeqHBlocks const & seqHBlocks;
    TSeqVBlocks const & seqVBlocks;
    TTileBuffer       & tileBuffer;
    TDPSettings const & dpSettings;
    TEvent            * ptrEvent{nullptr};

    //NOTE(rrahn) Bug in g++-4.9 prevents us from using as aggregate type.
    WavefrontAlignmentContext(size_t const _alignmentId,
                              TSeqHBlocks const & _seqHBlocks,
                              TSeqVBlocks const & _seqVBlocks,
                              TTileBuffer       & _tileBuffer,
                              TDPSettings const & _dpSettings) :
        alignmentId(_alignmentId),
        seqHBlocks(_seqHBlocks),
        seqVBlocks(_seqVBlocks),
        tileBuffer(_tileBuffer),
        dpSettings(_dpSettings)
    {}
};

// The abstract task that is executed as separat alignment instance.
template <typename TAlignmentContext>
class WavefrontTask
{
public:

    using TContext  = TAlignmentContext;

    TContext &                     context;
    std::array<WavefrontTask*, 2>  successor{{nullptr, nullptr}};
    size_t                         col{0};
    size_t                         row{0};
    std::atomic<size_t>            refCount{0};
    bool                           lastTileH{false};
    bool                           lastTileV{false};


    //-------------------------------------------------------------------------
    // Constructor
    WavefrontTask() = delete;

    WavefrontTask(TContext & context, std::array<WavefrontTask*, 2> successor,
                  size_t const col,
                  size_t const row,
                  size_t const refCount,
                  bool const lastTileH,
                  bool const lastTileV) :
        context(context),
        successor(std::move(successor)),
        col(col), row(row),
        refCount(refCount),
        lastTileH(lastTileH), lastTileV(lastTileV)
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

    using TSeqHBlocks       = typename std::decay<decltype(std::declval<TaskContext_>().seqHBlocks)>::type;
    using TSeqVBlocks       = typename std::decay<decltype(std::declval<TaskContext_>().seqVBlocks)>::type;
    using TWavefrontBuffer  = typename std::decay<decltype(std::declval<TaskContext_>().tileBuffer)>::type;
    using TDPSettings       = typename std::decay<decltype(std::declval<TaskContext_>().dpSettings)>::type;

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
    me.refCount.store(count, std::memory_order_relaxed);
}

template <typename ...TArgs>
inline unsigned
decrementRefCount(WavefrontTask<TArgs...> & me)
{
    return --me.refCount;
}

template <typename ...TArgs>
inline unsigned
incrementRefCount(WavefrontTask<TArgs...> & me)
{
    return ++me.refCount;
}

template <typename TTask>
inline auto
column(TTask const & task) -> decltype(task.col)
{
    return task.col;
}

template <typename TTask>
inline auto
row(TTask const & task) -> decltype(task.row)
{
    return task.row;
}

template <typename TTask>
inline bool
inLastColumn(TTask const & task)
{
    return task.lastTileH;
}

template <typename TTask>
inline bool
inLastRow(TTask const & task)
{
    return task.lastTileV;
}

template <typename TTask>
inline bool
isLastTask(TTask const & task)
{
    return inLastColumn(task) && inLastRow(task);
}

template <typename TTask>
inline auto
successor(TTask & task) -> std::add_lvalue_reference_t<decltype(task.successor)>
{
    return task.successor;
}

template <typename TTask>
inline auto
successor(TTask const & task) -> std::add_lvalue_reference_t<std::add_const_t<decltype(task.successor)>>
{
    return task.successor;
}

template <typename TTask>
inline auto
context(TTask & task) -> std::add_lvalue_reference_t<decltype(task.context)>
{
    return task.context;
}

template <typename TTask>
inline auto
context(TTask const & task) -> std::add_lvalue_reference_t<std::add_const_t<decltype(task.context)>>
{
    return task.context;
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

    auto & taskContext = context(task);
    // Load the cache from the local data.
    auto & dpCache = cache(dpLocal, taskContext.alignmentId);
    auto & buffer = taskContext.tileBuffer;

    // Capture the buffer.
    typename TExecTraits::TDPScoutState scoutState(buffer.horizontalBuffer[column(task)],
                                                   buffer.verticalBuffer[row(task)]);  // Task local

    typename TExecTraits::TDPScout scout(scoutState);

    impl::computeTile(dpCache, scout,
                      taskContext.seqHBlocks[column(task)],
                      taskContext.seqVBlocks[row(task)],
                      taskContext.dpSettings.scoringScheme,
                      taskContext.dpSettings);
    // We want to get the state here from the scout.
    if(impl::AlgorithmProperty<typename TExecTraits::TAlgorithmType>::isTrackingEnabled(task))
    {
        // TODO(rrahn): Implement the interface.
        // TODO(rrahn): Make it a member function of a policy so that we don't have to implement the specifics here
        updateMax(intermediate(dpLocal, taskContext.alignmentId),
                  {maxScore(scout), maxHostPosition(scout)},
                  column(task),
                  row(task));
    }
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
    using TTask = typename std::remove_pointer<typename Value<TTasks>::Type>::type;
    using TExecTraits = SimdTaskExecutionTraits<typename TTask::TContext>;

    auto offset = impl::computeOffset(tasks, TExecTraits{});
    // Has to be adapted to take the correct buffer from the corresponding task.
    auto simdBufferH = impl::gatherSimdBuffer(tasks,
                                              [] (auto & task)
                                              {
                                                  return &context(task).tileBuffer.horizontalBuffer[column(task)];
                                              },
                                              offset,
                                              TExecTraits{});
    auto simdBufferV = impl::gatherSimdBuffer(tasks,
                                              [] (auto & task)
                                              {
                                                  return &context(task).tileBuffer.verticalBuffer[row(task)];
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
                                return &context(task).tileBuffer.horizontalBuffer[column(task)];
                            },
                            offset,
                            TExecTraits{});
    impl::scatterSimdBuffer(tasks,
                            simdBufferV,
                            [](auto & task)
                            {
                                return &context(task).tileBuffer.verticalBuffer[row(task)];
                            },
                            offset,
                            TExecTraits{});
}
#endif  // SEQAN_SIMD_ENABLED

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_H_
