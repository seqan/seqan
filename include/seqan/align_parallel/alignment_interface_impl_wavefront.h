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
//  Implementation of the wavefront based alignment computation.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGNMENT_IMPL_WAVEFRONT_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGNMENT_IMPL_WAVEFRONT_H_

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

// ============================================================================
// Functions
// ============================================================================

template <typename TSpec,
          typename TSetH,
          typename TSetV,
          typename TSettings,
          typename TContinuator>
inline void
alignExecBatch(ExecutionPolicy<WavefrontAlignment<TSpec>, Serial> const & execPolicy,
               TSetH const & setH,
               TSetV const & setV,
               TSettings const & settings,
               TContinuator && callback)
{
    using TSeqH = typename Value<TSetH const>::Type;
    using TSeqV = typename Value<TSetV const>::Type;

    // Initialize instance of alignment scheduler.
    using TAlignTask    = WavefrontAlignmentTask<TSeqH, TSeqV, TSettings>;
    using TThreadLocal  = typename WavefrontAlignmentTaskConfig<TSettings>::TThreadLocal;
    using TExecutor     = WavefrontExecutorStd<WavefrontTaskScheduler, EnumerableThreadLocal<TThreadLocal, Limit>>;

    // Initialize the alignment scheduler.
    WavefrontAlignmentScheduler alignScheduler{parallelAlignments(execPolicy), numThreads(execPolicy)};

//    EnumerableThreadLocal<TThreadLocal> tls{TThreadLocal{parallelAlignments(execPolicy)}};
    EnumerableThreadLocal<TThreadLocal, Limit> tls{numThreads(execPolicy), TThreadLocal{parallelAlignments(execPolicy)}};
    TExecutor executor{&taskScheduler(alignScheduler), &tls};

    // Continuator for calling the alignment instance functor.
    std::function<std::function<void(uint16_t)>(TAlignTask)> aiCaller =
    [&executor, &callback](auto && func) mutable
    {
        return [&executor, &callback, func = std::move(func)](uint16_t id) mutable
        {
            func(id, executor, std::forward<TContinuator>(callback));
        };
    };

    // Iterate over pair of sequences.
//    auto zipCont = makeZipView(setH, setV);
    for(size_t i = 0u; i < length(setH); ++i)
    {
        scheduleTask(alignScheduler, aiCaller(TAlignTask{i, setH[i], setV[i], settings, blockSize(execPolicy)}));
    }
    notify(alignScheduler);
    wait(alignScheduler);
}

template <typename TSpec,
          typename TSetH,
          typename TSetV,
          typename TSettings,
          typename TContinuator>
inline void
alignExecBatch(ExecutionPolicy<WavefrontAlignment<TSpec>, Vectorial> const & execPolicy,
               TSetH const & setH,
               TSetV const & setV,
               TSettings const & settings,
               TContinuator && callback)
{
    using TSeqH = typename Value<TSetH const>::Type;
    using TSeqV = typename Value<TSetV const>::Type;

    // Translate dp settings into simd settings.
    using TSimdSettings = SimdDPSettings<TSettings>;
    TSimdSettings simdSettings{settings.mScoringScheme};

    // Initialize instance of alignment scheduler.
    using TAlignTask    = WavefrontAlignmentTask<TSeqH, TSeqV, TSimdSettings, WavefrontAlignmentSimdTaskConfig<TSimdSettings>>;
    using TThreadLocal  = typename WavefrontAlignmentSimdTaskConfig<TSimdSettings>::TThreadLocal;
    using TExecutor     = WavefrontExecutorStd<WavefrontTaskScheduler, EnumerableThreadLocal<TThreadLocal, Limit>>;

    // I need a resource here. And the resource should take the tasks.
    // But I don't have the task description here yet.
    using TWavefrontTask = WavefrontTask<typename TAlignTask::TTaskContext>;
    WavefrontSimdDPTasks<TWavefrontTask, LENGTH<typename TSimdSettings::TScoreValueSimd>::VALUE> simdTaskQueue;

    // Initialize the alignment scheduler.
    WavefrontAlignmentScheduler alignScheduler{parallelAlignments(execPolicy), numThreads(execPolicy)};

    //    EnumerableThreadLocal<TThreadLocal> tls{TThreadLocal{parallelAlignments(execPolicy)}};
    EnumerableThreadLocal<TThreadLocal, Limit> tls{numThreads(execPolicy), TThreadLocal{parallelAlignments(execPolicy)}};
    TExecutor executor{&taskScheduler(alignScheduler), &tls};

    // Continuator for calling the alignment instance functor.
    std::function<std::function<void(uint16_t)>(TAlignTask)> aiCaller =
    [&executor, &simdTaskQueue, &callback](auto && func) mutable
    {
        return [&executor, &simdTaskQueue, &callback, func = std::move(func)](uint16_t id) mutable
        {
                 func(id, executor, simdTaskQueue, std::forward<TContinuator>(callback));
        };
    };

    // Iterate over pair of sequences.
    //    auto zipCont = makeZipView(setH, setV);
    for(size_t i = 0u; i < length(setH); ++i)
    {
        scheduleTask(alignScheduler, aiCaller(TAlignTask{i, setH[i], setV[i], simdSettings, blockSize(execPolicy)}));
    }
    notify(alignScheduler);
    wait2(alignScheduler, simdTaskQueue);
}
}  // namespace impl
}  // namespace seqan
#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGNMENT_IMPL_WAVEFRONT_H_
