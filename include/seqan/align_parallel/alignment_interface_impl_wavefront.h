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

template <typename TSeqH, typename TSeqV, typename TSettings>
class AsyncWaveAlignExecutor
{
public:

    using TAlignmentTask = WavefrontAlignmentTask<TSeqH, TSeqV, TSettings>;
    using TThreadLocal   = typename WavefrontAlignmentTaskConfig<TSettings>::TThreadLocal;
    using TExecutor      = WavefrontExecutorStd<WavefrontTaskScheduler, EnumerableThreadLocal<TThreadLocal, Limit>>;

    TSettings                                   _settings;
    // Initialize the alignment scheduler.
    WavefrontAlignmentScheduler                 _alignScheduler;//{parallelAlignments(execPolicy), numThreads(execPolicy)};

    EnumerableThreadLocal<TThreadLocal, Limit>  _threadLocalStorage;//{numThreads(execPolicy), TThreadLocal{parallelAlignments(execPolicy)}};
    TExecutor                                   _executor{}; //{&taskScheduler(_alignScheduler), &_threadLocalStorage};
    unsigned                                    _alignCounter{0};
    unsigned                                    _blockSize{};

    template <typename TSpec>
    AsyncWaveAlignExecutor(TSettings settings,
                           ExecutionPolicy<WavefrontAlignment<TSpec>, Serial> const & execPolicy) :
        _settings(std::move(settings)),
        _alignScheduler(parallelAlignments(execPolicy), numThreads(execPolicy)),
        _threadLocalStorage(numThreads(execPolicy), TThreadLocal{parallelAlignments(execPolicy)}),
        _blockSize(blockSize(execPolicy))
    {
        _executor.mTaskSchedulerPtr = &taskScheduler(_alignScheduler);
        _executor.mThreadLocalPtr   = &_threadLocalStorage;
    }
};

template <typename ...TArgs,
          typename TSeqH,
          typename TSeqV,
          typename TCallable>
inline void
submit(AsyncWaveAlignExecutor<TArgs...> & me,
       TSeqH const & seqH,
       TSeqV const & seqV,
       TCallable && callback)
{
    using TAlignTask = typename AsyncWaveAlignExecutor<TArgs...>::TAlignmentTask;

    std::function<std::function<void(uint16_t)>(TAlignTask)> aiCaller =
    [&me, &callback](auto && func) mutable
    {
        return [&me, &callback, func = std::move(func)](uint16_t id) mutable
        {
            func(id, me._executor, std::forward<TCallable>(callback));
        };
    };
    scheduleTask(me._alignScheduler, aiCaller(TAlignTask{me._alignCounter++, seqH, seqV, me._settings, me._blockSize}));
}

template <typename ...TArgs>
inline void
wait(AsyncWaveAlignExecutor<TArgs...> & me)
{
    notify(me._alignScheduler);
    wait(me._alignScheduler);
}

#ifdef SEQAN_SIMD_ENABLED
template <typename TSeqH, typename TSeqV, typename TSettings, typename TWaveSpec>
class AsyncWaveAlignExecutorSimd
{
public:

    // Translate dp settings into simd settings.
    using TSimdSettings = SimdDPSettings<TSettings, TWaveSpec>;

    using TAlignmentTask = WavefrontAlignmentTask<TSeqH, TSeqV, TSimdSettings, WavefrontAlignmentSimdTaskConfig<TSimdSettings>>;
    using TWavefrontTask = WavefrontTask<typename TAlignmentTask::TTaskContext>;
    using TSimdTaskQueue = WavefrontSimdDPTasks<TWavefrontTask, LENGTH<typename TSimdSettings::TScoreValueSimd>::VALUE>;

    using TThreadLocal  = typename WavefrontAlignmentSimdTaskConfig<TSimdSettings>::TThreadLocal;
    using TExecutor     = WavefrontExecutorStd<WavefrontTaskScheduler, EnumerableThreadLocal<TThreadLocal, Limit>>;

    TSimdSettings                               _settings;
    // Initialize the alignment scheduler.
    WavefrontAlignmentScheduler                 _alignScheduler; //{parallelAlignments(execPolicy), numThreads(execPolicy)};

    //    EnumerableThreadLocal<TThreadLocal> tls{TThreadLocal{parallelAlignments(execPolicy)}};
    EnumerableThreadLocal<TThreadLocal, Limit>  _threadLocalStorage; //{numThreads(execPolicy), TThreadLocal{parallelAlignments(execPolicy)}};
    TExecutor                                   _executor{}; //{&taskScheduler(_alignScheduler), &_threadLocalStorage};
    TSimdTaskQueue                              _simdTaskQueue{};
    unsigned                                    _alignCounter{0};
    unsigned                                    _blockSize{};

    template <typename TSpec>
    AsyncWaveAlignExecutorSimd(TSettings const & settings,
                               ExecutionPolicy<WavefrontAlignment<TSpec>, Vectorial> const & execPolicy) :
        _settings(settings.mScoringScheme),
        _alignScheduler(parallelAlignments(execPolicy), numThreads(execPolicy)),
        _threadLocalStorage(numThreads(execPolicy), TThreadLocal{parallelAlignments(execPolicy)}),
        _blockSize(blockSize(execPolicy))
    {
        _executor.mTaskSchedulerPtr = &taskScheduler(_alignScheduler);
        _executor.mThreadLocalPtr   = &_threadLocalStorage;
    }
};

template <typename ...TArgs,
          typename TSeqH,
          typename TSeqV,
          typename TCallable>
inline void
submit(AsyncWaveAlignExecutorSimd<TArgs...> & me,
       TSeqH const & seqH,
       TSeqV const & seqV,
       TCallable && callback)
{
    using TAlignTask = typename AsyncWaveAlignExecutorSimd<TArgs...>::TAlignmentTask;

    // Continuator for calling the alignment instance functor.
    std::function<std::function<void(uint16_t)>(TAlignTask)> aiCaller =
    [&me, &callback](auto && func) mutable
    {
        return [&me, &callback, func = std::move(func)](uint16_t id) mutable
        {
            func(id, me._executor, me._simdTaskQueue, std::forward<TCallable>(callback));
        };
    };
    scheduleTask(me._alignScheduler, aiCaller(TAlignTask{me._alignCounter++, seqH, seqV, me._settings, me._blockSize}));
}

template <typename ...TArgs>
inline void
wait(AsyncWaveAlignExecutorSimd<TArgs...> & me)
{
    notify(me._alignScheduler);
    wait2(me._alignScheduler, me._simdTaskQueue);
}
#endif // SEQAN_SIMD_ENABLED

template <typename TSpec, typename TSimdSpec,
          typename TSetH,
          typename TSetV,
          typename TSettings,
          typename TCallable>
inline void
alignExecBatch(ExecutionPolicy<WavefrontAlignment<TSpec>, TSimdSpec> const & execPolicy,
               TSetH const & setH,
               TSetV const & setV,
               TSettings const & settings,
               TCallable && callback)
{
    using TSeqH = typename Value<TSetH const>::Type;
    using TSeqV = typename Value<TSetV const>::Type;

#ifdef SEQAN_SIMD_ENABLED
    using TExecutor = std::conditional_t<std::is_same<TSimdSpec, Vectorial>::value,
                                         AsyncWaveAlignExecutorSimd<TSeqH, TSeqV, TSettings, TSpec>,
                                         AsyncWaveAlignExecutor<TSeqH, TSeqV, TSettings>>;
#else
    using TExecutor = AsyncWaveAlignExecutor<TSeqH, TSeqV, TSettings>;
#endif
    TExecutor executor(settings, execPolicy);

    for(size_t i = 0u; i < length(setH); ++i)
    {
        submit(executor, setH[i], setV[i], std::forward<TCallable>(callback));
    }
    wait(executor);
}

}  // namespace impl
}  // namespace seqan
#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGNMENT_IMPL_WAVEFRONT_H_
