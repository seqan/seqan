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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_STD_TASK_CONTEXT_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_STD_TASK_CONTEXT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TScheduler, typename TThreadLocalStore>
struct WavefrontExecutorStd
{
    // Shared data in parallel context.
    TScheduler *                  mTaskSchedulerPtr{nullptr};
    TThreadLocalStore *           mThreadLocalPtr{nullptr};
};

//template <typename TGlobalExecutor>
//struct WavefrontExecutorEvent
//{
//    TGlobalExecutor             & mExecutor;
//    WavefrontAlignmentTaskEvent mThreadEvent;
//
//    WavefrontExecutorEvent(TGlobalExecutor & exec) : mExecutor(exec)
//    {}
//};
//
//template <typename TGlobalExecutor>
//struct WavefrontExecutorAtomic
//{
//    TGlobalExecutor &   mExecutor;
//    std::atomic_flag    mFlag;
//
//    WavefrontExecutorAtomic(TGlobalExecutor & exec) : mExecutor(exec)
//    {
//        mFlag.test_and_set(std::memory_order_relaxed);
//    }
//};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename ...TArgs,
          typename TTaskExecutor>
inline void
spawn(WavefrontExecutorStd<TArgs...> & executor,
      TTaskExecutor && taskExec)
{
    SEQAN_ASSERT(executor.mTaskSchedulerPtr != nullptr);
    scheduleTask(*executor.mTaskSchedulerPtr, std::forward<TTaskExecutor>(taskExec));
}

template <typename ...TArgs>
inline auto&
local(WavefrontExecutorStd<TArgs...> & executor)
{
    SEQAN_ASSERT(executor.mThreadLocalPtr != nullptr);
    return local(*executor.mThreadLocalPtr);
}

//template <typename ...TArgs>
//inline void
//notify(WavefrontExecutorStd<TArgs...> & executor)
//{
//    SEQAN_ASSERT(executor.mThreadEventPtr != nullptr);
//    notify(*executor.mThreadEventPtr);
//}
//
//template <typename ...TArgs>
//inline void
//wait(WavefrontExecutorStd<TArgs...> & executor)
//{
//    SEQAN_ASSERT(executor.mThreadEventPtr != nullptr);
//    wait(*executor.mThreadEventPtr);
//}

//template <typename ...TArgs,
//typename TTaskExecutor>
//inline void
//spawn(WavefrontExecutorEvent<TArgs...> & executor,
//      TTaskExecutor && taskExec)
//{
//    spawn(executor.mExecutor, std::forward<TTaskExecutor>(taskExec));
//}
//
//template <typename ...TArgs>
//inline auto&
//local(WavefrontExecutorEvent<TArgs...> & executor)
//{
//    return local(executor.mExecutor);
//}
//
//template <typename ...TArgs>
//inline void
//notify(WavefrontExecutorEvent<TArgs...> & executor)
//{
//    notify(executor.mThreadEvent);
//}
//
//template <typename ...TArgs>
//inline void
//wait(WavefrontExecutorEvent<TArgs...> & executor)
//{
//    wait(executor.mThreadEvent);
//}

//template <typename ...TArgs,
//typename TTaskExecutor>
//inline void
//spawn(WavefrontExecutorAtomic<TArgs...> & executor,
//      TTaskExecutor && taskExec)
//{
//    spawn(executor.mExecutor, std::forward<TTaskExecutor>(taskExec));
//}
//
//template <typename ...TArgs>
//inline auto&
//local(WavefrontExecutorAtomic<TArgs...> & executor)
//{
//    return local(executor.mExecutor);
//}

//template <typename ...TArgs>
//inline void
//notify(WavefrontExecutorAtomic<TArgs...> & executor)
//{
//    executor.mFlag.clear(std::memory_order_release);
//}
//
//template <typename ...TArgs>
//inline void
//wait(WavefrontExecutorAtomic<TArgs...> & executor)
//{
//    SpinDelay delay;
//    while (executor.mFlag.test_and_set(std::memory_order_acquire))
//        waitFor(delay);
//}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_STD_TASK_CONTEXT_H_
