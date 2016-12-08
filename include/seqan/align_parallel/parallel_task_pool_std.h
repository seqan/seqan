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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_TASK_POOL_STD_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_TASK_POOL_STD_H_

namespace seqan
{
    
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

//template <typename TWorkQueue, typename TTaskContext>
//struct DPThread<TWorkQueue, TTaskContext, Serial>
//{
//    using TTask = typename Value<TWorkQueue>::Type;
//    
//    TWorkQueue *  workQueuePtr = nullptr;
//    TTaskContext& mTaskContext;
//    uint8_t       mThreadId;
//    
//    inline void
//    operator()()
//    {
//        lockWriting(*workQueuePtr);
//        while (true)
//        {
//            TTask task = nullptr;
//            
//            if (!popFront(task, *workQueuePtr))
//                return;
//            
//            SEQAN_ASSERT(task != nullptr);
//            task->template execute(*workQueuePtr, Nothing(), mThreadId);
//        }
//    }
//};
//
//template <typename TWorkQueue, typename TTaskContext>
//struct DPThread<TWorkQueue, TTaskContext, Vectorial>
//{
//    using TTask = typename Value<TWorkQueue>::Type;
//    
//    TWorkQueue *  workQueuePtr = nullptr;
//    TTaskContext& mTaskContext;
//    uint8_t       mThreadId;
//    
//    inline void
//    operator()()
//    {
//        lockWriting(*workQueuePtr);
//        String<TTask> tasks;
//        while (true)
//        {
//            TTask task = nullptr;
//            clear(tasks);
//            {
//                std::lock_guard<decltype(mTaskContext.mLock)> scopedLock(mTaskContext.mLock);
//                if (!popFront(task, *workQueuePtr))
//                    return;
//                
//                if (length(*workQueuePtr) >= TTaskContext::VECTOR_SIZE - 1)
//                {
//                    for (unsigned i = 0; i < TTaskContext::VECTOR_SIZE - 1; ++i)
//                        appendValue(tasks, popFront(*workQueuePtr));
//                }
//            }
//            
//            SEQAN_ASSERT(task != nullptr);
//            task->template execute(*workQueuePtr, tasks, mThreadId);
//        }
//    }
//};

template <typename TTask, typename TVecSpec>
struct TaskPoolTrait<TTask, ExecutionPolicy<ThreadModelStd, TVecSpec>>
{
    using TTaskPoolType  = ConcurrentQueue<TTask>;
    using ThreadPoolType = ThreadPool<ThreadModelStd>;
};

template <typename TTask, typename TVecSpec, typename TTrait>
class TaskPool<TTask, ExecutionPolicy<ThreadModelStd, TVecSpec>, TTrait>
{
    using TTaskPool = typename TTrait::TPoolType;
    using TThreadPool = typename TTrait::ThreadPoolType;
    
    TThreadPool mThreadPool;
    TTaskPool   mTaskPool;
    
    TaskPool() : mThreadPool(/*worker*/)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
    
    
}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_TASK_POOL_STD_H_