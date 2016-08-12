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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_TBB_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_TBB_H_

#include <tbb/task.h>
#include "tbb/enumerable_thread_specific.h"

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TTaskConfig, typename TThreadLocalStorage, typename TVecExecPolicy>
class DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb> :
    public DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb> >,
    public tbb::task
{
public:

    using TBase = DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb> >;
    using TSimdQueue = ConcurrentQueue<DPTaskImpl*>;

    // ============================================================================
    // Member variables.

    std::shared_ptr<TThreadLocalStorage> _ptrTls;
    TSimdQueue*                          mSimdQueuePtr = nullptr;
    bool                                 mIsLastTask = false;
//    std::atomic<uint8_t>                 mRefCount;

    // ============================================================================
    // Constructor.

    DPTaskImpl() = default;

    DPTaskImpl(size_t pCol, size_t pRow, TTaskConfig & pConfig,
               std::shared_ptr<TThreadLocalStorage> & pPtrTls) :
        TBase(pCol, pRow, pConfig, *this),
        _ptrTls(pPtrTls)
    {}

    // ============================================================================
    // Member functions.

    void setRefCount(unsigned n)
    {
        set_ref_count(n);
    }

    auto decrementRefCount()
    {
        return decrement_ref_count();
    }

    auto incrementRefCount()
    {
        return increment_ref_count();
    }

    auto& getLocalDpContext()
    {
        return _ptrTls->local();
    }

    template <typename TVec>
    inline SEQAN_FUNC_ENABLE_IF(Not<IsVectorExecutionPolicy<TVec> >, void)
    updateAndSpawn()
    {
        for (auto& task : TBase::successor)
        {
            if (task != nullptr && task->decrementRefCount() == 0)
            {
                spawn(*task);
            }
        }
    }

    template <typename TVec>
    inline SEQAN_FUNC_ENABLE_IF(IsVectorExecutionPolicy<TVec> , void)
    updateAndSpawn()
    {
        for (auto& task : TBase::successor)
        {
            if (task != nullptr && task->decrementRefCount() == 0)
            {
                task->mSimdQueuePtr = mSimdQueuePtr;  // Forwarding the simd queue.
                appendValue(*mSimdQueuePtr, task);
                spawn(*task);
            }
        }
    }

    template <typename TVec>
    inline SEQAN_FUNC_ENABLE_IF(Not<IsVectorExecutionPolicy<TVec> >, void)
    executeImpl()
    {
        SEQAN_ASSERT(_ptrTls != nullptr);
        TBase::runScalar(_ptrTls->local());
        updateAndSpawn<TVecExecPolicy>();
    }

    template <typename TVec>
    inline SEQAN_FUNC_ENABLE_IF(IsVectorExecutionPolicy<TVec> , void)
    executeImpl()
    {
        SEQAN_ASSERT(localDpContext != nullptr);
        SEQAN_ASSERT(mSimdQueuePtr != nullptr);

        String<DPTaskImpl*> tasks;
        {  // Acquire scoped lock.
            std::lock_guard<decltype(TBase::_taskContext.mLock)> scopedLock(TBase::_taskContext.mLock);

            if (empty(*mSimdQueuePtr))
            {
                return;
            }

            if (length(*mSimdQueuePtr) < TBase::_taskContext.mSimdLength)
            {
                appendValue(tasks, popFront(*mSimdQueuePtr));
            }
            else
            {
                SEQAN_ASSERT_GEQ(length(*mSimdQueuePtr), TBase::_taskContext.mSimdLength);
                for (auto i = 0; i < TBase::_taskContext.mSimdLength; ++i)
                    appendValue(tasks, popFront(*mSimdQueuePtr));
            }
        }  // Release scoped lock.

        SEQAN_ASSERT_GT(length(tasks), 0u);

        if (length(tasks) == 1)
        {
            front(tasks)->TBase::runScalar(_ptrTls->local());
        }
        else
        {
            SEQAN_ASSERT_EQ(length(tasks), TBase::_taskContext.mSimdLength);
            // Make to simd version!
            TBase::runSimd(tasks, _ptrTls->local());
        }
        // Update and spawn sucessors.
        for (auto& task : tasks)
        {
            task->template updateAndSpawn<TVecExecPolicy>();
            if (task->mIsLastTask)
            {
                cancel_group_execution();  // Notify, that all queued tasks can be canceled.
                {
                    std::lock_guard<std::mutex> lk(TBase::_taskContext.mLockEvent);
                    TBase::_taskContext.mReady = true;
                }
                TBase::_taskContext.mReadyEvent.notify_all();
            }
        }
    }

    // Override of tbb::task::execute()
    task *
    execute()
    {
        executeImpl<TVecExecPolicy>();
        return nullptr;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy>
struct IsDPTask<DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb> > : True
{};

// ============================================================================
// Functions
// ============================================================================

template <typename TTaskContext, typename TVecExecPolicy>
inline auto
createGraph(TTaskContext & context,
            TVecExecPolicy const & /*vecExecPolicy*/,
            ParallelExecutionPolicyTbb const & /*taskImplTag*/)
{
    using TThreadLocalStorage = tbb::enumerable_thread_specific<typename TTaskContext::TDPContext>;
    using TDagTask = DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb>;

    DPTaskGraph<TDagTask> graph;

    auto tls = std::make_shared<TThreadLocalStorage>();  // Copy into shared ptr per task. So first if all tasks are destroyed this resources is freed.
    resize(graph.get(), length(context.getSeqH()));
    for (int i = length(context.getSeqH()); --i >= 0;)
    {
        resize(graph[i], length(context.getSeqV()));
        for (int j = length(context.getSeqV()); --j >= 0;)
        {
            graph[i][j] = new (tbb::task::allocate_root()) TDagTask(i, j, context, tls);
            graph[i][j]->successor[0] = (i + 1 < length(context.getSeqH())) ? graph[i+1][j] : nullptr;
            graph[i][j]->successor[1] = (j + 1 < length(context.getSeqV())) ? graph[i][j+1] : nullptr;
            graph[i][j]->setRefCount(((i > 0) ? 1 : 0) + ((j > 0) ? 1 : 0));
        }
    }
    lastTask(graph)->mIsLastTask = true;
    return graph;
}

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Not<IsVectorExecutionPolicy<TVecExecPolicy> >, void)
invoke(DPTaskGraph<DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb>, TSpec> & graph)
{
    lastTask(graph)->incrementRefCount();
    lastTask(graph)->spawn_and_wait_for_all(*firstTask(graph));
    lastTask(graph)->execute();
    tbb::task::destroy(*lastTask(graph));
}

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(IsVectorExecutionPolicy<TVecExecPolicy>, void)
invoke(DPTaskGraph<DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb>, TSpec> & graph)
{
    using DPDagTask = DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyTbb>;
    ConcurrentQueue<DPDagTask*> queue;
    firstTask(graph)->mSimdQueuePtr = &queue;
    appendValue(queue, firstTask(graph));
    lastTask(graph)->spawn(*firstTask(graph));

    {
        std::unique_lock<std::mutex> lk(lastTask(graph)->_taskContext.mLockEvent);
        lastTask(graph)->_taskContext.mReadyEvent.wait(lk, [&graph]{return lastTask(graph)->_taskContext.mReady;});
    }
    SEQAN_ASSERT(empty(queue));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_TBB_H_
