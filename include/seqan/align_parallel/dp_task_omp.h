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
// Author: Jonnhy Hancox <>
//         Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_OMP_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_OMP_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TTaskConfig, typename TThreadLocalStorage, typename TVecExecPolicy>
class DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp> :
    public DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp> >
{
public:

    using TSize = typename TTaskConfig::TSize;
    using TBase = DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp> >;

    // ============================================================================
    // Member variables.

    std::atomic<unsigned>   refCount;
    TThreadLocalStorage *   mTlsPtr;

    // ============================================================================
    // Constructor.

    DPTaskImpl(TSize pCol, TSize pRow, TTaskConfig & pConfig,
               TThreadLocalStorage & pTls) :
        TBase(pCol, pRow, pConfig, *this),
        mTlsPtr(&pTls)
    {}

    // ============================================================================
    // Member functions.

    inline void setRefCount(unsigned const n)
    {
        refCount.store(n, std::memory_order_relaxed);
    }

    inline unsigned decrementRefCount()
    {
        return --refCount;
    }

    inline unsigned incrementRefCount()
    {
        return ++refCount;
    }

    inline void
    updateAndSpawn()
    {
        if (DPTaskImpl* t = TBase::successor[0])
        {
            if (t->decrementRefCount() == 0)
            {
                if (TBase::_col % 4 == 0)  //only spawn new task every fourth item
                {
                    SEQAN_OMP_PRAGMA(task firstprivate(t) untied)
                    {
                        t->execute();
                    }
                }
                else
                { //use existing thread
                    t->execute();
                }
            }
        }
        if (DPTaskImpl* t = TBase::successor[1])
        {  //use existing thread.
            if (t->decrementRefCount() == 0)
            {
                t->execute();
            }
        }
    }

    template <typename TTaskQueue>
    inline void
    updateAndSpawn(TTaskQueue & pTaskQueue)
    {
        if (DPTaskImpl* t = TBase::successor[0])
        {
            if (t->decrementRefCount() == 0)
            {
                appendValue(pTaskQueue, t);
                if (TBase::_col % 4 == 0)  //only spawn new task every fourth item
                {
                    SEQAN_OMP_PRAGMA(task shared(pTaskQueue) firstprivate(t) untied)
                    {
                        t->execute(pTaskQueue);
                    }
                }
                else
                { //use existing thread
                    t->execute(pTaskQueue);
                }
            }
        }
        if (DPTaskImpl* t = TBase::successor[1])
        {  //use existing thread.
            if (t->decrementRefCount() == 0)
            {
                appendValue(pTaskQueue, t);
                t->execute(pTaskQueue);
            }
        }
    }

    DPTaskImpl* execute()
    {
        SEQAN_ASSERT(mTlsPtr != nullptr);
        TBase::runScalar((*mTlsPtr)[omp_get_thread_num()]);
        updateAndSpawn();
        return nullptr;
    }

    template <typename TTaskQueue>
    DPTaskImpl* execute(TTaskQueue& pTaskQueue)
    {
        SEQAN_ASSERT(mTlsPtr != nullptr);

        String<DPTaskImpl*> tasks;
        {  // Acquire scoped lock.
            std::lock_guard<decltype(TBase::_taskContext.mLock)> scopedLock(TBase::_taskContext.mLock);

            if (empty(pTaskQueue))
            {
                return nullptr;
            }

            lockReading(pTaskQueue);
            if (length(pTaskQueue) < TBase::_taskContext.mSimdLength)
            {
                appendValue(tasks, popFront(pTaskQueue));
            }
            else
            {
                SEQAN_ASSERT_GEQ(length(pTaskQueue), TBase::_taskContext.mSimdLength);
                for (auto i = 0; i < TBase::_taskContext.mSimdLength; ++i)
                    appendValue(tasks, popFront(pTaskQueue));
            }
            unlockReading(pTaskQueue);
        }  // Release scoped lock.

        SEQAN_ASSERT_GT(length(tasks), 0u);

        if (length(tasks) == 1)
        {
            front(tasks)->TBase::runScalar((*mTlsPtr)[omp_get_thread_num()]);
        }
        else
        {
            SEQAN_ASSERT_EQ(length(tasks), TBase::_taskContext.mSimdLength);
            TBase::runSimd(tasks, (*mTlsPtr)[omp_get_thread_num()]);
        }
        // Update and spawn sucessors.
        for (auto& task : tasks)
        {
            task->template updateAndSpawn(pTaskQueue);
//            if (task->mIsLastTask)
//            {
//                bool res = tbb::task::self().cancel_group_execution();  // Notify, that all queued tasks can be canceled.
//                SEQAN_ASSERT(res);
//                {
//                    std::lock_guard<std::mutex> lk(TBase::_taskContext.mLock);
//                    TBase::_taskContext.mReady = true;
//                }
//                TBase::_taskContext.mReadyEvent.notify_all();
//            }
        }
        return nullptr;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy>
struct IsDPTask<DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp> > : True
{};

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy>
struct Pointer_<DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp> >
{
    using TTask_ = DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp>;
    using Type  = std::unique_ptr<TTask_>;
};

namespace impl
{

namespace dp
{

namespace parallel
{

template <typename TLocalStore>
struct ThreadLocalStorage<TLocalStore, ParallelExecutionPolicyOmp>
{
    using Type = std::vector<TLocalStore>;
};
}  // namespace parallel
}  // namespace dp
}  // namespace impl

// ============================================================================
// Functions
// ============================================================================

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy>
inline auto
createGraph(TTaskContext & context,
            TThreadLocalStorage & tls,
            ExecutionPolicy<ParallelExecutionPolicyOmp, TVecExecPolicy> const & /*unused*/)
{
    using TDagTask = DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp>;

    DPTaskGraph<TDagTask> graph;

    resize(graph.get(), length(context.getSeqH()));
    for (int i = length(context.getSeqH()); --i >= 0;)
    {
        resize(graph[i], length(context.getSeqV()));
        for (int j = length(context.getSeqV()); --j >= 0;)
        {
            graph[i][j].reset(new TDagTask(i, j, context, tls));
            graph[i][j]->successor[0] = (i + 1 < length(context.getSeqH())) ? graph[i+1][j].get() : nullptr;
            graph[i][j]->successor[1] = (j + 1 < length(context.getSeqV())) ? graph[i][j+1].get() : nullptr;
            graph[i][j]->setRefCount(((i > 0) ? 1 : 0) + ((j > 0) ? 1 : 0));
        }
    }
    lastTask(graph)->incrementRefCount();
    return graph;
}

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Not<IsVectorExecutionPolicy<TVecExecPolicy> >, void)
invoke(DPTaskGraph<DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp>, TSpec> & graph,
       ExecutionPolicy<ParallelExecutionPolicyOmp, TVecExecPolicy> const & execPolicy)
{
    resize(*(firstTask(graph)->mTlsPtr), execPolicy.numThreads);
    SEQAN_OMP_PRAGMA(parallel num_threads(execPolicy.numThreads))
    {
        SEQAN_OMP_PRAGMA(master)
        {
            // Create context and pass to execute.
            //now kick the computaion off
            firstTask(graph)->execute();
        }
        SEQAN_OMP_PRAGMA(barrier)
    }
    SEQAN_ASSERT_EQ(lastTask(graph)->refCount, 1u);
    lastTask(graph)->execute();
}

template <typename TTaskContext, typename TThreadLocalStorage, typename TVecExecPolicy, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(IsVectorExecutionPolicy<TVecExecPolicy>, void)
invoke(DPTaskGraph<DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp>, TSpec> & graph,
       ExecutionPolicy<ParallelExecutionPolicyOmp, TVecExecPolicy> const & execPolicy)
{
    using TTask = DPTaskImpl<TTaskContext, TThreadLocalStorage, TVecExecPolicy, ParallelExecutionPolicyOmp>;

    resize(*(firstTask(graph)->mTlsPtr), execPolicy.numThreads);
    ConcurrentQueue<TTask*> queue;
    queue.writerCount = execPolicy.numThreads;
    appendValue(queue, firstTask(graph).get());

    SEQAN_OMP_PRAGMA(parallel num_threads(execPolicy.numThreads))
    {
        SEQAN_OMP_PRAGMA(master)
        {
            // Create context and pass to execute.
            //now kick the computaion off
            firstTask(graph)->execute(queue);
        }
        SEQAN_OMP_PRAGMA(barrier)
    }
    while (queue.writerCount != 0)
        unlockWriting(queue);

    SEQAN_ASSERT(empty(queue));
    appendValue(queue, lastTask(graph).get());
    lastTask(graph)->execute();
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_OMP_H_
