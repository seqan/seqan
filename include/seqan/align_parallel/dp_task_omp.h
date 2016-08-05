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

struct DPTaskOmp_;
typedef Tag<DPTaskOmp_> DPTaskOmp;

template <typename TTaskConfig, typename TThreadLocalStorage>
class DPTaskImpl<TTaskConfig, TThreadLocalStorage, DPTaskOmp> :
    public DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, DPTaskOmp> >
{
public:

    using TSize = typename TTaskConfig::TSize;
    using TBase = DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, DPTaskOmp> >;

    // ============================================================================
    // Member variables.

    std::atomic<unsigned>   refCount;

    // ============================================================================
    // Constructor.

    DPTaskImpl(TSize pCol, TSize pRow, TTaskConfig & pConfig) :
        TBase(pCol, pRow, pConfig, *this)
    {}

    // ============================================================================
    // Member functions.

    void setRefCount(unsigned const n)
    {
        refCount.store(n, std::memory_order_relaxed);
    }

    unsigned decrementRefCount()
    {
        return --refCount;
    }

    unsigned incrementRefCount()
    {
        return ++refCount;
    }

    template <typename TThreadStore>
    DPTaskImpl* execute(TThreadStore & tls)
    {
        TBase::execute(tls[omp_get_thread_num()]);

        // spawning right and downward neighbors
        if (DPTaskImpl* t = TBase::successor[0])
        {
            if (t->decrementRefCount() == 0)
            {
                if (TBase::_col % 4 == 0)  //only spawn new task every fourth item
                {
                    SEQAN_OMP_PRAGMA(task shared(tls) firstprivate(t) untied)
                    {
                        t->execute(tls);
                    }
                }
                else
                { //use existing thread
                    t->execute(tls);
                }
            }
        }
        if (DPTaskImpl* t = TBase::successor[1])
        {  //use existing thread.
            if (t->decrementRefCount() == 0)
            {
                t->execute(tls);
            }
        }
        return nullptr;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TTaskContext, typename TThreadLocalStorage>
struct IsDPTask<DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskOmp> > : True
{};

template <typename TTaskContext, typename TThreadLocalStorage>
struct Pointer_<DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskOmp> >
{
    using TTask_ = DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskOmp>;
    using Type  = std::unique_ptr<TTask_>;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TTaskContext>
inline auto
createGraph(TTaskContext & context, DPTaskOmp const & /*taskImplTag*/)
{
    using TThreadLocalStorage = typename TTaskContext::TDPContext;
    using TDagTask = DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskOmp>;

    DPTaskGraph<TDagTask> graph;

    resize(graph.get(), length(context.getSeqH()));
    for (int i = length(context.getSeqH()); --i >= 0;)
    {
        resize(graph[i], length(context.getSeqV()));
        for (int j = length(context.getSeqV()); --j >= 0;)
        {
            graph[i][j].reset(new TDagTask(i, j, context));
            graph[i][j]->successor[0] = (i + 1 < length(context.getSeqH())) ? graph[i+1][j].get() : nullptr;
            graph[i][j]->successor[1] = (j + 1 < length(context.getSeqV())) ? graph[i][j+1].get() : nullptr;
            graph[i][j]->setRefCount(((i > 0) ? 1 : 0) + ((j > 0) ? 1 : 0));
        }
    }
    lastTask(graph)->incrementRefCount();
    return graph;
}

template <typename TTaskContext, typename TThreadLocalStorage, typename TSpec>
inline void
invoke(DPTaskGraph<DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskOmp>, TSpec> & graph)
{
    std::vector<TThreadLocalStorage> tls(omp_get_max_threads());
    SEQAN_OMP_PRAGMA(parallel)
    {
        SEQAN_OMP_PRAGMA(master)
        {
            // Create context and pass to execute.
            //now kick the computaion off
            firstTask(graph)->execute(tls);
        }
        SEQAN_OMP_PRAGMA(barrier)
    }
    SEQAN_ASSERT_EQ(lastTask(graph)->refCount, 1u);
    lastTask(graph)->execute(tls);
}
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_OMP_H_
