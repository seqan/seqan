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

struct DPTaskTbb_;
typedef Tag<DPTaskTbb_> DPTaskTbb;

template <typename TTaskConfig, typename TThreadLocalStorage>
class DPTaskImpl<TTaskConfig, TThreadLocalStorage, DPTaskTbb> :
    public DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, DPTaskTbb> >,
    public tbb::task
{
public:

    using TBase = DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, DPTaskTbb> >;

    // ============================================================================
    // Member variables.

    std::shared_ptr<TThreadLocalStorage> _ptrTls;

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

    template <typename TArray>
    void updateAndSpawnSuccessors(TArray & successor)
    {
        for (auto& t : successor)
            if (t != nullptr && t->decrementRefCount() == 0)
                spawn(*t);
    }

    task* execute()
    {
        TBase::execute(_ptrTls->local());
        for (auto& t : TBase::successor)
            if (t != nullptr && t->decrementRefCount() == 0)
                spawn(*t);
        return nullptr;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TTaskContext, typename TThreadLocalStorage>
struct IsDPTask<DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskTbb> > : True
{};

// ============================================================================
// Functions
// ============================================================================

template <typename TTaskContext>
inline auto
createGraph(TTaskContext & context, DPTaskTbb const & /*taskImplTag*/)
{
    using TThreadLocalStorage = tbb::enumerable_thread_specific<typename TTaskContext::TDPContext>;
    using TDagTask = DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskTbb>;

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
    graph[length(context.getSeqH())-1][length(context.getSeqV())-1]->incrementRefCount();
    return graph;
}

template <typename TTaskContext, typename TThreadLocalStorage, typename TSpec>
inline void
invoke(DPTaskGraph<DPTaskImpl<TTaskContext, TThreadLocalStorage, DPTaskTbb>, TSpec> & graph)
{
    lastTask(graph)->spawn_and_wait_for_all(*firstTask(graph));
    lastTask(graph)->execute();
    tbb::task::destroy(*lastTask(graph));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_TBB_H_
