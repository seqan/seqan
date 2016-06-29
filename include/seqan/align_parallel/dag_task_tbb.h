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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DAG_TASK_TBB_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DAG_TASK_TBB_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TTraits>
class DpDagTask<TTraits, TbbTask> : public tbb::task
{
public:
    DagTask* successor[2];

private:
    const int i, j;
    bool isExecuted;
    TDpContext& dpContext;
    String<TTraceTile>& tracebackBlock;
    const TSeqH& seqH;
    const TSeqV& seqV;
    TTileBuffer& tileBuffer;
    const Score<TScoreValue, TScoreSpec>& score;
    unsigned simdVectorLength = 4;
    ConcurrentQueue<DagTask*> & queue;
    std::mutex& lock;

public:
    DagTask(const int i, const int j,
            TDpContext& dpContext,
            String<TTraceTile>& tracebackBlock,
            const TSeqH& seqH, const TSeqV& seqV,
            TTileBuffer& tileBuffer,
            const Score<TScoreValue, TScoreSpec>& score,
            ConcurrentQueue<DagTask*> & pQueue,
            std::mutex & pLock) :
    i(i), j(j), isExecuted(false),
    dpContext(dpContext),
    tracebackBlock(tracebackBlock),
    seqH(seqH), seqV(seqV),
    tileBuffer(tileBuffer),
    score(score),
    queue(pQueue),
    lock(pLock)
    {}

    task* execute()
    {
        String<DagTask*> tasks;
        {
            std::lock_guard<std::mutex> scopedLock(lock);

            if (empty(queue))  // Make sure task wasn't executed by another task.
                return nullptr;

            if (length(queue) < simdVectorLength)
                appendValue(tasks, popFront(queue));
            else
                for (auto i = 0; i < simdVectorLength; ++i)
                    appendValue(tasks, popFront(queue));
        }   // Release scoped lock.

        SEQAN_ASSERT_GT(length(tasks), 0u);

        if (length(tasks) == 1)
        {
            auto t = tasks[0];
            //                printf("%i %i\n", t->i, t->j);
            auto& localDpContext = dpContext.local();
            getDpTraceMatrix(localDpContext) = t->tracebackBlock[(t->i * length(seqV)) + t->j];
            TDPScoutState scoutState(tileBuffer.horizontalBuffer[t->i], tileBuffer.verticalBuffer[t->j]);
            String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.
            _computeAlignment(localDpContext, traceSegments, scoutState, seqH[t->i], seqV[t->j], score, DPBandConfig<BandOff>(), TDPProfile());
        }
        else
        {
            // Make to simd version!
            for (auto& task : tasks)
            {
                //                    printf("%i %i\n", task->i, task->j);
                auto& localDpContext = dpContext.local();
                getDpTraceMatrix(localDpContext) = task->tracebackBlock[(task->i * length(seqV)) + task->j];
                TDPScoutState scoutState(tileBuffer.horizontalBuffer[task->i], tileBuffer.verticalBuffer[task->j]);
                String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.
                _computeAlignment(localDpContext, traceSegments, scoutState, seqH[task->i], seqV[task->j], score, DPBandConfig<BandOff>(), TDPProfile());
            }
        }
        // spawning right and downward neighbors
        for (auto& task : tasks)
            for (int k = 0; k < 2; ++k)
                if (DagTask* t = task->successor[k])
                    if (t->decrement_ref_count() == 0)
                    {
                        SEQAN_ASSERT(t->i == 130 && t->j == 129);
                        appendValue(task->queue, t);
                        spawn(*t);
                    }
        return nullptr;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
}

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DAG_TASK_TBB_H_
