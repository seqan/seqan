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

template <typename TConfig_, typename TSpec = ParallelAlignWavefront<>>
struct WavefrontIncubator
{
    // ----------------------------------------------------------------------------
    // Typedefs.
    
    // DPTrait type forwarding.
    using TConfig           = TConfig_;
    using TScoreValue       = typename Value<typename TConfig::TScore>::Type;
    using TAlgorithmType    = typename TConfig::TAlgorithm;
    using TTracebackType    = typename TConfig::TTraceback;
    using TGapType          = typename TConfig::TGap;
    
    // Wavefront Context.
    using TDPCell           = DPCell_<TScoreValue, TGapType>;
    using TBufferValue      = Pair<TDPCell, typename TraceBitMap_<>::Type>;
    using TBuffer           = String<TBufferValue>;
    using TBlockBuffer      = DPTileBuffer<TBuffer>;
    
    // DP Execution Context.
    using TDPProfile        = DPProfile_<TAlgorithmType, TGapType, TTracebackType, Parallel>;
    using TDPCache          = DPContext<TDPCell, typename TraceBitMap_<>::Type>;
    using TDPScout          = DPScout_<TDPCell, typename ScoutSpecForAlignmentAlgorithm_<TAlgorithmType>::Type>;

    // Parallel Context.
    using TThreadLocal              = ThreadLocal<IntermediateDPResult<typename TIncubator::TScout>,
                                                  typename TIncubator::TDPCache>::Type;
    using TEnumerableThreadSpecific = EnumerableThreadSpecific<TThreadLocal>;

    // ----------------------------------------------------------------------------
    // Static member functions.

    template <typename TTaskContext>
    static auto createTaskGraph(TTaskContext & taskContext)
    {
        using TDagTask = DPTaskImpl<TTaskContext, ParallelStd>;

        DPTaskGraph<TDagTask> graph;

        resize(graph.get(), length(taskContext.getSeqH()));
        for (int i = length(taskContext.getSeqH()); --i >= 0;)
        {
            resize(graph[i], length(taskContext.getSeqV()));
            for (int j = length(taskContext.getSeqV()); --j >= 0;)
            {
                graph[i][j].reset(new TDagTask(i, j, taskContext));
                graph[i][j]->successor[0] = (i + 1 < length(context.getSeqH())) ? graph[i+1][j].get() : nullptr;
                graph[i][j]->successor[1] = (j + 1 < length(context.getSeqV())) ? graph[i][j+1].get() : nullptr;
                graph[i][j]->setRefCount(((i > 0) ? 1 : 0) + ((j > 0) ? 1 : 0));
                graph[i][j]->_lastHBlock = (i + 1 == length(context.getSeqH()));
                graph[i][j]->_lastVBlock = (j + 1 == length(context.getSeqV()));
            }
        }
        lastTask(graph)->mIsLastTask = true;
        return graph;
    }
    
    // ----------------------------------------------------------------------------
    // Function createBlockBuffer()
    // ----------------------------------------------------------------------------
    
    template <typename TSeqH, typename TSeqV>
    static auto createBlockBuffer(TSeqH const & seqH, TSeqV const & seqV, TScore const & score, size_t const blockSize)
    {
        TBlockBuffer buffer;
        resize(buffer.horizontalBuffer, length(seqH), Exact());
        resize(buffer.verticalBuffer, length(seqV), Exact());
        
        TBufferValue tmp;
        tmp.i2 = _doComputeScore(tmp.i1, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), score, RecursionDirectionZero(), TDPProfile());
        for (auto itH = begin(buffer.horizontalBuffer, Standard()); itH != end(buffer.horizontalBuffer, Standard()); ++itH)
        {
            resize(*itH, length(front(seqH)), Exact());
            for (auto it = begin(*itH, Standard()); it != end(*itH, Standard()); ++it)
            {
                it->i2 = _doComputeScore(it->i1, TDPCell(), tmp.i1, TDPCell(), Nothing(), Nothing(), score, RecursionDirectionHorizontal(), TDPProfile());
                tmp.i1 = it->i1;
            }
        }
        tmp.i1 = decltype(tmp.i1){};
        tmp.i2 = _doComputeScore(tmp.i1, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), score, RecursionDirectionZero(), TDPProfile());
        
        for (auto itV = begin(buffer.verticalBuffer, Standard()); itV != end(buffer.verticalBuffer, Standard()); ++itV)
        {
            resize(*itV, length(front(seqV)) + 1, Exact());
            auto it = begin(*itV, Standard());
            it->i2 = tmp.i2;
            it->i1 = tmp.i1;
            ++it;
            for (; it != end(*itV, Standard()); ++it)
            {
                it->i2 = _doComputeScore(it->i1, TDPCell(), TDPCell(), tmp.i1, Nothing(), Nothing(), score, RecursionDirectionVertical(), TDPProfile());
                tmp.i1 = it->i1;
                tmp.i2 = it->i2;  // TODO(rrahn): Move out of loop.
            }
        }
        return buffer;
    }
    
    // ----------------------------------------------------------------------------
    // Function createBlocks()
    // ----------------------------------------------------------------------------
    
    template <typename TSeq>
    static auto createBlocks(TSeq const & seq, size_t const blockSize)
    {
        using TIter = typename Iterator<typename Infix<TSeq const>::Type, Standard>::Type;
        String<Range<TIter>> blocks;
        resize(blocks, (length(seq) + blockSize - 1) / blockSize, Exact());
        
        for (unsigned id = 0; id < length(blocks); ++id)
            blocks[id] = toRange(infix(seq, id * blockSize, _min(length(seq),(id + 1) * blockSize)));
        return blocks;
    }
    
};

template <typename TExecutionPolicy>
struct BatchAlignmentExecutor;

template <typename TSpec>
struct BatchAlignmentExecutor<ExecutionPolicy<ParallelAlignWavefront<TSpec>, Serial>>
{
    // We use this thread here.
    struct AlignThread
    {
        template <typename TQueue>
        inline void
        operator()(TQueue & queue)
        {
            lockWriting(queue);
            while (true)
            {
                auto task = popFront(queue);

                if (task == nullptr)
                    return;

                SEQAN_ASSERT(task != nullptr);
                task->execute();
            }
        }
    };

    // We can only support score based alignment for now.
    template <typename TExecPolicy,
              typename TSeqBatchH,
              typename TSeqBatchV,
              typename TConfig,
              typename TCallback>
    inline static void run(TExecPolicy const & execPol,
                           TSeqBatchH const & seqBatchH,
                           TSeqBatchV const & seqBatchV,
                           TDPConfig const & config,
                           TCallback && callback)
    {
        // Sequence types.
        using TSeqH = typename Value<TSequenceH>::Type;
        using TSeqV = typename Value<TSequenceV>::Type;
        // Wavefront Alignment Incubator
        using TIncubator = WavefrontIncubator<TConfig, ParallelAlignWavefront<TSpec>>;
        // Definition of TaskContext
        using TBlocksSeqH = decltype(TIncubator::createBlocks(seqBatchH[0], 1));
        using TBlocksSeqV = decltype(TIncubator::createBlocks(seqBatchV[0], 1));
        using TTaskContext = DPTaskContext<TBlocksSeqH, TBlocksSeqV, TIncubator>;
        // Parallel data structures.
        using TTask = DPTaskImpl<TTaskContext, Parallel>;
        using TWorkQueue = ConcurrentQueue<TTask *>;
        using TEts = typename TIncubator::TEnumerableThreadSpecific;

        // Create instance of task-queue and thread pool.
        TWorkQueue queue;
        ThreadPool pool;
        TEts ets;

        // We need to initialize the ets data structure with the number of paralle alingment instances.

        // Initialiaze thread pool and task queue.
        for (unsigned i = 0; i < numThreads(execPol); ++i)
        {
            emplaceBack(pool, AlignThread, std::ref(queue));
        }
        waitForWriters(queueContext.mQueue, numThreads(execPol));

        // Initialize instance of alignment scheduler.
        using TAlignmentInstance = AlignmentInstance<TSeqH, TSeqV, TConfig, TIncubator>;

        std::function<std::function<void(uint16_t)>(TAlignmentInstance)> aiCaller =
            [&queue, &ets, &callback](auto && func) mutable {
                return [&queue, &ets, &callback, func = std::move(func)](uint16_t id) mutable {
                func(id, queue, ets, std::forward<TCallback>(callback));
            };
        };

        using TCallable = decltype(aiCaller(TAlignmentInstance{}));
        AlignmentScheduler<TCallable> scheduler{parallelInstances(execPol)};
        // Iterate over pair of sequences.
        auto zip_cont = make_zip(seqSetH, seqSetV);
        for (std::tie<seqH, seqV> tie : zip_cont)
        {
            schedule(scheduler, aiCaller(TAlignmentInstance{seqH, seqV, config, blockSize(execPol)}));
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace impl
}  // namespace seqan
#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGNMENT_IMPL_WAVEFRONT_H_
