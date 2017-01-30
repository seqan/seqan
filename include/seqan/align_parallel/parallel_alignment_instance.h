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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INSTANCE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INSTANCE_H_

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
    
template <typename T>
struct AICallableRaiiWrapper {
    
    unique_ptr<T> mData;
    
    template <typename ...TArgs>
    AICallableRaiiWrapper(TArgs && ...args) mData(std::make_unique<T>(std::forward<TArgs>(args)...))
    {}
    
    // Call operator() of the by mData pointed-to object.
    void operator()()
    {
        mData->operator();
    }
};

template <typename TSeqH,
          typename TSeqV,
          typename TConfig,
          typename TIncubator = WavefrontIncubator<TConfig>>
struct AlignmentInstance
{
public:

    // ----------------------------------------------------------------------------
    // Member Variables.
    // ----------------------------------------------------------------------------

    TSeqH const &   mSeqH;
    TSeqV const &   mSeqV;
    TConfig const & mDPConfig;
    size_t          mBlockSize;
    
    // ----------------------------------------------------------------------------
    // Constructors.
    // ----------------------------------------------------------------------------

    // ----------------------------------------------------------------------------
    // Member Functions.
    // ----------------------------------------------------------------------------
    
    // This function now run's in a separate thread.
    template <typename TConcurrentQueue,
              typename TEnumerableThreadSpecific,
              typename TCallback>
    inline void
    operator()(uint16_t const mAlignmentId,
               TConcurrentQueue & queue,
               TEnumerableThreadSpecific & ets,
               TCallback && callback)
    {
        // Initialize the strings.
        auto seqHBlocks = TIncubator::createBlocks(mSeqH, mBlockSize);
        auto seqVBlocks = TIncubator::createBlocks(mSeqV, mBlockSize);
        
        // Create the buffer for the matrix.
        auto buffer = TIncubator::createBlockBuffer(seqHBlocks, seqVBlocks, scoringScheme(mDPConfig), mBlockSize);

        // Create trace store -> where presumably the traceback is stored.
        // Now we can create the DPContext.
        
        // Setup DPTaskContext
        // task context with access to the blocks the score the band, and the buffer and the trace proxy
        DPTaskContext<TBlocksSeqH, TBlocksSeqV, TIncubator> dpContext{mAlignmentId,
                                                                      seqHBlocks, seqVBlocks,
                                                                      scoringScheme(mDPConfig),
                                                                      queue,
                                                                      ets};
        // What could we possibly need?
        // create DPLocalStorage as thread local storage.
        
        // create task graph -> nothing to do with the parallel framework
        auto taskGraph = TIncubator::createGraph(taskContext);

        appendValue(queue, firstTask(graph).get());  // Kick off execution.

        // Wait for threads to finish.
        {
            std::unique_lock<std::mutex> lck(lastTask(graph)->_taskContext.mLockEvent);
            lastTask(graph)->_taskContext.mReadyEvent.wait(lck, [&graph]{return lastTask(graph)->_taskContext.mReady;});
        }

        // invoke
            // std::threads -> vector of threads -> get executed with a shared queue.
            // per alignment instance
            // spawn(task) -> adds task* to the internal queue.
            // wait_for_all() clause to wait for finishing execution.
            // unlock writers and join threads. -> leave no threads dangling around. -> all tasks need to be computed.
        
            // tbb -> create task_group
            // spawn(first task of tg)
            // tg.run(with task->execution method)
            // tg.wait()  // barrier.
            // execute last task graph
        
            // omp -> more difficult to model.
            // global prallel section
            // master part triggers execution with queue.
            // sets barrier to end of parallel execution
        callback(score, trace);
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

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INSTANCE_H_
