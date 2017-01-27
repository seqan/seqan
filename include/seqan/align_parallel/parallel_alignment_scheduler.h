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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_ALIGNMENT_SCHEDULER_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_ALIGNMENT_SCHEDULER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TVal>
struct AlignmentSchedulerTraits
{
    using TTaskPool = ConcurrentQueue<TVal, Suspendable<Limit>>;
};

template <typename TCallable,
          typename TTraits = AlignmentSchedulerTraits<TCallable>>
class AlignmentScheduler
{
public:

    // Typedefs.
    // ----------------------------------------------------------------------------

    using TTaskPool = TTraits::TTaskPool;

    // Member Variables.
    // ----------------------------------------------------------------------------
    
    std::vector<std::thread> mPool;
    TTaskPool                mQueue;

    // Constructors.
    // ----------------------------------------------------------------------------

    AlignmentScheduler(unsigned const reader) :
    {
        SEQAN_ASSERT_GT(reader, 0u);  // Bad if reader is 0.
        
        setReaderWriterCount(mQueue, reader, 1);
        
        // Create the threads here, later we can try to make lazy thread creation.
        for (unsigned i = 0; i < reader; ++i)
            mPool.emplace_back([this](){
                for (;;) {
                    TCallable callable;
                    if (popFront(callable, mQueue))
                        callable();  // invokes the alingment.
                    else
                        return; // no elements in the queue and no writer registered.
                }
            });
    }
    
    // Default constructor.
    AlignmentScheduler() : AlignmentScheduler(1)
    {}
    
    // Implicitly deleted copy constructor (because of std::thread).
    // Default move constructor.
    
    // Destructor.
    ~AlignmentScheduler()
    {
        setReaderWriterCount(mQueue, 0, 0);  // We notify that there is
        for_each(mPool.begin(), mPool.end(), [](auto & t){ t.join(); });
    }

    // Member Functions.
    // ----------------------------------------------------------------------------
    
    // Implicitly deleted copy assignment operator (because of std::thread).
    // Default move assignment operator.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// TODO(rrahn): Remove!
//template <typename TDPScoreStore, typename TDPTraceStore, typename TExecutionPolicy>
//inline void
//init(ParallelAlignScheduler<TDPScoreStore, TDPTraceStore, TExecutionPolicy> & scheduler)
//{
//
//}
//
//template <typename TExecutionPolicy, typename TAlignInstances>
//execute(ParallelAlignScheduler<TExecutionPolicy> & execPolicy,
//        TAlignInstances & alignInstances)
//{
//    using TAlignInstance = typename Value<TAlignInstances>::Type;
//    using TInstanceTrait = typename Traits<TAlignInstances>::Type;
//
//    // Initialize scheduler:
//
//
//    while (!empty(alignInstances))
//    {
//        auto& instance = popFront(alignInstances);
//        SEQAN_ASSERT(instance != nullptr);
//        std::thread(*instance, std::ref()).detach();  // Invocation must run in detached thread!
//        ++sharedConfig.runningInstances;
//        auto& sharedConfig = getSharedConfig(instance);
//        {
//            std::unique_lock<decltype(sharedConfig.schedulerLock)> lk(sharedConfig.schedulerLock);
//            cond.wait_for(lk, []{ sharedConfig.runningInstances < TInstanceTrait::MAX_PARALLEL_INSTANCES; } );
//        }
//    }
//
//    // Finally wait for running jobs.
//    std::unique_lock<decltype(sharedConfig.schedulerLock)> lk(sharedConfig.schedulerLock);
//    cond.wait_for(lk, []{ sharedConfig.runningInstances < TInstanceTrait::MAX_PARALLEL_INSTANCES; } );
//}

// basic exception-safety guarantee.
// Throws if appendValue failed.
template <typename TCallable>
void schedule(AlignmentScheduler<TCallable> & alignScheduler,
              TCallable && callable)
{
    // Spins until there is enough space to add to the queue.
    if (!appendValue(alignScheduler.mAlignInstanceQueue, std::forward<TCallable>(callable)))
        throw std::runtime_exception("Called append value on invalid queue. Maybe reader count 0?");
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_ALIGNMENT_SCHEDULER_H_
