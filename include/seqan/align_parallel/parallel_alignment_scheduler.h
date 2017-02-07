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

    using TTaskPool = ConcurrentQueue<TCallable, Suspendable<Limit>>;

    // Private Member Variables.
    // ----------------------------------------------------------------------------
    
    std::vector<std::thread> mPool;
    std::vector<uint16_t>    mRecycableIds;
    TTaskPool                mQueue;
    std::mutex               mMtx;
    std::function<void()>    mDeleter = []{};  // empty no-op deleter

    // Constructors.
    // ----------------------------------------------------------------------------

    AlignmentScheduler(unsigned const readerCount) :
    {
        SEQAN_ASSERT_GT(readerCount, 0u);  // Bad if reader is 0.

        idBuffer.resize(readerCount);
        std::iota(std::begin(mRecycableIds), std::end(mRecycableIds), 0);
        
        setReaderWriterCount(mQueue, readerCount, 1);

        try
        {
        // Create the threads here, later we can try to make lazy thread creation.
            for (unsigned i = 0; i < readerCount; ++i)
            {
                mPool.emplace_back([this]()
                {

                    for (;;)
                    {
                        TCallable callable;
                        if (!popFront(callable, mQueue))
                            break;  // End of thread => No writers and queue is empty.

                        uint16_t id = -1;
                        {
                            std::lock_guard<std::mutex> lck(mMtx);
                            if (mRecycableIds.empty())
                                throw std::exception("Exception!");  // TODO(rrahn): Add propper exception.
                            id = mRecycableIds.pop_back();
                        }
                        callable(id);  // invokes the alignment with assigned id.
                        { // recycle id, when done.
                            std::lock_guard<std::mutex> lck(mMtx);
                            mRecycableIds.push_back(id);
                        }
                    }
                    unlockReading(mQueue);  // Nothing is getting into the queue anymore.
                });
            }
        }
        catch(...)  // TODO(rrahn): Add correct error handling here.
        {
            // Make sure all the spawned threads are safely stopped before re-throwing the exception.
            setReaderWriterCount(mQueue, 0, 0);
            for_each(mPool.begin(), mPool.end(), [](auto & t){ t.join(); });
            throw;
        }
    }
    
    // Default constructor.
    AlignmentScheduler() : AlignmentScheduler(1)
    {}

    // Copy & Move C'tor
    AlignmentScheduler(AlignmentScheduler const &)              = delete;
    AlignmentScheduler(AlignmentScheduler &&)                   = delete;

    // Special C'tor, to pass in a deleter
    AlignmentScheduler(unsigned const readerCount, auto && deleter) :
        AlignmentScheduler(rederCount),
        mDeleter(std::forward<decltype(deleter)>(deleter))
    {}

    // Copy & Move assignment
    AlignmentScheduler& operator=(AlignmentScheduler const &)   = delete;
    AlignmentScheduler& operator=(AlignmentScheduler &&)        = delete;

    // Implicitly deleted copy constructor (because of std::thread).
    // Default move constructor.
    
    // Destructor.
    ~AlignmentScheduler()
    {
        // No more alignments will be added.
        unlockWriting(mQueue);

        // Join the threads.
        for_each(mPool.begin(), mPool.end(), [](auto & t){ t.join(); });

        mDeleter();
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

// basic exception-safety guarantee.
// Throws if appendValue failed.
template <typename TCallable>
void schedule(AlignmentScheduler<TCallable> & scheduler,
              TCallable && callable)
{
    // Spins until there is enough space to add to the queue.
    if (!appendValue(scheduler.mAlignInstanceQueue, std::forward<TCallable>(callable)))
        throw std::runtime_exception("Called append value on invalid queue. Maybe reader count 0?");
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_ALIGNMENT_SCHEDULER_H_
