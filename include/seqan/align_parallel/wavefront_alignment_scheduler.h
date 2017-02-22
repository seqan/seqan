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

class WavefrontAlignmentScheduler
{
public:

    //-------------------------------------------------------------------------
    // Member Types.

    using TCallable         = std::function<void(uint16_t)>;
    using TAlignmentQueue   = ConcurrentQueue<TCallable, Suspendable<Limit>>;
    using TRecycleList      = std::list<uint16_t>;

    //-------------------------------------------------------------------------
    // Private Member Variables.

    WavefrontTaskScheduler  _mTaskScheduler;
    ThreadPool              _mPool;
    TRecycleList            _mRecycableIds;
    TAlignmentQueue         _mQueue;
    bool                    _mReceivedEndSignal;

    std::mutex              _mMutexRecycleId;
    unsigned                _mNumParallelAlignments;

    std::mutex                       _mMutexPushException;
    std::vector<std::exception_ptr>  _mExceptionPtrs;

    std::atomic<bool>               _mIsValid{true};


    std::function<void()> job = [this]()
    {
        for (;;)
        {
            TCallable callable;
            if (!popFront(callable, _mQueue))
                break;  // End of thread => No writers and queue is empty.

            uint16_t id = -1;
            { // Receive id.
                std::lock_guard<std::mutex> lck(_mMutexRecycleId);
                SEQAN_ASSERT_NOT(_mRecycableIds.empty());
                id = _mRecycableIds.front();
                _mRecycableIds.pop_front();
            }

            try
            {
                callable(id);  // invokes the alignment with assigned id.
            }
            catch(...)
            {  // Catch any exception thrown by callable. Store exception, and set this invalid.
               // We still keep running until the queue is empty. The thread is cleaned either by,
               // explicit wait or by destruction of this.
                _mIsValid.store(false, std::memory_order_release);
                {
                    std::lock_guard<std::mutex> lck(_mMutexPushException);
                    _mExceptionPtrs.push_back(std::current_exception());
                }
            }

            // Check if task scheduler is still valid.
            // If not, something went wrong, and we should not continue adding new tasks.
            // So we propagate the invalid state to this and break exceution chain.
            if (!isValid(_mTaskScheduler))
            {
                _mIsValid.store(false, std::memory_order_release);
            }

            { // recycle id, when done.
                std::lock_guard<std::mutex> lck(_mMutexRecycleId);
                _mRecycableIds.push_back(id);
            }
        }
        unlockReading(_mQueue);  // Notify that this reader is finished.
        unlockWriting(_mTaskScheduler);  // Notify that this writer is finished.
    };

    //-------------------------------------------------------------------------
    // Constructors.

    // implicitly deleted default constructor.

    WavefrontAlignmentScheduler(size_t const numParallelAlignments, size_t const numParallelTasks) :
        _mTaskScheduler(numParallelTasks),
        _mQueue(numParallelAlignments),
        _mReceivedEndSignal(false),
        _mNumParallelAlignments(numParallelAlignments)
    {
        SEQAN_ASSERT_GT(numParallelAlignments, 0u);  // Bad if reader is 0.

        // Setup recycable ids.
        _mRecycableIds.resize(numParallelAlignments);
        std::iota(std::begin(_mRecycableIds), std::end(_mRecycableIds), 0);

        setReaderWriterCount(_mQueue, numParallelAlignments, 1);

        _mExceptionPtrs.resize(numParallelAlignments, nullptr);

        try
        { // Create the threads here, later we can try to make lazy thread creation.
            for (unsigned i = 0; i < numParallelAlignments; ++i)
            {
                spawn(_mPool, job);
            }
        }
        catch(...)  // Make sure all the spawned threads are safely stopped before re-throwing the exception.
        {
            unlockWriting(_mQueue);
            waitForWriters(_mTaskScheduler);
            join(_mPool);
            throw;
        }

        setWriterCount(_mTaskScheduler, numParallelAlignments);
        // Notify task scheduler, that everything was setup correctly.
        for (unsigned i = 0; i < numParallelAlignments; ++i)
            lockWriting(_mTaskScheduler);
        waitForWriters(_mTaskScheduler);  // Invoke task scheduler.
    }
    
    // Default constructor.
    WavefrontAlignmentScheduler() : WavefrontAlignmentScheduler(16, 8)
    {}

    // Copy & Move C'tor
    WavefrontAlignmentScheduler(WavefrontAlignmentScheduler const &)              = delete;
    WavefrontAlignmentScheduler(WavefrontAlignmentScheduler &&)                   = delete;


    ///-------------------------------------------------------------------------
    // Destructor.

    ~WavefrontAlignmentScheduler()
    {
        // Signal that no more alignments will be added.
        if (!_mReceivedEndSignal)
            unlockWriting(_mQueue);

        SEQAN_ASSERT(_mQueue.writerCount == 0);

        // Wait until all remaining threads are finished with their execution.
        join(_mPool);

        // In destructor of thread pool we wait for the outstanding alignments to be finished
        // and then continue destruction of the remaining members and cleaning up the stack.
    }

    // ------------------------------------------------------------------------
    // Member Functions.

    // Copy & Move assignment
    WavefrontAlignmentScheduler& operator=(WavefrontAlignmentScheduler const &)   = delete;
    WavefrontAlignmentScheduler& operator=(WavefrontAlignmentScheduler &&)        = delete;
};

// ============================================================================
// Metafunctions
// ============================================================================

template<>
struct SchedulerTraits<WavefrontAlignmentScheduler>
{
    using TTask = typename WavefrontAlignmentScheduler::TCallable;
};

// ============================================================================
// Functions
// ============================================================================

inline bool
isValid(WavefrontAlignmentScheduler const & me)
{
    return me._mIsValid.load(std::memory_order_acquire);
}

// basic exception-safety guarantee.
// Throws if appendValue failed.
inline void
scheduleTask(WavefrontAlignmentScheduler & me,
             typename SchedulerTraits<WavefrontAlignmentScheduler>::TTask && callable)
{
    if (!isValid(me))
        throw std::runtime_error("Invalid alignment scheduler!");

    // Spins until there is enough space to add to the queue.
    if (!appendValue(me._mQueue, std::forward<decltype(callable)>(callable)))
        throw std::runtime_error("Invalid alignment scheduler 2!");
}

inline void
scheduleTask(WavefrontAlignmentScheduler & me,
             typename SchedulerTraits<WavefrontAlignmentScheduler>::TTask & callable)
{
    if (!isValid(me))
        throw std::runtime_error("Invalid alignment scheduler!");
    // Spins until there is enough space to add to the queue.
    if(!appendValue(me._mQueue, callable))
        throw std::runtime_error("Invalid alignment scheduler 2!");
}

// Only possible if some other thread is signaling the end of it.
inline void
signalAndWait(WavefrontAlignmentScheduler & me)
{
    unlockWriting(me._mQueue);
    me._mReceivedEndSignal = true;
    join(me._mPool);
    wait(me._mTaskScheduler);
}

inline auto
getExceptions(WavefrontAlignmentScheduler & me)
{
    auto vec = me._mExceptionPtrs;
    auto innerExceptions = getExceptions(me._mTaskScheduler);
    std::copy(std::begin(innerExceptions), std::end(innerExceptions), std::back_inserter(vec));
    return vec;
}

inline auto&
taskScheduler(WavefrontAlignmentScheduler & me)
{
    return me._mTaskScheduler;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_ALIGNMENT_SCHEDULER_H_
