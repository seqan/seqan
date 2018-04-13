// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_SCHEDULER_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_SCHEDULER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Scheduler for wavefront tasks.
class WavefrontTaskScheduler
{
public:

    //-------------------------------------------------------------------------
    // Memeber Types

    using TWrapper = std::function<void()>;
    using TTaskQueue = ConcurrentQueue<TWrapper>;

    //-------------------------------------------------------------------------
    // Member Variables

    ThreadPool  _threadPool;
    TTaskQueue  _taskQueue;

    unsigned    _writerCount;

    std::mutex                      _mutexPushException;
    std::vector<std::exception_ptr> _exceptionPointers;
    std::atomic<bool>               _isValid{true};

    std::function<void()> job = [this] ()
    {
        lockReading(_taskQueue);
        waitForFirstValue(_taskQueue);  // Wait for all writers to be setup.

        std::function<void()> _dummy = [] ()
        {  // TODO(rrahn): Could throw exception to signal something went terribly wrong.
            SEQAN_ASSERT_FAIL("Trying to exceute empty wavefront task in a thread");
        };
        TWrapper task{_dummy};

        while (true)
        {
            if (!popFront(task, _taskQueue))
                break;  // Empty queue and no writer registered.

            try
            {
                task();  // Execute the task;
            }
            catch (...)
            {  // Catch exception, and signal failure. Continue running until queue is empty.
                {
                    std::lock_guard<std::mutex> lck(_mutexPushException);
                    _exceptionPointers.push_back(std::current_exception());
                }
                _isValid.store(false, std::memory_order_release);
            }
        }
        unlockReading(_taskQueue);
    };

    //-------------------------------------------------------------------------
    // Constructor

    WavefrontTaskScheduler(size_t const threadCount, size_t const writerCount) :
        _writerCount(writerCount)
    {

        for (unsigned i = 0; i < threadCount; ++i)
        {
            spawn(_threadPool, job);
        }
        setCpuAffinity(_threadPool, 0, 1);
    }

    WavefrontTaskScheduler(size_t const threadCount) : WavefrontTaskScheduler(threadCount, 0)
    {}

    WavefrontTaskScheduler(WavefrontTaskScheduler const &) = delete;
    WavefrontTaskScheduler(WavefrontTaskScheduler &&) = delete;

    //-------------------------------------------------------------------------
    // Member Functions

    WavefrontTaskScheduler& operator=(WavefrontTaskScheduler const &) = delete;
    WavefrontTaskScheduler& operator=(WavefrontTaskScheduler &&) = delete;

    //-------------------------------------------------------------------------
    // Destructor

    ~WavefrontTaskScheduler()
    {}
    // In destructor of thread pool we wait for the outstanding alignments to be finished
    // and then continue destruction of the remaining members and cleaning up the stack.
    // Note the number of writers must be set to 0, for the queue to stop spinning.
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TScheduler>
struct SchedulerTraits;

template <>
struct SchedulerTraits<WavefrontTaskScheduler>
{
    using TWrapper_ = typename WavefrontTaskScheduler::TWrapper;
    using TTask  = TWrapper_;
};

// ============================================================================
// Functions
// ============================================================================

inline void
setWriterCount(WavefrontTaskScheduler & me, size_t const count) noexcept
{
    me._writerCount = count;
}

inline void
lockWriting(WavefrontTaskScheduler & me) noexcept
{
    lockWriting(me._taskQueue);
}

inline void
unlockWriting(WavefrontTaskScheduler & me) noexcept
{
    unlockWriting(me._taskQueue);
}

inline void
waitForWriters(WavefrontTaskScheduler & me) noexcept
{
    waitForWriters(me._taskQueue, me._writerCount);
}

inline bool
isValid(WavefrontTaskScheduler & me) noexcept
{
    return me._isValid.load(std::memory_order_acquire);
}

inline void
scheduleTask(WavefrontTaskScheduler & me,
             typename SchedulerTraits<WavefrontTaskScheduler>::TTask task)
{
    if (!isValid(me))
    {  // TODO(rrahn): Improve error handling.
        throw std::runtime_error("Invalid Task Scheduler");
    }
    appendValue(me._taskQueue, std::move(task));
}

inline void
wait(WavefrontTaskScheduler & me)
{
    SEQAN_ASSERT(me._taskQueue.writerCount == 0);

    join(me._threadPool);

    SEQAN_ASSERT(empty(me._taskQueue));
    SEQAN_ASSERT(me._taskQueue.readerCount == 0);
}

inline auto
getExceptions(WavefrontTaskScheduler & me)
{
    return me._exceptionPointers;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_SCHEDULER_H_
