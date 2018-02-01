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

#include <seqan/align_parallel.h>

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_alignment_scheduler_construct)
{
    using namespace seqan;

    // We need to be able to construct a thread pool.
    SEQAN_ASSERT(std::is_default_constructible<WavefrontAlignmentScheduler>::value);
    SEQAN_ASSERT(!std::is_copy_constructible<WavefrontAlignmentScheduler>::value);
    SEQAN_ASSERT(!std::is_move_constructible<WavefrontAlignmentScheduler>::value);
    SEQAN_ASSERT(!std::is_copy_assignable<WavefrontAlignmentScheduler>::value);
    SEQAN_ASSERT(!std::is_move_assignable<WavefrontAlignmentScheduler>::value);

    {  // Default construction and termination.
        try
        {
            WavefrontAlignmentScheduler scheduler{};
        }
        catch(...)
        {
            SEQAN_ASSERT_FAIL("Failed to default construct scheduler!");
        }
    }

    {  // Construction with some parameters and some values.

        try
        {
            WavefrontAlignmentScheduler scheduler{4, 2};
        }
        catch(...)
        {
            SEQAN_ASSERT_FAIL("Failed to default construct scheduler!");
        }
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_alignment_scheduler_async)
{
    using namespace seqan;

    using TTask = SchedulerTraits<WavefrontAlignmentScheduler>::TTask;

    WavefrontAlignmentScheduler scheduler{8,2};

    std::vector<uint16_t>  calledIds{0, 0, 0, 0, 0, 0, 0, 0};
    bool isRecycled{false};
    std::mutex mutexSetBool;

    TTask t = [&](uint16_t const id)
    {
        // Now how can we spawn a test to the underlying scheduler.
        auto & _taskScheduler = taskScheduler(scheduler);
        using TInnerTask = typename SchedulerTraits<typename std::decay<decltype(_taskScheduler)>::type>::TTask;

        bool eventState{false};
        std::mutex mutexEvent;
        std::condition_variable event;
        TInnerTask task = [&] ()
        {
            if (calledIds[id] > 0)
            {
                std::lock_guard<std::mutex> lck(mutexSetBool);
                isRecycled = true;
            }

            ++calledIds[id];

            {
                std::lock_guard<std::mutex> lck(mutexEvent);
                eventState = true;
                event.notify_one();
            }
        };
        scheduleTask(_taskScheduler, task);
        // need to wait for the internal task_scheduler.
        {
            std::unique_lock<std::mutex> lck(mutexEvent);
            event.wait(lck, [&]{ return eventState; });
        }
    };

    try
    {
        for (unsigned i = 0; i < 100; ++i)
        {
            scheduleTask(scheduler, t);
        }

        notify(scheduler);
        seqan::wait(scheduler);
    }
    catch (...)
    {
        SEQAN_ASSERT_FAIL("Unexpected exception!");
    }

    auto val = std::accumulate(std::begin(calledIds), std::end(calledIds), 0);
    SEQAN_ASSERT_EQ(val, 100);
    SEQAN_ASSERT(isValid(scheduler));
    SEQAN_ASSERT(isRecycled);
}

namespace test_align_parallel
{

struct test_error : public std::exception
{
    const char* msg;

    explicit test_error( const char* what_arg ) : std::exception(), msg(what_arg)
    {}

    explicit test_error( const std::string& what_arg ) : test_error(what_arg.c_str())
    {}

    virtual const char* what() const noexcept
    {
        return msg;
    }
};

struct RaiiEvent
{
    bool                    eventState{false};
    std::mutex              mutexEvent{};
    std::condition_variable event{};

    inline void wait(unsigned const /*id*/)
    {
        {
            std::unique_lock<std::mutex> lck(mutexEvent);
            event.wait(lck, [&]{ return eventState; });
        }
    }

    inline void notify(unsigned const /*id*/)
    {
        {
            std::lock_guard<std::mutex> lck(mutexEvent);
            eventState = true;
            event.notify_one();
        }
    }
};
}

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_alignment_scheduler_async_with_exception)
{
    using namespace seqan;

    using TTask = SchedulerTraits<WavefrontAlignmentScheduler>::TTask;

    WavefrontAlignmentScheduler scheduler{8,2};

    std::vector<uint16_t>  calledIds{0, 0, 0, 0, 0, 0, 0, 0};
    bool isRecycled{false};
    std::mutex mutexSetBool;

    TTask t = [&] (uint16_t const id)
    {
        // Now how can we spawn a test to the underlying scheduler.
        auto& _taskScheduler = taskScheduler(scheduler);
        using TInnerTask = typename SchedulerTraits<typename std::decay<decltype(_taskScheduler)>::type>::TTask;

        test_align_parallel::RaiiEvent event;
        TInnerTask task = [&]()
        {
            {
                std::lock_guard<std::mutex> lck(mutexSetBool);
                isRecycled = true;
                if (std::accumulate(std::begin(calledIds), std::end(calledIds), 0) == 50)
                {
                    event.notify(id);
                    throw test_align_parallel::test_error("Test");
                }
            }
            ++calledIds[id];
            event.notify(id);
        };

        scheduleTask(_taskScheduler, task);
        // need to wait for the internal task_scheduler.
        event.wait(id);
    };

    try
    {
        for (unsigned i = 0; i < 100; ++i)
        {
            scheduleTask(scheduler, t);
        }

        notify(scheduler);
        seqan::wait(scheduler);
    }
    catch (std::runtime_error & e)
    {
        std::string msg = e.what();
        SEQAN_ASSERT_EQ(msg, "Invalid alignment scheduler!");
    }

    auto exceptVec = getExceptions(scheduler);
    for (unsigned i = 0; i < exceptVec.size(); ++i)
    {
        try
        {
            if (exceptVec[i] != nullptr)
                std::rethrow_exception(exceptVec[i]);
        }
        catch (std::runtime_error const & e)
        {
            std::string msg = e.what();
            SEQAN_ASSERT_EQ(msg, "Invalid Task Scheduler");
        }
        catch (test_align_parallel::test_error const & e)
        {
            std::string msg = e.what();
            SEQAN_ASSERT_EQ(msg, "Test");
        }
        catch (...)
        {
            SEQAN_ASSERT_FAIL("Caught unknown exception!");
        }
    }

    notify(scheduler);
    seqan::wait(scheduler);
    auto val = std::accumulate(std::begin(calledIds), std::end(calledIds), 0);
    SEQAN_ASSERT_EQ(val, 50);
    SEQAN_ASSERT_NOT(isValid(scheduler));
}
