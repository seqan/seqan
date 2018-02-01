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

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_task_scheduler_construct)
{
    using namespace seqan;

    // We need to be able to construct a thread pool.
    SEQAN_ASSERT(!std::is_default_constructible<WavefrontTaskScheduler>::value);
    SEQAN_ASSERT(!std::is_copy_constructible<WavefrontTaskScheduler>::value);
    SEQAN_ASSERT(!std::is_move_constructible<WavefrontTaskScheduler>::value);
    SEQAN_ASSERT(!std::is_copy_assignable<WavefrontTaskScheduler>::value);
    SEQAN_ASSERT(!std::is_move_assignable<WavefrontTaskScheduler>::value);

    try
    {
        WavefrontTaskScheduler scheduler(4);
        waitForWriters(scheduler);
    }
    catch(...)
    {
        SEQAN_ASSERT_FAIL("Error during construction of scheduler!");
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_task_scheduler_async)
{
    using namespace seqan;

    using TTask = SchedulerTraits<WavefrontTaskScheduler>::TTask;

    WavefrontTaskScheduler scheduler{2, 1};

    bool t1Executed{false};
    TTask t1 = [&]()
    {
        t1Executed = true;
    };
    bool t2Executed{false};
    TTask t2 = [&]()
    {
        t2Executed = true;
    };

    // Register writer.
    lockWriting(scheduler);

    // Schedule jobs for async execution.
    scheduleTask(scheduler, t1);
    scheduleTask(scheduler, t2);

    SEQAN_ASSERT_NOT(t1Executed);
    SEQAN_ASSERT_NOT(t2Executed);

    // Trigger execution.
    waitForWriters(scheduler);

    // Unregister writer.
    unlockWriting(scheduler);

    // Wait until schduler is done with all jobs.
    seqan::wait(scheduler);

    SEQAN_ASSERT(t1Executed);
    SEQAN_ASSERT(t2Executed);
}
