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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_SCHEDULER_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_SCHEDULER_H_

namespace seqan
{
namespace impl
{
namespace dp
{
namespace parallel
{


// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TTaskQueue>
class SharedContext
{
public:
    // Typedefs.
    // ----------------------------------------------------------------------------

    // Member Variables.
    // ----------------------------------------------------------------------------

    uint8_t              mRunningInstances; // Used to syncronize alignment instances.
    std::vector<uint8_t> mIdleInstances;    // Used to syncronize alignment instances.
    TTaskQueue*          mTaskQueue;        // Shared by all alignment instances.

    // Constructors.
    // ----------------------------------------------------------------------------

    // Member Functions.
    // ----------------------------------------------------------------------------
};

template <typename TDPScoreStore, typename TDPTraceStore, typename TExecutionPolicy>
class ParallelAlignScheduler
{
public:

    // Typedefs.
    // ----------------------------------------------------------------------------
    using TDPContext             = ThreadSpecificDPContext<TDPScoreStore, TDPTraceStore>;
    using TThreadSpecificContext = typename ThreadSpecificContext<TDPContext>::Type;

    // Member Variables.
    // ----------------------------------------------------------------------------

    TExecutionPolicy mExecPolicy;

    // Constructors.
    // ----------------------------------------------------------------------------

    ParallelAlignScheduler() = default;

    ParallelAlignScheduler(TExecutionPolicy & execPolicy) : mExecPolicy(execPolicy)
    {}

    // Member Functions.
    // ----------------------------------------------------------------------------

    // Implements the task scheduler
    // Implements the thread local store
    // Implements the scheduling of multiple alignment jobs.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TDPScoreStore, typename TDPTraceStore, typename TExecutionPolicy>
inline void
init(ParallelAlignScheduler<TDPScoreStore, TDPTraceStore, TExecutionPolicy> & scheduler)
{

}

template <typename TExecutionPolicy, typename TAlignInstances>
execute(ParallelAlignScheduler<TExecutionPolicy> & execPolicy,
        TAlignInstances & alignInstances)
{
    using TAlignInstance = typename Value<TAlignInstances>::Type;
    using TInstanceTrait = typename Traits<TAlignInstances>::Type;

    // Initialize scheduler:


    while (!empty(alignInstances))
    {
        auto& instance = popFront(alignInstances);
        SEQAN_ASSERT(instance != nullptr);
        std::thread(*instance, std::ref()).detach();  // Invocation must run in detached thread!
        ++sharedConfig.runningInstances;
        auto& sharedConfig = getSharedConfig(instance);
        {
            std::unique_lock<decltype(sharedConfig.schedulerLock)> lk(sharedConfig.schedulerLock);
            cond.wait_for(lk, []{ sharedConfig.runningInstances < TInstanceTrait::MAX_PARALLEL_INSTANCES; } );
        }
    }

    // Finally wait for running jobs.
    std::unique_lock<decltype(sharedConfig.schedulerLock)> lk(sharedConfig.schedulerLock);
    cond.wait_for(lk, []{ sharedConfig.runningInstances < TInstanceTrait::MAX_PARALLEL_INSTANCES; } );
}

}  // namespace impl::dp::parallel
} // namespace impl::dp
} // namespace impl

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_SCHEDULER_H_
