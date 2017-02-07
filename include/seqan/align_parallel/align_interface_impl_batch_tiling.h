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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_IMPL_BATCH_TILING_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_IMPL_BATCH_TILING_H_

namespace seqan
{
    
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// We can define some interface here.
struct DPTilingStd_;
using DPTilingStd = Tag<DPTilingStd_>;

namespace impl
{

}  // namespace impl

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


namespace impl
{

template <>
struct BatchAlignmentExecutor<DPTilingStd>
{
    template <typename TScore, typename TDPTraits, typename TExecutionTraits,
    typename TSeqBatchH,
    typename TSeqBatchV,
    typename TDelegate>
    inline static void run(DPConfig<TScore, TDPTraits, TExecutionTraits> const & config,
                           TSeqBatchH const & seqBatchH,
                           TSeqBatchV const & seqBatchV,
                           TDelegate && delegate)
    {
    }
};

template <typename TSpec,
          typename TSetH,
          typename TSetV,
          typename TSettings,
          typename TContinuator>
inline void
impl::alignExecBatch(ExecutionPolicy<ParallelWavefront<TSpec>, Serial> const & execPolicy,
                     TSetH const & seqH,
                     TSetV const & seqV,
                     TSettings const & settings,
                     TContinuator && callback)
{
    // Now we need to have a parallel context:


    WavefrontIncubator<TSetH, TSetV, TSettings, TSpec> incubator;
//    Can we now check for the appropriate types that are necessary to define the functionaliy?
//    Actually we don't know until it is to late.

    // Implement me!
    using TDPThreadLocalContext = typename SelectScorePolicy<TDPTraits, TExecutionTraits>::Type;

    impl::ThreadEnumerabelSpecific<DPThreadLocalContext>
    // config must fulfill a certain interface.
    std::vector<std::thread> thread_pool;  // resize to getNumThreads(config);

    // We need the local_thread data.


    ParallelDPAlignmentContext<Parallel<TSpec> >           // Defines the thread pool, the tasks, buffer and thread_local storage and so on ...

    // Normally we would create a single AlingmentInstance and call it.
    // We get the traits for the alingment instance and we need to add a policy for the alignment instance.
    // All these parameters depend on certain trait values.

    AlignmentInstance<ParallelDPAlignmentContext> parallelDpContext;

    // We need to set the number of concurrent Alignments running
    // We need to propagate the dpTaskPool to the underlying structures.


    using TSchedulerTraits = typename Traits<Dynamic<TSchedulingPolicy>>::Type;
    AlignmentScheduler<AICallableRAIIWrapper<TAlignmentInstance>, TSchedulerTraits> scheduler();

    using TSeqH = typename Value<TSequenceH>::Type;
    using TSeqV = typename Value<TSequenceV>::Type;

    auto zip_cont = make_zip(seqSetH, seqSetV);

    for (std::tie<seqH, seqV> tie : zip_cont)
    {
        // Here we call asynchronsly into a concurrent queue.
        // TODO(rrahn): Lazy thread creation mechanism.
        // Is supposed to be a callable.
        // Recycling, becuase we don't allocate more than currently consumed by the list
        schedule(scheduler, impl::AICallableRaiiWrapper<TAlignmentInstance>{seqH, seqV, config, callback});
    }
}

}  // namespace impl
}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_IMPL_BATCH_TILING_H_
