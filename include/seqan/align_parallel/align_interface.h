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
//  Implements the new interface for calling alingment algorithms.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_

namespace seqan {
    
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
    
//struct DPStandardExecution_;
//using DPStandardExecution = Tag<DPStandardExecution_>;
//    
//struct DPBlockedExecution_;
//using DPBlockedExecution = Tag<DPBlockedExecution_>;
//
//struct DPExecutionConfig
//{
//    size_t blockSize;
//};
//    
//template <typename TDPExecutionMode>
//using DPExecutionConfigType = typename std::conditional<std::is_same<TDPExecutionMode, DPBlockedExecution>::value, DPExecutionConfig, Nothing>::type;
//
//struct DPTraits
//{
//    using TAlgorithmType    = GlobalAlignment_<>;
//    using TGapsType         = AffineGaps;
//    using TTracebackType    = TracebackOn<TracebackConfig_<SingleTrace, GapsLeft>>;
//    using TExecutionType    = DPStandardExecution;
//    using TBandType         = BandOff;
//    using TDPScoutStateType = Default;
//    using TOutputType       = ArrayGaps;
//};
//    
//struct DPTraitsBlocked : public DPTraits
//{
//    using TExecutionType = DPBlockedExecution;
//    using TTracebackType = TracebackOff;
//};
//
//// Runtime options
//template <typename TDPTraits>
//struct DPConfig : public DPExecutionConfigType<typename TDPTraits::TExecutionType>
//{
//    // Runtime information, like band parameter, dpState and block size.
//    DPScoutState_<typename TDPTraits::TDPScoutStateType> dpState;
//};

// ============================================================================
// Metafunctions
// ============================================================================
    
//template <typename TObject>
//struct Traits;
//    
//template <typename TTraits>
//struct Traits<DPConfig<TTraits>>
//{
//    using Type = TTraits;
//};
//    
template <typename TCont>
struct IsContainerOfContainer : IsContainerConcept<typename Value<TCont>::Type>
{};

// ============================================================================
// Functions
// ============================================================================
    
template <typename TExecutionPolicy,
          typename TSeqH,
          typename TSeqV,
          typename TScore,
          typename TDelegate>
std::enable_if_t<And<IsContainerOfContainer<TSeqH>, IsContainerOfContainer<TSeqV> >::VALUE>
alignParallel(TExecutionPolicy const & execPolicy,
              TSeqH const & seqH,
              TSeqV const & seqV,
              TScore const & score,  // replaced by the config object later on.
              TDelegate && delegate)
{
    // We need to create the alignment instance here.
    
    // We need to setup the ThreadContext here.
    auto thread_pool = createThreadPool<TTask, TExecutionPolicy>();   // depends on the exec_policy. -> We need the DPTask Type for this.
    // For c++ this also creates the threads that read from the queue. In all other cases it is just the queue.
    
    // We need to generate the Task Scheduler here.
    
    
    // We need to create the Alignment Instances here.
    auto zipSeqs = makeZippedView(seqH, seqV);
    for (auto tuple : zipSeqs)
    {
        std::get<0>(tuple)
        AlignInstanz()
    }
}
    
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_
