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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_EXECUTOR_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_EXECUTOR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Executor class for an alignment task in the wave-front model.
// Stores the scheduler and the thread local storage.
template <typename TScheduler, typename TThreadLocalStore>
struct WavefrontAlignmentExecutor
{
    // Shared data in parallel context.
    TScheduler *                  ptrTaskScheduler{nullptr};
    TThreadLocalStore *           ptrThreadLocal{nullptr};

    //NOTE(rrahn) Bug in g++-4.9 prevents us from using as aggregate type.
    WavefrontAlignmentExecutor() = default;

    WavefrontAlignmentExecutor(TScheduler * _ptrScheduler,
                               TThreadLocalStore * _ptrTls) :
            ptrTaskScheduler{_ptrScheduler},
            ptrThreadLocal(_ptrTls)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Asynchronosly schedule a new alignment job.
template <typename ...TArgs,
          typename TTaskExecutor>
inline void
spawn(WavefrontAlignmentExecutor<TArgs...> & executor,
      TTaskExecutor && taskExec)
{
    SEQAN_ASSERT(executor.ptrTaskScheduler != nullptr);
    scheduleTask(*executor.ptrTaskScheduler, std::forward<TTaskExecutor>(taskExec));
}

// Access thread local storage.
template <typename ...TArgs>
inline auto &
local(WavefrontAlignmentExecutor<TArgs...> & executor)
{
    SEQAN_ASSERT(executor.ptrThreadLocal != nullptr);
    return local(*executor.ptrThreadLocal);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_EXECUTOR_H_
