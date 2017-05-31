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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_SIMD_DP_TASKS_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_SIMD_DP_TASKS_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TValue,
          size_t VECTOR_SIZE_>
class WavefrontSimdDPTasks
{
public:

    //-------------------------------------------------------------------------
    // Member Types.

    using TQueue = ConcurrentQueue<TValue*>;
    using ResultType = std::vector<TValue*>;
    using ValueType = TValue;

    //-------------------------------------------------------------------------
    // Members.
    static constexpr size_t VECTOR_SIZE{VECTOR_SIZE_};

    TQueue      mQueue;
    std::mutex  mMutexPopQueue;
    bool        mHasNotified{false};

    //-------------------------------------------------------------------------
    // Constructors.

    WavefrontSimdDPTasks()
    {
        lockWriting(mQueue);
        lockReading(mQueue);
    }

    WavefrontSimdDPTasks(WavefrontSimdDPTasks const&) = delete;
    WavefrontSimdDPTasks(WavefrontSimdDPTasks &&) = delete;

    WavefrontSimdDPTasks& operator=(WavefrontSimdDPTasks const &) = delete;
    WavefrontSimdDPTasks& operator=(WavefrontSimdDPTasks &&) = delete;

    ~WavefrontSimdDPTasks()
    {
        if (!mHasNotified)
            unlockWriting(mQueue);
        unlockReading(mQueue);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TValue, size_t VECTOR_SIZE>
inline bool
tryPopTasks(typename WavefrontSimdDPTasks<TValue, VECTOR_SIZE>::ResultType & tasks,
            WavefrontSimdDPTasks<TValue, VECTOR_SIZE> & me)
{
    clear(tasks);
    std::lock_guard<std::mutex> lck(me.mMutexPopQueue);
    if (length(me.mQueue) < WavefrontSimdDPTasks<TValue, VECTOR_SIZE>::VECTOR_SIZE)
    {
        resize(tasks, 1);
        if (!popFront(tasks[0], me.mQueue, Serial()))
        {
            return false;
        }
    }
    else
    {
        for (size_t lane = 0u; lane < VECTOR_SIZE; ++lane)
            tasks.push_back(popFront(me.mQueue, Serial()));
    }
    return true;
}

template <typename TValue, size_t VECTOR_SIZE>
inline void
appendValue(WavefrontSimdDPTasks<TValue, VECTOR_SIZE> & me,
            TValue & newTask)
{
    appendValue(me.mQueue, &newTask);
}

template <typename TValue, size_t VECTOR_SIZE>
inline void
notify(WavefrontSimdDPTasks<TValue, VECTOR_SIZE> & me)
{
    me.mHasNotified = true;
    unlockWriting(me.mQueue);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_SIMD_DP_TASKS_H_
