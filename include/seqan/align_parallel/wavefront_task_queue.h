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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_WAVEFRONT_TASK_QUEUE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_WAVEFRONT_TASK_QUEUE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Central task queue used in simd mode to gather multiple blocks to gather full simd registers.
template <typename TValue,
          size_t VECTOR_SIZE_>
class WavefrontTaskQueue
{
public:


    // Member Types.
    using TQueue = ConcurrentQueue<TValue*>;
    using ResultType = std::vector<TValue*>;
    using ValueType = TValue;

    // Members.
    static constexpr size_t VECTOR_SIZE{VECTOR_SIZE_};

    TQueue      queue;
    std::mutex  mutexPopQueue;
    bool        hasNotified{false};

    // Constructors.
    WavefrontTaskQueue()
    {
        lockWriting(queue);
        lockReading(queue);
    }

    WavefrontTaskQueue(WavefrontTaskQueue const&) = delete;
    WavefrontTaskQueue(WavefrontTaskQueue &&) = delete;

    WavefrontTaskQueue& operator=(WavefrontTaskQueue const &) = delete;
    WavefrontTaskQueue& operator=(WavefrontTaskQueue &&) = delete;

    ~WavefrontTaskQueue()
    {
        if (!hasNotified)
            unlockWriting(queue);
        unlockReading(queue);
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
tryPopTasks(typename WavefrontTaskQueue<TValue, VECTOR_SIZE>::ResultType & tasks,
            WavefrontTaskQueue<TValue, VECTOR_SIZE> & me)
{
    clear(tasks);
    std::lock_guard<std::mutex> lck(me.mutexPopQueue);
    if (length(me.queue) < WavefrontTaskQueue<TValue, VECTOR_SIZE>::VECTOR_SIZE)
    {
        resize(tasks, 1);
        if (!popFront(tasks[0], me.queue, Serial()))
        {
            return false;
        }
    }
    else
    {
        for (size_t lane = 0u; lane < VECTOR_SIZE; ++lane)
            tasks.push_back(popFront(me.queue, Serial()));
    }
    return true;
}

template <typename TValue, size_t VECTOR_SIZE>
inline void
appendValue(WavefrontTaskQueue<TValue, VECTOR_SIZE> & me,
            TValue & newTask)
{
    appendValue(me.queue, &newTask);
}

template <typename TValue, size_t VECTOR_SIZE>
inline void
notify(WavefrontTaskQueue<TValue, VECTOR_SIZE> & me)
{
    me.hasNotified = true;
    unlockWriting(me.queue);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_WAVEFRONT_TASK_QUEUE_H_
