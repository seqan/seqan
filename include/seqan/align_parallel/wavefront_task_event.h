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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_EVENT_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_EVENT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Event to signal end of one alignment instance.
class WavefrontTaskEvent
{
public:
    std::mutex                  mutexLastTask{};
    std::condition_variable     conditionLastTask{};
    bool                        readyLastTask{false};

    WavefrontTaskEvent() = default;

    WavefrontTaskEvent(WavefrontTaskEvent const &) = delete;
    WavefrontTaskEvent(WavefrontTaskEvent &&) = delete;

    WavefrontTaskEvent& operator=(WavefrontTaskEvent const &) = delete;
    WavefrontTaskEvent& operator=(WavefrontTaskEvent &&) = delete;

    ~WavefrontTaskEvent()
    {
        if (!readyLastTask)
        {
            {
                std::lock_guard<decltype(mutexLastTask)> lck(mutexLastTask);
                readyLastTask = true;
            }
            conditionLastTask.notify_one();
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

inline void
notify(WavefrontTaskEvent & event)
{
    std::lock_guard<decltype(event.mutexLastTask)> lck(event.mutexLastTask);
    event.readyLastTask = true;
    event.conditionLastTask.notify_one();  // We require a strict synchronization between waiting and notifying thread.
}

inline void
wait(WavefrontTaskEvent & event)
{
    std::unique_lock<decltype(event.mutexLastTask)> lck(event.mutexLastTask);
    if (!event.readyLastTask)
        event.conditionLastTask.wait(lck, [&] { return event.readyLastTask; });
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_EVENT_H_
