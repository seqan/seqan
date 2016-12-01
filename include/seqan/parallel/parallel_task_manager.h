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
// Manager class for tasks.
// ==========================================================================

#ifndef INCLUDE_SEQAN_PARALLEL_PARALLEL_TASK_MANAGER_H_
#define INCLUDE_SEQAN_PARALLEL_PARALLEL_TASK_MANAGER_H_

namespace seqan {
    
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
    
// ----------------------------------------------------------------------------
// Traits TaskTraits
// ----------------------------------------------------------------------------

template <typename TTask>
struct AlignmentTaskTraits : public TaskTraits<TTask>
{

    static void invoke(TTask && task, uint16_t const managerId)
    {
        task(managerId);  // That calls the alignment kernel.
        // In this kernel, we want to use thread local storage through
        // all the threads. How can we do this?
        // In static function we can only access
    }
    
    // We need to implement a round-robin thread pinning strategy
    template <typename TNativeThreadHandle>
    static void pinThread(TNativeThreadHandle && threadHandle)
    {
        // no-op!
    }
};

    
// ----------------------------------------------------------------------------
// Traits TaskTraits
// ----------------------------------------------------------------------------
    
template <typename TTask>
struct TaskTraits
{
    using TTaskQueue = ConcurrentQueue<TTask, Suspendable<>>;
    
    static void invoke(TTask && task, uint16_t const managerId)
    {
        task(managerId);
    }
    
    // We need to implement a round-robin thread pinning strategy
    template <typename TNativeThreadHandle>
    static void pinThread(TNativeThreadHandle && threadHandle)
    {
        // no-op!
    }
};

// ----------------------------------------------------------------------------
// Class TaskManager
// ----------------------------------------------------------------------------
    
template <typename TTask, typename TTaskTraits = TaskTraits<TTask> >
class TaskManager
{
public:
    
    // ----------------------------------------------------------------------------
    // Members
    // ----------------------------------------------------------------------------
    
    typename TTaskTraits::TTaskQueue mTodoQueue;
    std::vector<std::thread>         mHelperThreads;
    
    // ----------------------------------------------------------------------------
    // ManagerTask
    // ----------------------------------------------------------------------------
    
    struct ManagerTask
    {
        uint16_t mId;
        
        operator()()
        {
            while (true)
            {
                TTask task;
                if (!popFront(task, mTodoQueue))
                {
                    return;
                }
                TTaskTraits::invoke(task, mId);
            }
        }
    };
    
    // ----------------------------------------------------------------------------
    // Constructor
    // ----------------------------------------------------------------------------
    
    TaskManager(uint16_t const helperThreadsCount)
    {
        setReaderWriterCount(mTodoQueue, helperThreadsCount, 1);
        for (uint16_t i = 0; i < helperThreadsCount; ++i)
        {
            mHelperThreads.emplace_back(ManagerTask, i);
            // todo: Implement me!
//            TTaskTraits::pinThread(mHelperThreads[i].native_handle());
        }
    }
    
    TaskManager(TaskManager const & /*other*/) = delete;
    TaskManager(TaskManager && /*other*/) = delete;
    
    // ----------------------------------------------------------------------------
    // Destructor
    // ----------------------------------------------------------------------------
    
    ~TaskManager()
    {
        unlockWriting(mTodoQueue);
        for (auto& thread : mHelperThreads)
        {
            if (thread.joinable())
                thread.join();
        }
    }
    
    // ----------------------------------------------------------------------------
    // Member functions
    // ----------------------------------------------------------------------------
    
    TaskManager& operator=(TaskManager /*other*/) = delete;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
    
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_PARALLEL_PARALLEL_TASK_MANAGER_H_
