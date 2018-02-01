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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_EXECUTOR_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_EXECUTOR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename TResource>
struct WavefrontTaskExecutionPolicy;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Task executor. Manages the execution of single alignment blocks.
template <typename TResource, typename TWavefrontExecutor>
struct WavefrontTaskExecutor
{
    TResource*            _ptrResource{nullptr};
    TWavefrontExecutor *  _ptrWavefrontExecutor{nullptr};

    //NOTE(rrahn) Bug in g++-4.9 prevents us from using as aggregate type.
    WavefrontTaskExecutor() = default;

    WavefrontTaskExecutor(TResource * _resource,
                          TWavefrontExecutor * _wavefrontExecutor) :
            _ptrResource{_resource},
            _ptrWavefrontExecutor(_wavefrontExecutor)
    {}

    inline void operator()()
    {
        WavefrontTaskExecutionPolicy<TResource>::execute(*_ptrResource, *_ptrWavefrontExecutor);
    }
};

// Policy for no SIMD execution.
template <typename ...TArgs>
struct WavefrontTaskExecutionPolicy<WavefrontTask<TArgs...>>
{

    template <typename TResource, typename TWavefrontExecutor>
    inline static void
    execute(TResource & task, TWavefrontExecutor & wavefrontExec)
    {
        using TWaveTaskExec = WavefrontTaskExecutor<TResource, TWavefrontExecutor>;

        executeScalar(task, local(wavefrontExec));
        for (auto succ : successor(task))
        {
            if (succ && decrementRefCount(*succ) == 0)
                spawn(wavefrontExec, TWaveTaskExec{succ, &wavefrontExec});
        }
        if (isLastTask(task))
        {
            notify(*(context(task).ptrEvent));
        }
    }
};

// Policy for SIMD execution.
template <typename TValue, size_t VECTOR_SIZE>
struct WavefrontTaskExecutionPolicy<WavefrontTaskQueue<TValue, VECTOR_SIZE>>
{
    template <typename TResource, typename TWavefrontExecutor>
    inline static void
    execute(TResource & resource, TWavefrontExecutor & wavefrontExec)
    {
        using TWaveTaskExec = WavefrontTaskExecutor<TResource, TWavefrontExecutor>;

        typename TResource::ResultType tasks;
        if (!tryPopTasks(tasks, resource))
            return;

        SEQAN_ASSERT(!empty(tasks));
        if (tasks.size()  == 1)
            executeScalar(*front(tasks), local(wavefrontExec));
        else
            executeSimd(tasks, local(wavefrontExec));

        for (auto task : tasks)
        {
            for (auto succ : successor(*task))
            {
                if (succ && decrementRefCount(*succ) == 0)
                {
                    appendValue(resource, *succ);
                    spawn(wavefrontExec, TWaveTaskExec{&resource, &wavefrontExec});
                }
            }
            if (isLastTask(*task))
            {
                notify(*(context(*task).ptrEvent));
            }
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_EXECUTOR_H_
