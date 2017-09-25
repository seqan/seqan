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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_THREAD_POOL_STD_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_THREAD_POOL_STD_H_

#if defined(__linux__)
#include <sched.h>
#endif  // defined(__linux__)

namespace seqan
{
    
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

std::mutex _globalMutexCout;

class TestThread
{
public:

    std::string name;
    unsigned id;
    std::thread thread;

    TestThread(TestThread && moveable) = default;

    template <typename TFunc, typename ...TArgs>
    TestThread(std::string _name, unsigned const _id, TFunc && f, TArgs && ...args) :
        name(std::move(_name)),
        id(_id),
        thread(std::forward<TFunc>(f), std::forward<TArgs>(args)...)
    {
        std::lock_guard<std::mutex> lck(_globalMutexCout);
        std::cout << "Created Thread: " << name << " with id: " << id << '\n';
    }

    bool joinable()
    {
        return thread.joinable();
    }

    void join()
    {
        {
            std::lock_guard<std::mutex> lck(_globalMutexCout);
            std::cout << "Join Thread: " << name << " with id: " << id << '\n';
        }
        thread.join();
    }

    ~TestThread()
    {
        std::lock_guard<std::mutex> lck(_globalMutexCout);
        std::cout << "Killed Thread: " << name << " with id: " << id << '\n';
    }
};

// Simple Thread Pool.
class ThreadPool
{
public:

    //-------------------------------------------------------------------------
    // Constructor.

    ThreadPool() = default;
    ThreadPool(ThreadPool const &) = delete;
    ThreadPool(ThreadPool &&) = delete;

    ThreadPool& operator=(ThreadPool const &) = delete;
    ThreadPool& operator=(ThreadPool &&) = delete;

    ~ThreadPool()
    {
        for_each(std::begin(_mPool), std::end(_mPool),
        [](auto & t)
        {
            if (t.joinable())
                t.join();
        });
    }

    //-------------------------------------------------------------------------
    // Private Member Variables.

    std::vector<std::thread> _mPool;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

//template <typename TCallable, typename ...TArgs>
//inline void
//spawn(ThreadPool & pool,
//      std::string const & str,
//      unsigned id,
//      TCallable callable, TArgs && ...args)
//{
//    pool._mPool.emplace_back(str, id, callable, std::forward<TArgs>(args)...);
//}

//template <typename TCallable, typename ...TArgs>
//inline void
//spawn(ThreadPool & pool,
//      size_t const cpuId,
//      TCallable callable, TArgs && ...args)
//{
//    pool._mPool.emplace_back(callable, std::forward<TArgs>(args)...);
//    cpu_set_t cpuset;
//    CPU_ZERO(&cpuset);
//    CPU_SET(cpuId, &cpuset);
//    int rc = pthread_setaffinity_np(back(pool).native_handle(),
//                                    sizeof(cpu_set_t), &cpuset);
//    if (rc != 0) {
//        std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
//    }
//}

template <typename TCallable, typename ...TArgs>
inline void
spawn(ThreadPool & pool,
      TCallable callable, TArgs && ...args)
{
    pool._mPool.emplace_back(callable, std::forward<TArgs>(args)...);
}

inline void
join(ThreadPool & me)
{
    for_each(std::begin(me._mPool), std::end(me._mPool),
    [](auto & t)
    {
        if (t.joinable())
            t.join();
    });
}

#if defined(__linux__)
inline bool
setCpuAffinity(ThreadPool & me, size_t firstCpu = 0, size_t const scale = 1)
{
    SEQAN_ASSERT_GEQ(scale, 1);
    bool success{true};
    for (auto & t : me._mPool)
    {
        // Create a cpu_set_t object representing a set of CPUs. Clear it and mark
        // only CPU i as set.
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET((firstCpu * scale) % std::thread::hardware_concurrency(), &cpuset);
        int rc = pthread_setaffinity_np(t.native_handle(),
                                        sizeof(cpu_set_t), &cpuset);
        success &= (rc == 0);
        ++firstCpu;
    }
    return success;
}
#else
inline bool
setCpuAffinity(ThreadPool & /*me*/, size_t /*firstCpu = 0*/, size_t const /*scale = 0*/)
{
    return false;
}
#endif  // defined(__linux__)
    
}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_PARALLEL_THREAD_POOL_STD_H_
