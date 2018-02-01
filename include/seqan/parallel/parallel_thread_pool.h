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

#ifndef INCLUDE_SEQAN_PARALLEL_PARALLEL_THREAD_POOL_H_
#define INCLUDE_SEQAN_PARALLEL_PARALLEL_THREAD_POOL_H_

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

/*!
 * @class ThreadPool
 * @headerfile <seqan/parallel.h>
 * @brief Handles a pool of concurrent <a href="http://en.cppreference.com/w/cpp/thread/thread">std::threads</a>.
 * @signature class ThreadPool;
 *
 * This is a simple raii-wrapper class to manage a set of <tt>std::threads</tt>.
 */
class ThreadPool
{
public:

    //-------------------------------------------------------------------------
    // Constructor.

    /*!
     * @fn ThreadPool::ThreadPool
     * @brief The constructor.
     * @signature ThreadPool::TThreadPool() = default;
     * @signature ThreadPool::TThreadPool(ThreadPool const &) = delete;
     * @signature ThreadPool::TThreadPool(ThreadPool &&) = delete;
     *
     * Creates a new instance with an empty pool.
     */
    ThreadPool() = default;
    ThreadPool(ThreadPool const &) = delete;
    ThreadPool(ThreadPool &&) = delete;

    ThreadPool& operator=(ThreadPool const &) = delete;
    ThreadPool& operator=(ThreadPool &&) = delete;

    /*!
     * @fn ThreadPool::~ThreadPool
     * @brief The destructor.
     * @signature ThreadPool::~ThreadPool()
     *
     * Safely destroys the thread pool instance by joining all registered threads.
     * This is an implicit barrier for the owning thread.
     *
     * @warning If threads cannot be joined (e.g. dead lock) the destructor will wait forever.
     */
    ~ThreadPool()
    {
        for_each(std::begin(_mPool), std::end(_mPool),
                 [] (auto & t)
                 {
                     if (t.joinable())
                        t.join();
                 });
    }

    //-------------------------------------------------------------------------
    // Private Member Variables.

    std::deque<std::thread> _mPool;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*!
 * @fn ThreadPool#spawn
 * @brief Spawns a new thread and registers it in the thread pool.
 * @headerfile <seqan/parallel.h>
 *
 * @signature void spawn(pool, callable, ...args);
 * @param[in,out] pool The @link ThreadPool @endlink to register the new spawned thread in.
 * @param[in] callable A callable object (e.g. functor or function) that should be executed by the spawned thread.
 * @param     args Any number of arguments passed to the callable object.
 *
 * Emplaces the thread in the <tt>pool</tt> and associates the thread with the callable by passing the <tt>args</tt>
 * to the callable instance.
 *
 * @datarace This function is not thread safe. Concurrent invocations of this function for the same <tt>pool</tt> might
 * result in undefined behavior.
 */
template <typename TCallable, typename ...TArgs>
inline void
spawn(ThreadPool & pool,
      TCallable callable, TArgs && ...args)
{
    pool._mPool.emplace_back(callable, std::forward<TArgs>(args)...);
}

/*!
 * @fn ThreadPool#join
 * @brief Explicit barrier over the pool.
 * @headerfile <seqan/parallel.h>
 *
 * @signature void join(pool);
 * @param[in,out] pool The @link ThreadPool @endlink to be joined.
 *
 * Allows the user to wait for all registered threads to be finished before the calling thread continues its execution.
 *
 * @datarace This function is not thread safe. Concurrent invocations of this function for the same <tt>pool</tt> might
 * result in undefined behavior.
 */
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

/*!
 * @fn ThreadPool#setCpuAffinity
 * @brief Pins the spawned threads in a round-robin fashion.
 * @headerfile <seqan/parallel.h>
 *
 * @signature void setCpuAffinity(pool, cpu, scale);
 * @param[in,out] pool The @link ThreadPool @endlink to pin the threads for.
 * @param[in] cpu The number of the first cpu to be pinned.
 * @param[in] scale A scaling factor to
 *
 * Iterates over all threads registered in the pool and pins each of them to a cpu in a round-robin fashion.
 * Using <tt>cpu</tt> and <tt>scale</tt> the cpu ids for pinning the threads can be configured dynamically.
 * @section Possible implementation:
 * @code{.cpp}
 * cpu * scale % std::thread::hardware_concurrency();
 * @endcode
 *
 * @note This function uses pthread functions on the native thread handle and is only available for linux platforms.
 *       On all other platforms this function is a no-op.
 *
 * @datarace This function is not thread safe. Concurrent invocations of this function for the same <tt>pool</tt> might
 * result in undefined behavior.
 */
#if defined(__linux__)
inline bool
setCpuAffinity(ThreadPool & me, size_t firstCpu = 0, size_t const scale = 1)
{
    SEQAN_ASSERT_GEQ(scale, static_cast<size_t>(1));
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

#endif  // INCLUDE_SEQAN_PARALLEL_PARALLEL_THREAD_POOL_H_
