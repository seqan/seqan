// ==========================================================================
//                                 BASIL
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_BASIL_THREAD_SAFE_QUEUE_H_
#define SANDBOX_HERBARIUM_APPS_BASIL_THREAD_SAFE_QUEUE_H_

// A simple thread-safe queue implementation based on std::list<>::splice
// after a tip in a talk by Sean Parent of Adobe.
//
// Uses standard library threading and synchronization primitives together
// with a std::list<> container for implementing a thread-safe queue.  The
// only special thing is that the queue uses std::list<>::splice in the
// locked region to minimize locked time.
//
// Also implements a maximal size and can thus be used as a buffer between
// the elements in a pipeline with limited buffer bloat.
//
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// License: 3-clause BSD

#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <stdexcept>
#include <iostream>
#include <list>
#include <vector>

// Simple thread-safe queue with maximal size based on std::list<>::splice().
//
// Implemented after a tip by Sean Parent at Adobe.

template <typename T>
class ThreadSafeQueue
{
public:
    // Returns whether could push/pop or queue was closed.  Currently, we only implemented the blocking
    // operations.
    enum QueueResult
    {
        OK,
        CLOSED
    };

    // Initialize the queue with a maximal size.
    explicit ThreadSafeQueue(size_t maxSize = 0) : state(State::OPEN), currentSize(0), maxSize(maxSize)
    {}

    // Push v to queue.  Blocks if queue is full.
    void push(T const & v)
    {
        // Create temporary queue.
        decltype(list) tmpList;
        tmpList.push_back(v);

        // Pushing with lock, only using list<>::splice().
        {
            std::unique_lock<std::mutex> lock(mutex);

            // Wait until there is space in the queue.
            while (currentSize == maxSize)
                cvPush.wait(lock);

            // Check that the queue is not closed and thus pushing is allowed.
            if (state == State::CLOSED)
                throw std::runtime_error("Trying to push to a closed queue.");

            // Push to queue.
            currentSize += 1;
            list.splice(list.end(), tmpList, tmpList.begin());

            // Wake up one popping thread.
            if (currentSize == 1u)
                cvPop.notify_one();
        }
    }

    // Push to queue with rvalue reference.
    void push(T && v)
    {
        // Create temporary queue.
        decltype(list) tmpList;
        tmpList.push_back(v);

        // Pushing with lock, only using list<>::splice().
        {
            std::unique_lock<std::mutex> lock(mutex);

            // Wait until there is space in the queue.
            while (currentSize == maxSize)
                cvPush.wait(lock);

            // Check that the queue is not closed and thus pushing is allowed.
            if (state == State::CLOSED)
                throw std::runtime_error("Trying to push to a closed queue.");

            // Push to queue.
            currentSize += 1;
            list.splice(list.end(), tmpList, tmpList.begin());

            // Wake up one popping thread.
            cvPop.notify_one();
        }
    }

    // Pop value from queue and write to v.
    //
    // If this succeeds, OK is returned.  CLOSED is returned if the queue is empty and was closed.
    QueueResult pop(T & v)
    {
        decltype(list) tmpList;

        // Pop into tmpList which is finally written out.
        {
            std::unique_lock<std::mutex> lock(mutex);

            // If there is no item then we wait until there is one.
            while (list.empty() && state != State::CLOSED)
                cvPop.wait(lock);

            // If the queue was closed and the list is empty then we cannot return anything.
            if (list.empty() && state == State::CLOSED)
                return CLOSED;

            // If we reach here then there is an item, get it.
            currentSize -= 1;
            tmpList.splice(tmpList.begin(), list, list.begin());
            // Wake up one pushing thread.
            cvPush.notify_one();
        }

        // Write out value.
        v = tmpList.front();

        return OK;
    }

    // Pushing data is not allowed any more, when the queue is
    void close() noexcept
    {
        std::unique_lock<std::mutex> lock(mutex);
        state = State::CLOSED;

        // Notify all consumers.
        cvPop.notify_all();
    }

private:

    // Whether the queue is running or closed.
    enum class State
    {
        OPEN,
        CLOSED
    };

    // The current state.
    State state;
    // The current size.
    size_t currentSize;
    // The maximal size.
    size_t maxSize;
    // The condition variables to use for pushing/popping.
    std::condition_variable cvPush, cvPop;
    // The mutex for locking the queue.
    std::mutex mutex;
    // The list that the queue is implemented with.
    std::list<T> list;
};

// A simple pipeline.
//
// Limitations:
//  * Can only work on one fixed type.
//  * Creates one thread for each pipeline step for now.
//  * Cannot be copied after begin started.

template <typename T>
class Pipeline
{
public:
    typedef ThreadSafeQueue<T> TQueue;

    explicit Pipeline(int bufferSize = 10) : bufferSize(bufferSize)
    {
        // Create initial queue.
        queues.push_back(std::make_shared<TQueue>(bufferSize));
    }

    // Append a pipeline step.
    //
    // Each step consists of a function that takes a T & as the value and returns a T.
    Pipeline & append(std::function<T(T&, bool)> fun, int bufferSize)
    {
        queues.push_back(std::make_shared<TQueue>(bufferSize));
        TQueue * inQ = queues[queues.size() - 2].get();
        TQueue * outQ = queues[queues.size() - 1].get();
        steps.push_back(
                [fun, inQ, outQ]() {
                    typename TQueue::QueueResult res;
                    T buffer;
                    while ((res = inQ->pop(buffer)) == TQueue::OK)
                        outQ->push(fun(buffer, false));
                    outQ->push(fun(buffer, true));
                    // Error checking and queue closing.
                    if (res == TQueue::CLOSED)
                        outQ->close();
                    else
                        throw std::runtime_error("Problem reading from in queue in pipeline step.");
                });
        return *this;
    }

    // Start to run the pipeline.
    void start()
    {
        if (!threads.empty())
            return;
        for (auto step : steps)
            threads.push_back(std::thread(step));
    }

    // Run pipeline threads to completion.
    void join()
    {
        for (auto & t : threads)
            t.join();
    }

    TQueue & frontQueue()
    {
        return *queues.front();
    }

    TQueue & backQueue()
    {
        return *queues.back();
    }

private:
    // The size the buffers in the queue.
    int bufferSize;
    // The threads to use for the processing.
    std::vector<std::thread> threads;
    // The functors to use for the steps.
    std::vector<std::function<void()>> steps;
    // The queues to use in this pipeline.
    std::vector<std::shared_ptr<TQueue>> queues;
};

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_BASIL_THREAD_SAFE_QUEUE_H_
