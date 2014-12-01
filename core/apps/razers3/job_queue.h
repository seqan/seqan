// Simple lock-based queuing system.

#ifndef SEQAN_PARALLEL_JOB_QUEUE_H_
#define SEQAN_PARALLEL_JOB_QUEUE_H_

#include <omp.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// TODO(holtgrew): Specializations?

// Simple job deque with atomic access through OpenMP locks.
template <typename TJob>
class JobQueue
{
public:
    std::deque<TJob> queue_;
    omp_lock_t lock_;

    JobQueue()
    {
        omp_init_lock(&lock_);
    }

    ~JobQueue()
    {
        omp_destroy_lock(&lock_);
    }

};

template <typename TJob, typename TSpec>
class ThreadLocalStorage;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

// Class JobQueue.

template <typename TJob>
inline
size_t
length(JobQueue<TJob> & queue)
{
    omp_set_lock(&queue.lock_);
    return queue.queue_.size();

    omp_unset_lock(&queue.lock_);
}

template <typename TJob>
inline
void
pushFront(JobQueue<TJob> & queue, TJob const & job)
{
    omp_set_lock(&queue.lock_);
    queue.queue_.push_front(job);
    omp_unset_lock(&queue.lock_);
}

template <typename TJob>
inline
bool
popFront(TJob & job, JobQueue<TJob> & queue)
{
    omp_set_lock(&queue.lock_);
    if (queue.queue_.size() == 0u)
    {
        omp_unset_lock(&queue.lock_);
        return false;
    }
    job = queue.queue_.front();
    queue.queue_.pop_front();
    omp_unset_lock(&queue.lock_);
    return true;
}

template <typename TJob>
inline
void
pushBack(JobQueue<TJob> & queue, TJob const & job)
{
    omp_set_lock(&queue.lock_);
    queue.queue_.push_back(job);
    omp_unset_lock(&queue.lock_);
}

template <typename TJob>
inline
bool
popBack(TJob & job, JobQueue<TJob> & queue)
{
    omp_set_lock(&queue.lock_);
    if (queue.queue_.size() == 0u)
    {
        omp_unset_lock(&queue.lock_);
        return false;
    }
    job = queue.queue_.back();
    queue.queue_.pop_back();
    omp_unset_lock(&queue.lock_);
    return true;
}

// ThreadLocalStorage interface.

// template <typename TJob, typename TSpec>
// JobQueue<TJob> & jobQueue(ThreadLocalStorage<TJob, TSpec> & tls);

// template <typename TJob>
// bool isWorking(ThreadLocalStorage<TJob> & tls);

// template <typename TJob>
// void
// setWorking(ThreadLocalStorage<TJob> & tls, bool b);

template <typename TJob, typename TSpec>
inline
bool
work(ThreadLocalStorage<TJob, TSpec> & tls)
{
    // fprintf(stderr, "Thread %d tries to work...\n", omp_get_thread_num());
    TJob job;
    if (!popFront(job, jobQueue(tls)))
        return false;

    work(tls, job);
    return true;
}

template <typename TJob, typename TSpec>
inline
bool
stealWork(TJob & job, ThreadLocalStorage<TJob, TSpec> & tls)
{
    return popBack(job, jobQueue(tls));
}

// Load balancing / work stealing algorithm.

struct StealOne_;
typedef Tag<StealOne_> StealOne;

template <typename TJob, typename TSpec>
void work(String<ThreadLocalStorage<TJob, TSpec> > & threadLocalStorages, StealOne const &)
{
    // Initialization.
    int maxThreads = omp_get_max_threads();
    volatile int workingCount = maxThreads;
    SEQAN_OMP_PRAGMA(parallel)
    {
        setWorking(threadLocalStorages[omp_get_thread_num()], true);
    }

    // Work loop.
    SEQAN_OMP_PRAGMA(parallel)
    {
        int p = omp_get_thread_num();
        unsigned x = 73 * p;  // LCG pseudorandomess :-O
        while (workingCount > 0)
        {
            // Try to steal work if not working.
            if (!isWorking(threadLocalStorages[p]) && maxThreads > 1)
            {
                // Select thread to steal from.
                x = 1664525 * x + 1013904223;
                int targetId = x % (maxThreads - 1);
                targetId += targetId >= p;

                // Try to steal a job and process it if successful.
                TJob stolenJob;
                if (stealWork(stolenJob, threadLocalStorages[targetId]))
                {
                    // Stealing was successful, become active again.
                    SEQAN_OMP_PRAGMA(atomic)
                    workingCount += 1;
                    setWorking(threadLocalStorages[p], true);
                    SEQAN_OMP_PRAGMA(flush(workingCount))

                    // fprintf(stderr, "%d could steal\n", p);
                    work(threadLocalStorages[p], stolenJob);
                }
                continue;  // Next iteration.
            }

            bool res = work(threadLocalStorages[p]);
            if (!res)
            {
                // No more work, could steal some in next
                // iteration and become active again, though.
                // fprintf(stderr, "Thread %d is done for now.\n", omp_get_thread_num());
                SEQAN_OMP_PRAGMA(atomic)
                workingCount -= 1;
                setWorking(threadLocalStorages[p], false);
                SEQAN_OMP_PRAGMA(flush(workingCount))
            }
        }
    }

    // Sanity check: All queues empty? No thread working?
    for (int i = 0; i < omp_get_max_threads(); ++i)
    {
        if (isWorking(threadLocalStorages[i]))
        {
            std::cerr << "ERROR: Thread " << i << " still working!" << std::endl;
        }
        if (length(jobQueue(threadLocalStorages[i])) > 0u)
        {
            std::cerr << "ERROR: Queue of thread " << i << " not empty!" << std::endl;
        }
    }
}

#endif  // #ifndef SEQAN_PARALLEL_JOB_QUEUE_H_
