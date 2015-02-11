// Simple lock-based queuing system.

#ifndef RAZERS_PARALLEL_JOB_QUEUE_H_
#define RAZERS_PARALLEL_JOB_QUEUE_H_

#define SEQAN_PARALLEL_DEBUG 1
#ifndef SEQAN_PARALLEL_DEBUG
#define SEQAN_PARALLEL_DEBUG 0
#endif  // #ifndef SEQAN_PARALLEL_DEBUG


#include <omp.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// TODO(holtgrew): Specializations?

template <typename TJob, typename TSpec>
class TaskQueue;

struct OmpLock_;
typedef Tag<OmpLock_> OmpLock;

// Simple job deque with atomic access through OpenMP locks.
template <typename TJob>
class TaskQueue<TJob, OmpLock>
{
public:
    std::deque<TJob> queue_;
    omp_lock_t lock_;

    TaskQueue()
    {
        omp_init_lock(&lock_);
    }

    ~TaskQueue()
    {
        omp_destroy_lock(&lock_);
    }

};

// template <typename TJob, typename TSpec>
// class ThreadLocalStorage;

template <typename TSpec>
class Job;

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename T>
struct JobQueue;

// ===========================================================================
// Functions
// ===========================================================================

template <typename T>
inline
bool
predicateTrue(T const &)
{
    return true;
}

// Class TaskQueue.

template <typename TJob>
inline
size_t
length(TaskQueue<TJob, OmpLock> & queue)
{
    omp_set_lock(&queue.lock_);
    size_t res = queue.queue_.size();
    omp_unset_lock(&queue.lock_);
    return res;
}

template <typename TJob>
inline
void
pushFront(TaskQueue<TJob, OmpLock> & queue, TJob const & job)
{
    omp_set_lock(&queue.lock_);
    queue.queue_.push_front(job);
    omp_unset_lock(&queue.lock_);
}

template <typename TJob>
inline
void
pushFront(TaskQueue<TJob, OmpLock> & queue, String<TJob> const & jobs)
{
    typedef typename Iterator<String<TJob> const>::Type TIterator;

    omp_set_lock(&queue.lock_);
    for (TIterator it = begin(jobs); it != end(jobs); ++it)
        queue.queue_.push_front(*it);
    omp_unset_lock(&queue.lock_);
}

template <typename TJob>
inline
bool
popFront(TJob & job, TaskQueue<TJob, OmpLock> & queue)
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

template <typename TJob, typename TPredicate>
inline
bool
popFront(TJob & job, TaskQueue<TJob, OmpLock> & queue, TPredicate const & predicate)
{
    omp_set_lock(&queue.lock_);
    if (queue.queue_.size() == 0u)
    {
        omp_unset_lock(&queue.lock_);
        return false;
    }
    job = queue.queue_.front();
    if (!predicate(job))
    {
        omp_unset_lock(&queue.lock_);
        return false;
    }
    queue.queue_.pop_front();
    omp_unset_lock(&queue.lock_);
    return true;
}

template <typename TJob>
inline
void
pushBack(TaskQueue<TJob, OmpLock> & queue, TJob const & job)
{
    omp_set_lock(&queue.lock_);
    queue.queue_.push_back(job);
    omp_unset_lock(&queue.lock_);
}

template <typename TJob>
inline
void
pushBack(TaskQueue<TJob, OmpLock> & queue, String<TJob> const & jobs)
{
    typedef typename Iterator<String<TJob> const>::Type TIterator;

    omp_set_lock(&queue.lock_);
    for (TIterator it = begin(jobs); it != end(jobs); ++it)
        queue.queue_.push_back(*it);
    omp_unset_lock(&queue.lock_);
}

template <typename TJob>
inline
bool
popBack(TJob & job, TaskQueue<TJob, OmpLock> & queue)
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

template <typename TJob, typename TPredicate>
inline
bool
popBack(TJob & job, TaskQueue<TJob, OmpLock> & queue, TPredicate const & predicate)
{
    omp_set_lock(&queue.lock_);
    if (queue.queue_.size() == 0u)
    {
        omp_unset_lock(&queue.lock_);
        return false;
    }
    job = queue.queue_.back();
    if (!predicate(job))
    {
        omp_unset_lock(&queue.lock_);
        return false;
    }
    queue.queue_.pop_back();
    omp_unset_lock(&queue.lock_);
    return true;
}

// ThreadLocalStorage interface.

// template <typename TJob, typename TSpec>
// typename JobQueue<ThreadLocalStorage<TJob, TSpec> >::Type & jobQueue(ThreadLocalStorage<TJob, TSpec> & tls);

// template <typename TJob>
// bool isWorking(ThreadLocalStorage<TJob> & tls);

// template <typename TJob>
// void
// setWorking(ThreadLocalStorage<TJob> & tls, bool b);

// template <typename TJob, typename TSpec, typename TPredicate, typename TThreadId>
// inline
// bool
// work(ThreadLocalStorage<TJob, TSpec> & tls, TPredicate & predicate, TThreadId threadId)
// {
//     // fprintf(stderr, "Thread %d tries to work...\n", omp_get_thread_num());
//     TJob job;
//     if (!popFront(job, jobQueue(tls), predicate))
//         return false;
//     work(tls, job, threadId, threadId);
//     return true;
// }

// template <typename TJob, typename TSpec, typename TThreadId>
// inline
// bool
// work(ThreadLocalStorage<TJob, TSpec> & tls, TThreadId threadId)
// {
//     return work(tls, predicateTrue<TJob>, threadId);
// }

// template <typename TJob, typename TSpec, typename TPredicate>
// inline
// bool
// stealWork(TJob & job, ThreadLocalStorage<TJob, TSpec> & tls, TPredicate const & predicate)
// {
//     return popBack(job, jobQueue(tls), predicate);
// }

// // Load balancing / work stealing algorithm.

// struct StealOne_;
// typedef Tag<StealOne_> StealOne;

// // Main loop for parallel work stealing with multiple jobs.
// template <typename TJob, typename TSpec, typename TPredicate>
// void workInParallel(String<ThreadLocalStorage<TJob, TSpec> > & threadLocalStorages, TPredicate const & predicate, StealOne const &)
// {
//     // Initialization.
//     int maxThreads = omp_get_max_threads();
//     volatile int workingCount = maxThreads;
//     SEQAN_OMP_PRAGMA(parallel)
//     {
//         setWorking(threadLocalStorages[omp_get_thread_num()], true);
//     }

//     // Work loop.
//     SEQAN_OMP_PRAGMA(parallel)
//     {
//         int p = omp_get_thread_num();
//         unsigned x = 73 * p;  // LCG pseudorandomess :-O
//         while (workingCount > 0) {
//             // Try to steal work if not working.
//             if (!isWorking(threadLocalStorages[p]) && maxThreads > 1) {
//                 // Select thread to steal from.
//                 x = 1664525 * x + 1013904223;
//                 int targetId = x % (maxThreads - 1);
//                 targetId += targetId >= p;

//                 // Try to steal a job and process it if successful.
//                 TJob stolenJob;
//                 if (stealWork(stolenJob, threadLocalStorages[targetId], predicate)) {
//                     //fprintf(stderr, "Thread %d stole from %d\n", p, targetId);
//                     // Stealing was successful, become active again.
//                     SEQAN_OMP_PRAGMA(atomic)
//                     workingCount += 1;
//                     setWorking(threadLocalStorages[p], true);
//                     SEQAN_OMP_PRAGMA(flush(workingCount))

//                     // fprintf(stderr, "%d could steal\n", p);
//                     work(threadLocalStorages[p], stolenJob, p, targetId);
//                 }
//                 continue;  // Next iteration.
//             }

//             bool res = work(threadLocalStorages[p], p);
//             if (!res) {
//                 // No more work, could steal some in next
//                 // iteration and become active again, though.
//                 // fprintf(stderr, "Thread %d is done for now.\n", omp_get_thread_num());
//                 SEQAN_OMP_PRAGMA(atomic)
//                 workingCount -= 1;
//                 setWorking(threadLocalStorages[p], false);
//                 SEQAN_OMP_PRAGMA(flush(workingCount))
//             }
//         }
//     }

//     // Sanity check: All queues empty? No thread working?
//     for (int i = 0; i < omp_get_max_threads(); ++i) {
//         if (isWorking(threadLocalStorages[i])) {
//             std::cerr << "ERROR: Thread " << i << " still working!" << std::endl;
//         }
//         if (length(jobQueue(threadLocalStorages[i])) > 0u) {
//             std::cerr << "ERROR: Queue of thread " << i << " not empty!" << std::endl;
//         }
//     }
// }

// template <typename TJob, typename TSpec>
// void workInParallel(String<ThreadLocalStorage<TJob, TSpec> > & threadLocalStorages, StealOne const & tag)
// {
//     work(threadLocalStorages, predicateTrue<TJob>, tag);
// }

}  // namespace seqan

#endif  // #ifndef RAZERS_PARALLEL_JOB_QUEUE_H_
