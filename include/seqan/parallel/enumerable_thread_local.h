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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

struct EnumerableThreadLocalIterSpec_;
using EnumerableThreadLocalIterSpec = Tag<EnumerableThreadLocalIterSpec_>;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class SimpleThreadLocalManager
 * @headerfile <seqan/parallel.h>
 * @brief Internal policy used to manage thread specific storage.
 * @signature struct SimpleThreadLocalManager;
 *
 * Uses a shared mutex to synchronize concurrent access to the storage buffer (multiple read, single write).
 * This policy checks if the thread-id was already registered and if not creates a new storage instance for the
 * given thread-id.
 *
 * @see CountingThreadLocalManager
 */
struct SimpleThreadLocalManager
{
    std::shared_timed_mutex  _mutex{};        // TODO(rrahn): Replace by shared_mutex when c++17 is available.

    /*!
     * @fn SimpleThreadLocalManager::local
     * @brief Implements the specific <tt>local</tt> policy.
     * @headerfile <seqan/parallel.h>
     *
     * @signature auto local(map, init);
     * @param map The buffer provided by the @link EnumerableThreadLocal @endlink instance to the hold the thread
     *            local storage. Must satisfy the <tt>std::map</tt> interface.
     * @param init The value used to initialize created storage.
     *
     * @datarace thread-safe.
     * @return auto A tuple containing a lvalue reference to the associated storage and a boolean indicating first time
     * access of the given thread-id.
     */
    template <typename TResourceMap, typename TValue>
    inline auto &
    local(TResourceMap & map, TValue const & initValue, bool & exists)
    {
        decltype(map.find(std::this_thread::get_id())) elemIt;

        { // try to read
            std::shared_lock<decltype(_mutex)> read_lck(_mutex);
            elemIt = map.find(std::this_thread::get_id());
            exists = elemIt != map.end();
        }

        if (!exists)
        {
            {  // Create new entry.
                std::unique_lock<decltype(_mutex)> write_lck(_mutex);
                std::tie(elemIt, exists) = map.emplace(std::this_thread::get_id(), initValue);
            }

            SEQAN_ASSERT(exists);
            exists = false;  // Notify that element was added for the first time.
        }
        return elemIt->second;
    }
};

/*!
 * @class CountingThreadLocalManager
 * @headerfile <seqan/parallel.h>
 * @brief Internal policy used to manage thread specific storage for a limited number of threads.
 * @signature struct SimpleThreadLocalManager;
 *
 * Uses a shared mutex to synchronize concurrent access to the storage buffer (multiple read, single write).
 * This policy checks if the thread-id was already registered and if not creates a new storage instance for the
 * given thread-id.
 * In addition maintains an atomic counter to only check for existing values as long as the counter is not 0.
 * This can be used if the number of threads registering at the @link EnumerableThreadLocal @endlink instance,
 * is known beforehand. In this case, as soon as all threads have been registered and obtained their local storage,
 * no further synchronization is necessary.
 *
 * @see SimpleThreadLocalManager
 */
struct CountingThreadLocalManager
{
    std::atomic<size_t>      _count{0};
    std::shared_timed_mutex  _mutex{};        // TODO(rrahn): Replace by shared_mutex when c++17 is available.

    /*!
     * @fn CountingThreadLocalManager::local
     * @brief Implements the specific <tt>local</tt> policy.
     * @headerfile <seqan/parallel.h>
     *
     * @signature auto local(map, init);
     * @param map The buffer provided by the @link EnumerableThreadLocal @endlink instance to the hold the thread
     *            local storage. Must satisfy the <tt>std::map</tt> interface.
     * @param init The value used to initialize created storage.
     *
     * @datarace thread-safe.
     * @return auto A tuple containing a lvalue reference to the associated storage and a boolean indicating first time
     * access of the given thread-id.
     */
    template <typename TResourceMap, typename TValue>
    inline auto &
    local(TResourceMap & map, TValue const & initValue, bool & exists)
    {
        if (_count.load(std::memory_order_relaxed) == 0)
            return map.find(std::this_thread::get_id())->second;

        decltype(map.find(std::this_thread::get_id())) elemIt;

        { // try to read
            std::shared_lock<decltype(_mutex)> read_lck(_mutex);
            elemIt = map.find(std::this_thread::get_id());
            exists = elemIt != map.end();
        }

        if (!exists)
        {
            {  // Create new entry.
                std::unique_lock<decltype(_mutex)> write_lck(_mutex);
                std::tie(elemIt, exists) = map.emplace(std::this_thread::get_id(), initValue);
            }

            --_count;
            SEQAN_ASSERT(exists);
            exists = false;  // Notify that element was added for the first time.
        }
        return elemIt->second;
    }
};

// ----------------------------------------------------------------------------
//  Function setCount()
// ----------------------------------------------------------------------------

inline void
setCount(CountingThreadLocalManager & mngr, size_t const count)
{
    mngr._count.store(count, std::memory_order_relaxed);
}

/*!
 * @class EnumerableThreadLocal
 * @headerfile <seqan/parallel.h>
 * @brief Manages thread local storage.
 * @signature template <typename TValue, typename TManager, typename TSpec>
 *            class ThreadPool;
 * @tparam TValue The type of the stored value.
 * @tparam TManager A policy used to manage the thread local storage. Defaults to @link SimpleThreadLocalManager @endlink.
 * @tparam TSpec Specialization of the <tt>EnumerableThreadLocal</tt> class. Defaults to <tt>void</tt>.
 *
 * The enumerable thread local class can be used to manage thread local storage using a map as
 * internal buffer. The class offers an iterator interface, such that the thread specific values can be
 * enumerated allowing to apply reduce operations at the end of the parallel execution.
 * Creating thread local storage happens via a lazy evaluation using the <tt>local</tt> function.
 * If a thread, identified by its <a href="http://en.cppreference.com/w/cpp/thread/get_id">thread id</a>,
 * requests storage for the first time a new thread specific storage will be created and a lvalue reference pointing
 * to this storage is returned.
 * If the thread id was already registered, then a lvalue reference to the associated storage will be returned.
 * The access to the <tt>local</tt> function is thread safe and can be called concurrently.
 *
 * A thread local manager can be selected via template argument. This manager ensures safe access from concurrent
 * invocations. The manager can further be replaced to change the behavior of storing the thread local data.
 *
 * @see SimpleThreadLocalManager
 * @see CountingThreadLocalManager
 */
template <typename TValue, typename TManager = SimpleThreadLocalManager, typename TSpec = void>
class EnumerableThreadLocal
{
public:

    //-------------------------------------------------------------------------
    // Member types.

    using TMap           = std::unordered_map<std::thread::id, TValue>;
    using TIterator      = Iter<EnumerableThreadLocal, EnumerableThreadLocalIterSpec>;
    using TConstIterator = Iter<EnumerableThreadLocal const, EnumerableThreadLocalIterSpec>;

    //-------------------------------------------------------------------------
    // Private Members.

    TMap                     _map{};
    TValue                   _initValue{};
    TManager                 _manager{};

    //-------------------------------------------------------------------------
    // Constructor.

    /*!
     * @fn EnumerableThreadLocal::EnumerableThreadLocal
     * @brief Constructing an instance of this class.
     * @signature EnumerableThreadLocal::EnumerableThreadLocal() = default;
     * @signature EnumerableThreadLocal::EnumerableThreadLocal(TValue init);
     *
     * @param[in] init An optional value used to initialize the newly created storage.
     * @note The class is not @link CopyConstructibleConcept copy constructible @endlink and not
     * <a href="http://en.cppreference.com/w/cpp/concept/MoveConstructible">move constructible</a>.
     */
    EnumerableThreadLocal() = default;

    explicit EnumerableThreadLocal(TValue initValue) : _initValue(std::move(initValue))
    {}

    EnumerableThreadLocal(EnumerableThreadLocal const &) = delete;
    EnumerableThreadLocal(EnumerableThreadLocal &&) = delete;

    //-------------------------------------------------------------------------
    // Member Functions.

    EnumerableThreadLocal & operator=(EnumerableThreadLocal const &) = delete;
    EnumerableThreadLocal & operator=(EnumerableThreadLocal &&) = delete;
};

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn EnumerableThreadLocal#Iterator
 * @brief Constructing an instance of this class.
 * @headerfile <seqan/parallel.h>
 * @signature Iterator<TEnumerableThreadLocal>::Type;
 * @tparam  TEnumerableThreadLocal The type of the enumerable thread local class.
 * @return Type The type of the iterator.
 */
template <typename TValue, typename TManager, typename TSpec,
          typename TIterTag>
struct Iterator<EnumerableThreadLocal<TValue, TManager, TSpec>, TIterTag>
{
    using Type = typename EnumerableThreadLocal<TValue, TManager, TSpec>::TIterator;
};

template <typename TValue, typename TManager, typename TSpec,
          typename TIterTag>
struct Iterator<EnumerableThreadLocal<TValue, TManager, TSpec> const, TIterTag>
{
    using Type = typename EnumerableThreadLocal<TValue, TManager, TSpec>::TConstIterator;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
//  Function storageManager()
// ----------------------------------------------------------------------------

/*!
 * @fn EnumerableThreadLocal#storageManager
 * @brief Constructing an instance of this class.
 * @headerfile <seqan/parallel.h>
 * @signature TManager & storageManager(etl);
 * @param[in] etl An instance of @link EnumerableThreadLocal @endlink.
 * @return TManager& A lvalue reference to the associated storage manager.
 * @datarace not thread-safe.
 */
template <typename TValue, typename TManager, typename TSpec>
inline TManager &
storageManager(EnumerableThreadLocal<TValue, TManager, TSpec> & me)
{
    return me._manager;
}

template <typename TValue, typename TManager, typename TSpec>
inline TManager const &
storageManager(EnumerableThreadLocal<TValue, TManager, TSpec> const & me)
{
    return me._manager;
}

// ----------------------------------------------------------------------------
//  Function local()
// ----------------------------------------------------------------------------

/*!
 * @fn EnumerableThreadLocal#local
 * @brief Constructing an instance of this class.
 * @headerfile <seqan/parallel.h>
 *
 * @signature auto & local(etl);
 * @signature auto & local(etl, b);
 * @param[in,out] etl An instance of @link EnumerableThreadLocal @endlink.
 * @param[in,out] b A boolean to indicate successful creation.
 *
 * @return auto& A lvalue reference to the associated thread specific storage.
 * @datarace thread-safe. Concurrent invocations of this function are synchronized via the storage manager.
 *
 * Calls the internal <tt>local</tt> member function of the associated storage manager. If the thread-id was used
 * for the first time <tt>b</tt> will be set to <tt>true</tt> to indicate successful creation of the storage.
 * If the storage was already created for the given thread-id, then <tt>b</tt> is set to <tt>false</tt>.
 */
template <typename TValue, typename TManager, typename TSpec>
inline auto&
local(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
      bool & exists)
{
    return  me._manager.local(me._map, me._initValue, exists);
}

template <typename TValue, typename TManager, typename TSpec>
inline auto&
local(EnumerableThreadLocal<TValue, TManager, TSpec> & me)
{
    bool exists{true};
    return local(me, exists); // Double indirection to to iterator and pointer to therad_local storage.
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

/*!
 * @fn EnumerableThreadLocal#begin
 * @brief Returns a bidirectional iterator to the thread specific storage.
 * @headerfile <seqan/parallel.h>
 *
 * @signature TIteator begin(etl, tag);
 * @param[in] etl An instance of @link EnumerableThreadLocal @endlink.
 * @param[in] tag A tag to choose the type of the iterator. One of @link ContainerIteratorTags @endlink.
 *
 * @return TIteator Iterator to the begin of the thread stores.
 * @datarace thread-safe.
 */
template <typename TValue, typename TManager, typename TSpec,
          typename TIterSpec>
inline auto
begin(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
      Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TManager, TSpec>::TIterator{me};
}

template <typename TValue, typename TManager, typename TSpec,
          typename TIterSpec>
inline auto
begin(EnumerableThreadLocal<TValue, TManager, TSpec> const & me,
      Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TManager, TSpec>::TConstIterator{me};
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

/*!
 * @fn EnumerableThreadLocal#end
 * @brief Returns a bidirectional iterator to the thread specific storage.
 * @headerfile <seqan/parallel.h>
 *
 * @signature TIteator end(etl, tag);
 * @param[in] etl An instance of @link EnumerableThreadLocal @endlink.
 * @param[in] tag A tag to choose the type of the iterator. One of @link ContainerIteratorTags @endlink.
 *
 * @return TIteator Iterator to the end of the thread stores.
 * @datarace thread-safe.
 */
template <typename TValue, typename TManager, typename TSpec,
          typename TIterSpec>
inline auto
end(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
    Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TManager, TSpec>::TIterator{me._map.end()};
}

template <typename TValue, typename TManager, typename TSpec,
          typename TIterSpec>
inline auto
end(EnumerableThreadLocal<TValue, TManager, TSpec> const & me,
    Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TManager, TSpec>::TConstIterator{me._map.cend()};
}

// ----------------------------------------------------------------------------
// Function combine()
// ----------------------------------------------------------------------------

/*!
 * @fn EnumerableThreadLocal#combineEach
 * @brief Enumerates thread local stores and applies an unary functor for each store.
 * @headerfile <seqan/parallel.h>
 *
 * @signature void combineEach(etl, f);
 * @param[in] etl An instance of @link EnumerableThreadLocal @endlink.
 * @param[in] f An unary functor called on each thread local storage.
 *
 * @datarace not thread-safe.
 * @see EnumerableThreadLocal#combine
 *  @section Possible implementation:
 * @code{.cpp}
 * std::for_each(begin(etl), end(etl), f);
 * @endcode
 */
template <typename TValue, typename TManager, typename TSpec,
          typename TUnaryCombine>
inline void
combineEach(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
            TUnaryCombine && fUnaryCombine)
{
    std::shared_lock<decltype(storageManager(me)._mutex)> read_lck(storageManager(me)._mutex);
    std::for_each(begin(me), end(me), std::forward<TUnaryCombine>(fUnaryCombine));
}

/*!
 * @fn EnumerableThreadLocal#combine
 * @brief Enumerates thread local stores and applies a binary functor for each store.
 * @headerfile <seqan/parallel.h>
 *
 * @signature void combine(etl, f);
 * @param[in] etl An instance of @link EnumerableThreadLocal @endlink.
 * @param[in] f A binary combinator called on each thread local storage.
 *
 * @datarace not thread-safe.
 * @see EnumerableThreadLocal#combineEach
 * @section Possible implementation:
 * @code{.cpp}
 * std::accumulate(begin(etl), end(etl), TValue{}, f);
 * @endcode
 */
template <typename TValue, typename TManager, typename TSpec,
          typename TBinaryCombine>
inline auto
combine(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
        TBinaryCombine && fBinaryCombine)
{
    std::shared_lock<decltype(storageManager(me)._mutex)> read_lck(storageManager(me)._mutex);
    return std::accumulate(begin(me), end(me), TValue{}, std::forward<TBinaryCombine>(fBinaryCombine));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_H_
