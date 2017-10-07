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

struct SimpleThreadLocalManager
{
    std::shared_timed_mutex  _mutex{};        // TODO(rrahn): Replace by shared_mutex when c++17 is available.

    template <typename TResourceMap, typename TValue>
    inline auto local(TResourceMap & map, TValue const & initValue)
    {
        decltype(map.find(std::this_thread::get_id())) elemIt;
        bool exists;
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
        return std::forward_as_tuple(elemIt->second, exists);
    }
};

struct CountingThreadLocalManager
{
    std::atomic<size_t>      _count{0};
    std::shared_timed_mutex  _mutex{};        // TODO(rrahn): Replace by shared_mutex when c++17 is available.

    template <typename TResourceMap, typename TValue>
    inline auto
    local(TResourceMap & map, TValue const & initValue)
    {
        bool exists{true};
        if (_count.load(std::memory_order_relaxed) == 0)
            return std::forward_as_tuple(map.find(std::this_thread::get_id())->second, exists);


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
        return std::forward_as_tuple(elemIt->second, exists);
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

// TODO(rrahn): Legacy metafunction to work with SeqAn 2.x
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

template <typename TValue, typename TManager, typename TSpec>
inline auto&
local(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
      bool & exists)
{
    auto res = me._manager.local(me._map, me._initValue);
    exists = std::get<1>(res);
    return std::get<0>(res);
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

template <typename TValue, typename TManager, typename TSpec,
          typename TUnaryCombine>
inline void
combineEach(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
             TUnaryCombine && fUnaryCombine)
{
    std::for_each(begin(me), end(me), std::forward<TUnaryCombine>(fUnaryCombine));
}

template <typename TValue, typename TManager, typename TSpec,
          typename TBinaryCombine>
inline auto
combine(EnumerableThreadLocal<TValue, TManager, TSpec> & me,
        TBinaryCombine && fBinaryCombine)
{
    TValue init{};
    return std::accumulate(begin(me), end(me), TValue{}, std::forward<TBinaryCombine>(fBinaryCombine));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_H_
