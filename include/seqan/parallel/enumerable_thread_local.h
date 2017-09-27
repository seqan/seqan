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

// We leave it open but we do not use it yet.
template <typename TValue>
struct EnumerableThreadLocalTraits
{};

template <typename TValue, typename TTraits = EnumerableThreadLocalTraits<TValue>>
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

    TMap                     _mMap{};
    std::shared_timed_mutex  _mMutex{};        // TODO(rrahn): Replace by shared_mutex when c++17 is available.
    TValue                   _mInitValue{};

    //-------------------------------------------------------------------------
    // Constructor.

    EnumerableThreadLocal() = default;

    EnumerableThreadLocal(TValue initValue) : _mInitValue(std::move(initValue))
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
template <typename TValue, typename TTraits,
          typename TIterTag>
struct Iterator<EnumerableThreadLocal<TValue, TTraits>, TIterTag>
{
    using Type = typename EnumerableThreadLocal<TValue, TTraits>::TIterator;
};

template <typename TValue, typename TTraits,
          typename TIterTag>
struct Iterator<EnumerableThreadLocal<TValue, TTraits> const, TIterTag>
{
    using Type = typename EnumerableThreadLocal<TValue, TTraits>::TConstIterator;
};

// ============================================================================
// Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
//  Function createNewEntry()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline auto
createNewEntry(EnumerableThreadLocal<TValue, TTraits> & me)
{
    // TODO(rrahn): Make cache aligned allocation of value and only store pointer in map.
    return me._mMap.emplace(std::this_thread::get_id(), me._mInitValue);
}
}  // namespace impl

// ----------------------------------------------------------------------------
//  Function local()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
inline auto&
local(EnumerableThreadLocal<TValue, TTraits> & me,
      bool & exists)
{
    decltype(me._mMap.find(std::this_thread::get_id())) elemIt;
    { // try to read
        std::shared_lock<decltype(me._mMutex)> read_lck(me._mMutex);
        elemIt = me._mMap.find(std::this_thread::get_id());
        exists = elemIt != me._mMap.end();
    }
    if (!exists)
    {
        {  // Create new entry.
            std::unique_lock<decltype(me._mMutex)> write_lck(me._mMutex);
            std::tie(elemIt, exists) = impl::createNewEntry(me);
        }
        SEQAN_ASSERT(exists);
        exists = false;  // Notify that element was added for the first time.
    }
    return elemIt->second;
}

template <typename TValue, typename TTraits>
inline auto&
local(EnumerableThreadLocal<TValue, TTraits> & me)
{
    bool exists{true};
    return local(me, exists); // Double indirection to to iterator and pointer to therad_local storage.
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits,
          typename TIterSpec>
inline auto
begin(EnumerableThreadLocal<TValue, TTraits> & me,
      Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TTraits>::TIterator{me};
}

template <typename TValue, typename TTraits,
          typename TIterSpec>
inline auto
begin(EnumerableThreadLocal<TValue, TTraits> const & me,
      Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TTraits>::TConstIterator{me};
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits,
          typename TIterSpec>
inline auto
end(EnumerableThreadLocal<TValue, TTraits> & me,
    Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TTraits>::TIterator{me._mMap.end()};
}

template <typename TValue, typename TTraits,
          typename TIterSpec>
inline auto
end(EnumerableThreadLocal<TValue, TTraits> const & me,
    Tag<TIterSpec> const & /*tag*/)
{
    return typename EnumerableThreadLocal<TValue, TTraits>::TConstIterator{me._mMap.cend()};
}

// ----------------------------------------------------------------------------
// Function combine()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits,
          typename TUnaryCombine>
inline void
combineEach(EnumerableThreadLocal<TValue, TTraits> & me,
             TUnaryCombine && fUnaryCombine)
{
    std::for_each(begin(me), end(me), std::forward<TUnaryCombine>(fUnaryCombine));
}

template <typename TValue, typename TTraits,
          typename TBinaryCombine>
inline auto
combine(EnumerableThreadLocal<TValue, TTraits> & me,
        TBinaryCombine && fBinaryCombine)
{
    TValue init{};
    return std::accumulate(begin(me), end(me), TValue{}, std::forward<TBinaryCombine>(fBinaryCombine));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_H_
