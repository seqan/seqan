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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_LIMIT_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_LIMIT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TValue>
class EnumerableThreadLocal<TValue, Limit>
{
public:

    //-------------------------------------------------------------------------
    // Member types.

    using TMap           = std::map<std::thread::id, TValue>;
    using TIterator      = Iter<EnumerableThreadLocal, EnumerableThreadLocalIterSpec>;
    using TConstIterator = Iter<EnumerableThreadLocal const, EnumerableThreadLocalIterSpec>;

    //-------------------------------------------------------------------------
    // Private Members.

    TMap                     _mMap{};
    std::shared_timed_mutex  _mMutex{};        // TODO(rrahn): Replace by shared_mutex when c++17 is available.
    TValue                   _mInitValue{};
    std::atomic<size_t>      _mNumThreads{0};

    //-------------------------------------------------------------------------
    // Constructor.

    EnumerableThreadLocal() = delete;

    EnumerableThreadLocal(size_t const numThreads) : _mNumThreads(numThreads)
    {}

    EnumerableThreadLocal(size_t const numThreads, TValue initValue) : _mInitValue(std::move(initValue)), _mNumThreads(numThreads)
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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
//  Function local()
// ----------------------------------------------------------------------------

template <typename TValue>
inline auto&
local(EnumerableThreadLocal<TValue, Limit> & me,
      bool & exists)
{
    if (me._mNumThreads.load(std::memory_order_relaxed) == 0)
    {
        exists = true;
        return me._mMap.find(std::this_thread::get_id())->second;
    }

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
        --me._mNumThreads;
        SEQAN_ASSERT(exists);
        exists = false;  // Notify that element was added for the first time.
    }
    return elemIt->second;
}


}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_LIMIT_H_
