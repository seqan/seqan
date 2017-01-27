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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_THREAD_ENUMERABLE_SPECIFIC_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_THREAD_ENUMERABLE_SPECIFIC_H_

namespace seqan
{
namespace impl
{
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// We leave it open but we do not use it yet.
template <typename TSpec>
struct EnumerableThreadSpecificTraits
{};

template <typename TValue, typename TTraits = EnumerableThreadSpecificTraits<TValue>>
class EnumerableThreadSpecific
{
public:
    
    std::map<std::thread::id, TValue*> _mMap;
    std::shared_mutex                  _mMutex;

    // Constructor with default value.
    // So whenever we create one for the first time, we choose the default value
    // We can move construct it as well.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
//  Function createThreadLocal()
// ----------------------------------------------------------------------------

template <typename TValue>
auto createThreadLocal()
{
    static thread_local alignas(128) TValue storage{};
    return &storage;
}

// ----------------------------------------------------------------------------
//  Function createNewEntry()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
auto createNewEntry(EnumerableThreadSpecific<TValue, TTraits> & me)
{
    std::uniqe_lock<decltype(me._mMutex)> write_lck(me._mMutex);
    return me._mMap.insert({std::thread::id, createThreadLocal()});
}

// ----------------------------------------------------------------------------
//  Function local()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits>
auto& local(EnumerableThreadSpecific<TValue, TTraits> & me, bool & exists)
{
    std::shared_lock<decltype(me._mMutex)> read_lck(me._mMutex);
    exists = true;
    auto elemIt = me._mMap.find(std::thread::id);
    if (elemIt == me._mMap.end())
    {
        std::tie<elemIt, exists> = createNewEntry(me);
        if (!exists) {
            // We need to throw an exception here, becuase something went terribly wrong.
        }
        exists = false;
    }
    return *(*elemIt); // Double indirection to to iterator and pointer to therad_local storage.
}

template <typename TValue, typename TTraits>
auto& local(EnumerableThreadSpecific<TValue, TTraits> & me)
{
    bool exists{true};
    return local(me, exists); // Double indirection to to iterator and pointer to therad_local storage.
}

// ----------------------------------------------------------------------------
// Function combine()
// ----------------------------------------------------------------------------

template <typename TValue, typename TTraits, typename TUnaryCombine>
void combine(EnumerableThreadSpecific<TValue, TTraits> const & me,
             TUnaryCombine && fUnaryCombine)
{
    std::shared_lock<decltype(me._mMutex)> read_lck(me._mMutex);
    std::for_each(me._mMap.begin(), me._mMap.end(), fUnaryCombine);
}

template <typename TValue, typename TTraits>
auto combine(EnumerableThreadSpecific<TValue, TTraits> const & me,
             std::function<TValue(TValue const &, TValue const &)> && fBinaryCombine)
{
    std::shared_lock<decltype(me._mMutex)> read_lck(me._mMutex);
    TValue init{};
    std::accumulate(me._mMap.begin(), me._mMap.end(), init, fBinaryCombine);
    return init;
}

}  // namespace impl
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_THREAD_SPECIFIC_DP_CONTEXT_H_
