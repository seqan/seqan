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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Thread-safe / lock-free sequence operations.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_SEQUENCE_H_
#define SEQAN_PARALLEL_PARALLEL_SEQUENCE_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

// Concurrent lock-free append for alloc stirngs.
template <typename TSpec = void>
struct ConcurrentAppend{};

// Subclassed alloc string.
template <typename TValue, typename TSpec>
class String<TValue, Alloc<ConcurrentAppend<TSpec> > >
{
public:
    using TPtr = typename Value<String>::Type *;

    std::atomic<TPtr>           data_begin;
    std::atomic<TPtr>           data_end;
    typename Size<String>::Type data_capacity;
    ReadWriteLock               lock;

    String() : data_capacity(0)
    {
        data_begin.store(nullptr);
        data_end.store(nullptr);
    }

    // Custom constructor to create string from source.
    template <typename TSource>
    String(TSource && source) : String()
    {
        if (length(source) > 0u)
            assign(*this, std::forward<TSource>(source));
        // NOTE(rrahn): The memory order is relaxed, since concurrent construction is not allowed.
        SEQAN_ASSERT_LEQ_MSG(data_begin, data_end.load(std::memory_order_relaxed), "String end is before begin!");
    }

    template <typename TSource, typename TSize>
    String(TSource && source, TSize limit) : String()
    {
        if (length(source) > 0u)
            assign(*this, std::forward<TSource>(source), limit);
        SEQAN_ASSERT_LEQ_MSG(data_begin, data_end, "String end is before begin!");
    }

    // Copy ctor.
    String(String const & other) : String()
    {
        reserve(*this, std::min(capacity(other), computeGenerousCapacity(other, length(other))), Exact());
        if (length(other) > 0u)
            assign(*this, other);
        SEQAN_ASSERT_LEQ_MSG(data_begin, data_end, "String end is before begin!");
    }

    // Move ctor.
    String(String && other) : String()
    {
        swap(*this, other);
    }

    // Assignment
    String & operator=(String other)
    {
        swap(*this, other);
        return *this;
    }

    // Custom Assignment.
    template <typename TSource>
    inline
    String & operator=(TSource && source)
    {
        assign(*this, std::forward<TSource>(source));
        SEQAN_ASSERT_LEQ_MSG(data_begin, data_end, "String end is before begin!");
        return *this;
    }

    ~String()
    {
        arrayDestruct(this->data_begin.load(), this->data_end.load());
        _deallocateStorage(*this, this->data_begin.load(), data_capacity);
    }


    // Subscript operator.
    template <typename TPos>
    inline typename Reference<String>::Type
    operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<String const>::Type
    operator[] (TPos pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function swap()
// ----------------------------------------------------------------------------

// NOTE(rrahn): No conccurrent operation.
template <typename TTargetValue, typename TSpec>
inline void
swap(String<TTargetValue, Alloc<ConcurrentAppend<TSpec> > > & lhs,
     String<TTargetValue, Alloc<ConcurrentAppend<TSpec> > > & rhs)
{
    auto tmp = lhs.data_begin.load();
    lhs.data_begin.store(rhs.data_begin.load());
    rhs.data_begin.store(tmp);
    tmp = lhs.data_end.load();
    lhs.data_end.store(rhs.data_end.load());
    rhs.data_end.store(tmp);
    std::swap(lhs.data_capacity, rhs.data_capacity);
    // NOTE(rrahn): The lock should not be swapped.
}

// ----------------------------------------------------------------------------
// Function _incLength()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<String<TValue, Alloc<ConcurrentAppend<TSpec> > > >::Type
_incLength(String<TValue, Alloc<ConcurrentAppend<TSpec> > > & me)
{
    return ++me.data_end - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function _decLength()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<String<TValue, Alloc<ConcurrentAppend<TSpec> > > >::Type
_decLength(String<TValue, Alloc<ConcurrentAppend<TSpec> > > & me)
{
    return --me.data_end - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TTargetValue, typename TSpec, typename TValue>
inline void
appendValue(String<TTargetValue, Alloc<ConcurrentAppend<TSpec> > > & me,
            TValue && _value,
            Insist const &,
            Parallel const &)
{
    valueConstruct(begin(me, Standard()) + _incLength(me) - 1, std::forward<TValue>(_value));
}

template <typename TTargetValue, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(String<TTargetValue, Alloc<ConcurrentAppend<TSpec> > > & me,
            TValue && val,
            Tag<TExpand> const & expandTag,
            Parallel const &)
{
    while (true)
    {
        // try to append the value
        {
            ScopedReadLock<> readLock(me.lock);
            decltype(capacity(me)) newLen = _incLength(me);
            if (newLen <= capacity(me))
            {
                valueConstruct(begin(me, Standard()) + newLen - 1, std::forward<TValue>(val));
                break;
            }
            _decLength(me);
        }
        
        // try to extend capacity
        {
            ScopedWriteLock<> writeLock(me.lock);
            decltype(length(me)) cap = capacity(me);
            if (cap == length(me))
                reserve(me, cap + 1, expandTag);
        }
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_SEQUENCE_H_
