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

#ifndef INCLDUE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_ITERATOR_H_
#define INCLDUE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_ITERATOR_H_

namespace std
{

template<typename TContainer>
struct iterator_traits<seqan::Iter<TContainer, seqan::EnumerableThreadLocalIterSpec> > // nolint
{
    typedef seqan::Iter<TContainer, seqan::EnumerableThreadLocalIterSpec> TIter; // nolint

    typedef forward_iterator_tag iterator_category; // nolint
    typedef typename seqan::Value<TIter>::Type value_type; // nolint
    typedef typename seqan::Difference<TIter>::Type difference_type; // nolint
    typedef typename seqan::Value<TIter>::Type * pointer; // nolint
    typedef typename seqan::Reference<TIter>::Type reference; // nolint
};

}

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TEnumerableThreadLocal>
class Iter<TEnumerableThreadLocal, EnumerableThreadLocalIterSpec>
{
public:

    //-------------------------------------------------------------------------
    // Member Types.

    using TMap = typename TEnumerableThreadLocal::TMap;
    using TMapIterator = typename IfC<std::is_const<TEnumerableThreadLocal>::value,
                                      typename TMap::const_iterator,
                                      typename TMap::iterator>::Type;

    //-------------------------------------------------------------------------
    // Member Variables.

    TMapIterator _mInnerIter{};

    //-------------------------------------------------------------------------
    // Constructors.

    Iter() = default;

    Iter(typename TMap::iterator const & it) : _mInnerIter(it)
    {}

    Iter(typename TMap::const_iterator const & cit) : _mInnerIter(cit)
    {}

    Iter(TEnumerableThreadLocal & iterable) : Iter(iterable._map.begin())
    {}

    //-------------------------------------------------------------------------
    // Member Functions.

    inline auto *
    operator->()
    {
        return &*(*this);
    }

    inline auto *
    operator->() const
    {
        return &*(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TEnumerableThreadLocal>
inline auto &
operator*(Iter<TEnumerableThreadLocal, EnumerableThreadLocalIterSpec> & me)
{
    return me._mInnerIter->second;
}

template <typename TEnumerableThreadLocal>
inline auto const &
operator*(Iter<TEnumerableThreadLocal, EnumerableThreadLocalIterSpec> const & me)
{
    return me._mInnerIter->second;
}

template <typename TEnumerableThreadLocal>
inline void
operator++(Iter<TEnumerableThreadLocal, EnumerableThreadLocalIterSpec> & me)
{
    ++me._mInnerIter;
}

template <typename TEnumerableThreadLocal>
inline void
operator++(Iter<TEnumerableThreadLocal, EnumerableThreadLocalIterSpec> & me, int /*postfix*/)
{
    me._mInnerIter++;
}

template <typename TEnumerableThreadLocal>
inline void
operator--(Iter<TEnumerableThreadLocal, EnumerableThreadLocalIterSpec> & me)
{
    --me._mInnerIter;
}

template <typename TEnumerableThreadLocal>
inline void
operator--(Iter<TEnumerableThreadLocal, EnumerableThreadLocalIterSpec> & me, int /*postfix*/)
{
    me._mInnerIter--;
}

template <typename TEnumerableThreadLocalLeft,
          typename TEnumerableThreadLocalRight,
          typename = std::enable_if<std::is_same<typename std::decay<TEnumerableThreadLocalLeft>::type,
                                                 typename std::decay<TEnumerableThreadLocalRight>::type>::value>>
inline bool
operator==(Iter<TEnumerableThreadLocalLeft, EnumerableThreadLocalIterSpec> const & lhs,
           Iter<TEnumerableThreadLocalRight, EnumerableThreadLocalIterSpec> const & rhs)
{
    return lhs._mInnerIter == rhs._mInnerIter;
}

template <typename TEnumerableThreadLocalLeft,
          typename TEnumerableThreadLocalRight>
inline bool
operator!=(Iter<TEnumerableThreadLocalLeft, EnumerableThreadLocalIterSpec> const & lhs,
           Iter<TEnumerableThreadLocalRight, EnumerableThreadLocalIterSpec> const & rhs)
{
    return !(lhs._mInnerIter == rhs._mInnerIter);
}

}  // namespace seqan

#endif  // #ifndef INCLDUE_SEQAN_ALIGN_PARALLEL_ENUMERABLE_THREAD_LOCAL_ITERATOR_H_
