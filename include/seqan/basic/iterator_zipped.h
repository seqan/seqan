// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_BASIC_ITERATOR_ZIPPED_H_
#define INCLUDE_SEQAN_BASIC_ITERATOR_ZIPPED_H_

#include <tuple>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// This function enables us to use the initializer-list trick to expand a parameter pack,
// in order to call a function on each element of the parameter pack.
template <typename T> inline void
_seqanUnpackFunc(std::initializer_list<T> const /*unusued*/)
{}

// I use this macro to hide wrapping the function call as an initializer-list.
#define SEQAN_UNPACK_FUNC(f) _seqanUnpackFunc({(f, 0)...})

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Traits to hold the indices of the corresponding tuple.
template <unsigned... >
struct TupleIndices
{};

// Tag to determine the zip iterator.
struct ZipIterator_;
typedef Tag<ZipIterator_> ZipIterator;

// The ZipIterator wraps a pack of Iterator Types.
// We use a std::tuple to store the different iterators.
// The tuple also helps to use variadic templates, while being compliant with
// the SeqAn class definitions.
template <typename... TIteratorPack>
class Iter<std::tuple<TIteratorPack...>, ZipIterator>
{
public:
    std::tuple<TIteratorPack...> dataIter;  // tuple stores the different iterators.

    // Default c'tor.
    Iter()
    {
        static_assert(sizeof...(TIteratorPack) > 0, "Requires at least one argument!");
    }

    // Custom c'tor to create from tuple.
    Iter(std::tuple<TIteratorPack...> && iterTuple) : dataIter(iterTuple)
    {
        static_assert(sizeof...(TIteratorPack) > 0, "Requires at least one argument!");
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction MakeTupleIndices
// ----------------------------------------------------------------------------

// This metafunction recursively calls itself, while adding a new index to the growing
// Indices pack.
template<unsigned ID, unsigned... Indices>
struct MakeTupleIndices :
MakeTupleIndices<ID-1, ID-1, Indices...>
{};

// Base case to stop the recursion and to define the TupleIndices with the indices pack.
template<unsigned... Indices>
struct MakeTupleIndices<0, Indices...>
{
    typedef TupleIndices<Indices...> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// We expand each iterator type in the param pack to call its Value metafunction
// and define a tuple with the value types.
template <typename... TIteratorPack>
struct Value<Iter<std::tuple<TIteratorPack...>, ZipIterator> >
{
    typedef std::tuple<typename Value<TIteratorPack>::Type...>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

// See above.
template <typename... TIteratorPack>
struct GetValue<Iter<std::tuple<TIteratorPack...>, ZipIterator> >
{
    typedef std::tuple<typename GetValue<TIteratorPack>::Type...>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

// See above.
template <typename... TIteratorPack>
struct Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> >
{
    typedef std::tuple<typename Reference<TIteratorPack>::Type...>  Type;
};

// Unclear why we need that! It should be possible to let the generic implementation
// in iterator_base.h to take care of it. Apparently, it is the same semantic.
template <typename... TIteratorPack>
struct Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> const> :
    Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> > {};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

// The size type is only defined by the first type as only the first iterator
// is used to determine the difference between two iterators.
template <typename TFirst, typename... TIteratorPack>
struct Size<Iter<std::tuple<TFirst, TIteratorPack...>, ZipIterator> >
{
    typedef typename Size<TFirst>::Type Type;  // The position and size is determined only by the first argument in the tuple.
};

// Use the default Position-MF

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

// We use only the first iterator to determine the difference.
template <typename TFirst, typename... TIteratorPack>
struct Difference<Iter<std::tuple<TFirst, TIteratorPack...>, ZipIterator> >
{
    typedef typename Difference<TFirst>::Type Type;  // The difference is determined by only the first argument in the tuple.
};

// ============================================================================
// Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::increment
// ----------------------------------------------------------------------------

// private increment interface.
// Uses the initializer-trick to apply the ++operator to all elements in the tuple
// using pack expansion.
template <typename... TIteratorPack, unsigned... TUPLE_INDICES>
inline void
increment(std::tuple<TIteratorPack...> & me,
          TupleIndices<TUPLE_INDICES...> const)
{
    SEQAN_UNPACK_FUNC(++std::get<TUPLE_INDICES>(me));
}

// ----------------------------------------------------------------------------
// Function impl::decrement
// ----------------------------------------------------------------------------

// Same here for decrement.
template <typename... TIteratorPack, unsigned... TUPLE_INDICES>
inline void
decrement(std::tuple<TIteratorPack...> & me,
          TupleIndices<TUPLE_INDICES...> const)
{
    SEQAN_UNPACK_FUNC(--std::get<TUPLE_INDICES>(me));
}

// ----------------------------------------------------------------------------
// Function impl::advance
// ----------------------------------------------------------------------------

// ... and for advancing the iterator.
template <typename... TIteratorPack, unsigned... INDICES, typename TIntegral>
inline void
advance(std::tuple<TIteratorPack...> & me,
        TupleIndices<INDICES...> const,
        TIntegral const steps)
{
    SEQAN_UNPACK_FUNC(std::advance(std::get<INDICES>(me), steps));
}

// ----------------------------------------------------------------------------
// Function impl::dereference
// ----------------------------------------------------------------------------

// Dereference each iterator in the tuple and forward return values as tuple.
template <typename... TIteratorPack, unsigned... TUPLE_INDICES>
inline auto
dereference(std::tuple<TIteratorPack...> & me,
            TupleIndices<TUPLE_INDICES...> const) -> decltype(std::forward_as_tuple(*std::get<TUPLE_INDICES>(me)...))
{
    return std::forward_as_tuple(*std::get<TUPLE_INDICES>(me)...);
}

template <typename... TIteratorPack, unsigned... TUPLE_INDICES>
inline auto
dereference(std::tuple<TIteratorPack...> const & me,
            TupleIndices<TUPLE_INDICES...> const) -> decltype(std::forward_as_tuple(*std::get<TUPLE_INDICES>(me)...))
{
    return std::forward_as_tuple(*std::get<TUPLE_INDICES>(me)...);
}

}

// ----------------------------------------------------------------------------
// Function operator*
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline typename Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> >::Type
operator*(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me)
{
    return impl::dereference(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type());
}

template <typename... TIteratorPack>
inline typename Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> const>::Type
operator*(Iter<std::tuple<TIteratorPack...>, ZipIterator> const & me)
{
    return impl::dereference(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type());
}

// ----------------------------------------------------------------------------
// Function operator++; pre-increment
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator> &
operator++(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me /*pre-increment*/)
{
    impl::increment(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type());
    return me;
}

// ----------------------------------------------------------------------------
// Function operator++; post-increment
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator>
operator++(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me, int const /*post-increment*/)
{
    typename std::remove_reference<decltype(me)>::type tmp(me);
    impl::increment(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type());
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator+
// ----------------------------------------------------------------------------

template <typename... TIteratorPack, typename TIntegral>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator>
operator+(Iter<std::tuple<TIteratorPack...>, ZipIterator> me,
          TIntegral const steps)
{
    impl::advance(me.dataIter,
                  typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type(),
                  steps);
    return me;
}

template <typename... TIteratorPack, typename TIntegral>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator>
operator+(TIntegral const steps,
          Iter<std::tuple<TIteratorPack...>, ZipIterator> me)
{
    impl::advance(me.dataIter,
                  typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type(),
                  steps);
    return me;
}

// ----------------------------------------------------------------------------
// Function operator+=
// ----------------------------------------------------------------------------

template <typename... TIteratorPack, typename TIntegral>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator> &
operator+=(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me,
           TIntegral const steps)
{
    impl::advance(me.dataIter,
                  typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type(),
                  steps);
    return me;
}

// ----------------------------------------------------------------------------
// Function operator--; pre-increment
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator>&
operator--(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me /*pre-increment*/)
{
    impl::decrement(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type());
    return me;
}

// ----------------------------------------------------------------------------
// Function operator--; post-increment
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator>
operator--(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me, int const /*post-increment*/)
{
    typename std::remove_reference<decltype(me)>::type tmp(me);
    impl::decrement(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type());
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator-
// ----------------------------------------------------------------------------

template <typename... TIteratorPack, typename TIntegral>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator>
operator-(Iter<std::tuple<TIteratorPack...>, ZipIterator> me,
          TIntegral const steps)
{
    typedef typename MakeSigned<TIntegral>::Type TSigned;
    impl::advance(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type(),
                  -static_cast<TSigned>(steps));
    return me;
}

template <typename... TIteratorPack>
inline typename Difference<Iter<std::tuple<TIteratorPack...>, ZipIterator> >::Type
operator-(Iter<std::tuple<TIteratorPack...>, ZipIterator> const & lhs,
          Iter<std::tuple<TIteratorPack...>, ZipIterator> const & rhs)
{
    return std::get<0>(lhs.dataIter) - std::get<0>(rhs.dataIter);
}

// ----------------------------------------------------------------------------
// Function operator-=
// ----------------------------------------------------------------------------

template <typename... TIteratorPack, typename TIntegral>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator> &
operator-=(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me,
           TIntegral const steps)
{
    typedef typename MakeSigned<TIntegral>::Type TSigned;
    impl::advance(me.dataIter, typename MakeTupleIndices<sizeof...(TIteratorPack)>::Type(),
                  -static_cast<TSigned>(steps));
    return me;
}

// ----------------------------------------------------------------------------
// Function operator==
// ----------------------------------------------------------------------------

template <typename... TIteratorPackL, typename... TIteratorPackR>
inline bool
operator==(Iter<std::tuple<TIteratorPackL...>, ZipIterator> const & lhs,
           Iter<std::tuple<TIteratorPackR...>, ZipIterator> const & rhs)
{
    return lhs.dataIter == rhs.dataIter;
}

// ----------------------------------------------------------------------------
// Function operator!=
// ----------------------------------------------------------------------------

template <typename... TIteratorPackL, typename... TIteratorPackR>
inline bool
operator!=(Iter<std::tuple<TIteratorPackL...>, ZipIterator> const & lhs,
           Iter<std::tuple<TIteratorPackR...>, ZipIterator> const & rhs)
{
    return lhs.dataIter != rhs.dataIter;
}

// ----------------------------------------------------------------------------
// Function operator>=
// ----------------------------------------------------------------------------

template <typename... TIteratorPackL, typename... TIteratorPackR>
inline bool
operator>=(Iter<std::tuple<TIteratorPackL...>, ZipIterator> const & lhs,
           Iter<std::tuple<TIteratorPackR...>, ZipIterator> const & rhs)
{
    return lhs.dataIter >= rhs.dataIter;
}

// ----------------------------------------------------------------------------
// Function operator>
// ----------------------------------------------------------------------------

template <typename... TIteratorPackL, typename... TIteratorPackR>
inline bool
operator>(Iter<std::tuple<TIteratorPackL...>, ZipIterator> const & lhs,
          Iter<std::tuple<TIteratorPackR...>, ZipIterator> const & rhs)
{
    return lhs.dataIter > rhs.dataIter;
}

// ----------------------------------------------------------------------------
// Function operator<=
// ----------------------------------------------------------------------------

template <typename... TIteratorPackL, typename... TIteratorPackR>
inline bool
operator<=(Iter<std::tuple<TIteratorPackL...>, ZipIterator> const & lhs,
           Iter<std::tuple<TIteratorPackR...>, ZipIterator> const & rhs)
{
    return lhs.dataIter <= rhs.dataIter;
}

// ----------------------------------------------------------------------------
// Function operator<
// ----------------------------------------------------------------------------

template <typename... TIteratorPackL, typename... TIteratorPackR>
inline bool
operator<(Iter<std::tuple<TIteratorPackL...>, ZipIterator> const & lhs,
          Iter<std::tuple<TIteratorPackR...>, ZipIterator> const & rhs)
{
    return lhs.dataIter < rhs.dataIter;
}

// ----------------------------------------------------------------------------
// Function swap()
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline void
swap(Iter<std::tuple<TIteratorPack...>, ZipIterator> & lhs,
     Iter<std::tuple<TIteratorPack...>, ZipIterator> & rhs)
{
    std::swap(lhs.dataIter, rhs.dataIter);
}

// ----------------------------------------------------------------------------
// Function makeZippedIterator()
// ----------------------------------------------------------------------------

//  helper function
template <typename... TIteratorPack>
inline Iter<std::tuple<TIteratorPack... >, ZipIterator>
makeZippedIterator(TIteratorPack... iterPack)
{
    return Iter<std::tuple<TIteratorPack... >, ZipIterator>(std::make_tuple(iterPack...));
}

}

#endif  // #ifndef INCLUDE_SEQAN_BASIC_ITERATOR_ZIPPED_H_
