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

#ifndef INCLUDE_SEQAN_BASIC_ITERATOR_ZIP_H_
#define INCLUDE_SEQAN_BASIC_ITERATOR_ZIP_H_

#include <tuple>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// This function enables us to use the initializer-list trick to expand a parameter pack,
// in order to call a function on each element of the parameter pack.
template <typename T>
inline void
_seqanUnpackFunc(std::initializer_list<T> const /*unused*/)
{}

// This macro is used to hide wrapping the function call as an initializer-list.
#define SEQAN_UNPACK_FUNC(f) _seqanUnpackFunc({(f, 0)...})

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Traits to hold the indices of the corresponding tuple.
template <unsigned... >
struct IndexSequence
{};

// Tag to determine the zip iterator.
struct ZipIterator_;
typedef Tag<ZipIterator_> ZipIterator;

/*!
 * @class ZipIterator
 *
 * @extends Iter
 * @implements IteratorAssociatedTypesConcept
 *
 * @headerfile <seqan/basic.h>
 *
 * @brief Zips multiple iterators over different containers into a single iterator.
 *
 * @signature class Iter<std::tuple<TIteratorTypes...>, ZipIterator>;
 *
 * @tparam TIteratorTypes A template parameter pack with one or more @link ContainerConcept#Iterator @endlink types.
 *
 * This iterator ties together different iterator types for different containers of the same size.
 * It allows one to operate on a single iterator, if multiple containers need to be traversed simultaneously.
 * Note, that all operations are still executed in a serial fashion.
 * If the zip iterator is dereferenced it returns a <a href="http://en.cppreference.com/w/cpp/utility/tuple">std::tuple</a>
 * containing the dereferenced values of all embedded iterators. 
 * The metafunctions @link Value @endlink, @link GetValue @endlink and @link Reference @endlink
 * are overloaded accordingly.
 *
 * To easily create a zip iterator one can use the helper function @link makeZipIterator @endlink.
 *
 * @section Example
 *
 * The following example demonstrates the function of the zip iterator:
 *
 * @include demos/dox/basic/zip_iterator.cpp
 *
 * This outputs the following to the console:
 * @include demos/dox/basic/zip_iterator.cpp.stdout
 *
 * @note Throws an assertion if <tt>sizeof...(TIteratorTypes) == 0</tt> is true.
 */

template <typename... TIteratorPack>
class Iter<std::tuple<TIteratorPack...>, ZipIterator>
{
public:
    std::tuple<TIteratorPack...> dataIter;  // tuple stores the different iterators.

    /*!
     * @fn ZipIterator#ZipIterator
     * @brief Constructor.
     * @signature Iter()
     *            Iter(TIteratorTypes... args)
     *
     * @param args The iterator instances to create the @link ZipIterator @endlink from. Default constructors are not listed.
     */

    // Default c'tor.
    Iter()
    {
        static_assert(sizeof...(TIteratorPack) > 0, "Requires at least one argument!");
    }

    // Custom c'tor to create from iterator pack
    Iter(TIteratorPack... iterPack) : dataIter(std::make_tuple(iterPack...))
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction MakeIndexSequence
// ----------------------------------------------------------------------------

// This metafunction recursively calls itself, while adding a new index to the growing
// Indices pack.
template<unsigned ID, unsigned... Indices>
struct MakeIndexSequence :
MakeIndexSequence<ID-1, ID-1, Indices...>
{};

// Base case to stop the recursion and to define the IndexSequence with the indices pack.
template<unsigned... Indices>
struct MakeIndexSequence<0, Indices...>
{
    typedef IndexSequence<Indices...> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

/*!
 * @mfn ZipIterator#Container
 * @headerfile <seqan/basic.h>
 * @brief Retruns the Container type of the @link ZipIterator @endlink.
 *
 * @signature typename Container<TZipIterator>::Type;
 *
 * @tparam TZipIterator The type of the @link ZipIterator @endlink.
 *
 * @return TRes A <a href="http://en.cppreference.com/w/cpp/utility/tuple">std::tuple</a> containing the container
 *               types of all iterator types embedded in the <tt>TZipIterator</tt> type.
 */
template <typename... TIteratorPack>
struct Container<Iter<std::tuple<TIteratorPack...>, ZipIterator> >
{
    typedef std::tuple<typename Container<TIteratorPack>::Type...> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// We expand each iterator type in the param pack to call its Value metafunction
// and define a tuple with the value types.

/*!
 * @mfn ZipIterator#Value
 * @headerfile <seqan/basic.h>
 * @brief Retruns the Value type of the @link ZipIterator @endlink.
 *
 * @signature typename Value<TZipIterator>::Type;
 *
 * @tparam TZipIterator The type of the @link ZipIterator @endlink.
 *
 * @return TRes A <a href="http://en.cppreference.com/w/cpp/utility/tuple">std::tuple</a> containing the value
 *               types of all iterator types embedded in the <tt>TZipIterator</tt> type.
 */
template <typename... TIteratorPack>
struct Value<Iter<std::tuple<TIteratorPack...>, ZipIterator> >
{
    typedef std::tuple<typename Value<TIteratorPack>::Type...> Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

/*!
 * @mfn ZipIterator#GetValue
 * @headerfile <seqan/basic.h>
 * @brief Retruns the GetValue type of the @link ZipIterator @endlink.
 *
 * @signature typename GetValue<TZipIterator>::Type;
 *
 * @tparam TZipIterator The type of the @link ZipIterator @endlink.
 *
 * @return TRes A <a href="http://en.cppreference.com/w/cpp/utility/tuple">std::tuple</a> containing the value
 *               types of all iterator types embedded in the <tt>TZipIterator</tt> type.
 */
template <typename... TIteratorPack>
struct GetValue<Iter<std::tuple<TIteratorPack...>, ZipIterator> >
{
    typedef std::tuple<typename GetValue<TIteratorPack>::Type...> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

/*!
 * @mfn ZipIterator#Reference
 * @headerfile <seqan/basic.h>
 * @brief Retruns the Reference type of the @link ZipIterator @endlink.
 *
 * @signature typename GetValue<TZipIterator>::Type;
 *
 * @tparam TZipIterator The type of the @link ZipIterator @endlink.
 *
 * @return TRes A <a href="http://en.cppreference.com/w/cpp/utility/tuple">std::tuple</a> containing the value
 *               types of all iterator types embedded in the <tt>TZipIterator</tt> type.
 */
template <typename... TIteratorPack>
struct Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> >
{
    typedef std::tuple<typename Reference<TIteratorPack>::Type...> Type;
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
template <typename... TIteratorPack, unsigned... INDEX_SEQUENCE>
inline void
increment(std::tuple<TIteratorPack...> & me,
          IndexSequence<INDEX_SEQUENCE...> const)
{
    SEQAN_UNPACK_FUNC(++std::get<INDEX_SEQUENCE>(me));
}

// ----------------------------------------------------------------------------
// Function impl::decrement
// ----------------------------------------------------------------------------

// Same here for decrement.
template <typename... TIteratorPack, unsigned... INDEX_SEQUENCE>
inline void
decrement(std::tuple<TIteratorPack...> & me,
          IndexSequence<INDEX_SEQUENCE...> const)
{
    SEQAN_UNPACK_FUNC(--std::get<INDEX_SEQUENCE>(me));
}

// ----------------------------------------------------------------------------
// Function impl::advance
// ----------------------------------------------------------------------------

// ... and for advancing the iterator.
template <typename... TIteratorPack, unsigned... INDICES, typename TIntegral>
inline void
advance(std::tuple<TIteratorPack...> & me,
        IndexSequence<INDICES...> const,
        TIntegral const steps)
{
    SEQAN_UNPACK_FUNC(std::advance(std::get<INDICES>(me), steps));
}

// ----------------------------------------------------------------------------
// Function impl::dereference
// ----------------------------------------------------------------------------

// Dereference each iterator in the tuple and forward return values as tuple.
template <typename... TIteratorPack, unsigned... INDEX_SEQUENCE>
inline auto
dereference(std::tuple<TIteratorPack...> & me,
            IndexSequence<INDEX_SEQUENCE...> const) -> decltype(std::forward_as_tuple(*std::get<INDEX_SEQUENCE>(me)...))
{
    return std::forward_as_tuple(*std::get<INDEX_SEQUENCE>(me)...);
}

template <typename... TIteratorPack, unsigned... INDEX_SEQUENCE>
inline auto
dereference(std::tuple<TIteratorPack...> const & me,
            IndexSequence<INDEX_SEQUENCE...> const) -> decltype(std::forward_as_tuple(*std::get<INDEX_SEQUENCE>(me)...))
{
    return std::forward_as_tuple(*std::get<INDEX_SEQUENCE>(me)...);
}

// ----------------------------------------------------------------------------
// Function impl::assignValue
// ----------------------------------------------------------------------------

// Dereference each iterator in the tuple and forward return values as tuple.
template <typename... TRefPack, typename... TValuePack, unsigned... INDEX_SEQUENCE>
inline void
assignValue(std::tuple<TRefPack...> targetPack,
            std::tuple<TValuePack...> const & sourcePack,
            IndexSequence<INDEX_SEQUENCE...> const)
{
    SEQAN_UNPACK_FUNC(assign(std::get<INDEX_SEQUENCE>(targetPack), std::get<INDEX_SEQUENCE>(sourcePack)));
}

// ----------------------------------------------------------------------------
// Function impl::moveValue
// ----------------------------------------------------------------------------

// Dereference each iterator in the tuple and forward return values as tuple.
template <typename... TRefPack, typename... TValuePack, unsigned... INDEX_SEQUENCE>
inline void
moveValue(std::tuple<TRefPack...> targetPack,
          std::tuple<TValuePack...> const & sourcePack,
          IndexSequence<INDEX_SEQUENCE...> const)
{
    SEQAN_UNPACK_FUNC(move(std::get<INDEX_SEQUENCE>(targetPack), std::get<INDEX_SEQUENCE>(sourcePack)));
}

}  // namespace impl

// ----------------------------------------------------------------------------
// Function operator*
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline typename Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> >::Type
operator*(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me)
{
    return impl::dereference(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
}

template <typename... TIteratorPack>
inline typename Reference<Iter<std::tuple<TIteratorPack...>, ZipIterator> const>::Type
operator*(Iter<std::tuple<TIteratorPack...>, ZipIterator> const & me)
{
    return impl::dereference(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
}

// ----------------------------------------------------------------------------
// Function operator++; pre-increment
// ----------------------------------------------------------------------------

template <typename... TIteratorPack>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator> &
operator++(Iter<std::tuple<TIteratorPack...>, ZipIterator> & me /*pre-increment*/)
{
    impl::increment(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
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
    impl::increment(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
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
                  typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type(),
                  steps);
    return me;
}

template <typename... TIteratorPack, typename TIntegral>
inline Iter<std::tuple<TIteratorPack...>, ZipIterator>
operator+(TIntegral const steps,
          Iter<std::tuple<TIteratorPack...>, ZipIterator> me)
{
    impl::advance(me.dataIter,
                  typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type(),
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
                  typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type(),
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
    impl::decrement(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
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
    impl::decrement(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
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
    impl::advance(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type(),
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
    impl::advance(me.dataIter, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type(),
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
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename... TIteratorPack, typename... TValuePack>
inline void
assignValue(Iter<std::tuple<TIteratorPack...>, ZipIterator> & target,
            std::tuple<TValuePack...> const & value)
{
    static_assert(sizeof...(TIteratorPack) == sizeof...(TValuePack), "Wrong number of arguments.");
    impl::assignValue(*target, value, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

template <typename... TIteratorPack, typename... TValuePack>
inline void
moveValue(Iter<std::tuple<TIteratorPack...>, ZipIterator> & target,
          std::tuple<TValuePack...> const & value)
{
    static_assert(sizeof...(TIteratorPack) == sizeof...(TValuePack), "Wrong number of arguments.");
    impl::moveValue(*target, value, typename MakeIndexSequence<sizeof...(TIteratorPack)>::Type());
}

// ----------------------------------------------------------------------------
// Function makeZipIterator()
// ----------------------------------------------------------------------------

/*!
 * @fn makeZipIterator
 * @headerfile <seqan/basic_iterator.h>
 * @brief Creates a @link ZipIterator @endlink, deducing the iterator types from the arguments.
 *
 * @signature iter makeZipIterator(TIteratorTypes... args)
 *
 * @param [in] args One or more iterator instances to construct the @link ZipIterator @endlink from.
 *
 * @return iter A @link ZipIterator @endlink containing the given iterator instances.
 */
template <typename... TIteratorPack>
inline Iter<std::tuple<TIteratorPack... >, ZipIterator>
makeZipIterator(TIteratorPack... iterPack)
{
    return Iter<std::tuple<TIteratorPack... >, ZipIterator>(iterPack...);
}

}

#endif  // #ifndef INCLUDE_SEQAN_BASIC_ITERATOR_ZIP_H_
