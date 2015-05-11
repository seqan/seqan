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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Adaptions for STL containers to SeqAn sequences.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_ADAPT_STL_CONTAINER_H_
#define SEQAN_SEQUENCE_ADAPT_STL_CONTAINER_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Concepts
// ===========================================================================

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::vector<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::vector<TChar, TAlloc> const), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::deque<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::deque<TChar, TAlloc> const), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::list<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::list<TChar, TAlloc> const), (StlContainerConcept));

template<typename TChar, typename TTraits, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::basic_string<TChar, TTraits, TAlloc>), (StlContainerConcept));

template<typename TChar, typename TTraits, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::basic_string<TChar, TTraits, TAlloc> const), (StlContainerConcept));

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::forward_list<TChar, TAlloc>), (StlContainerConcept));

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::forward_list<TChar, TAlloc> const), (StlContainerConcept));

template <typename TChar, std::size_t N>
SEQAN_CONCEPT_IMPL((std::array<TChar, N>), (StlContainerConcept));

template <typename TChar, std::size_t N>
SEQAN_CONCEPT_IMPL((std::array<TChar, N> const), (StlContainerConcept));
#endif

// ===========================================================================
// Metafunctions
// ===========================================================================

// ----------------------------------------------------------------------------
// Mfn IsContiguous (default is False)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct IsContiguous<std::vector<TChar, TAlloc> > :
    public True
{};

template <typename TChar, typename TAlloc>
struct IsContiguous<std::vector<TChar, TAlloc> const> :
    public True
{};

template <typename TChar, typename TTraits, typename TAlloc>
struct IsContiguous<std::basic_string<TChar, TTraits, TAlloc> > :
    public True
{};

template <typename TChar, typename TTraits, typename TAlloc>
struct IsContiguous<std::basic_string<TChar, TTraits, TAlloc> const> :
    public True
{};

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, size_t N>
struct IsContiguous<std::array<TChar, N> > :
    public True
{};

template <typename TChar, size_t N>
struct IsContiguous<std::array<TChar, N> const> :
    public True
{};
#endif

// ----------------------------------------------------------------------------
// Mfn Value
// ----------------------------------------------------------------------------

#define COMMA ,

#define SUPERMACRO__(MTFN, CONT, CONST, TMPL, MEMB) \
template <TMPL> \
struct MTFN<CONT CONST> \
{ \
    typedef typename CONT::MEMB Type; \
};

SUPERMACRO__(Value, std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc, value_type)
SUPERMACRO__(Value, std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc, value_type)
SUPERMACRO__(Value, std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc, value_type)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Value, std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc, value_type)
SUPERMACRO__(Value, std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N,   value_type)
#endif
SUPERMACRO__(Value, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc, value_type)

SUPERMACRO__(Value, std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc, value_type)
SUPERMACRO__(Value, std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc, value_type)
SUPERMACRO__(Value, std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc, value_type)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Value, std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc, value_type)
SUPERMACRO__(Value, std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N,   value_type)
#endif
SUPERMACRO__(Value, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc, value_type)

// ----------------------------------------------------------------------------
// Mfn Reference
// ----------------------------------------------------------------------------

SUPERMACRO__(Reference, std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc, reference)
SUPERMACRO__(Reference, std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc, reference)
SUPERMACRO__(Reference, std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc, reference)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Reference, std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc, reference)
SUPERMACRO__(Reference, std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N,   reference)
#endif
SUPERMACRO__(Reference, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc, reference)

SUPERMACRO__(Reference, std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(Reference, std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(Reference, std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc, const_reference)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Reference, std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(Reference, std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N,   const_reference)
#endif
SUPERMACRO__(Reference, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc, const_reference)

// ----------------------------------------------------------------------------
// Mfn GetValue
// ----------------------------------------------------------------------------

SUPERMACRO__(GetValue, std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(GetValue, std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(GetValue, std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc, const_reference)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(GetValue, std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(GetValue, std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N,   const_reference)
#endif
SUPERMACRO__(GetValue, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc, const_reference)

SUPERMACRO__(GetValue, std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(GetValue, std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(GetValue, std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc, const_reference)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(GetValue, std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc, const_reference)
SUPERMACRO__(GetValue, std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N,   const_reference)
#endif
SUPERMACRO__(GetValue, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc, const_reference)

// ----------------------------------------------------------------------------
// Mfn Size
// ----------------------------------------------------------------------------

SUPERMACRO__(Size, std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Size, std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Size, std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc, size_type)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Size, std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Size, std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N,   size_type)
#endif
SUPERMACRO__(Size, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc, size_type)

SUPERMACRO__(Size, std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Size, std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Size, std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc, size_type)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Size, std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Size, std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N,   size_type)
#endif
SUPERMACRO__(Size, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc, size_type)

// ----------------------------------------------------------------------------
// Mfn Position
// ----------------------------------------------------------------------------

SUPERMACRO__(Position, std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Position, std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Position, std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc, size_type)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Position, std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Position, std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N,   size_type)
#endif
SUPERMACRO__(Position, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc, size_type)

SUPERMACRO__(Position, std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Position, std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Position, std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc, size_type)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(Position, std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc, size_type)
SUPERMACRO__(Position, std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N,   size_type)
#endif
SUPERMACRO__(Position, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc, size_type)

// ----------------------------------------------------------------------------
// Mfn StdContainerIterator
// ----------------------------------------------------------------------------

SUPERMACRO__(StdContainerIterator, std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc, iterator)
SUPERMACRO__(StdContainerIterator, std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc, iterator)
SUPERMACRO__(StdContainerIterator, std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc, iterator)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(StdContainerIterator, std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc, iterator)
SUPERMACRO__(StdContainerIterator, std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N,   iterator)
#endif
SUPERMACRO__(StdContainerIterator, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc, iterator)

SUPERMACRO__(StdContainerIterator, std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc, const_iterator)
SUPERMACRO__(StdContainerIterator, std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc, const_iterator)
SUPERMACRO__(StdContainerIterator, std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc, const_iterator)
#ifdef SEQAN_CXX11_STANDARD
SUPERMACRO__(StdContainerIterator, std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc, const_iterator)
SUPERMACRO__(StdContainerIterator, std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N,   const_iterator)
#endif
SUPERMACRO__(StdContainerIterator, std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc, const_iterator)

// ----------------------------------------------------------------------------
// Mfn Iterator (Standard)
// ----------------------------------------------------------------------------

#define ITMACRO__(CONT, CONST, TMPL) \
template <TMPL> \
struct Iterator<CONT CONST, Standard> \
{ \
    typedef Iter<CONT CONST, StdIteratorAdaptor> Type;\
};

ITMACRO__(std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc)
ITMACRO__(std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc)
ITMACRO__(std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc)
#ifdef SEQAN_CXX11_STANDARD
ITMACRO__(std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc)
ITMACRO__(std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N)
#endif
ITMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc)

ITMACRO__(std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc)
ITMACRO__(std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc)
ITMACRO__(std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc)
#ifdef SEQAN_CXX11_STANDARD
ITMACRO__(std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc)
ITMACRO__(std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N)
#endif
ITMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc)

// TODO(h-2): remove this crap overload and resolve to stl's RandomAcccesIterator like for std::vector
// template <typename TChar, typename TCharTraits, typename TAlloc>
// struct Iterator<std::basic_string<TChar, TCharTraits, TAlloc>, Standard >
// {
//     typedef TChar * Type;
// };
//
// template <typename TChar, typename TCharTraits, typename TAlloc>
// struct Iterator<std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>
// {
//     typedef TChar const * Type;
// };

// ----------------------------------------------------------------------------
// Mfn Iterator (Rooted)
// ----------------------------------------------------------------------------

#define ITRMACRO__(CONT, CONST, TMPL) \
template <TMPL> \
struct Iterator<CONT CONST, Rooted> \
{ \
    typedef Iter<CONT CONST, AdaptorIterator<Iter<CONT CONST, StdIteratorAdaptor> > > Type;\
};

ITRMACRO__(std::vector<TChar COMMA TAlloc>,             , typename TChar COMMA typename TAlloc)
ITRMACRO__(std::deque<TChar COMMA TAlloc>,              , typename TChar COMMA typename TAlloc)
ITRMACRO__(std::list<TChar COMMA TAlloc>,               , typename TChar COMMA typename TAlloc)
#ifdef SEQAN_CXX11_STANDARD
ITRMACRO__(std::forward_list<TChar COMMA TAlloc>,       , typename TChar COMMA typename TAlloc)
ITRMACRO__(std::array<TChar COMMA N>,                   , typename TChar COMMA std::size_t N)
#endif

ITRMACRO__(std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc)
ITRMACRO__(std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc)
ITRMACRO__(std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc)
#ifdef SEQAN_CXX11_STANDARD
ITRMACRO__(std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc)
ITRMACRO__(std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N)
#endif

//NOTE(h-2): adding these breaks stuff; apperently Rooted is default for begin(), why?
// ITRMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc)
// ITRMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc)

// ----------------------------------------------------------------------------
// Mfn IsSequence
// ----------------------------------------------------------------------------

//NOTE(h-2): shouldnt we be using Is<StringConcept> or IsContiguous<> or ...?

// template <typename TChar, typename TAlloc>
// struct IsSequence<std::vector<TChar, TAlloc> > : True {};
//
// template <typename TChar, typename TCharTraits, typename TAlloc>
// struct IsSequence<std::basic_string<TChar, TCharTraits, TAlloc> > : True {};

// ----------------------------------------------------------------------------
// Mfn Chunk
// ----------------------------------------------------------------------------

// #define ITRMACRO__(CONT, CONST, TMPL) \
// template <TMPL> \
// struct Chunk<Iter<CONT CONST, AdaptorIterator<TChar*, TSpec> > > \
// { \
//     typedef Iter<CONT CONST, AdaptorIterator<Iter<CONT CONST, StdIteratorAdaptor> > > Type;\
// };

// ----------------------------------------------------------------------------
// Mfn HasCapacity_ (private)
// ----------------------------------------------------------------------------

// template <typename TContainer>
// struct HasCapacity_ : public False {};
//
// template <typename TChar, typename TAlloc>
// struct HasCapacity_<std::vector<TChar, TAlloc> > : public True {};
//
// template <typename TChar, typename TAlloc>
// struct HasCapacity_<std::vector<TChar, TAlloc> const> : public True {};
//
// template <typename TChar, typename TTraits, typename TAlloc>
// struct HasCapacity_<std::basic_string<TChar, TTraits, TAlloc> > : public True {};
//
// template <typename TChar, typename TTraits, typename TAlloc>
// struct HasCapacity_<std::basic_string<TChar, TTraits, TAlloc> const> : public True {};

// ----------------------------------------------------------------------------
// Mfn FixedSize_ (private)
// ----------------------------------------------------------------------------

template <typename TContainer>
struct FixedSize_ : public False {};

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, size_t N>
struct FixedSize_<std::array<TChar, N> > : public True {};
#endif

// ===========================================================================
// Functions
// ===========================================================================

// ----------------------------------------------------------------------------
// Function getObjectId
// ----------------------------------------------------------------------------

// default value for SFINAE type in sequence_forwards.h
template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type>
inline void const *
getObjectId(TContainer SEQAN_FORWARD_CARG me)
{
    if (me.empty())
        return NULL;
    else
        return (& *(me.end() - 1)) + 1;
}

// ----------------------------------------------------------------------------
// Function begin (standard)
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Iterator<TContainer, Standard>::Type
begin(TContainer & me, Standard const &)
{
//     return typename Iterator<TContainer, Standard>::Type(me.begin());
    return me.begin();
}

template <typename TContainer, typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Iterator<TContainer const, Standard>::Type
begin(TContainer const & me, Standard const &)
{
//     return typename Iterator<TContainer const, Standard>::Type(me.begin());
    return me.begin();
}

// ----------------------------------------------------------------------------
// Function begin (rooted)
// ----------------------------------------------------------------------------

// template <typename TContainer,typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
// inline typename Iterator<TContainer, Rooted>::Type
// begin(TContainer & me, Rooted const &)
// {
// //     static_assert(std::is_same<typename Iterator<TContainer, Rooted>::Type, Iter<TContainer, StdIteratorAdaptor> >::value, "AHH");
//     return typename Iterator<TContainer, Rooted>::Type(begin(me, Standard()));
// }
//
// template <typename TContainer,typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
// inline typename Iterator<TContainer const, Rooted>::Type
// begin(TContainer const & me, Rooted const &)
// {
// //     static_assert(std::is_same<typename Iterator<TContainer const, Rooted>::Type, Iter<TContainer const, StdIteratorAdaptor> >::value, "AHH");
//     return typename Iterator<TContainer const, Rooted>::Type(begin(me, Standard()));
// }

// ----------------------------------------------------------------------------
// Function end (standard)
// ----------------------------------------------------------------------------

template <typename TContainer, typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Iterator<TContainer, Standard>::Type
end(TContainer & me, Standard const &)
{
//     return typename Iterator<TContainer, Standard>::Type(me.end());
    return me.end();
}

template <typename TContainer, typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Iterator<TContainer const, Standard>::Type
end(TContainer const & me, Standard const &)
{
//     return typename Iterator<TContainer const, Standard>::Type(me.end());
    return me.end();
}

// // ----------------------------------------------------------------------------
// // Function end (rooted)
// // ----------------------------------------------------------------------------
//
// template <typename TContainer, typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
// inline typename Iterator<TContainer, Rooted>::Type
// end(TContainer & me, Rooted const &)
// {
//     return typename Iterator<TContainer, Rooted>::Type(end(me, Standard()));
// }
//
// template <typename TContainer, typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
// inline typename Iterator<TContainer const, Rooted>::Type
// end(TContainer const & me, Rooted const &)
// {
//     return typename Iterator<TContainer const, Rooted>::Type(end(me, Standard()));
// }

// ----------------------------------------------------------------------------
// Function length
// ----------------------------------------------------------------------------

// default value for SFINAE type in sequence_forwards.h
template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type>
inline typename Size<typename RemoveReference<TContainer>::Type>::Type
length(TContainer const & me)
{
    return me.size();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, std::size_t N>
constexpr std::size_t
length(std::array<TChar, N> const & me)
{
    return me.size();
}
#endif

// ----------------------------------------------------------------------------
// Function capacity
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline typename Size<std::vector<TChar, TAlloc> >::Type
capacity(std::vector<TChar, TAlloc> const & me)
{
    return me.capacity();
}

template <typename TChar, typename TTraits, typename TAlloc>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
capacity(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return me.capacity();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, std::size_t N>
constexpr std::size_t
capacity(std::array<TChar, N> const & me)
{
    return me.size();
}
#endif

// for other types the default overload to length holds

// ----------------------------------------------------------------------------
// Function empty
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline bool
empty(TContainer const & me)
{
    return me.empty();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, std::size_t N>
constexpr bool
empty(std::array<TChar, N> const & me)
{
    return me.empty();
}
#endif

// ----------------------------------------------------------------------------
// Function reserve
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TSize>
inline typename Size<std::vector<TChar, TAlloc> >::Type
reserve(std::vector<TChar, TAlloc> & me,
        TSize const s)
{
    me.reserve(s);
    return capacity(me);
}

template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::vector<TChar, TAlloc> >::Type
reserve(std::vector<TChar, TAlloc> & me,
        TSize const s,
        Tag<TExpand> const &)
{
    return reserve(me, s);
}

template <typename TChar, typename TTraits, typename TSize, typename TAlloc>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
reserve(std::basic_string<TChar, TTraits, TAlloc> & me,
        TSize const s)
{
    me.reserve(s);
    return capacity(me);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
reserve(std::basic_string<TChar, TTraits, TAlloc> & me,
        TSize const s,
        Tag<TExpand> const &)
{
    return reserve(me, s);
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc, typename TSize>
inline typename Size<std::vector<TChar, TAlloc> >::Type
reserve(std::vector<TChar, TAlloc> && me,
        TSize const s)
{
    me.reserve(s);
    return capacity(me);
}

template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::vector<TChar, TAlloc> >::Type
reserve(std::vector<TChar, TAlloc> && me,
        TSize const s,
        Tag<TExpand> const &)
{
    return reserve(std::move(me), s);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSize>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
reserve(std::basic_string<TChar, TTraits, TAlloc> && me,
        TSize const s)
{
    me.reserve(s);
    return capacity(me);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<std::basic_string<TChar, TTraits, TAlloc> >::Type
reserve(std::basic_string<TChar, TTraits, TAlloc> && me,
        TSize const s,
        Tag<TExpand> const &)
{
    return reserve(std::move(me), s);
}
#endif

// for other types the default overload holds

// ----------------------------------------------------------------------------
// Function resize
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSize,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Size<typename RemoveReference<TContainer>::Type>::Type
resize(TContainer SEQAN_FORWARD_ARG me,
        TSize const s)
{
    me.resize(s);
    return length(me);
}

template <typename TContainer,
          typename TSize,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Size<typename RemoveReference<TContainer>::Type>::Type
resize(TContainer SEQAN_FORWARD_ARG me,
       TSize const s,
       Tag<TExpand> const &)
{
    return resize(SEQAN_FORWARD(TContainer, me), s);
}

template <typename TContainer,
          typename TSize,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Size<typename RemoveReference<TContainer>::Type>::Type
resize(TContainer SEQAN_FORWARD_ARG me,
       TSize const s,
       typename Value<typename RemoveReference<TContainer>::Type>::Type const & val)
{
    me.resize(s, val);
    return length(me);
}

template <typename TContainer,
          typename TSize,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Size<typename RemoveReference<TContainer>::Type>::Type
resize(TContainer SEQAN_FORWARD_ARG me,
       TSize const s,
       typename Value<typename RemoveReference<TContainer>::Type>::Type const & val,
       Tag<TExpand> const &)
{
    return resize(SEQAN_FORWARD(TContainer, me), s, val);
}

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
clear(TContainer SEQAN_FORWARD_ARG me)
{
    me.clear();
}

// ----------------------------------------------------------------------------
// Function value
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                IsContiguous<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Reference<typename RemoveReference<TContainer>::Type>::Type
value(TContainer & me, TPos const pos)
{
    return me[pos];
}

template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                IsContiguous<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Reference<TContainer const>::Type
value(TContainer const & me, TPos const pos)
{
    return me[pos];
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                IsContiguous<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Value<typename RemoveReference<TContainer>::Type>::Type
value(TContainer && me, TPos const pos)
{
    return me[pos];
}
#endif

// linear complexity for list and fwd list
template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<IsContiguous<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Reference<typename RemoveReference<TContainer>::Type>::Type
value(TContainer & me, TPos const pos)
{
    return *(begin(me) + pos);
}

template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<IsContiguous<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Reference<TContainer const>::Type
value(TContainer const & me, TPos const pos)
{
    return *(begin(me) + pos);
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<IsContiguous<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Value<typename RemoveReference<TContainer>::Type>::Type
value(TContainer && me, TPos const pos)
{
    return *(begin(me) + pos);
}
#endif

// --------------------------------------------------------------------------
// Function getValue
// --------------------------------------------------------------------------

template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Reference<TContainer const>::Type
getValue(TContainer const & me, TPos const pos)
{
    return me[pos];
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Value<typename RemoveReference<TContainer>::Type>::Type
getValue(TContainer && me, TPos const pos)
{
    return me[pos];
}
#endif

// linear complexity for list and fwd list
template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<IsContiguous<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Reference<TContainer const>::Type
getValue(TContainer const & me, TPos const pos)
{
    return *(begin(me) + pos);
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<IsContiguous<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline typename Value<typename RemoveReference<TContainer>::Type>::Type
getValue(TContainer && me, TPos const pos)
{
    return *(begin(me) + pos);
}
#endif

// ----------------------------------------------------------------------------
// Function front
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Reference<typename RemoveReference<TContainer>::Type>::Type
front(TContainer & me)
{
    return me.front();
}

template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Reference<TContainer const>::Type
front(TContainer const & me)
{
    return me.front();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Value<typename RemoveReference<TContainer>::Type>::Type
front(TContainer && me)
{
    return me.front();
}
#endif

// pre-c++11 stdstring doesnt have front and back
#ifndef SEQAN_CXX11_STANDARD
template <typename TChar, typename TTraits, typename TAlloc, typename TExpand>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> >::Type
front(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    return *(end(me) - 1);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TExpand>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> const>::Type
front(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return *(end(me) - 1);
}
#endif

// ----------------------------------------------------------------------------
// Function back
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Reference<typename RemoveReference<TContainer>::Type>::Type
back(TContainer & me)
{
    return me.back();
}

template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Reference<TContainer const>::Type
back(TContainer const & me)
{
    return me.back();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline typename Value<typename RemoveReference<TContainer>::Type>::Type
back(TContainer && me)
{
    return me.back();
}

// forward_list doesnt have back, we achieve it in linear time
template <typename TChar, typename TAlloc>
inline TChar &
back(std::forward_list<TChar, TAlloc> & me)
{
    return *(begin(me) + (length(me) - 1));
}

template <typename TChar, typename TAlloc>
inline TChar const &
back(std::forward_list<TChar, TAlloc> const & me)
{
    return *(begin(me) + (length(me) - 1));
}

template <typename TChar, typename TAlloc>
inline TChar
back(std::forward_list<TChar, TAlloc> && me)
{
    return *(begin(me) + (length(me) - 1));
}
#endif

// pre-c++11 stdstring doesnt have back and back
#ifndef SEQAN_CXX11_STANDARD
template <typename TChar, typename TTraits, typename TAlloc, typename TExpand>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> >::Type
back(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    return *(end(me) - 1);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TExpand>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> const>::Type
back(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return *(end(me) - 1);
}
#endif

// ----------------------------------------------------------------------------
// Function assign
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline void
assign(TContainer SEQAN_FORWARD_ARG me,
       TContainer SEQAN_FORWARD_CARG source)
{
    me = SEQAN_FORWARD(TContainer, source);
}

template <typename TContainer,
          typename TSource,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline void
assign(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source)
{
    me.assign(begin(source, Standard()), end(source, Standard()));
}

template <typename TContainer,
          typename TSource,
          typename EnableIf<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, int>::Type = 0>
inline void
assign(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source,
       typename Size<TSource>::Type limit)
{
    me.assign(begin(source, Standard()), begin(source, Standard()) + limit);
}

// ----------------------------------------------------------------------------
// Function insert
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insert(TContainer SEQAN_FORWARD_ARG me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source)
{
    me.insert(me.begin() + pos, begin(source, Standard()), end(source, Standard()));
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insert(TContainer SEQAN_FORWARD_ARG me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source,
       Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), pos, source);
}

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insert(TContainer SEQAN_FORWARD_ARG me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source,
       typename Size<TSource>::Type const limit)
{
    me.insert(me.begin() + pos, begin(source, Standard()), begin(source, Standard()) + limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insert(TContainer SEQAN_FORWARD_ARG me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source,
       typename Size<TSource>::Type const limit,
       Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), pos, source, limit);
}

#ifdef SEQAN_CXX11_STANDARD
// forward_list doesnt have insert, we achieve it slower
template <typename TChar, typename TAlloc, typename TSource>
inline void
insert(std::forward_list<TChar, TAlloc> & me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source)
{
    me.insert_after(me.before_begin() + pos, begin(source), end(source));
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
insert(std::forward_list<TChar, TAlloc> && me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source)
{
    insert(me, pos, source);
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
insert(std::forward_list<TChar, TAlloc> & me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const limit)
{
    me.insert_after(me.before_begin() + pos, begin(source), begin(source) + limit);
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
insert(std::forward_list<TChar, TAlloc> && me,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
       TSource const & source,
       typename Size<std::forward_list<TChar, TAlloc> >::Type const limit)
{
    insert(me, pos, source, limit);
}
#endif

// ----------------------------------------------------------------------------
// Function append
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
append(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source)
{
    insert(SEQAN_FORWARD(TContainer, me), length(me), source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
append(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source,
       Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), length(me), source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
append(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const limit)
{
    insert(SEQAN_FORWARD(TContainer, me), length(me), source, limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
append(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const limit,
       Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), length(me), source, limit);
}

// forward list is handled downstream in insert

// ----------------------------------------------------------------------------
// Function prepend
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prepend(TContainer SEQAN_FORWARD_ARG me,
        TSource const & source)
{
    insert(SEQAN_FORWARD(TContainer, me), 0, source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prepend(TContainer SEQAN_FORWARD_ARG me,
        TSource const & source,
        Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), 0, source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prepend(TContainer SEQAN_FORWARD_ARG me,
        TSource const & source,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const limit)
{
    insert(SEQAN_FORWARD(TContainer, me), 0, source, limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prepend(TContainer SEQAN_FORWARD_ARG me,
        TSource const & source,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const limit,
        Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), 0, source, limit);
}

// forward list is handled downstream in insert

// ----------------------------------------------------------------------------
// Function insertValue
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insertValue(TContainer SEQAN_FORWARD_ARG me,
            typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
            TSource const & source)
{
    me.insert(begin(me) + pos, source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insertValue(TContainer SEQAN_FORWARD_ARG me,
            typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
            TSource const & source,
            Tag<TExpand> const &)
{
    insertValue(SEQAN_FORWARD(TContainer, me), pos, source);
}


#ifdef SEQAN_CXX11_STANDARD
// move semantics with double universal references, yay!
template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insertValue(TContainer && me,
            typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
            TSource && source)
{
    me.insert(begin(me) + pos, std::forward<TSource>(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
insertValue(TContainer && me,
            typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
            TSource && source,
            Tag<TExpand> const &)
{
    insertValue(SEQAN_FORWARD(TContainer, me), pos, SEQAN_FORWARD(TSource, source));
}

// forward_list doesnt have push_back, we achieve it in linear time
template <typename TChar, typename TAlloc, typename TSource>
inline void
insertValue(std::forward_list<TChar, TAlloc> & me,
            typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
            TSource const & source)
{
    me.insert_after(me.before_begin() + pos, source);
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
insertValue(std::forward_list<TChar, TAlloc> & me,
            typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
            TSource && source)
{
    me.insert_after(me.before_begin() + pos, std::forward<TSource>(source));
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
insertValue(std::forward_list<TChar, TAlloc> && me,
            typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
            TSource const & source)
{
    insertValue(me, pos, source);
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
insertValue(std::forward_list<TChar, TAlloc> && me,
            typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
            TSource && source)
{
    insertValue(me, pos, source);
}
#endif

// ----------------------------------------------------------------------------
// Function appendValue
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
appendValue(TContainer SEQAN_FORWARD_ARG me,
            TSource const & source)
{
    me.push_back(source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
appendValue(TContainer SEQAN_FORWARD_ARG me,
            TSource const & source,
            Tag<TExpand> const &)
{
    appendValue(SEQAN_FORWARD(TContainer, me), source);
}

#ifdef SEQAN_CXX11_STANDARD
// move semantics with double universal references, yay!
template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
appendValue(TContainer && me,
            TSource && source)
{
    me.push_back(std::forward<TSource>(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
appendValue(TContainer && me,
            TSource && source,
            Tag<TExpand> const &)
{
    me.push_back(std::forward<TSource>(source));
}

// forward_list doesnt have push_back, we achieve it in linear time
template <typename TChar, typename TAlloc, typename TSource>
inline void
appendValue(std::forward_list<TChar, TAlloc> & me,
            TSource const & source)
{
    insertValue(me, length(me), source);
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
appendValue(std::forward_list<TChar, TAlloc> & me,
            TSource && source)
{
    insertValue(me, length(me), std::forward<TSource>(source));
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
appendValue(std::forward_list<TChar, TAlloc> && me,
            TSource const & source)
{
    appendValue(me, source);
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
appendValue(std::forward_list<TChar, TAlloc> && me,
            TSource && source)
{
    appendValue(me, source);
}
#endif

// ----------------------------------------------------------------------------
// Function prependValue
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prependValue(TContainer SEQAN_FORWARD_ARG me,
             TSource const & source)
{
    me.push_front(source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prependValue(TContainer SEQAN_FORWARD_ARG me,
             TSource const & source,
             Tag<TExpand> const &)
{
    prependValue(SEQAN_FORWARD(TContainer, me), source);
}

// vector and string dont have push_front, we achieve it in linear time
template <typename TChar, typename TAlloc, typename TSource>
inline void
prependValue(std::vector<TChar, TAlloc> & me,
             TSource const & source)
{
    insertValue(me, 0, source);
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSource>
inline void
prependValue(std::basic_string<TChar, TTraits, TAlloc> & me,
             TSource const & source)
{
    insertValue(me, 0, source);
}

#ifdef SEQAN_CXX11_STANDARD
// move semantics with double universal references, yay!
template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prependValue(TContainer && me,
             TSource && source)
{
    me.push_front(std::forward<TSource>(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
prependValue(TContainer && me,
            TSource && source,
            Tag<TExpand> const &)
{
    prependValue(SEQAN_FORWARD(TContainer, me), SEQAN_FORWARD(TSource, source));
}

// vector and string dont have push_front, we achieve it in linear time
template <typename TChar, typename TAlloc, typename TSource>
inline void
prependValue(std::vector<TChar, TAlloc> & me,
             TSource && source)
{
    insertValue(me, 0, std::forward<TSource>(source));
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSource>
inline void
prependValue(std::basic_string<TChar, TTraits, TAlloc> & me,
             TSource && source)
{
    insertValue(me, 0, std::forward<TSource>(source));
}

// && to &
template <typename TChar, typename TAlloc, typename TSource>
inline void
prependValue(std::vector<TChar, TAlloc> && me,
             TSource && source)
{
    prependValue(me, std::forward<TSource>(source));
}

template <typename TChar, typename TTraits, typename TAlloc, typename TSource>
inline void
prependValue(std::basic_string<TChar, TTraits, TAlloc> && me,
             TSource && source)
{
    prependValue(me, std::forward<TSource>(source));
}
#endif

// ----------------------------------------------------------------------------
// Function erase
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
erase(TContainer SEQAN_FORWARD_ARG me,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd)
{
    me.erase(me.begin() + pos, me.begin() + posEnd);
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
inline void
erase(std::forward_list<TChar, TAlloc> & me,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const posEnd)
{
    me.erase_after(me.begin() + (pos - 1), me.begin() + (posEnd - 1));
}

template <typename TChar, typename TAlloc>
inline void
erase(std::forward_list<TChar, TAlloc> && me,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const posEnd)
{
    erase(me, pos, posEnd);
}
#endif

// ----------------------------------------------------------------------------
// Function eraseFront
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
eraseFront(TContainer SEQAN_FORWARD_ARG me)
{
    me.pop_front();
}

template <typename TChar, typename TAlloc>
inline void
eraseFront(std::vector<TChar, TAlloc> & me)
{
    me.erase(me.begin());
}

template <typename TChar, typename TTraits, typename TAlloc>
inline void
eraseFront(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    me.erase(me.begin());
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
inline void
eraseFront(std::vector<TChar, TAlloc> && me)
{
    me.erase(me.begin());
}

template <typename TChar, typename TTraits, typename TAlloc>
inline void
eraseFront(std::basic_string<TChar, TTraits, TAlloc> && me)
{
    me.erase(me.begin());
}
#endif

// ----------------------------------------------------------------------------
// Function eraseBack
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
eraseBack(TContainer SEQAN_FORWARD_ARG me)
{
    me.pop_back();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
inline void
eraseBack(std::forward_list<TChar, TAlloc> & me)
{
    if (length(me) > 0)
        me.erase_after(me.begin() + (length(me) - 1));
}

template <typename TChar, typename TAlloc>
inline void
eraseBack(std::forward_list<TChar, TAlloc> && me)
{
    if (length(me) > 0)
        me.erase_after(me.begin() + (length(me) - 1));
}
#endif

#ifndef SEQAN_CXX11_STANDARD
template <typename TChar, typename TTraits, typename TAlloc>
inline void
eraseBack(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    if (length(me) > 0)
        me.erase(me.end() - 1);
}

template <typename TChar, typename TTraits, typename TAlloc>
inline void
eraseBack(std::basic_string<TChar, TTraits, TAlloc> && me)
{
    if (length(me) > 0)
        me.erase(me.end() - 1);
}
#endif


// ----------------------------------------------------------------------------
// Function replace
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
replace(TContainer SEQAN_FORWARD_ARG me,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd,
        TSource const & source,
        typename Size<TSource>::Type const limit)
{
    erase(me, pos, posEnd);
    insert(me, pos, source, limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
replace(TContainer SEQAN_FORWARD_ARG me,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd,
        TSource const & source,
        typename Size<TSource>::Type const limit,
        Tag<TExpand> const &)
{
    replace(SEQAN_FORWARD(TContainer, me), pos, posEnd, source, limit);
}

template <typename TContainer,
          typename TSource,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
replace(TContainer SEQAN_FORWARD_ARG me,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd,
        TSource const & source)
{
    erase(me, pos, posEnd);
    insert(me, pos, source, length(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand,
          typename EnableIf<And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, int>::Type = 0>
inline void
replace(TContainer SEQAN_FORWARD_ARG me,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd,
        TSource const & source,
        Tag<TExpand> const &)
{
    replace(SEQAN_FORWARD(TContainer, me), pos, posEnd, source);
}

// ----------------------------------------------------------------------------
// Function toCString
// ----------------------------------------------------------------------------

template <typename TChar, typename TTraits, typename TAlloc>
inline const TChar*
toCString(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    return me.c_str();
}

}

#undef SUPERMACRO__
#undef ITMACRO__
#undef ITRMACRO__
#undef COMMA

#endif
