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

/* NOTE(h-2) on ConceptChecking and universal references
 *
 * When type deduction takes place before && then && can act as both && or &
 * [https://isocpp.org/blog/2012/11/universal-references-in-c11-scott-meyers]
 *
 * This has the effect that in
 *      foobar(T &&)
 * T actually can be "int &". Consequently the above becomes "int & &&" which
 * collapses to "int &". The important thing is that in
 *      template<typename EnableIf<Is<SuperConcept<T> >
 * this will result in the ConceptCheck always failing (because the Concept is
 * not defined for reference types).
 *
 * So in the following you will see a RemoveReference< > around the type
 * whenever the function takes a && argument or SEQAN_FORWARD_ARG argument AND
 * it shall act as univeral reference (e.g. length()). But in some cases it
 * shall not (because there is a different overload for the l-value-reference
 * arg) and it is  intentionally omitted, e.g. value().
 */

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
// Mfn HasSubscriptOperator (different for e.g. std::deque)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct HasSubscriptOperator<std::deque<TChar, TAlloc> > :
    public True
{};

template <typename TChar, typename TAlloc>
struct HasSubscriptOperator<std::deque<TChar, TAlloc> const> :
    public True
{};

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
ITRMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc)

ITRMACRO__(std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc)
ITRMACRO__(std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc)
ITRMACRO__(std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc)
#ifdef SEQAN_CXX11_STANDARD
ITRMACRO__(std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc)
ITRMACRO__(std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N)
#endif
ITRMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc)

// ----------------------------------------------------------------------------
// Mfn IsSequence
// ----------------------------------------------------------------------------

//NOTE(h-2): shouldnt we be using Is<StringConcept> or IsContiguous<> or ...?

template <typename TChar, typename TAlloc>
struct IsSequence<std::vector<TChar, TAlloc> > : True {};

template <typename TChar, typename TAlloc>
struct IsSequence<std::deque<TChar, TAlloc> > : True {};

template <typename TChar, typename TAlloc>
struct IsSequence<std::list<TChar, TAlloc> > : True {};

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
struct IsSequence<std::forward_list<TChar, TAlloc> > : True {};

//NOTE(h-2): array isn't by SeqAn definition, right?
#endif

template <typename TChar, typename TCharTraits, typename TAlloc>
struct IsSequence<std::basic_string<TChar, TCharTraits, TAlloc> > : True {};

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
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void const *)
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

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Standard>::Type)
begin(TContainer & me, Standard const &)
{
    return me.begin();
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Standard>::Type)
begin(TContainer const & me, Standard const &)
{
    return me.begin();
}

// ----------------------------------------------------------------------------
// Function begin (rooted)
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Rooted>::Type)
begin(TContainer & me, Rooted const &)
{
    return _beginDefault(me, Rooted());
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Rooted>::Type)
begin(TContainer const & me, Rooted const &)
{
    return _beginDefault(me, Rooted());
}

// ----------------------------------------------------------------------------
// Function end (standard)
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Standard>::Type)
end(TContainer & me, Standard const &)
{
    return me.end();
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Standard>::Type)
end(TContainer const & me, Standard const &)
{
    return me.end();
}

// ----------------------------------------------------------------------------
// Function end (rooted)
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer, Rooted>::Type)
end(TContainer & me, Rooted const &)
{
    return _endDefault(me, Rooted());
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Iterator<TContainer const, Rooted>::Type)
end(TContainer const & me, Rooted const &)
{
    return _endDefault(me, Rooted());
}

// ----------------------------------------------------------------------------
// Function _iterStl (like seqan's iter() but for stl-iterators)
// ----------------------------------------------------------------------------

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos>
inline typename TContainer::iterator
_iterStl(TContainer & me,
         TPos const pos)
{
    SEQAN_ASSERT_LEQ_MSG(pos,
                         static_cast<TPos>(length(me)),
                         "Trying to get an iterator behind a container through _iterStl().");
    return std::next(me.begin(), pos);
}

template <typename TContainer,
          typename TPos>
inline typename TContainer::const_iterator
_iterStl(TContainer const & me,
         TPos const pos)
{
    SEQAN_ASSERT_LEQ_MSG(pos,
                         static_cast<TPos>(length(me)),
                         "Trying to get an iterator behind a container through _iterStl().");
    return std::next(me.begin(), pos);
}
#else
template <typename TContainer,
          typename TPos>
inline typename TContainer::iterator
_iterStl(TContainer & me,
         TPos const pos)
{
    SEQAN_ASSERT_LEQ_MSG(pos,
                         static_cast<TPos>(length(me)),
                         "Trying to get an iterator behind a container through _iterStl().");
    typename TContainer::iterator it = me.begin();
    std::advance(it, pos);
    return it;
}

template <typename TContainer,
          typename TPos>
inline typename TContainer::const_iterator
_iterStl(TContainer const & me,
         TPos const pos)
{
    SEQAN_ASSERT_LEQ_MSG(pos,
                         static_cast<TPos>(length(me)),
                         "Trying to get an iterator behind a container through _iterStl().");
    typename TContainer::const_iterator it = me.begin();
    std::advance(it, pos);
    return it;
}
#endif

// ----------------------------------------------------------------------------
// Function iter (Standard) (only overloaded for non-contiguous)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc>, Standard>::Type
iter(std::list<TChar, TAlloc> & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc> const, Standard>::Type
iter(std::list<TChar, TAlloc> const & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc>, Standard const>::Type
iter(std::forward_list<TChar, TAlloc> & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc> const, Standard>::Type
iter(std::forward_list<TChar, TAlloc> const & me,
     TPos const pos,
     Standard const &)
{
    return _iterStl(me, pos);
}
#endif

// ----------------------------------------------------------------------------
// Function iter (Rooted) (only overloaded for non-contiguous)
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc>, Rooted>::Type
iter(std::list<TChar, TAlloc> & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::list<TChar, TAlloc>, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::list<TChar, TAlloc> const, Rooted>::Type
iter(std::list<TChar, TAlloc> const & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::list<TChar, TAlloc> const, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc>, Rooted>::Type
iter(std::forward_list<TChar, TAlloc> & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::forward_list<TChar, TAlloc>, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}

template <typename TChar, typename TAlloc, typename TPos>
inline typename Iterator<std::forward_list<TChar, TAlloc> const, Rooted>::Type
iter(std::forward_list<TChar, TAlloc> const & me,
     TPos const pos,
     Rooted const &)
{
    typedef typename Iterator<std::forward_list<TChar, TAlloc> const, Rooted>::Type TIterator;
    return TIterator(me, iter(me, pos, Standard()));
}
#endif

// ----------------------------------------------------------------------------
// Function length
// ----------------------------------------------------------------------------

// default value for SFINAE type in sequence_forwards.h
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Size<TContainer>::Type)
length(TContainer const & me)
{
    return me.size();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
inline typename Size<std::forward_list<TChar, TAlloc> >::Type
length(std::forward_list<TChar, TAlloc> const & me)
{
    typename Size<std::forward_list<TChar, TAlloc> >::Type l = 0;

    for (auto it = me.begin(), itEnd = me.end(); it != itEnd; ++it, ++l);

    return l;
}

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

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, bool)
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
          typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >,
                            typename Size<typename RemoveReference<TContainer>::Type>::Type)
resize(TContainer SEQAN_FORWARD_ARG me,
        TSize const s)
{
    me.resize(s);
    return length(me);
}

template <typename TContainer,
          typename TSize,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >,
                            typename Size<typename RemoveReference<TContainer>::Type>::Type)
resize(TContainer SEQAN_FORWARD_ARG me,
       TSize const s,
       Tag<TExpand> const &)
{
    return resize(SEQAN_FORWARD(TContainer, me), s);
}

template <typename TContainer,
          typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >,
                            typename Size<typename RemoveReference<TContainer>::Type>::Type)
resize(TContainer SEQAN_FORWARD_ARG me,
       TSize const s,
       typename Value<typename RemoveReference<TContainer>::Type>::Type const & val)
{
    me.resize(s, val);
    return length(me);
}

template <typename TContainer,
          typename TSize,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >,
                            typename Size<typename RemoveReference<TContainer>::Type>::Type)
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

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
clear(TContainer SEQAN_FORWARD_ARG me)
{
    me.clear();
}

// ----------------------------------------------------------------------------
// Function value
// ----------------------------------------------------------------------------

// default value for SFINAE type in sequence_forwards.h
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Reference<TContainer>::Type)
value(TContainer & me, TPos const pos)
{
    return me[pos];
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Reference<TContainer const>::Type)
value(TContainer const & me, TPos const pos)
{
    return me[pos];
}

#ifdef SEQAN_CXX11_STANDARD
// RemoveReference explicitly not set on TContainer in SFINAE-clause, because this function shall not
// pick up on TContainer & &&, i.e. the && shall act as strict R-Value reference, not universal
// [same below for some overloads]
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Value<TContainer>::Type)
value(TContainer && me, TPos const pos)
{
    return me[pos];
}
#endif

// linear complexity for list and fwd list
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer>::Type)
value(TContainer & me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer const>::Type)
value(TContainer const & me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Value<TContainer>::Type)
value(TContainer && me, TPos const pos)
{
    return *(_iterStl(me, pos));
}
#endif

// --------------------------------------------------------------------------
// Function getValue
// --------------------------------------------------------------------------

template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Reference<TContainer const>::Type)
getValue(TContainer const & me, TPos const pos)
{
    return me[pos];
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                HasSubscriptOperator<TContainer> >, typename Value<TContainer>::Type)
getValue(TContainer && me, TPos const pos)
{
    return me[pos];
}
#endif

// linear complexity for list and fwd list
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer const>::Type)
getValue(TContainer const & me, TPos const pos)
{
    return *(_iterStl(me, pos));
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename TPos>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >,
                                Not<HasSubscriptOperator<TContainer> > >, typename Value<TContainer>::Type)
getValue(TContainer && me, TPos const pos)
{
    return *(_iterStl(me, pos));
}
#endif

// ----------------------------------------------------------------------------
// Function front
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Reference<TContainer>::Type)
front(TContainer & me)
{
    return me.front();
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Reference<TContainer const>::Type)
front(TContainer const & me)
{
    return me.front();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer,
          typename EnableIf<Is<StlContainerConcept<TContainer> >, int>::Type = 0>
inline typename Value<TContainer>::Type
front(TContainer && me)
{
    return me.front();
}
#endif

// pre-c++11 stdstring doesnt have front and back
#ifndef SEQAN_CXX11_STANDARD
template <typename TChar, typename TTraits, typename TAlloc>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> >::Type
front(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    return *(me.begin());
}

template <typename TChar, typename TTraits, typename TAlloc>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> const>::Type
front(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return *(me.begin());
}
#endif

// ----------------------------------------------------------------------------
// Function back
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Reference<TContainer>::Type)
back(TContainer & me)
{
    return me.back();
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Reference<TContainer const>::Type)
back(TContainer const & me)
{
    return me.back();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Value<TContainer>::Type)
back(TContainer && me)
{
    return me.back();
}

// forward_list doesnt have back, we achieve it in linear time
template <typename TChar, typename TAlloc>
inline TChar &
back(std::forward_list<TChar, TAlloc> & me)
{
    return *(std::next(me.before_begin(), length(me)));
}

template <typename TChar, typename TAlloc>
inline TChar const &
back(std::forward_list<TChar, TAlloc> const & me)
{
    return *(std::next(me.before_begin(), length(me)));
}

template <typename TChar, typename TAlloc>
inline TChar
back(std::forward_list<TChar, TAlloc> && me)
{
    return *(std::next(me.before_begin(), length(me)));
}
#endif

// pre-c++11 stdstring doesnt have back and back
#ifndef SEQAN_CXX11_STANDARD
template <typename TChar, typename TTraits, typename TAlloc>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> >::Type
back(std::basic_string<TChar, TTraits, TAlloc> & me)
{
    return *(me.end() - 1);
}

template <typename TChar, typename TTraits, typename TAlloc>
inline typename Reference<std::basic_string<TChar, TTraits, TAlloc> const>::Type
back(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return *(me.end() - 1);
}
#endif

// ----------------------------------------------------------------------------
// Function assign
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, void)
assign(TContainer & me,
       TContainer const & source)
{
    me = source;
}
#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, void)
assign(TContainer & me,
       TContainer && source)
{
    me = std::move(source);
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, void)
assign(TContainer && me,
       TContainer const & source)
{
    assign(me, source);
}

// double universal reference doesn't work above, because it would be ambiguous
// with the fundamental_assign implementation on gcc

#endif
template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void)
assign(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source)
{
    me.assign(begin(source, Standard()), end(source, Standard()));
}

template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void)
assign(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source,
       typename Size<TSource>::Type limit)
{
    me.assign(begin(source, Standard()),
              iter(source, std::min(length(source), limit), Standard()));
}

// ----------------------------------------------------------------------------
// Function insert
// ----------------------------------------------------------------------------

template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                            Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insert(TContainer SEQAN_FORWARD_ARG me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source)
{
    me.insert(_iterStl(SEQAN_FORWARD(TContainer, me), pos),
              begin(source, Standard()),
              end(source, Standard()));
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insert(TContainer SEQAN_FORWARD_ARG me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source,
       Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), pos, source);
}

template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insert(TContainer SEQAN_FORWARD_ARG me,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
       TSource const & source,
       typename Size<TSource>::Type const limit)
{
    me.insert(_iterStl(SEQAN_FORWARD(TContainer, me), pos),
              begin(source, Standard()),
              iter(source, std::min(length(source), limit), Standard()));
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
    me.insert_after(std::next(me.before_begin(), pos), begin(source, Standard()), end(source, Standard()));
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
    me.insert_after(std::next(me.before_begin(), pos),
                    begin(source, Standard()),
                    iter(source, std::min(length(source), limit), Standard()));
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                            Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
append(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source)
{
    insert(SEQAN_FORWARD(TContainer, me), length(me), source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
append(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source,
       Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), length(me), source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
append(TContainer SEQAN_FORWARD_ARG me,
       TSource const & source,
       typename Size<typename RemoveReference<TContainer>::Type>::Type const limit)
{
    insert(SEQAN_FORWARD(TContainer, me), length(me), source, limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prepend(TContainer SEQAN_FORWARD_ARG me,
        TSource const & source)
{
    insert(SEQAN_FORWARD(TContainer, me), 0, source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prepend(TContainer SEQAN_FORWARD_ARG me,
        TSource const & source,
        Tag<TExpand> const &)
{
    insert(SEQAN_FORWARD(TContainer, me), 0, source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prepend(TContainer SEQAN_FORWARD_ARG me,
        TSource const & source,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const limit)
{
    insert(SEQAN_FORWARD(TContainer, me), 0, source, limit);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insertValue(TContainer SEQAN_FORWARD_ARG me,
            typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
            TSource const & source)
{
    me.insert(_iterStl(SEQAN_FORWARD(TContainer, me), pos), source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
insertValue(TContainer && me,
            typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
            TSource && source)
{
    me.insert(_iterStl(SEQAN_FORWARD(TContainer, me), pos), std::forward<TSource>(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
    me.insert_after(std::next(me.before_begin(), pos), source);
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
insertValue(std::forward_list<TChar, TAlloc> & me,
            typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
            TSource && source)
{
    me.insert_after(std::next(me.before_begin(), pos), std::forward<TSource>(source));
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
appendValue(TContainer SEQAN_FORWARD_ARG me,
            TSource const & source)
{
    me.push_back(source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
appendValue(TContainer SEQAN_FORWARD_ARG me,
            TSource const & source,
            Tag<TExpand> const &)
{
    appendValue(SEQAN_FORWARD(TContainer, me), source);
}

#ifdef SEQAN_CXX11_STANDARD
// move semantics with double universal references, yay!
template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
appendValue(TContainer && me,
            TSource && source)
{
    me.push_back(std::forward<TSource>(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prependValue(TContainer SEQAN_FORWARD_ARG me,
             TSource const & source)
{
    me.push_front(source);
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
prependValue(TContainer && me,
             TSource && source)
{
    me.push_front(std::forward<TSource>(source));
}

template <typename TContainer,
          typename TSource,
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
erase(TContainer SEQAN_FORWARD_ARG me,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd)
{
    me.erase(_iterStl(me, pos), _iterStl(me, posEnd));
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
inline void
erase(std::forward_list<TChar, TAlloc> & me,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const pos,
      typename Size<std::forward_list<TChar, TAlloc> >::Type const posEnd)
{
    me.erase_after(std::next(me.before_begin(), pos), std::next(me.begin(), posEnd));
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

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
erase(TContainer SEQAN_FORWARD_ARG me,
      typename Size<typename RemoveReference<TContainer>::Type>::Type const pos)
{
    erase(SEQAN_FORWARD(TContainer, me), pos, pos + 1);
}

// ----------------------------------------------------------------------------
// Function eraseFront
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
eraseBack(TContainer SEQAN_FORWARD_ARG me)
{
    me.pop_back();
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
inline void
eraseBack(std::forward_list<TChar, TAlloc> & me)
{
    if (!empty(me))
    {
        auto && it = me.before_begin(), itN = std::next(me.begin()), itEnd = me.end();
        while (itN != itEnd)
        {
            ++it;
            ++itN;
        }
        me.erase_after(it);
    }
}

template <typename TChar, typename TAlloc>
inline void
eraseBack(std::forward_list<TChar, TAlloc> && me)
{
    eraseBack(me);
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
#endif

// ----------------------------------------------------------------------------
// Function replace
// ----------------------------------------------------------------------------

//TODO(h-2): could be done a little more efficient?
template <typename TContainer,
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TSource>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
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
          typename TExpand>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, void)
replace(TContainer SEQAN_FORWARD_ARG me,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const pos,
        typename Size<typename RemoveReference<TContainer>::Type>::Type const posEnd,
        TSource const & source,
        Tag<TExpand> const &)
{
    replace(SEQAN_FORWARD(TContainer, me), pos, posEnd, source);
}

// ----------------------------------------------------------------------------
// Function reverse
// ----------------------------------------------------------------------------

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, HasSubscriptOperator<TContainer> >)
reverse(TContainer SEQAN_FORWARD_ARG me)
{
    std::reverse(me.begin(), me.end());
}

template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, Not<HasSubscriptOperator<TContainer> > >)
reverse(TContainer SEQAN_FORWARD_ARG me)
{
    me.reverse();
}

// ----------------------------------------------------------------------------
// Function toCString
// ----------------------------------------------------------------------------

template <typename TChar, typename TTraits, typename TAlloc>
inline const TChar*
toCString(std::basic_string<TChar, TTraits, TAlloc> const & me)
{
    return me.c_str();
}

template <typename TStream, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >,
                                Not<FixedSize_<typename RemoveReference<TContainer>::Type> > >, TStream &)
operator<<(TStream & stream, TContainer const & string)
{
    for (unsigned i = 0; i < length(string); ++i)
        stream << string[i] << ' ';

    return stream;
}

// ----------------------------------------------------------------------------
// Function assign (from CString)
// ----------------------------------------------------------------------------

// NOTE(h-2): is this optimal if container is not Contiguous?
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, void)
assign(TContainer & me,
       char const * source)
{
    resize(me, length(source));
    std::copy(source, source + std::strlen(source), me.begin());
}

template <typename TTraits, typename TAlloc>
inline void
assign(std::basic_string<char, TTraits, TAlloc> & me,
       char const * source)
{
    me = source;
}

#ifdef SEQAN_CXX11_STANDARD
template <typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, void)
assign(TContainer && me,
       char const * source)
{
    assign(me, source);
}

template <typename TTraits, typename TAlloc>
inline void
assign(std::basic_string<char, TTraits, TAlloc> && me,
       char const * source)
{
    me = source;
}
#endif

}

#undef SUPERMACRO__
#undef ITMACRO__
#undef ITRMACRO__
#undef COMMA

#endif
