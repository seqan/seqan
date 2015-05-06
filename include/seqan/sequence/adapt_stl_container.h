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
// ITMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, , typename TChar COMMA typename TTraits COMMA typename TAlloc)

ITMACRO__(std::vector<TChar COMMA TAlloc>,        const, typename TChar COMMA typename TAlloc)
ITMACRO__(std::deque<TChar COMMA TAlloc>,         const, typename TChar COMMA typename TAlloc)
ITMACRO__(std::list<TChar COMMA TAlloc>,          const, typename TChar COMMA typename TAlloc)
#ifdef SEQAN_CXX11_STANDARD
ITMACRO__(std::forward_list<TChar COMMA TAlloc>,  const, typename TChar COMMA typename TAlloc)
ITMACRO__(std::array<TChar COMMA N>,              const, typename TChar COMMA std::size_t N)
#endif
// ITMACRO__(std::basic_string<TChar COMMA TTraits COMMA TAlloc>, const, typename TChar COMMA typename TTraits COMMA typename TAlloc)

// TODO(h-2): remove this crap overload and resolve to stl's RandomAcccesIterator like for std::vector
template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator<std::basic_string<TChar, TCharTraits, TAlloc>, Standard >
{
    typedef TChar * Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator<std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>
{
    typedef TChar const * Type;
};

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


}

#undef SUPERMACRO__
#undef ITMACRO__
#undef ITRMACRO__
#undef COMMA

#endif
