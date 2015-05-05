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

// template <typename TContainer>
// struct StlContainerConcept :
//     ContainerConcept<TContainer>
// {
//     SEQAN_CONCEPT_USAGE(StlContainerConcept)
//     {
//     }
// };

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
// Mfn IsStlContiguous_
// ----------------------------------------------------------------------------

template <typename T>
struct IsStlContiguous_ :
    public True
{};

template <typename TChar, typename TAlloc>
struct IsStlContiguous_<std::list<TChar, TAlloc> > :
    public False
{};

template <typename TChar, typename TAlloc>
struct IsStlContiguous_<std::list<TChar, TAlloc> const> :
    public False
{};

#ifdef SEQAN_CXX11_STANDARD
template <typename TChar, typename TAlloc>
struct IsStlContiguous_<std::forward_list<TChar, TAlloc> > :
    public False
{};

template <typename TChar, typename TAlloc>
struct IsStlContiguous_<std::forward_list<TChar, TAlloc> const> :
    public False
{};
#endif

// ----------------------------------------------------------------------------
// Mfn IsContiguous
// ----------------------------------------------------------------------------

// general case
template <typename TContainer>
struct IsContiguous<TContainer, typename EnableIf<Is<StlContainerConcept<TContainer> > >::Type> :
    public IsStlContiguous_<TContainer>
{};

//TODO change, only overload positively for string, vector, array

// ----------------------------------------------------------------------------
// Mfn Value
// ----------------------------------------------------------------------------

// template <typename TContainer>
// struct Value<TContainer, typename EnableIf<Is<StlContainerConcept<TContainer> > >::Type>
// {
//     typedef typename TContainer::value_type Type;
// };

}

#endif
