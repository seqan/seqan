// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Concept definitions for containers.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_CONTAINER_CONCEPT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_CONTAINER_CONCEPT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T> struct Infix;
template <typename T> struct Prefix;
template <typename T> struct Suffix;
    
// ============================================================================
// Concepts
// ============================================================================

// Forwards.
struct Standard_;
typedef Tag<Standard_> const Standard;
template <typename TContainer, typename TSpec> struct Iterator;

/**
.Concept.ContainerConcept
..baseconcept:Concept.AssignableConcept
..baseconcept:Concept.DefaultConstructibleConcept
..baseconcept:Concept.CopyConstructibleConcept
..signature:ContainerConcept
..summary:Concept for mutable containers.
..include:seqan/basic.h
*/

// mutable container concept
template <typename TContainer>
struct ContainerConcept :
    Assignable<TContainer>,
    DefaultConstructible<TContainer>,
    CopyConstructible<TContainer>
{
    typedef typename Value<TContainer>::Type                TValue;
    typedef typename GetValue<TContainer>::Type             TGetValue;
    typedef typename Reference<TContainer>::Type            TReference;
    typedef typename Size<TContainer>::Type                 TSize;
    typedef typename Position<TContainer>::Type             TPosition;
    typedef typename Difference<TContainer>::Type           TDifference;
    typedef typename Iterator<TContainer, Standard>::Type   TIterator;

    TContainer  c, c2;
    TValue      val;
    TSize       size;
    TPosition   pos;
    TDifference diff;
    TIterator   iter;
    
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<TValue>));
    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TDifference>));
    SEQAN_CONCEPT_ASSERT((UnsignedIntegerConcept<TSize>));

    SEQAN_CONCEPT_USAGE(ContainerConcept)
    {
        // test return of const values
        sameType(getValue(c, 0), val);

        // TODO(holtgrew): Index based access requires random-access and *linear* container.
        // test whether returned references/proxies
        // can be assigned to val and vice versa
        val = value(c, 0);
        value(c, 0) = val;
        moveValue(c, pos, val);

        // sameType(value(str, 0), val); would not work
        // for Strings returning proxies, e.g. String<.., Packed>

        // TODO(holtgrew): Too strong assumption about random access iterators IMO.
        // test iterators
        sameType(iter, begin(c, Standard()));
        sameType(iter, end(c, Standard()));
        sameType(diff, end(c, Standard()) - begin(c, Standard()));
        
        // length and empty
        sameType(size, length(c));
        sameType(true, empty(c));

        // clear
        clear(c);
        
        // TODO: infix/suffix/prefix 
        // maybe we need a SequenceConcept between Container and String
        
        // swap containers
//        swap(c, c2);          // swap is not yet supported by every string
    }
};

/**
.Concept.SequenceConcept
..baseconcept:Concept.ContainerConcept
..summary:Concept for sequences.
..include:seqan/basic.h
*/

SEQAN_CONCEPT_REFINE(SequenceConcept, (TString), (ContainerConcept))
{
    typedef typename Value<TString>::Type                 TValue;
    typedef typename Size<TString>::Type                  TSize;
    typedef typename Position<TString>::Type              TPos;
    typedef typename Difference<TString>::Type            TDifference;
    typedef typename Iterator<TString, Standard>::Type    TIterator;

    TValue      val;
    TSize       size;
    TPos        pos;

    TString     str, str2;
    
    SEQAN_CONCEPT_USAGE(SequenceConcept)
    {
        pos = 0u;

        // append
        append(str, str2);
        appendValue(str, val);

        // capacity
        sameType(size, capacity(str));
    }
};

//void testStringConcepts()
//{
//    SEQAN_CONCEPT_ASSERT((StringConcept<String<char, Alloc<> > >));
//    SEQAN_CONCEPT_ASSERT((StringConcept<String<Pair<int, double>, Alloc<> > >));
////    SEQAN_CONCEPT_ASSERT((StringConcept<String<bool, Packed<> > >));  // doesn't compile yet
////    SEQAN_CONCEPT_ASSERT((StringConcept<String<Dna5, Packed<> > >));
//    SEQAN_CONCEPT_ASSERT((StringConcept<String<int, Array<50> > >));
//}
    
}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_CONTAINER_CONCEPT_H_
