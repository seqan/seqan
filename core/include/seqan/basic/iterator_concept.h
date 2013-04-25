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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef CORE_INCLUDE_SEQAN_BASIC_ITERATOR_CONCEPT_H_
#define CORE_INCLUDE_SEQAN_BASIC_ITERATOR_CONCEPT_H_

namespace seqan {

/**
.Metafunction.Pointer
..cat:Basic
..summary:Returns pointer to an object, required for @Function.operator->@, for example.
..signature:Pointer<T>::Type
..param.T:The type to query.
..returns:Pointer type.
..include:seqan/basic.h
 */

// Forward Declaration / Prototype.
template <typename T> struct Pointer;

/**
.Concept.IteratorAssociatedTypesConcept
..cat:Iterators
..summary:Requires metafunctions for the associated types used in the iterator concepts.
..signature:IteratorAssociatedTypesConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h

.Metafunction.Value.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.GetValue.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.Difference.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.Reference.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.Pointer.concept:Concept.IteratorAssociatedTypesConcept
 */

SEQAN_CONCEPT(IteratorAssociatedTypesConcept, (T))
{
    typedef typename Value<T>::Type      TValue;
    typedef typename GetValue<T>::Type   TGetValue;
    typedef typename Difference<T>::Type TDifference;
    typedef typename Reference<T>::Type  TReference;
    typedef typename Pointer<T>::Type    TPointer;
    
    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TDifference>));

    SEQAN_CONCEPT_USAGE(IteratorAssociatedTypesConcept)
    {
    }
};

/**
.Concept.InputIteratorConcept
..cat:Iterators
..summary:Iterator that allows dereferenced reading.
..baseconcept:Concept.IteratorAssociatedTypesConcept
..baseconcept:Concept.CopyConstructibleConcept
..baseconcept:Concept.EqualityComparableConcept
..signature:InputIteratorConcept<T>
..see:Concept.BasicOutputIteratorConcept
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h

.Function.operator->.concept:Concept.InputIteratorConcept
.Function.operator++ (prefix).concept:Concept.InputIteratorConcept
.Function.operator++ (suffix).concept:Concept.InputIteratorConcept
.Function.goNext.concept:Concept.InputIteratorConcept
.Function.operator*.concept:Concept.InputIteratorConcept
.Function.operator!=.concept.Concept.InputIteratorConcept
 */

SEQAN_CONCEPT_REFINE(InputIteratorConcept, (T), (IteratorAssociatedTypesConcept)(CopyConstructible)(EqualityComparable))
{
    typedef typename Value<T>::Type      TValue;
    typedef typename GetValue<T>::Type   TGetValue;
    typedef typename Difference<T>::Type TDifference;
    typedef typename Reference<T>::Type  TReference;
    typedef typename Pointer<T>::Type    TPointer;

    TValue v;
    T      x, y;

    SEQAN_CONCEPT_USAGE(InputIteratorConcept)
    {
        TReference & rv = v;
        T & rx =          x;

        SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TDifference>));
        // TODO(holtgrew): requires Convertible<reference, value_type>;
        // TODO(holtgrew): requires Convertible<pointer, cont value_type*>;

        // TODO(holtgrew): requires Dereferenceable<postincrement_result>;
        // TODO(holtgrew): requires requires Dereferenceable<postincrement_result>;

        // operator->: Cannot check this, need to know member for this.
        
        sameType(++x, rx);
        sameType(x++, y);
        goNext(x);

        sameType(rv, *x);

        x != x;
    }
};

/**
.Concept.BasicOutputIteratorConcept
..cat:Iterators
..summary:Iterator that allows dereferenced writing.
..baseconcept:Concept.IteratorAssociatedTypesConcept
..baseconcept:Concept.CopyConstructibleConcept
..signature:OutputIteratorConcept<T>
..see:Concept.InputIteratorConcept
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$, $t$ is a valid rvalue of type $X$.
..example.text:The following expressions must be valid.
..example.code:
*x = t     // Dereference assignment.
++x        // Preincrement.
(void)x++  // Postincrement.
*x++ = t   // Postincrement and assign.

assignValue(x, t);

.Function.assignValue.concept:Concept.BasicOutputIteratorConcept
.Function.operator++ (prefix).concept:Concept.BasicOutputIteratorConcept
.Function.operator++ (suffix).concept:Concept.BasicOutputIteratorConcept
.Function.goNext.concept:Concept.BasicOutputIteratorConcept
 */

SEQAN_CONCEPT_REFINE(BasicOutputIteratorConcept, (T), (CopyConstructible))
{
    typedef typename Value<T>::Type TValue;

    SEQAN_CONCEPT_ASSERT((Is<Assignable<TValue> >));
    
    T      x;
    TValue v;

    SEQAN_CONCEPT_USAGE(BasicOutputIteratorConcept)
    {
        *x = v;
        assignValue(x, v);
        value(x) = v;

        ++x;
        ignoreUnusedVariableWarning(x++);
        goNext(x);
        *x++ = v;

        ignoreUnusedVariableWarning(x);
    }
};

/**
.Concept.ForwardIteratorConcept
..cat:Iterators
..summary:Iterator that allows passing over a linear sequence multiple times.
..baseconcept:Concept.InputIteratorConcept
..baseconcept:Concept.DefaultConstructibleConcept
..signature:ForwardIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$.
..example.text:The following expressions must be valid.
..example.code:
++x  // Preincrement.
x++  // Postincrement.
 */

SEQAN_CONCEPT_REFINE(ForwardIteratorConcept, (T), (InputIteratorConcept)(DefaultConstructible))
{
    typedef typename Value<T>::Type TValue;

    T x;

    SEQAN_CONCEPT_USAGE(ForwardIteratorConcept)
    {
        ++x;
        ignoreUnusedVariableWarning(*x++);
    }
};

/**
.Concept.MutableForwardIteratorConcept
..cat:Iterators
..summary:A @Concept.ForwardIteratorConcept|Forward Iterator@ that allows dereferenced assignment.
..baseconcept:Concept.ForwardIteratorConcept
..baseconcept:Concept.BasicOutputIteratorConcept
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableForwardIteratorConcept, (T), (ForwardIteratorConcept)(BasicOutputIteratorConcept))
{
    typedef typename Value<T>::Type      TValue;

    T      x;
    TValue v;

    SEQAN_CONCEPT_USAGE(MutableForwardIteratorConcept)
    {
        *x = v;
    }
};

/**
.Concept.BidirectionalIteratorConcept
..cat:Iterators
..summary:Iterator that can be both incremented and decremented.
..baseconcept:Concept.ForwardIteratorConcept
..signature:BidirectionalIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$.
..example.text:The following expressions must be valid.
..example.code:
--x  // Predecrement.
x--  // Postdecrement.

.Function.operator-- (prefix).concept:Concept.BidirectionalIteratorConcept
.Function.operator-- (suffix).concept:Concept.BidirectionalIteratorConcept
.Function.goPrevious.concept:Concept.BidirectionalIteratorConcept
 */

SEQAN_CONCEPT_REFINE(BidirectionalIteratorConcept, (T), (ForwardIteratorConcept))
{
    T x;

    SEQAN_CONCEPT_USAGE(BidirectionalIteratorConcept)
    {
        --x;
        x--;
        goPrevious(x);
    }
};

/**
.Concept.MutableBidirectionalIteratorConcept
..cat:Iterators
..summary:A @Concept.BidirectionalIteratorConcept|Bidirectional Iterator@ that allows dereferenced assignment
..baseconcept:Concept.ForwardIteratorConcept
..signature:MutableBidirectionalIteratorConcept<T>
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableBidirectionalIteratorConcept, (T), (BidirectionalIteratorConcept)(MutableForwardIteratorConcept))
{
    typedef typename Value<T>::Type TValue;

    T      x;
    TValue v;

    SEQAN_CONCEPT_USAGE(MutableBidirectionalIteratorConcept)
    {
        *x = v;
    }
};

/**
.Concept.RandomAccessIteratorConcept
..cat:Iterators
..summary:An iterator allowing random access.
..baseconcept:Concept.BidirectionalIteratorConcept
..baseconcept:Concept.LessThanComparableConcept
..signature:RandomAccessIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$, $t$ is a valid rvalue of type $X$, $n$ is a distance type.
..example.text:The following expressions must be valid.
..example.code:
x += n    // Iterator addition assignment.
x + n     // Iterator addition.
n + i     // Iterator addition.
x -= n    // Iterator subtraction assignment.
x - n     // Iterator subtraction.
x - a     // Difference.
x[n]      // Element operator.

.Metafunction.Difference.concept:Concept.RandomAccessIteratorConcept
.Function.operator+=.concept:Concept.RandomAccessIteratorConcept
.Function.operator+.concept:Concept.RandomAccessIteratorConcept
.Function.operator-=.concept:Concept.RandomAccessIteratorConcept
.Function.operator-.concept:Concept.RandomAccessIteratorConcept
.Function.difference.concept:Concept.RandomAccessIteratorConcept
.Function.operator[].concept:Concept.RandomAccessIteratorConcept
.Function.goFurther.concept:Concept.RandomAccessIteratorConcept

.Function.operator>=.concept:Concept.RandomAccessIteratorConcept
.Function.operator>.concept:Concept.RandomAccessIteratorConcept
.Function.operator<=.concept:Concept.RandomAccessIteratorConcept
 */

SEQAN_CONCEPT_REFINE(RandomAccessIteratorConcept, (T), (BidirectionalIteratorConcept)(LessThanComparable))
{
    typedef typename Difference<T>::Type TDifference;
    typedef typename Value<T>::Type      TValue;

    T           x, y, a;
    TValue      t;
    TDifference n;

    SEQAN_CONCEPT_USAGE(RandomAccessIteratorConcept)
    {

        x += n;
        goFurther(x, n);        
        ignoreUnusedVariableWarning(x + n);
        ignoreUnusedVariableWarning(n + x);
        x -= n;
        y = x - n;
        n = x - a;
        n = difference(x, a);
        ignoreUnusedVariableWarning(x[n]);

        ignoreUnusedVariableWarning(x);
        ignoreUnusedVariableWarning(y);
    }
};

/**
.Concept.MutableRandomAccessIteratorConcept
..cat:Iterators
..summary:A @Concept.RandomAccessIteratorConcept@ that allows assignable derefentiation.
..baseconcept:Concept.BidirectionalIteratorConcept
..baseconcept:Concept.LessThanComparableConcept
..signature:RandomAccessIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..example.text:The following expressions should be valid.
..example.code:
value(x, n) = t
x[n] = t
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableRandomAccessIteratorConcept, (T), (RandomAccessIteratorConcept)(MutableBidirectionalIteratorConcept))
{
    typedef typename Difference<T>::Type TDifference;
    typedef typename Value<T>::Type      TValue;

    T           x;
    TValue      t;
    TDifference n;

    SEQAN_CONCEPT_USAGE(MutableRandomAccessIteratorConcept)
    {
        value(x, n) = t;  // TODO(holtgrew): Not supported?
        x[n] = t;
    }
};

/**
.Concept.RootedIteratorConcept
..cat:Iterators
..summary:Iterator that knows its container.
..baseconcept:Concept.ForwardIteratorConcept
..signature:RootedIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h

.Metafunction.Container.concept:Concept.RootedIteratorConcept
.Function.container.concept:Concept.RootedIteratorConcept
.Function.atBegin.concept:Concept.RootedIteratorConcept
.Function.atEnd.concept:Concept.RootedIteratorConcept
 */

SEQAN_CONCEPT_REFINE(RootedIteratorConcept, (T), (IteratorAssociatedTypesConcept))
{
    typedef typename Container<T>::Type TContainer;

    T x;

    SEQAN_CONCEPT_USAGE(RootedIteratorConcept)
    {
        T xs;
        ignoreUnusedVariableWarning(xs);

        TContainer & c = container(x);
        atBegin(x);
        atEnd(x);
        ignoreUnusedVariableWarning(c);
    }
};

/**
.Concept.MutableRootedIteratorConcept
..cat:Iterators
..summary:A @Concept.RootedIteratorConcept|Rooted Iterator@ that allows dereferenced assignment.
..baseconcept:Concept.ForwardIteratorConcept
..signature:RootedIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableRootedIteratorConcept, (T), (RootedIteratorConcept)(MutableForwardIteratorConcept))
{
    SEQAN_CONCEPT_USAGE(MutableRootedIteratorConcept)
    {
    }
};

/**
.Concept.RootedRandomAccessIteratorConcept
..cat:Iterators
..summary:An iterator that is both rooted and random access, allowing to implement @Function.position@.
..signature:RootedRandomAccessIteratorConcept<T>
..baseconcept:Concept.RandomAccessIteratorConcept
..baseconcept:Concept.RootedIteratorConcept

.Metafunction.Position.concept:Concept.RootedRandomAccessIteratorConcept
.Function.position.concept:Concept.RootedRandomAccessIteratorConcept
.Function.setPosition.concept:Concept.RootedRandomAccessIteratorConcept
.Function.goBegin.concept:Concept.RootedRandomAccessIteratorConcept
.Function.goEnd.concept:Concept.RootedRandomAccessIteratorConcept
 */

SEQAN_CONCEPT_REFINE(RootedRandomAccessIteratorConcept, (T), (RootedIteratorConcept)(RandomAccessIteratorConcept))
{
    typedef typename Position<T>::Type TPosition;

    SEQAN_CONCEPT_USAGE(RootedRandomAccessIteratorConcept)
    {
        T x;

        TPosition p = position(x);
        setPosition(x, p);
        goBegin(x);
        goEnd(x);
    }
};

/**
.Concept.MutableRootedRandomAccessIteratorConcept
..cat:Iterators
..baseconcept:Concept.RootedRandomAccessIteratorConcept
..baseconcept:Concept.MutableBidirectionalIteratorConcept
..summary:A @Concept.RootedIteratorConcept|Rooted Iterator@ that allows dereferenced assignment.
..signature:MutableRootedRandomAccessIteratorConcept<T>
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableRootedRandomAccessIteratorConcept, (T), (RootedRandomAccessIteratorConcept)(MutableBidirectionalIteratorConcept))
{
    SEQAN_CONCEPT_USAGE(MutableRootedRandomAccessIteratorConcept)
    {
    }
};

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_ITERATOR_CONCEPT_H_
