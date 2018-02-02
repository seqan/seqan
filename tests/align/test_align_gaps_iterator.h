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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GAPS_ITERATOR_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_GAPS_ITERATOR_H_

#include <sstream>

#include <seqan/align.h>

// ==========================================================================
// Generic Tests For Gaps Specializations
// ==========================================================================

template <typename TGapsSpec>
void testAlignGapsIteratorMetafunctions(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                            TString;
    typedef Gaps<TString, TGapsSpec>              TGaps;
    typedef Iter<TGaps, GapsIterator<TGapsSpec> > TIter;

    // We only check that we can call these metafunctions.  We should also test
    // for expected resulting types.  This could/should be done with concept
    // checking.

    // Sequence metafunctions.
    typedef typename Difference<TIter>::Type         TDifference SEQAN_UNUSED_TYPEDEF;
    typedef typename GetValue<TIter>::Type           TGetValue   SEQAN_UNUSED_TYPEDEF;
    typedef typename Position<TIter>::Type           TPosition   SEQAN_UNUSED_TYPEDEF;
    typedef typename Reference<TIter>::Type          TReference  SEQAN_UNUSED_TYPEDEF;
    typedef typename Size<TIter>::Type               TSize       SEQAN_UNUSED_TYPEDEF;
    typedef typename Value<TIter>::Type              TValue      SEQAN_UNUSED_TYPEDEF;
}

// Check that dereferencing works.
template <typename TGapsSpec>
void testAlignGapsIteratorTrivialIteratorFunctions(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    TIter it = begin(gaps);

    SEQAN_ASSERT_EQ(*it, 'C');
}

// Test the functions for rooted container (container(), atBegin(),
// atEnd()) and rooted andom acess (position(), setPosition()).

template <typename TGapsSpec>
void testAlignGapsIteratorRootedRandomAccessIteratorFunctions(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    insertGaps(gaps, 2, 2);

    // 012345
    // CG--AT

    // Test container().
    {
        TIter it = begin(gaps);

        SEQAN_ASSERT_EQ(&container(it), &gaps);
    }

    // Test position().
    {
        TIter it = begin(gaps);
        TIter itEnd = end(gaps);

        SEQAN_ASSERT_EQ(position(it), 0);
        SEQAN_ASSERT_EQ(position(itEnd), 6);
    }

    // Test atBegin().
    {
        TIter it = begin(gaps);
        TIter itEnd = end(gaps);

        SEQAN_ASSERT(atBegin(it));
        SEQAN_ASSERT_NOT(atBegin(itEnd));
    }

    // Test atEnd().
    {
        TIter it = begin(gaps);
        TIter itEnd = end(gaps);

        SEQAN_ASSERT_NOT(atEnd(it));
        SEQAN_ASSERT(atEnd(itEnd));
    }

    // Test setPosition().
    {
        TIter it = begin(gaps);
        TIter itEnd = end(gaps);

        SEQAN_ASSERT_EQ(*it, 'C');
        --itEnd;
        SEQAN_ASSERT_EQ(*itEnd, 'T');
    }
}

// Test movement functions of the iterators.

template <typename TGapsSpec>
void testAlignGapsIteratorMovement(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    insertGaps(gaps, 2, 2);

    // 012345
    // CG--AT

    // Test prefix increment.
    {
        TIter it = iter(gaps, 1, Standard());
        SEQAN_ASSERT_EQ(*it, 'G');
        ++it;
        SEQAN_ASSERT(isGap(it));
        ++it;
        SEQAN_ASSERT(isGap(it));
        ++it;
        SEQAN_ASSERT_EQ(*it, 'A');
    }

    // Test postfix increment.
    {
        TIter it = iter(gaps, 1, Standard());
        SEQAN_ASSERT_EQ(*it, 'G');
        it++;
        SEQAN_ASSERT(isGap(it));
        it++;
        SEQAN_ASSERT(isGap(it));
        it++;
        SEQAN_ASSERT_EQ(*it, 'A');
    }

    // Test goNext().
    {
        TIter it = iter(gaps, 1, Standard());
        SEQAN_ASSERT_EQ(*it, 'G');
        goNext(it);
        SEQAN_ASSERT(isGap(it));
        goNext(it);
        SEQAN_ASSERT(isGap(it));
        goNext(it);
        SEQAN_ASSERT_EQ(*it, 'A');
    }

    // Test prefix decrement.
    {
        TIter it = iter(gaps, 4, Standard());
        SEQAN_ASSERT_EQ(*it, 'A');
        --it;
        SEQAN_ASSERT(isGap(it));
        --it;
        SEQAN_ASSERT(isGap(it));
        --it;
        SEQAN_ASSERT_EQ(*it, 'G');
    }

    // Test postfix decrement.
    {
        TIter it = iter(gaps, 4, Standard());
        SEQAN_ASSERT_EQ(*it, 'A');
        it--;
        SEQAN_ASSERT(isGap(it));
        it--;
        SEQAN_ASSERT(isGap(it));
        it--;
        SEQAN_ASSERT_EQ(*it, 'G');
    }

    // Test goPrevious().
    {
        TIter it = iter(gaps, 4, Standard());
        SEQAN_ASSERT_EQ(*it, 'A');
        goPrevious(it);
        SEQAN_ASSERT(isGap(it));
        goPrevious(it);
        SEQAN_ASSERT(isGap(it));
        goPrevious(it);
        SEQAN_ASSERT_EQ(*it, 'G');
    }

    // Test operator+=.
    {
        TIter it = iter(gaps, 1, Standard());
        SEQAN_ASSERT_EQ(*it, 'G');
        it += 3;
        SEQAN_ASSERT_EQ(*it, 'A');
    }

    // Test goFurther() forward.
    {
        TIter it = iter(gaps, 1, Standard());
        SEQAN_ASSERT_EQ(*it, 'G');
        goFurther(it, 3);
        SEQAN_ASSERT_EQ(*it, 'A');
    }

    // Test operator-=.
    {
        TIter it = iter(gaps, 4, Standard());
        SEQAN_ASSERT_EQ(*it, 'A');
        it -= 3;
        SEQAN_ASSERT_EQ(*it, 'G');
    }

    // Test goFurther() backwards.
    {
        TIter it = iter(gaps, 4, Standard());
        SEQAN_ASSERT_EQ(*it, 'A');
        goFurther(it, -3);
        SEQAN_ASSERT_EQ(*it, 'G');
    }

    // Test moving to end.
    {
        TIter it = iter(gaps, 4, Standard());
        goEnd(it);
        SEQAN_ASSERT(atEnd(it));
    }

    // Test moving to begin.
    {
        TIter it = iter(gaps, 4, Standard());
        goBegin(it);
        SEQAN_ASSERT(atBegin(it));
    }
}

template <typename TGapsSpec>
void testAlignGapsIteratorRelations(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    insertGaps(gaps, 2, 2);

    TIter it1 = iter(gaps, 1, Standard());
    TIter it2 = iter(gaps, 2, Standard());

    SEQAN_ASSERT(it1 < it2);
    SEQAN_ASSERT_NOT(it2 < it1);

    SEQAN_ASSERT(it2 > it1);
    SEQAN_ASSERT_NOT(it1 > it2);

    SEQAN_ASSERT(it1 == it1);
    SEQAN_ASSERT_NOT(it1 == it2);

    SEQAN_ASSERT_NOT(it1 != it1);
    SEQAN_ASSERT(it1 != it2);

    SEQAN_ASSERT(it1 <= it1);
    SEQAN_ASSERT(it1 <= it2);
    SEQAN_ASSERT_NOT(it2 <= it1);

    SEQAN_ASSERT(it1 >= it1);
    SEQAN_ASSERT(it2 >= it1);
    SEQAN_ASSERT_NOT(it1 >= it2);
}

template <typename TGapsSpec>
void testAlignGapsIteratorPointerArithmetic(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    insertGaps(gaps, 2, 2);

    TIter it1 = iter(gaps, 1, Standard());
    TIter it2 = iter(gaps, 4, Standard());

    SEQAN_ASSERT_EQ(it2 - it1, 3);
    SEQAN_ASSERT_EQ(it1 - it2, -3);

    SEQAN_ASSERT_EQ(difference(it2, it1), 3);
    SEQAN_ASSERT_EQ(difference(it1, it2), -3);

    SEQAN_ASSERT(it1 + 3 == it2);
    SEQAN_ASSERT(it1 == it2 - 3);
}

// Iterate once over the gaps and collect result.
template <typename TGapsSpec>
void testAlignGapsIteratorForwardIteration(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    insertGaps(gaps, 4, 2);
    insertGaps(gaps, 2, 2);
    insertGaps(gaps, 0, 2);

    // 0123456789
    // --CG--AT--

    std::stringstream ss;
    for (TIter it = begin(gaps); !atEnd(it); goNext(it))
    {
        if (isGap(it))
            ss << '-';
        else
            ss << *it;
    }

    SEQAN_ASSERT_EQ(ss.str(), "--CG--AT--");
}

template <typename TGapsSpec>
void testAlignGapsIteratorReverseIteration(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    insertGaps(gaps, 4, 2);
    insertGaps(gaps, 2, 2);
    insertGaps(gaps, 0, 2);

    // 0123456789
    // --CG--AT--

    std::stringstream ss;
    TIter it = end(gaps);
    do
    {
        goPrevious(it);
        if (isGap(it))
            ss << '-';
        else
            ss << *it;
    } while (!atBegin(it));

    SEQAN_ASSERT_EQ(ss.str(), "--TA--GC--");
}

template <typename TGapsSpec>
void testAlignGapsIteratorCountGapsCountCharactersIsGap(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    TString seq = "CGAT";
    TGaps gaps(seq);
    insertGaps(gaps, 4, 2);
    insertGaps(gaps, 2, 2);
    insertGaps(gaps, 0, 2);

    // 0123456789
    // --CG--AT--

    TIter it = begin(gaps);

    SEQAN_ASSERT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 2u);
    SEQAN_ASSERT_EQ(countGaps(it, RightOfViewPos()), 2u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, RightOfViewPos()), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
    ++it;

    SEQAN_ASSERT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 1u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 1u);
    SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
    ++it;

    SEQAN_ASSERT_NOT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 0u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 2u);
    SEQAN_ASSERT_EQ(countCharacters(it), 2u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
    ++it;

    SEQAN_ASSERT_NOT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 0u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it), 1u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 1u);
    ++it;

    SEQAN_ASSERT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 2u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 2u);
    ++it;

    SEQAN_ASSERT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 1u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 1u);
    SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
    ++it;

    SEQAN_ASSERT_NOT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 0u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 2u);
    SEQAN_ASSERT_EQ(countCharacters(it), 2u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
    ++it;

    SEQAN_ASSERT_NOT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 0u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it), 1u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 1u);
    ++it;

    SEQAN_ASSERT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 2u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 2u);
    ++it;

    SEQAN_ASSERT(isGap(it));
    SEQAN_ASSERT_EQ(countGaps(it), 1u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 1u);
    SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
    ++it;

    SEQAN_ASSERT(atEnd(it));
    SEQAN_ASSERT_EQ(countGaps(it), 0u);
    SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 2u);
    SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
}

template <typename TGapsSpec>
void testAlignGapsIteratorClippedCountGapsCountCharactersIsGap(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    // Clip into leading/trailing gaps.
    {
        TString seq = "CGAT";
        TGaps gaps(seq);
        insertGaps(gaps, 4, 2);
        insertGaps(gaps, 2, 2);
        insertGaps(gaps, 0, 2);

        // 0123456789
        // --CG--AT--
        //  XXXXXXXX

        setClippedEndPosition(gaps, 9);
        setClippedBeginPosition(gaps, 1);

        std::stringstream ss;
        ss << gaps;
        SEQAN_ASSERT_EQ(ss.str(), "-CG--AT-");


        TIter it = begin(gaps);

        SEQAN_ASSERT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 1u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
        ++it;

        SEQAN_ASSERT_NOT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 1u);
        SEQAN_ASSERT_EQ(countCharacters(it), 2u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
        ++it;

        SEQAN_ASSERT_NOT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 1u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 1u);
        ++it;

        SEQAN_ASSERT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 2u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 2u);
        ++it;

        SEQAN_ASSERT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 1u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 1u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
        ++it;

        SEQAN_ASSERT_NOT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 2u);
        SEQAN_ASSERT_EQ(countCharacters(it), 2u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
        ++it;

        SEQAN_ASSERT_NOT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 1u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 1u);
        ++it;

        SEQAN_ASSERT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 1u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 2u);
        ++it;

        SEQAN_ASSERT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countGaps(it, LeftOfViewPos()), 1u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it, LeftOfViewPos()), 0u);
    }

    // Clip into leading/trailing characters.
    {
        TString seq = "CGAT";
        TGaps gaps(seq);
        insertGaps(gaps, 4, 2);
        insertGaps(gaps, 2, 2);
        insertGaps(gaps, 0, 2);

        // 0123456789
        // --CG--AT--
        //    XXXX

        setClippedEndPosition(gaps, 7);
        setClippedBeginPosition(gaps, 3);

        std::stringstream ss;
        ss << gaps;
        SEQAN_ASSERT_EQ(ss.str(), "G--A");

        TIter it = begin(gaps);

        SEQAN_ASSERT_NOT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 1u);
        ++it;

        SEQAN_ASSERT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 2u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
        ++it;

        SEQAN_ASSERT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 1u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
        ++it;

        SEQAN_ASSERT_NOT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 1u);
        ++it;

        SEQAN_ASSERT_NOT(isGap(it));
        SEQAN_ASSERT_EQ(countGaps(it), 0u);
        SEQAN_ASSERT_EQ(countCharacters(it), 0u);
    }
}

template <typename TGapsSpec>
void testAlignGapsIteratorGapOperationsCenter(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    // Note: The string's characters are chosen such that it is clear from which
    // part of the char is from by its value, i.e. A = left, C = left adjacent
    // to gaps etc.

    // Insert one gap in the center and do queries on the result.
    {
        TString seq("AAACGTTT");
        TGaps gaps(seq);

        TIter it = iter(gaps, 4);
        insertGap(it);

        SEQAN_ASSERT_EQ(length(seq), 8u);
        SEQAN_ASSERT_EQ(length(gaps), 9u);

        // Query sequence.
        SEQAN_ASSERT(gaps[3] == 'C');
        // SEQAN_ASSERT(gaps[4] == '-');
        SEQAN_ASSERT(gaps[5] == 'G');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 3));
        SEQAN_ASSERT(isGap(gaps, 4));
        SEQAN_ASSERT_NOT(isGap(gaps, 5));
    }

    // Insert two gaps in the center and do queries on the result.
    {
        TString seq("AAACGTTT");
        TGaps gaps(seq);

        TIter it = iter(gaps, 4);
        insertGaps(it, 2);

        SEQAN_ASSERT_EQ(length(seq), 8u);
        SEQAN_ASSERT_EQ(length(gaps), 10u);

        // Query sequence.
        SEQAN_ASSERT(gaps[3] == 'C');
        // SEQAN_ASSERT(gaps[4] == '-');
        // SEQAN_ASSERT(gaps[5] == '-');
        SEQAN_ASSERT(gaps[6] == 'G');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 3));
        SEQAN_ASSERT(isGap(gaps, 4));
        SEQAN_ASSERT(isGap(gaps, 5));
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
    }

    // Insert three gaps in the center, remove one, then perform queries.
    {
        TString seq("AAACGTTT");
        TGaps gaps(seq);

        TIter it = iter(gaps, 4);
        insertGaps(it, 3);
        ++it;
        SEQAN_ASSERT_EQ(removeGap(it), 1u);

        SEQAN_ASSERT_EQ(length(seq), 8u);
        SEQAN_ASSERT_EQ(length(gaps), 10u);

        // Query sequence.
        SEQAN_ASSERT(gaps[3] == 'C');
        // SEQAN_ASSERT(gaps[4] == '-');
        // SEQAN_ASSERT(gaps[5] == '-');
        SEQAN_ASSERT(gaps[6] == 'G');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 3));
        SEQAN_ASSERT(isGap(gaps, 4));
        SEQAN_ASSERT(isGap(gaps, 5));
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
    }

    // Insert four gaps in the center, remove one, then perform queries.
    {
        TString seq("AAACGTTT");
        TGaps gaps(seq);

        TIter it = iter(gaps, 4);
        insertGaps(it, 4);
        ++it;
        SEQAN_ASSERT_EQ(removeGaps(it, 2), 2u);

        SEQAN_ASSERT_EQ(length(seq), 8u);
        SEQAN_ASSERT_EQ(length(gaps), 10u);

        // Query sequence.
        SEQAN_ASSERT(gaps[3] == 'C');
        // SEQAN_ASSERT(gaps[4] == '-');
        // SEQAN_ASSERT(gaps[5] == '-');
        SEQAN_ASSERT(gaps[6] == 'G');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 3));
        SEQAN_ASSERT(isGap(gaps, 4));
        SEQAN_ASSERT(isGap(gaps, 5));
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
    }

    // Insert gaps in the center, then clear.
    {
        TString seq("AAACGTTT");
        TGaps gaps(seq);

        TIter it = iter(gaps, 4);
        insertGaps(it, 4);
        clearGaps(gaps);

        SEQAN_ASSERT_EQ(length(seq), 8u);
        SEQAN_ASSERT_EQ(length(gaps), 8u);

        // Query sequence.
        SEQAN_ASSERT(gaps[2] == 'A');
        SEQAN_ASSERT(gaps[3] == 'C');
        SEQAN_ASSERT(gaps[4] == 'G');
        SEQAN_ASSERT(gaps[5] == 'T');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 2));
        SEQAN_ASSERT_NOT(isGap(gaps, 3));
        SEQAN_ASSERT_NOT(isGap(gaps, 4));
        SEQAN_ASSERT_NOT(isGap(gaps, 5));
    }
}

template <typename TGapsSpec>
void testAlignGapsIteratorGapOperationsLeading(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    // Note: The string's characters are chosen such that it is clear from which
    // part of the char is from by its value, i.e. G = right adjacent to gap,
    // T = right of gap.

    // Insert one gap at the beginning and do queries on the result.
    {
        TString seq("GTTTTTT");
        TGaps gaps(seq);

        TIter it = begin(gaps);
        insertGap(it);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 8u);

        // Query sequence.
        // SEQAN_ASSERT(gaps[0] == '-');
        SEQAN_ASSERT(gaps[1] == 'G');

        // Query gaps.
        SEQAN_ASSERT(isGap(gaps, 0));
        SEQAN_ASSERT_NOT(isGap(gaps, 1));
    }

    // Insert two gaps at the beginning and do queries on the result.
    {
        TString seq("GTTTTTT");
        TGaps gaps(seq);

        TIter it = begin(gaps);
        insertGaps(it, 2);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 9u);

        // Query sequence.
        // SEQAN_ASSERT(gaps[0] == '-');
        // SEQAN_ASSERT(gaps[1] == '-');
        SEQAN_ASSERT(gaps[2] == 'G');

        // Query gaps.
        SEQAN_ASSERT(isGap(gaps, 0));
        SEQAN_ASSERT(isGap(gaps, 1));
        SEQAN_ASSERT_NOT(isGap(gaps, 2));
    }

    // Insert three gaps at the beginning, remove one, then perform queries.
    {
        TString seq("GTTTTTT");
        TGaps gaps(seq);

        TIter it = begin(gaps);
        insertGaps(it, 3);
        ++it;
        removeGap(it);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 9u);

        // Query sequence.
        // SEQAN_ASSERT(gaps[0] == '-');
        // SEQAN_ASSERT(gaps[1] == '-');
        SEQAN_ASSERT(gaps[2] == 'G');

        // Query gaps.
        SEQAN_ASSERT(isGap(gaps, 0));
        SEQAN_ASSERT(isGap(gaps, 1));
        SEQAN_ASSERT_NOT(isGap(gaps, 2));
    }

    // Insert four gaps at the beginning, remove one, then perform queries.
    {
        TString seq("GTTTTTT");
        TGaps gaps(seq);

        TIter it = begin(gaps);
        insertGaps(it, 4);
        it += 2;
        SEQAN_ASSERT_EQ(removeGaps(it, 2), 2u);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 9u);

        // Query sequence.
        // SEQAN_ASSERT(gaps[0] == '-');
        // SEQAN_ASSERT(gaps[1] == '-');
        SEQAN_ASSERT(gaps[2] == 'G');

        // Query gaps.
        SEQAN_ASSERT(isGap(gaps, 0));
        SEQAN_ASSERT(isGap(gaps, 1));
        SEQAN_ASSERT_NOT(isGap(gaps, 2));
    }

    // Insert gaps at the beginning, then clear.
    {
        TString seq("GTTTTTT");
        TGaps gaps(seq);

        TIter it = begin(gaps);
        insertGaps(it, 4);
        clearGaps(gaps);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 7u);

        // Query sequence.
        SEQAN_ASSERT(gaps[0] == 'G');
        SEQAN_ASSERT(gaps[1] == 'T');
        SEQAN_ASSERT(gaps[2] == 'T');
        SEQAN_ASSERT(gaps[3] == 'T');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 3));
        SEQAN_ASSERT_NOT(isGap(gaps, 4));
        SEQAN_ASSERT_NOT(isGap(gaps, 5));
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
    }
}

template <typename TGapsSpec>
void testAlignGapsIteratorGapOperationsTrailing(TGapsSpec const & /*spec*/)
{
    using namespace seqan;

    typedef Dna5String                               TString;
    typedef Gaps<TString, TGapsSpec>                 TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TIter;

    // Note: The string's characters are chosen such that it is clear from which
    // part of the char is from by its value, i.e. A = left of gap, C = left
    // adjacent with gap.

    // Insert one gap at the end and do queries on the result.
    {
        TString seq("AAAAAAC");
        TGaps gaps(seq);

        TIter it = end(gaps);
        insertGap(it);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 8u);

        // Query sequence.
        SEQAN_ASSERT(gaps[6] == 'C');
        // SEQAN_ASSERT(gaps[7] == '-');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
        SEQAN_ASSERT(isGap(gaps, 7));
    }

    // Insert two gaps at the end and do queries on the result.
    {
        TString seq("AAAAAAC");
        TGaps gaps(seq);

        TIter it = end(gaps);
        insertGaps(it, 2);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 9u);

        // Query sequence.
        SEQAN_ASSERT(gaps[6] == 'C');
        // SEQAN_ASSERT(gaps[7] == '-');
        // SEQAN_ASSERT(gaps[8] == '-');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
        SEQAN_ASSERT(isGap(gaps, 7));
        SEQAN_ASSERT(isGap(gaps, 8));
    }

    // Insert three gaps at the end, remove one, then perform queries.
    {
        TString seq("AAAAAAC");
        TGaps gaps(seq);

        TIter it = end(gaps);
        insertGaps(it, 3);
        SEQAN_ASSERT_EQ(removeGap(it), 1u);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 9u);

        // Query sequence.
        SEQAN_ASSERT(gaps[6] == 'C');
        // SEQAN_ASSERT(gaps[7] == '-');
        // SEQAN_ASSERT(gaps[8] == '-');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
        SEQAN_ASSERT(isGap(gaps, 7));
        SEQAN_ASSERT(isGap(gaps, 8));
    }

    // Insert four gaps at the end, remove two, then perform queries.
    {
        TString seq("AAAAAAC");
        TGaps gaps(seq);

        TIter it = end(gaps);
        insertGaps(it, 4);
        SEQAN_ASSERT_EQ(removeGaps(it, 2), 2u);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 9u);

        // Query sequence.
        SEQAN_ASSERT(gaps[6] == 'C');
        // SEQAN_ASSERT(gaps[7] == '-');
        // SEQAN_ASSERT(gaps[8] == '-');

        // Query gaps.
        SEQAN_ASSERT_NOT(isGap(gaps, 6));
        SEQAN_ASSERT(isGap(gaps, 7));
        SEQAN_ASSERT(isGap(gaps, 8));
    }

    // Insert gaps at the end, then clear.
    {
        TString seq("AAAAAAC");
        TGaps gaps(seq);

        TIter it = end(gaps);
        insertGaps(it, 4);
        clearGaps(gaps);

        SEQAN_ASSERT_EQ(length(seq), 7u);
        SEQAN_ASSERT_EQ(length(gaps), 7u);
    }
}

// ==========================================================================
// Tests for Array Gaps Iter
// ==========================================================================

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_metafunctions)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorMetafunctions(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_trivial_iterator_array_functions)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorTrivialIteratorFunctions(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_rooted_random_access_iterator_array_functions)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorRootedRandomAccessIteratorFunctions(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_movement)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorMovement(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_relations)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorRelations(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_pointer_arithmetic)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorPointerArithmetic(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_forward_iteration)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorForwardIteration(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_reverse_iteration)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorReverseIteration(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_count_gaps_count_characters_is_gap)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorCountGapsCountCharactersIsGap(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_clipped_count_gaps_count_characters_is_gap)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorClippedCountGapsCountCharactersIsGap(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_gap_operations_center)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorGapOperationsCenter(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_gap_operations_leading)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorGapOperationsLeading(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_array_gap_operations_trailing)
{
    using namespace seqan;
    typedef ArrayGaps TTag;
    testAlignGapsIteratorGapOperationsTrailing(TTag());
}

// ==========================================================================
// Tests for Anchor Gaps Iter
// ==========================================================================

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_metafunctions)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorMetafunctions(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_trivial_iterator_anchor_functions)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorTrivialIteratorFunctions(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_rooted_random_access_iterator_anchor_functions)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorRootedRandomAccessIteratorFunctions(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_movement)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorMovement(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_relations)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorRelations(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_pointer_arithmetic)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorPointerArithmetic(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_forward_iteration)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorForwardIteration(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_reverse_iteration)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorReverseIteration(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_count_gaps_count_characters_is_gap)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorCountGapsCountCharactersIsGap(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_clipped_count_gaps_count_characters_is_gap)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorClippedCountGapsCountCharactersIsGap(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_gap_operations_center)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorGapOperationsCenter(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_gap_operations_leading)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorGapOperationsLeading(TTag());
}

SEQAN_DEFINE_TEST(test_align_gaps_iterator_anchor_gap_operations_trailing)
{
    using namespace seqan;
    typedef AnchorGaps<> TTag;
    testAlignGapsIteratorGapOperationsTrailing(TTag());
}

#endif  // #ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GAPS_ITERATOR_H_
