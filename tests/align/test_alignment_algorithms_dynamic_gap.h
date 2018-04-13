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
// Unit tests for the global interfaces to test dynamic gap costs.
// ==========================================================================

#ifndef TESTS_ALIGN_TEST_ALIGNMENT_ALGORITHMS_DYNAMIC_GAP_H_
#define TESTS_ALIGN_TEST_ALIGNMENT_ALGORITHMS_DYNAMIC_GAP_H_

#include <seqan/basic.h>
#include <seqan/align.h>

template <typename TAlignment, bool FreeTop, bool FreeLeft, bool FreeRight, bool FreeBottom, typename TACSpec,
          typename TScore>
bool _validateAlignment(TAlignment const & align,
                        seqan::AlignConfig<FreeTop, FreeLeft, FreeRight, FreeBottom, TACSpec> const & /*conf*/,
                        TScore const & score,
                        seqan::Score<TScore, seqan::Simple> const & scoreScheme)
{
    using namespace seqan;

    typedef typename Row<TAlignment const>::Type TRow;
    typedef typename Iterator<TRow, Standard>::Type TRowIterator;

    SEQAN_ASSERT_EQ(length(row(align,0)), length(row(align, 1)));

    if (length(row(align, 0)) == 0u && length(row(align, 1)) == 0u)
        return true;

    TRowIterator itH = begin(row(align, 0), Standard());
    TRowIterator itV = begin(row(align, 1), Standard());

    bool inGapH = false;
    bool inGapV = false;
    int eScore = 0;
    while (itH != end(row(align, 0), Standard()) && itV != end(row(align, 1), Standard()))
    {
        if (isGap(itH))
        {
            inGapV = false;
            if (!inGapH)
            {
                eScore += scoreGapOpen(scoreScheme);
                inGapH = true;
            }
            else
            {
                eScore += scoreGapExtend(scoreScheme);
            }
        }
        else if (isGap(itV))
        {
            inGapH = false;
            if (!inGapV)
            {
                eScore += scoreGapOpen(scoreScheme);
                inGapV = true;
            }
            else
            {
                eScore += scoreGapExtend(scoreScheme);
            }
        }
        else
        {
            inGapH = false;
            inGapV = false;
            if (*itH == *itV)
                eScore += scoreMatch(scoreScheme);
            else
                eScore += scoreMismatch(scoreScheme);
        }
        ++itH; ++itV;
    }

    // Resolve the begin gaps.
    if (FreeLeft)
    {
        itH = begin(row(align, 0), Standard());
        inGapH = false;
        while (isGap(itH))
        {
            if (inGapH)
                eScore -= scoreGapExtend(scoreScheme);
            else
            {
                inGapH = true;
                eScore -= scoreGapOpen(scoreScheme);
            }
            ++itH;
        }
    }
    if (FreeTop)
    {
        itV = begin(row(align, 1), Standard());
        inGapV = false;
        while (isGap(itV))
        {
            if (inGapV)
                eScore -= scoreGapExtend(scoreScheme);
            else
            {
                inGapV = true;
                eScore -= scoreGapOpen(scoreScheme);
            }
            ++itV;
        }
    }
    if (FreeRight)
    {
        itH = end(row(align, 0), Standard()) - 1;
        inGapH = false;
        while (isGap(itH))
        {
            if (inGapH)
                eScore -= scoreGapExtend(scoreScheme);
            else
            {
                inGapH = true;
                eScore -= scoreGapOpen(scoreScheme);
            }
            --itH;
        }
    }
    if (FreeBottom)
    {
        itV = end(row(align, 1), Standard()) - 1;
        inGapV = false;
        while (isGap(itV))
        {
            if (inGapV)
                eScore -= scoreGapExtend(scoreScheme);
            else
            {
                inGapV = true;
                eScore -= scoreGapOpen(scoreScheme);
            }
            --itV;
        }
    }

    return eScore == score;
}

template <typename TAlign, typename TAlgorithm, typename TAlignConfig>
void testDynamicGapInterfaces(TAlign & alignObj,
                              TAlgorithm const & /*algorithm*/,
                              TAlignConfig const & alignConfig,
                              int lDiag = 0,
                              int uDiag = 0)
{
    using namespace seqan;

    Score<int, Simple> scoreScheme(2, -2, -1, -4);

    if (IsSameType<TAlgorithm, SmithWaterman>::VALUE)
    {
        ignoreUnusedVariableWarning(alignConfig);
        int score;
        if (lDiag == 0 && uDiag == 0)
            score = localAlignment(alignObj, scoreScheme, DynamicGaps());
        else
            score = localAlignment(alignObj, scoreScheme, lDiag, uDiag, DynamicGaps());
        SEQAN_ASSERT(_validateAlignment(alignObj, alignConfig, score, scoreScheme));

    }
    else
    {
        // search global alignment
        int score;
        if (lDiag == 0 && uDiag == 0)
            score = globalAlignment(alignObj, scoreScheme, alignConfig, DynamicGaps());
        else
            score = globalAlignment(alignObj, scoreScheme, alignConfig, lDiag, uDiag, DynamicGaps());
        SEQAN_ASSERT(_validateAlignment(alignObj, alignConfig, score, scoreScheme));
    }
}

template <typename TAlgorithm, typename TAlignConfig>
void testDynamicGapInterfaces(TAlgorithm const & algo,
                              TAlignConfig const & alignConfig,
                              int lDiag = 0,
                              int uDiag = 0)
{
    using namespace seqan;

    {
        DnaString str1 = "AAAAACTACGTACGTTTCTGGCCCCC";  // A G T A G C T A C G T A C G T T T C T G G A T G A C
        DnaString str2 = "CCCCCCACGTGTTACGTACGTAAAAA";   // G G T G A C A C G T G T T A C G T A C G T A A
        Align<DnaString> alignObj;
        resize(rows(alignObj), 2);
        assignSource(row(alignObj, 0), str1);
        assignSource(row(alignObj, 1), str2);
        testDynamicGapInterfaces(alignObj, algo, alignConfig, lDiag, uDiag);
    }

    {
        DnaString str1 = "AAAAACTACGTACGTTTCTGGCCCCC";
        DnaString str2 = "CCCCCCACGTGTTACGTACGTAAAAA";
        Align<DnaString> alignObj;
        resize(rows(alignObj), 2);
        assignSource(row(alignObj, 0), str1);
        assignSource(row(alignObj, 1), str2);
        testDynamicGapInterfaces(alignObj, algo, alignConfig, lDiag, uDiag);
    }

    {
        DnaString str1 = "AAAAAGGGGTTTT";
        DnaString str2 = "AAAGTT";
        Align<DnaString> alignObj;
        resize(rows(alignObj), 2);
        assignSource(row(alignObj, 0), str1);
        assignSource(row(alignObj, 1), str2);
        testDynamicGapInterfaces(alignObj, algo, alignConfig, lDiag, uDiag);
    }

    {
        DnaString str1 = "AAAAAATTTTTGGG";
        DnaString str2 = "TTTTTTTTGGGGGGGG";
        Align<DnaString> alignObj;
        resize(rows(alignObj), 2);
        assignSource(row(alignObj, 0), str1);
        assignSource(row(alignObj, 1), str2);
        testDynamicGapInterfaces(alignObj, algo, alignConfig, lDiag, uDiag);
    }


    {
        DnaString str1 = "GGGGCTTTTTTAAAGAGCGCCCTTTTTTGGGG";
        DnaString str2 = "AAAACTTTTTTGGTTTTTTAAAA";
        Align<DnaString> alignObj;
        resize(rows(alignObj), 2);
        assignSource(row(alignObj, 0), str1);
        assignSource(row(alignObj, 1), str2);
        testDynamicGapInterfaces(alignObj, algo, alignConfig, lDiag, uDiag);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_global_dynamic_cost)
{
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<>());
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<false, true, true, false>());
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<true, false, false, true>());
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<true, true, true, true>());
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_global_dynamic_cost_banded)
{
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<>(), -3, 3);
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<false, true, true, false>(), -3, 3);
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<true, false, false, true>(), -3, 3);
    testDynamicGapInterfaces(seqan::NeedlemanWunsch(), seqan::AlignConfig<true, true, true, true>(), -3, 3);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_local_dynamic_cost)
{
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>());
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>());
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>());
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>());
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_local_dynamic_cost_banded)
{
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>(), -3, 3);
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>(), -3, 3);
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>(), -3, 3);
    testDynamicGapInterfaces(seqan::SmithWaterman(), seqan::AlignConfig<>(), -3, 3);
}

#endif // TESTS_ALIGN_TEST_ALIGNMENT_ALGORITHMS_DYNAMIC_GAP_H_
