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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for align_extend
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_EXTEND_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_EXTEND_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/align_extend.h>

SEQAN_DEFINE_TEST(test_align_extend_simple)
{
    using namespace seqan;
    Score<int> sc(2, -1, -2);

    Align<typename Infix<CharString const>::Type, ArrayGaps> align;

    resize(rows(align), 2);

    // plain string, no gaps, stop extension before ends of both strings
    {
        CharString const s1("NNNNNNNNNNTTCCGGGAC" "GGTA""CACACACGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA""CAGGCACGGTTTTTTTT");

        assignSource(row(align, 0), infix(s1, 19, 23));
        assignSource(row(align, 1), infix(s2, 8, 12));

        int score = globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {19u,8u,23u,12u} };
        score = extendAlignment(align,
                                score,
                                s1,
                                s2,
                                positions,
                                EXTEND_BOTH,
                                sc);

        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""CAGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);

        SEQAN_ASSERT_EQ(score, 12 + 8 +  12);
    }

    // gaps in to-be-extended-regions, stop extension before ends of both strings
    {
        CharString const s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTTTT");

        assignSource(row(align, 0), infix(s1, 18, 22));
        assignSource(row(align, 1), infix(s2, 8, 12));

        int score = globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        score = extendAlignment(align,
                                score,
                                s1,
                                s2,
                                positions,
                                EXTEND_BOTH,
                                sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);

        SEQAN_ASSERT_EQ(score, 8 + 8 + 8);
    }

    // gap in center region, stop extension before ends of both strings
    {
        CharString const s1("NNNNNNNNNNTTCCGGGAC" "GGTA""CACACACGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC"  "GTA""CAGGCACGGTTTTTTTT");

        assignSource(row(align, 0), infix(s1, 19, 23));
        assignSource(row(align, 1), infix(s2, 8, 11));

        int score = globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {19u,8u,23u,11u} };
        score = extendAlignment(align,
                                score,
                                s1,
                                s2,
                                positions,
                                EXTEND_BOTH,
                                sc);

        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""-GTA""CAGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);

        SEQAN_ASSERT_EQ(score, 12 + 4 + 12);
    }
}

SEQAN_DEFINE_TEST(test_align_extend_simple_infixes_recomp)
{
    // same as above with the exception that the source strings
    // are picked from a concat direct string set (thereby being infixes)
    // this also extensively tests https://reviewboard.seqan.de/r/325/

    // in addition this test set tests the recompute capability of alignExtend
    // by omitting origScore

    using namespace seqan;
    Score<int> sc(2, -1, -2);

    Align<typename Infix<CharString const>::Type, ArrayGaps> align;

    resize(rows(align), 2);

    // plain string, no gaps, stop extension before ends of both strings
    {
        CharString s1("NNNNNNNNNNTTCCGGGAC" "GGTA""CACACACGGGGGGGGGG");
        CharString s2(           "CTCGGGAC" "GGTA""CAGGCACGGTTTTTTTT");

        StringSet<CharString, Owner<ConcatDirect<> > > set;
        appendValue(set, CharString("XXXXXXXXXXX"));
        appendValue(set, s1);
        appendValue(set, s2);
        appendValue(set, CharString("XXXXXXXXXXX"));

        assignSource(row(align, 0), infix(value(set,1), 19, 23));
        assignSource(row(align, 1), infix(value(set,2), 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {19u,8u,23u,12u} };
        int score = extendAlignment(align,
                                    value(set,1),
                                    value(set,2),
                                    positions,
                                    EXTEND_BOTH,
                                    sc);

        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""CAGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);

        SEQAN_ASSERT_EQ(score, 12 + 8 +  12);
    }

    // gaps in to-be-extended-regions, stop extension before ends of both strings
    {
        CharString s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGG");
        CharString s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTTTT");

        StringSet<CharString, Owner<ConcatDirect<> > > set;
        appendValue(set, CharString("XXXXXXXXXXX"));
        appendValue(set, s1);
        appendValue(set, s2);
        appendValue(set, CharString("XXXXXXXXXXX"));

        assignSource(row(align, 0), infix(value(set, 1), 18, 22));
        assignSource(row(align, 1), infix(value(set, 2), 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        int score = extendAlignment(align,
                                    value(set,1),
                                    value(set,2),
                                    positions,
                                    EXTEND_BOTH,
                                    sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);

        SEQAN_ASSERT_EQ(score, 8 + 8 + 8);
    }

    // gap in center region, stop extension before ends of both strings
    {
        CharString s1("NNNNNNNNNNTTCCGGGAC" "GGTA""CACACACGGGGGGGGGG");
        CharString s2(           "CTCGGGAC"  "GTA""CAGGCACGGTTTTTTTT");

        StringSet<CharString, Owner<ConcatDirect<> > > set;
        appendValue(set, CharString("XXXXXXXXXXX"));
        appendValue(set, s1);
        appendValue(set, s2);
        appendValue(set, CharString("XXXXXXXXXXX"));

        assignSource(row(align, 0), infix(value(set, 1), 19, 23));
        assignSource(row(align, 1), infix(value(set, 2), 8, 11));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {19u,8u,23u,11u} };
        int score = extendAlignment(align,
                                    value(set,1),
                                    value(set,2),
                                    positions,
                                    EXTEND_BOTH,
                                    sc);

        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""-GTA""CAGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);

        SEQAN_ASSERT_EQ(score, 12 + 4 + 12);
    }
}

SEQAN_DEFINE_TEST(test_align_extend_banded)
{
    using namespace seqan;
    Score<int> sc(2, -1, -2);

    Align<typename Infix<CharString const>::Type, ArrayGaps> align;

    resize(rows(align), 2);

    // gaps in to-be-extended-regions, inside band
    {
        CharString const s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTTTT");

        assignSource(row(align, 0), infix(s1, 18, 22));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        -1,
                        +1,
                        sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);
    }

    // gaps in to-be-extended-regions, inside band
    {
        CharString const s1("NNNNNNNNNNTTCCGGG"   "GGTA" "CACACACGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA"  "AGGCACGGTTTTTTTT");

        assignSource(row(align, 0), infix(s1, 17, 21));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {17u,8u,21u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        -2,
                        +2,
                        sc);

        SEQAN_ASSERT_EQ(CharString("CGGG--""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);
    }

    // gaps in to-be-extended-regions, outside of band (resulting in "shift")
    {
        CharString const s1("NNNNNNNNNNTTCCGGG"   "GGTA" "CACACACGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA"  "AGGCACGGTTTTTTTT");

        assignSource(row(align, 0), infix(s1, 17, 21));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {17u,8u,21u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        -1,
                        +1,
                        sc);

        SEQAN_ASSERT_EQ(CharString("TCCGGG-""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("TCGGGAC""GGTA""-AGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 11);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 1);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 31);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);
    }

}


SEQAN_DEFINE_TEST(test_align_extend_xdrop)
{
    using namespace seqan;
    Score<int> sc(2, -1, -2);

    Align<typename Infix<CharString const>::Type, ArrayGaps> align;

    resize(rows(align), 2);

    // no XDrop -> alignment passes through local minimum of scores
    {
        CharString const s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTGGGG");

        assignSource(row(align, 0), infix(s1, 18, 22));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGGGGGGGGGGG"),
                        row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGGTTTTTGGGG"),
                        row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 41);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 30);
    }

    // XDrop -> alignment stops before local minimum of scores
    {
        CharString const s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTGGGG");

        assignSource(row(align, 0), infix(s1, 18, 22));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        3,
                        sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);
    }
}

SEQAN_DEFINE_TEST(test_align_extend_xdrop_banded)
{
    using namespace seqan;
    Score<int> sc(2, -1, -2);

    Align<typename Infix<CharString const>::Type, ArrayGaps> align;

    resize(rows(align), 2);

    // no XDrop -> alignment spans local minimum of scores
    {
        CharString const s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTGGGG");

        assignSource(row(align, 0), infix(s1, 18, 22));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGGGGGGGGGGG"),
                        row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGGTTTTTGGGG"),
                        row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 41);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 30);
    }

    // XDrop -> alignment spans local minimum of scores,
    // because there is no band and XDropOff not reached
    {
        CharString const s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTGGGG");

        assignSource(row(align, 0), infix(s1, 18, 22));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        4,
                        sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGGGGGGGGGGG"),
                        row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGGTTTTTGGGG"),
                        row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 41);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 30);
    }

    // XDrop -> alignment doesnt span local minimum of scores,
    // because the same XDropOff as above is reached, because of the band
    {
        CharString const s1("NNNNNNNNNNTTCCGGGA"  "GGTA""CACACACGGGGGGGGGGG");
        CharString const s2(           "CTCGGGAC" "GGTA" "AGGCACGGTTTTTGGGG");

        assignSource(row(align, 0), infix(s1, 18, 22));
        assignSource(row(align, 1), infix(s2, 8, 12));

        globalAlignment(align, sc);

        Tuple<unsigned, 4> const positions = { {18u,8u,22u,12u} };
        extendAlignment(align,
                        s1,
                        s2,
                        positions,
                        EXTEND_BOTH,
                        -2,
                        +2,
                        4,
                        sc);

        SEQAN_ASSERT_EQ(CharString("CGGGA-""GGTA""CACACACGG"), row(align, 0));
        SEQAN_ASSERT_EQ(CharString("CGGGAC""GGTA""-AGGCACGG"), row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 13);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 2);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 32);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 21);
    }

}

SEQAN_DEFINE_TEST(test_align_extend_semiglobal)
{
    using namespace seqan;
    typedef Align<typename Infix<CharString const>::Type, ArrayGaps> TAlign;
    Score<int> sc(1, -1, -1);

    //                                                                                |---------- INFIX ---------|
    //                0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    CharString s1 = "TGACGATCCAGCGCGCACAGCAGGAGGACTCGGCCGTGTATCTCTGTGCCAGCAGCTTAGGGGACACGTACGAGCAGTACTTCGGGCCAGGCACGCTTCT";
    //                                                                                 | |||||||||||||||||||| |||||    ||
    CharString s2 =                                                                 "CTCCTACGAGCAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAG";
    //                                                                                01234567890123456789012345678901234567890123456
    //                                                                                                |-- INFIX -|

    TAlign alignOrig;
    resize(rows(alignOrig), 2);
    Segment<CharString> seg1(s1, 64, 92), seg2(s2, 16, 28);
    assignSource(row(alignOrig, 0), seg1);
    assignSource(row(alignOrig, 1), seg2);
    globalAlignment(alignOrig, sc, AlignConfig<true,false,false,true>());

    // Only right extension, left gaps not clipped
    {
        TAlign align(alignOrig);

        Tuple<unsigned, 4> positions = { { (unsigned)(beginPosition(seg1)+beginPosition(row(align, 0))),
                                           (unsigned)(beginPosition(seg2)+beginPosition(row(align, 1))),
                                           92u, 28u } };

        extendAlignment(align,
                s1,
                s2,
                positions,
                EXTEND_RIGHT,
                -2,
                2,
                sc);

        SEQAN_ASSERT_EQ(CharString("CACGTACGAGCAGTACTTCGGGCCAGGCAC"),
                        row(align, 0));
        SEQAN_ASSERT_EQ(CharString("----------------TTCGGGCCGGGCAC"),
                        row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 64);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 16);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 94);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 46);
    }

    // Only right extension, left gaps partially clipped
    {
        TAlign align(alignOrig);

        setClippedBeginPosition(row(align, 0), 14);
        setClippedBeginPosition(row(align, 1), 14);

        Tuple<unsigned, 4> positions = { { (unsigned)(beginPosition(seg1)+beginPosition(row(align, 0))),
                                           (unsigned)(beginPosition(seg2)+beginPosition(row(align, 1))),
                                           92u, 28u } };

        extendAlignment(align,
                s1,
                s2,
                positions,
                EXTEND_RIGHT,
                -2,
                2,
                sc);

        SEQAN_ASSERT_EQ(CharString("ACTTCGGGCCAGGCAC"),
                        row(align, 0));
        SEQAN_ASSERT_EQ(CharString("--TTCGGGCCGGGCAC"),
                        row(align, 1));


        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 78);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 16);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 94);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 32);
    }

    // Extension in both directions, left gaps partially clipped,
    // generates "poor" alignment at the joint between the center
    // and the left alignment
    {
        TAlign align(alignOrig);

        setClippedBeginPosition(row(align, 0), 14);
        setClippedBeginPosition(row(align, 1), 14);

        Tuple<unsigned, 4> positions = { { (unsigned)(beginPosition(seg1)+beginPosition(row(align, 0))),
                                           (unsigned)(beginPosition(seg2)+beginPosition(row(align, 1))),
                                           92u, 28u } };

        extendAlignment(align,
                s1,
                s2,
                positions,
                EXTEND_BOTH,
                -2,
                2,
                sc);

        SEQAN_ASSERT_EQ(CharString("TACGAGCAGT--ACTTCGGGCCAGGCAC"),
                        row(align, 0));
        SEQAN_ASSERT_EQ(CharString("TACGAGCAGTAC--TTCGGGCCGGGCAC"),
                        row(align, 1));

        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 0)), 68);
        SEQAN_ASSERT_EQ(clippedBeginPosition(row(align, 1)), 4);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 0)), 96);
        SEQAN_ASSERT_EQ(clippedEndPosition(row(align, 1)), 32);
    }

}

#endif  // SEQAN_TESTS_ALIGN_SPLIT_TEST_ALIGN_SPLIT_H_
