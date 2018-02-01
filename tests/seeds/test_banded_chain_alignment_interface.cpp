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
// Tests the different global interfaces of the banded chain alignment.
// ==========================================================================

#include <vector>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/stream.h>  // for printing seqan::String<>

#include <seqan/seeds.h>

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_linear_global_one_score)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGG--CACA");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGG--CACA");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_linear_global_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -47);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGG--CACA");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -47);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGG--CACA");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_linear_semi_one_score)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 9);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGGCACA--");
    }

    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2, 5, 6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<false, true, true, false>(),  2);
        SEQAN_ASSERT_EQ(result, 9);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGG--CACA");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_linear_semi_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);

        SEQAN_ASSERT_EQ(result, -11);
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "CGAATCCAT--CC--CACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "---GGCGATNNNCATGGCACA");
    }

    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(),  2);
        SEQAN_ASSERT_EQ(result, -147);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "------CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGGCAC------A");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_linear_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 13);
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGGCACA--");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_linear_overlap_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);

        SEQAN_ASSERT_EQ(result, -42);
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "------CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGGCACA------");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_affine_global_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 8);
    
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGGCACA--");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 8);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGGCACA--");
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_affine_global_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -98);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGA-----ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGGCAC------A");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -98);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGA-----ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGGCAC------A");
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_affine_semi_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 11);
    
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGGCACA--");
    }

    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, 11);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGGCACA--");
    }

}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_affine_semi_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, -8);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "CGAATCCAT--CC----CACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "---GGCGATNNNCATGGCACA--");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, -97);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGGCAC------A");
        SEQAN_ASSERT_EQ(row(alignment, 0), "------CGA-ATCCATCCCACACA");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_affine_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 14);
    
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 0), "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCG-ATNNNCATGGCACA--");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_align_affine_overlap_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceH = "CGAATCCATCCCACACA";
        CharString sequenceV = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);

        SEQAN_ASSERT_EQ(result, -42);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(row(alignment, 1), "GGCGATNNNCATGGCACA------");
        SEQAN_ASSERT_EQ(row(alignment, 0), "------CGA-ATCCATCCCACACA");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_linear_global_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGG--CACA");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGG--CACA");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_linear_global_two_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -47);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical,   "GGCGATNNNCATGG--CACA");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -47);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCGATNNNCATGG--CACA");
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_linear_semi_one_score)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 9);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }

    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<false, true, true, false>(),  2);
        SEQAN_ASSERT_EQ(result, 9);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal,  "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGG--CACA");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_linear_semi_two_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);

        SEQAN_ASSERT_EQ(result, -11);
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "CGAATCCAT--CC--CACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "---GGCGATNNNCATGGCACA");
    }

    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(),  2);
        SEQAN_ASSERT_EQ(result, -147);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "------CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCGATNNNCATGGCAC------A");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_linear_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 13);
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_linear_overlap_two_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);

        SEQAN_ASSERT_EQ(result, -42);
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "------CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCGATNNNCATGGCACA------");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_affine_global_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 8);
    
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 8);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_affine_global_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -98);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGA-----ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCGATNNNCATGGCAC------A");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -98);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGA-----ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCGATNNNCATGGCAC------A");
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_affine_semi_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 11);
    
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }

    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, 11);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }

}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_affine_semi_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, -8);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "CGAATCCAT--CC----CACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "---GGCGATNNNCATGGCACA--");
    }

    // Test on whole strings without AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, -97);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "------CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCGATNNNCATGGCAC------A");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_affine_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 14);
    
        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_gaps_affine_overlap_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        CharString sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);

        SEQAN_ASSERT_EQ(result, -42);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "------CGA-ATCCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCGATNNNCATGGCACA------");
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_linear_global_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 5);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||    ||||\n"
                   << "        GGCG-ATNNNCATGG--CACA\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 5);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||    ||||\n"
                   << "        GGCG-ATNNNCATGG--CACA\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_linear_global_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -47);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    : \n"
                   << "        --CGA-ATCCATCCCACACA\n"
                   << "          |||    |||    ||||\n"
                   << "        GGCGATNNNCATGG--CACA\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Test on whole strings without AlignConfig.
    {

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -47);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    : \n"
                   << "        --CGA-ATCCATCCCACACA\n"
                   << "          |||    |||    ||||\n"
                   << "        GGCGATNNNCATGG--CACA\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_linear_semi_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 9);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||  ||||  \n"
                   << "        GGCG-ATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);
        
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<false, true, true, false>(),  2);
        SEQAN_ASSERT_EQ(result, 9);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||    ||||\n"
                   << "        GGCG-ATNNNCATGG--CACA\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_linear_semi_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, -11);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        CGAATCCAT--CC--CACACA\n"
                   << "             | ||   |    ||||\n"
                   << "        ---GGCGATNNNCATGGCACA\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Test on whole strings without AlignConfig.
    {

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, -147);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :     \n"
                   << "        ------CGA-ATCCATCCCACACA\n"
                   << "                  ||    |      |\n"
                   << "        GGCGATNNNCATGGCAC------A\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_linear_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 13);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||  ||||  \n"
                   << "        GGCG-ATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_linear_overlap_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));


        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, -42);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :     \n"
                   << "        ------CGA-ATCCATCCCACACA\n"
                   << "                  ||    |       \n"
                   << "        GGCGATNNNCATGGCACA------\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}


SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_affine_global_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));
        
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);
        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 8);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||  ||||  \n"
                   << "        GGCG-ATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Test on whole strings without AlignConfig.
    {
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 8);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||  ||||  \n"
                   << "        GGCG-ATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_affine_global_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -98);
        
        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :     \n"
                   << "        --CGA-----ATCCATCCCACACA\n"
                   << "          |||     ||    |      |\n"
                   << "        GGCGATNNNCATGGCAC------A\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Test on whole strings without AlignConfig.
    {

        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -98);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :     \n"
                   << "        --CGA-----ATCCATCCCACACA\n"
                   << "          |||     ||    |      |\n"
                   << "        GGCGATNNNCATGGCAC------A\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_affine_semi_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    {
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 11);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||  ||||  \n"
                   << "        GGCG-ATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    {
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, 11);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||  ||||  \n"
                   << "        GGCG-ATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_affine_semi_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, -8);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :    \n"
                   << "        CGAATCCAT--CC----CACACA\n"
                   << "             | ||   |    ||||  \n"
                   << "        ---GGCGATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Test on whole strings without AlignConfig.
    {
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, -97);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :     \n"
                   << "        ------CGA-ATCCATCCCACACA\n"
                   << "                  ||    |      |\n"
                   << "        GGCGATNNNCATGGCAC------A\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_affine_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    {
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 14);
    
        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :  \n"
                   << "        --CGAAT--CCATCCCACACA\n"
                   << "          || ||   |||  ||||  \n"
                   << "        GGCG-ATNNNCATGGCACA--\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_alignmentgraph_affine_overlap_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef StringSet<CharString> TStringSet;
    typedef StringSet<CharString, Dependent<> > TDependentStringSet;
    typedef Graph<Alignment<TDependentStringSet, void> > TAlignmentGraph;
    // Test on whole strings with AlignConfig.
    {
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        TAlignmentGraph alignGraph(strings);

        int result = bandedChainAlignment(alignGraph, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);

        SEQAN_ASSERT_EQ(result, -42);

        std::stringstream ss;
        ss << alignGraph;

        // Compare alignment graph.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :     \n"
                   << "        ------CGA-ATCCATCCCACACA\n"
                   << "                  ||    |       \n"
                   << "        GGCGATNNNCATGGCACA------\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_linear_global_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 5);

        SEQAN_ASSERT(length(fragments) == 4u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 13, 1, 14, 4));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 5, 1, 8, 6));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[3] == TFragment(0, 0, 1, 2, 2));
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 5);

        SEQAN_ASSERT(length(fragments) == 4u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 13, 1, 14, 4));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 5, 1, 8, 6));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[3] == TFragment(0, 0, 1, 2, 2));
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_linear_global_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -47);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 13, 1, 14, 4));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 6, 8));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 3));
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -47);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 13, 1, 14, 4));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 6, 8));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 3));
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_linear_semi_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 9);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 5, 1, 8, 10));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 2));
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, 9);

        SEQAN_ASSERT(length(fragments) == 4u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 13, 1, 14, 4));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 5, 1, 8, 6));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[3] == TFragment(0, 0, 1, 2, 2));
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_linear_semi_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, -11);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 11, 1, 12, 6));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 9, 1, 8, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 3, 1, 0, 6));
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, -147);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 16, 1, 17, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 10, 7));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 6, 3));
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_linear_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 13);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 5, 1, 8, 10));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 2));
    }

}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_linear_overlap_two_scores)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");

        Score<int, Simple> scoringSchemeAnchor(5, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, -42);

        SEQAN_ASSERT(length(fragments) == 2u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 10, 8));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 0, 1, 6, 3));
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_affine_global_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, 8);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 5, 1, 8, 10));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 2));
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, 2);
        SEQAN_ASSERT_EQ(result, 8);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 5, 1, 8, 10));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 2));
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_affine_global_two_scores)
{
    using namespace seqan;
    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, false, false, false>(), 2);
        SEQAN_ASSERT_EQ(result, -98);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 16, 1, 17, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 10, 7));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 3));
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, 2);
        SEQAN_ASSERT_EQ(result, -98);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 16, 1, 17, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 10, 7));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 3));
    }
}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_affine_semi_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 11);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 5, 1, 8, 10));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 2));
    }

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, 11);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 5, 1, 8, 10));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 2));
    }

}
SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_affine_semi_two_scores)
{
    using namespace seqan;
    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    typedef Seed<Simple> TSeed;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, -8);
        
        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 11, 1, 14, 4));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 9, 1, 8, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 3, 1, 0, 6));
    }

    // Test on whole strings without AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<false, true, true, false>(), 2);
        SEQAN_ASSERT_EQ(result, -97);

        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 16, 1, 17, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 10, 7));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 6, 3));
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_affine_overlap_one_score)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;
    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringScheme(2, -1, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringScheme, AlignConfig<true, true, true, true>(), 2);
        SEQAN_ASSERT_EQ(result, 14);
        
        SEQAN_ASSERT(length(fragments) == 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 5, 1, 8, 10));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 4, 2));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 2, 2));
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_fragments_affine_overlap_two_scores)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;
    // Test on whole strings with AlignConfig.
    {
        TStringSet strings;
        appendValue(strings, "CGAATCCATCCCACACA");
        appendValue(strings, "GGCGATNNNCATGGCACA");
        Score<int, Simple> scoringSchemeAnchor(5, -10, -10, -20);
        Score<int, Simple> scoringSchemeGap(0, -2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        TFragmentString fragments;
        int result = bandedChainAlignment(fragments, strings, seedChain, scoringSchemeAnchor, scoringSchemeGap, AlignConfig<true, true, true, true>(), 2);

        SEQAN_ASSERT_EQ(result, -42);
        SEQAN_ASSERT(length(fragments) == 2u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 10, 8));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 0, 1, 6, 3));
    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_stl_vector_adaption)
{
	using namespace seqan;

    typedef Seed<Simple> TSeed;
    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        std::vector<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<true, false, false, true>(), 2);
        SEQAN_ASSERT_EQ(result, 9);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal, "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGGCACA--");
    }

    {
        CharString sequenceV = "CGAATCCATCCCACACA";
        Dna5String sequenceH = "GGCGATNNNCATGGCACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        std::vector<TSeed> seedChain;
        appendValue(seedChain, TSeed(0,2,5,6));
        appendValue(seedChain, TSeed(6,9,9,12));
        appendValue(seedChain, TSeed(11,14,17,16));

        Gaps<CharString, ArrayGaps> gapsHorizontal;
        Gaps<Dna5String, ArrayGaps> gapsVertical;
        assignSource(gapsHorizontal, sequenceV);
        assignSource(gapsVertical, sequenceH);

        int result = bandedChainAlignment(gapsHorizontal, gapsVertical, seedChain, scoringScheme, AlignConfig<false, true, true, false>(),  2);
        SEQAN_ASSERT_EQ(result, 9);

        // Compare alignment rows.
        SEQAN_ASSERT_EQ(gapsHorizontal,  "--CGAAT--CCATCCCACACA");
        SEQAN_ASSERT_EQ(gapsVertical, "GGCG-ATNNNCATGG--CACA");
    }
}

SEQAN_BEGIN_TESTSUITE(test_banded_chain_alignment_interface)
{
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_linear_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_linear_global_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_linear_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_linear_semi_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_linear_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_linear_overlap_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_affine_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_affine_global_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_affine_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_affine_semi_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_affine_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_align_affine_overlap_two_scores);

    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_linear_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_linear_global_two_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_linear_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_linear_semi_two_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_linear_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_linear_overlap_two_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_affine_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_affine_global_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_affine_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_affine_semi_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_affine_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_gaps_affine_overlap_two_scores);

    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_linear_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_linear_global_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_linear_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_linear_semi_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_linear_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_linear_overlap_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_affine_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_affine_global_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_affine_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_affine_semi_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_affine_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_alignmentgraph_affine_overlap_two_scores);

    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_linear_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_linear_global_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_linear_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_linear_semi_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_linear_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_linear_overlap_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_affine_global_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_affine_global_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_affine_semi_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_affine_semi_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_affine_overlap_one_score);
    SEQAN_CALL_TEST(test_banded_chain_alignment_fragments_affine_overlap_two_scores);
    SEQAN_CALL_TEST(test_banded_chain_alignment_stl_vector_adaption);
}
SEQAN_END_TESTSUITE
