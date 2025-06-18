// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2025, Knut Reinert, FU Berlin
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

#include <seqan/basic.h>
#include <seqan/stream.h>  // for printing seqan2::String<>

#include <seqan/seeds.h>

#include <string>

using TSeed = seqan2::Seed<seqan2::Simple>;

SEQAN_DEFINE_TEST(test_issue_2540_no_glue_point_with_default_band_extension)
{
    using namespace std::literals;
                                     // 0               1               2         3
                                     // 012345      6789012      345678901234567890123
    std::string const query =          "GTGCCG"s + "TTCTTTT"s + "TTTTTTTTTTTTTTTTTTTTT"s;
    std::string const ref =   "CCTTTTTACGTGCTT"s + "TTCTTTT"s + "CTTGCATGTAAATCTCTTTATTTTATTTTATTTTA"s;
                            // 0         1               2               3         4         5
                            // 012345678901234      5678901      23456789012345678901234567890123456

    seqan2::Align<seqan2::DnaString, seqan2::ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), query);
    assignSource(row(alignment, 1), ref);

    seqan2::String<TSeed> chain{};
    appendValue(chain, TSeed{6, 15, 13, 22});

    seqan2::AlignConfig<false, true, false, false> alignConfig;

    seqan2::Score<int, seqan2::Simple> scoringScheme(2, -1, -2);
    try {
        seqan2::bandedChainAlignment(alignment, chain, scoringScheme, scoringScheme, alignConfig, 15);
    } catch (std::logic_error const & e){
        SEQAN_ASSERT(!std::string{e.what()}.empty());
    }
}

SEQAN_BEGIN_TESTSUITE(test_banded_align_chain_issue_2540)
{
    SEQAN_CALL_TEST(test_issue_2540_no_glue_point_with_default_band_extension);
}
SEQAN_END_TESTSUITE
