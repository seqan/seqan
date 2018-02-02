// ==========================================================================
//                              test_consensus.h
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

#ifndef TESTS_CONSENSUS_TEST_CONSENSUS_REALIGN_H_
#define TESTS_CONSENSUS_TEST_CONSENSUS_REALIGN_H_

#include <fstream>

#include <seqan/basic.h>
#include <seqan/consensus.h>

#include <seqan/misc/svg.h>

SEQAN_DEFINE_TEST(test_consensus_realign_one_contig_small)
{
    // Load example SAM from file.
    //
    // There are many superfluous gaps in the SAM file that we will get rid of below.
    seqan::FragmentStore<> store;
    seqan::CharString samPath = seqan::getAbsolutePath("/tests/consensus/small_example.sam");
    seqan::BamFileIn fSamIn(toCString(samPath));
    readRecords(store, fSamIn);

    seqan::AlignedReadLayout layout;
    // layoutAlignment(layout, store);
    // printAlignment(std::cout, seqan::Raw(), layout, store, 0, 0, 160, 0, 1000);

    // Call Realignment method.
    seqan::Score<int, seqan::WeightedConsensusScore<
                          seqan::Score<int, seqan::FractionalScore>,
                          seqan::Score<int, seqan::ConsensusScore> > > combinedScore;
    reAlignment(store, 0, 1, 30, false);

    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 160, 0, 1000);

    // Check Result.
    char const * expected =
            "TTCTATCTCCTATAGTCTGATATTACTGTAGGTACAGTAGCTTTTCTTCATTAATGTTTGCATAATATAGCTTCTTCCATGCTTTTACTTCCAATTATTTTGGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTATATATATATA-------------\n"
            ".....................................................................................................\n"
            "  .....................................................................................................\n"
            "                                   ......................................................................**...............................\n"
            "                                              .....................................................................................................\n";

    SEQAN_ASSERT_EQ(ss.str(), expected);
    // SEQAN_ASSERT_EQ(length(store.alignedReadStore), 4u);

    // SEQAN_ASSERT_EQ(store.alignedReadStore[0].beginPos, 0);
    // SEQAN_ASSERT_EQ(store.alignedReadStore[0].endPos, 101);
    // // SEQAN_ASSERT_EQ(length(store.alignedReadStore[0].cigar), 1u);
    // // SEQAN_ASSERT_EQ(store.alignedReadStore[0].cigar[0].count, 1u);

    // SEQAN_ASSERT_EQ(store.alignedReadStore[1].beginPos, 2);
    // SEQAN_ASSERT_EQ(store.alignedReadStore[1].endPos, 103);

    // SEQAN_ASSERT_EQ(store.alignedReadStore[2].beginPos, 33);
    // SEQAN_ASSERT_EQ(store.alignedReadStore[2].endPos, 136);

    // SEQAN_ASSERT_EQ(store.alignedReadStore[3].beginPos, 45);
    // SEQAN_ASSERT_EQ(store.alignedReadStore[3].endPos, 146);
}

#endif  // #ifndef TESTS_CONSENSUS_TEST_CONSENSUS_REALIGN_H_
