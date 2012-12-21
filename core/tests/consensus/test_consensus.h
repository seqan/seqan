// ==========================================================================
//                              test_consensus.h
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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

#ifndef CORE_TESTS_CONSENSUS_TEST_CONSENSUS_H_
#define CORE_TESTS_CONSENSUS_TEST_CONSENSUS_H_

#include <seqan/basic.h>
#include <seqan/consensus.h>

template <typename TScoringScheme>
void testConsensusConsensusScoreSequenceEntry()
{
    using namespace seqan;

    {
        ConsensusScoreSequenceEntry<DnaString> consensusScoreSequenceEntry;
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._seq, (DnaString *) 0);
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._pos, 0);
    }

    {
        DnaString seq = "ACGTACGACT";
        ConsensusScoreSequenceEntry<DnaString> consensusScoreSequenceEntry(seq, 3);
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._seq, &seq);
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._pos, 3);
    }
}

template <typename TScoringScheme>
void testConsensusSequenceEntryForScoreMetafunction()
{
    using namespace seqan;

    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TSequenceEntry;

    bool result = IsSameType<TSequenceEntry, ConsensusScoreSequenceEntry<DnaString> >::VALUE;
    SEQAN_ASSERT_EQ(result, true);
}

template <typename TScoringScheme>
void testConsensusSequenceEntryForScore()
{
    using namespace seqan;

    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TSequenceEntry;

    DnaString seq = "ACGTACGTAC";

    TSequenceEntry entry = sequenceEntryForScore(TScoringScheme(), seq, 0);
    SEQAN_ASSERT_EQ(entry._pos, 0);
    SEQAN_ASSERT_EQ(*entry._seq, seq);

    entry = sequenceEntryForScore(TScoringScheme(), seq, 4);
    SEQAN_ASSERT_EQ(entry._pos, 4);
    SEQAN_ASSERT_EQ(*entry._seq, seq);
}

template <typename TScoringScheme>
void testConsensusValue()
{
    using namespace seqan;

    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TSequenceEntry;
    typedef typename Value<DnaString>::Type TDnaStringValue;

    DnaString seq = "ACGTACGTAC";

    TSequenceEntry entry = sequenceEntryForScore(TScoringScheme(), seq, 0);
    SEQAN_ASSERT_EQ(value(entry), 'A');

    entry = sequenceEntryForScore(TScoringScheme(), seq, 3);
    SEQAN_ASSERT_EQ(value(entry), 'T');
}

template <typename TScoringScheme>
void testConsensusPosition()
{
    using namespace seqan;

    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TSequenceEntry;
    typedef typename Value<DnaString>::Type TDnaStringValue;

    DnaString seq = "ACGTACGTAC";

    TSequenceEntry entry = sequenceEntryForScore(TScoringScheme(), seq, 0);
    SEQAN_ASSERT_EQ(position(entry), 0);

    entry = sequenceEntryForScore(TScoringScheme(), seq, 3);
    SEQAN_ASSERT_EQ(position(entry), 3);
}

template <typename TScoringScheme>
void testConsensusHost()
{
    using namespace seqan;

    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TSequenceEntry;
    typedef typename Value<DnaString>::Type TDnaStringValue;

    DnaString seq = "ACGTACGTAC";

    TSequenceEntry entry = sequenceEntryForScore(TScoringScheme(), seq, 0);
    SEQAN_ASSERT_EQ(host(entry), seq);
}

SEQAN_DEFINE_TEST(test_consensus_consensus_score_sequence_entry_consensus_score)
{
    typedef seqan::Score<int, seqan::ConsensusScore> TScore;
    testConsensusConsensusScoreSequenceEntry<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_consensus_score_sequence_entry_fractional_score)
{
    typedef seqan::Score<int, seqan::FractionalScore> TScore;
    testConsensusConsensusScoreSequenceEntry<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_consensus_score_sequence_entry_weightedconsensus_score)
{
    typedef seqan::Score<int, seqan::WeightedConsensusScore<seqan::ConsensusScore, seqan::FractionalScore> > TScore;
    testConsensusConsensusScoreSequenceEntry<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_sequence_entry_for_score_metafunction_consensus_score)
{
    typedef seqan::Score<int, seqan::ConsensusScore> TScore;
    testConsensusSequenceEntryForScoreMetafunction<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_sequence_entry_for_score_metafunction_fractional_score)
{
    typedef seqan::Score<int, seqan::FractionalScore> TScore;
    testConsensusSequenceEntryForScoreMetafunction<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_sequence_entry_for_score_metafunction_weightedconsensus_score)
{
    typedef seqan::Score<int, seqan::WeightedConsensusScore<seqan::ConsensusScore, seqan::FractionalScore> > TScore;
    testConsensusSequenceEntryForScoreMetafunction<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_sequence_entry_for_score_consensus_score)
{
    typedef seqan::Score<int, seqan::ConsensusScore> TScore;
    testConsensusSequenceEntryForScore<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_sequence_entry_for_score_fractional_score)
{
    typedef seqan::Score<int, seqan::FractionalScore> TScore;
    testConsensusSequenceEntryForScore<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_sequence_entry_for_score_weightedconsensus_score)
{
    typedef seqan::Score<int, seqan::WeightedConsensusScore<seqan::ConsensusScore, seqan::FractionalScore> > TScore;
    testConsensusSequenceEntryForScore<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_value_consensus_score)
{
    typedef seqan::Score<int, seqan::ConsensusScore> TScore;
    testConsensusValue<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_value_fractional_score)
{
    typedef seqan::Score<int, seqan::FractionalScore> TScore;
    testConsensusValue<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_value_weightedconsensus_score)
{
    typedef seqan::Score<int, seqan::WeightedConsensusScore<seqan::ConsensusScore, seqan::FractionalScore> > TScore;
    testConsensusValue<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_position_consensus_score)
{
    typedef seqan::Score<int, seqan::ConsensusScore> TScore;
    testConsensusPosition<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_position_fractional_score)
{
    typedef seqan::Score<int, seqan::FractionalScore> TScore;
    testConsensusPosition<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_position_weightedconsensus_score)
{
    typedef seqan::Score<int, seqan::WeightedConsensusScore<seqan::ConsensusScore, seqan::FractionalScore> > TScore;
    testConsensusPosition<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_host_consensus_score)
{
    typedef seqan::Score<int, seqan::ConsensusScore> TScore;
    testConsensusHost<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_host_fractional_score)
{
    typedef seqan::Score<int, seqan::FractionalScore> TScore;
    testConsensusHost<TScore>();
}

SEQAN_DEFINE_TEST(test_consensus_host_weightedconsensus_score)
{
    typedef seqan::Score<int, seqan::WeightedConsensusScore<seqan::ConsensusScore, seqan::FractionalScore> > TScore;
    testConsensusHost<TScore>();
}


#endif  // #ifndef CORE_TESTS_CONSENSUS_TEST_CONSENSUS_H_
