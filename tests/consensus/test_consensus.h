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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_CONSENSUS_TEST_CONSENSUS_H_
#define TESTS_CONSENSUS_TEST_CONSENSUS_H_

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/consensus.h>

template <typename TScoringScheme>
void testConsensusConsensusScoreSequenceEntry()
{
    using namespace seqan;

    {
        typedef typename Position<DnaString>::Type TPosition;
        ConsensusScoreSequenceEntry<DnaString> consensusScoreSequenceEntry;
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._seq, (DnaString *) 0);
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._pos, static_cast<TPosition>(0));
    }

    {
        typedef typename Position<DnaString>::Type TPosition;
        DnaString seq = "ACGTACGACT";
        ConsensusScoreSequenceEntry<DnaString> consensusScoreSequenceEntry(seq, 3);
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._seq, &seq);
        SEQAN_ASSERT_EQ(consensusScoreSequenceEntry._pos, static_cast<TPosition>(3));
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
    typedef typename Position<DnaString>::Type TPosition;

    DnaString seq = "ACGTACGTAC";

    TSequenceEntry entry = sequenceEntryForScore(TScoringScheme(), seq, 0);
    SEQAN_ASSERT_EQ(entry._pos, static_cast<TPosition>(0));
    SEQAN_ASSERT_EQ(*entry._seq, seq);

    entry = sequenceEntryForScore(TScoringScheme(), seq, 4);
    SEQAN_ASSERT_EQ(entry._pos, static_cast<TPosition>(4));
    SEQAN_ASSERT_EQ(*entry._seq, seq);
}

template <typename TScoringScheme>
void testConsensusValue()
{
    using namespace seqan;

    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TSequenceEntry;
    //typedef typename Value<DnaString>::Type TDnaStringValue;

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
    //typedef typename Value<DnaString>::Type TDnaStringValue;
    typedef typename Position<DnaString>::Type TPosition;

    DnaString seq = "ACGTACGTAC";

    TSequenceEntry entry = sequenceEntryForScore(TScoringScheme(), seq, 0);
    SEQAN_ASSERT_EQ(position(entry), static_cast<TPosition>(0));

    entry = sequenceEntryForScore(TScoringScheme(), seq, 3);
    SEQAN_ASSERT_EQ(position(entry), static_cast<TPosition>(3));
}

template <typename TScoringScheme>
void testConsensusHost()
{
    using namespace seqan;

    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TSequenceEntry;
    //typedef typename Value<DnaString>::Type TDnaStringValue;

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

SEQAN_DEFINE_TEST(test_consensus_write_celera_cgb)
{
    // Get path to input files.
    seqan::CharString inPathSam = seqan::getAbsolutePath("/tests/consensus/toy.sam");
    // Get path to temporary file.
    seqan::CharString outPathCgb = SEQAN_TEMP_FILENAME();

    // Read in SAM and FASTA.
    seqan::FragmentStore<> store;
    seqan::BamFileIn fSamIn(toCString(inPathSam));
    readRecords(store, fSamIn);

    // Write out CGB file.
    std::fstream fCgbOut(toCString(outPathCgb), std::ios::binary | std::ios::out);
    _writeCeleraCgb(fCgbOut, store);
    fCgbOut.close();

    // Compare result.
    seqan::CharString goldPathCgb = seqan::getAbsolutePath("/tests/consensus/sam_to_cgb_result.cgb");
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(toCString(outPathCgb), toCString(goldPathCgb)));
}

SEQAN_DEFINE_TEST(test_consensus_write_celera_frg)
{
    // Get path to input files.
    seqan::CharString inPathSam = seqan::getAbsolutePath("/tests/consensus/toy.sam");
    // Get path to temporary file.
    seqan::CharString outPathFrg = SEQAN_TEMP_FILENAME();

    // Read in SAM and FASTA.
    seqan::FragmentStore<> store;
    seqan::BamFileIn fSamIn(toCString(inPathSam));
    readRecords(store, fSamIn);

    // Write out FRG file.
    std::fstream fFrgOut(toCString(outPathFrg), std::ios::binary | std::ios::out);
    _writeCeleraFrg(fFrgOut, store);
    fFrgOut.close();

    // Compare result.
    seqan::CharString goldPathFrg = seqan::getAbsolutePath("/tests/consensus/sam_to_frg_result.frg");
    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(outPathFrg), toCString(goldPathFrg)));
}

SEQAN_DEFINE_TEST(test_consensus_write_fasta_read_format)
{
    // Get path to input files.
    seqan::CharString inPathSam = seqan::getAbsolutePath("/tests/consensus/toy.sam");
    seqan::CharString inPathFasta = seqan::getAbsolutePath( "/tests/consensus/toy.fa");
    // Get path to temporary file.
    seqan::CharString outPathFasta = SEQAN_TEMP_FILENAME();

    // Read in SAM and FASTA.
    seqan::FragmentStore<> store;
    SEQAN_ASSERT(loadContigs(store, toCString(inPathFasta)));
    seqan::BamFileIn fSamIn(toCString(inPathSam));
    readRecords(store, fSamIn);

    // Write out FASTA file.
    std::fstream fFastaOut(toCString(outPathFasta), std::ios::binary | std::ios::out);
    write(fFastaOut, store, seqan::FastaReadFormat());
    fFastaOut.close();

    // Compare result.
    seqan::CharString goldPathFasta = seqan::getAbsolutePath("/tests/consensus/sam_to_fasta_read_result.fa");
    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(outPathFasta), toCString(goldPathFasta)));
}

SEQAN_DEFINE_TEST(test_consensus_convert_simple_read_file)
{
    // Get path to input files.
    seqan::CharString inPathFasta = seqan::getAbsolutePath("/tests/consensus/simulated_reads.fasta");
    std::string filePath(toCString(inPathFasta));
    // Get path to temporary file.
    std::string outPathSam = (std::string)SEQAN_TEMP_FILENAME() + ".sam";

    // Read in FASTA file.
    seqan::FragmentStore<> store;
    std::fstream fFastaIn(toCString(inPathFasta), std::ios::binary | std::ios::in);
    SEQAN_ASSERT(fFastaIn.good());
    _convertSimpleReadFile(fFastaIn, store, filePath, false);

    // Write out as SAM.
    seqan::BamFileOut fSamOut(outPathSam.c_str());
    writeRecords(fSamOut, store);
    close(fSamOut);

    // Compare result.
    seqan::CharString goldPathSam = seqan::getAbsolutePath("/tests/consensus/reads_to_sam_read_result.sam");
    SEQAN_ASSERT(seqan::_compareTextFiles(outPathSam.c_str(), toCString(goldPathSam)));
}

#endif  // #ifndef TESTS_CONSENSUS_TEST_CONSENSUS_H_
