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
// Tests for the consensus module extensions.
// ==========================================================================

#include <sstream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/consensus.h>

// A test for consensusAlignment() with full coordinates.

SEQAN_DEFINE_TEST(test_consensus_consensus_alignment_coordinates)
{
    // -----------------------------------------------------------------------
    // Prepare Input
    // -----------------------------------------------------------------------

    // Build a FragmentStore with simulated reads from a reference sequence.
    seqan::Dna5String ref = 
            "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTT"
            "CAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGTGGCAT";

    // Read length and step width for reads.
    int const READ_LENGTH = 50;
    int const STEP = 5;

    // Compute reads and append to FragmentStore.
    seqan::FragmentStore<> store;
    for (unsigned pos = 0, i = 0; pos + READ_LENGTH < length(ref); pos += STEP, ++i)
    {
        // Append a new read sequence.
        unsigned readID = appendRead(store, infix(ref, pos, pos + READ_LENGTH));
        // Create small perturbation of the position but not left of position 0.
        int pos2 = std::max(0, (int)pos + ((int)i % 6 - 3));
        // Append a new read alignment for the just added read.
        appendAlignedRead(store, readID, 0, pos2, pos2 + READ_LENGTH);
    }

    // -----------------------------------------------------------------------
    // Compute Consensus Alignment
    // -----------------------------------------------------------------------

    // Compute consensus alignment.
    seqan::ConsensusAlignmentOptions options;
    consensusAlignment(store, options);

    // -----------------------------------------------------------------------
    // Check Results
    // -----------------------------------------------------------------------

    SEQAN_ASSERT_EQ(length(store.contigStore), 1u);

    // Print final consensus alignment into buffer.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, /*contigID=*/0, /*beginPos=*/0, /*endPos=*/(int)length(ref), 0, 30);
    
    char const * EXPECTED =
            "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGT-----\n"
            "..................................................     ..................................................\n"
            "     ..................................................     ..................................................\n"
            "          ..................................................     ..................................................\n"
            "               ..................................................     ..................................................\n"
            "                    ..................................................     ..................................................\n"
            "                         ..................................................     ..................................................\n"
            "                              ..................................................     ..................................................\n"
            "                                   ..................................................\n"
            "                                        ..................................................\n"
            "                                             ..................................................\n"
            "                                                  ..................................................\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// A test for consensusAlignment() with contig IDs only.

SEQAN_DEFINE_TEST(test_consensus_consensus_alignment_contig_ids)
{
    // -----------------------------------------------------------------------
    // Prepare Input
    // -----------------------------------------------------------------------

    // Build a FragmentStore with simulated reads from a reference sequence.  We will assign them to two contigs, however.
    seqan::Dna5String ref = 
            "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTT"
            "CAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGTGGCAT";

    // Read length and step width for reads.
    int const READ_LENGTH = 50;
    int const STEP = 5;

    // Compute reads and append to FragmentStore.
    seqan::FragmentStore<> store;
    for (unsigned pos = 0, i = 0; pos + READ_LENGTH < length(ref); pos += STEP, ++i)
    {
        // Append a new read sequence.
        unsigned readID = appendRead(store, infix(ref, pos, pos + READ_LENGTH));
        // Create small perturbation of the position but not left of position 0.
        int pos2 = std::max(0, (int)pos + ((int)i % 6 - 3));
        // Round-robin of contig sequences.
        unsigned contigID = (pos / STEP) % 2u;
        // Append a new read alignment for the just added read.
        appendAlignedRead(store, readID, contigID, pos2, pos2 + READ_LENGTH);
    }

    // -----------------------------------------------------------------------
    // Compute Consensus Alignment
    // -----------------------------------------------------------------------

    // Compute consensus alignment.
    seqan::ConsensusAlignmentOptions options;
    options.usePositions = false;
    consensusAlignment(store, options);

    // -----------------------------------------------------------------------
    // Check Results
    // -----------------------------------------------------------------------

    SEQAN_ASSERT_EQ(length(store.contigStore), 2u);

    seqan::AlignedReadLayout layout;

    // Check first contig.
    {
        layoutAlignment(layout, store);
        std::stringstream ss;
        printAlignment(ss, layout, store, /*contigID=*/0, /*beginPos=*/0, /*endPos=*/(int)length(ref), 0, 30);

        char const * EXPECTED = 
                "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATAT----------\n"
                "..................................................          ..................................................\n"
                "          ..................................................          ..................................................\n"
                "                    ..................................................          ..................................................\n"
                "                              ..................................................\n"
                "                                        ..................................................\n"
                "                                                  ..................................................\n";
        SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
    }
    // Check second contig.
    {
        layoutAlignment(layout, store);
        std::stringstream ss;
        printAlignment(ss, layout, store, /*contigID=*/1, /*beginPos=*/0, /*endPos=*/(int)length(ref), 0, 30);


        char const * EXPECTED = 
                "ATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGT----------\n"
                "..................................................          ..................................................\n"
                "          ..................................................          ..................................................\n"
                "                    ..................................................          ..................................................\n"
                "                              ..................................................\n"
                "                                        ..................................................\n"
                "                                                  ..................................................\n";
        SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
    }
}

// A test for consensusAlignment() without using contig IDs.

SEQAN_DEFINE_TEST(test_consensus_consensus_alignment_no_contig_ids)
{
    // -----------------------------------------------------------------------
    // Prepare Input
    // -----------------------------------------------------------------------

    // Build a FragmentStore with simulated reads from a reference sequence.
    seqan::Dna5String ref = 
            "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTT"
            "CAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGTGGCAT";

    // Read length and step width for reads.
    int const READ_LENGTH = 50;
    int const STEP = 5;

    // Compute reads and append to FragmentStore.
    seqan::FragmentStore<> store;
    for (unsigned pos = 0, i = 0; pos + READ_LENGTH < length(ref); pos += STEP, ++i)
        appendRead(store, infix(ref, pos, pos + READ_LENGTH));

    // -----------------------------------------------------------------------
    // Compute Consensus Alignment
    // -----------------------------------------------------------------------

    // Compute consensus alignment.
    seqan::ConsensusAlignmentOptions options;
    options.useContigID = false;
    consensusAlignment(store, options);

    // -----------------------------------------------------------------------
    // Check Results
    // -----------------------------------------------------------------------

    SEQAN_ASSERT_EQ(length(store.contigStore), 1u);

    // Print final consensus alignment into buffer.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, /*contigID=*/0, /*beginPos=*/0, /*endPos=*/(int)length(ref), 0, 30);
    
    char const * EXPECTED =
            "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGT-----\n"
            "..................................................     ..................................................\n"
            "     ..................................................     ..................................................\n"
            "          ..................................................     ..................................................\n"
            "               ..................................................     ..................................................\n"
            "                    ..................................................     ..................................................\n"
            "                         ..................................................     ..................................................\n"
            "                              ..................................................     ..................................................\n"
            "                                   ..................................................\n"
            "                                        ..................................................\n"
            "                                             ..................................................\n"
            "                                                  ..................................................\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// A test for consensusAlignment() without using global alignments.

SEQAN_DEFINE_TEST(test_consensus_consensus_alignment_global_alignment)
{
    // -----------------------------------------------------------------------
    // Prepare Input
    // -----------------------------------------------------------------------

    // Build a FragmentStore with infixes from a reference sequence.
    seqan::Dna5String ref = 
            "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGTGGCAT";

    seqan::FragmentStore<> store;
    appendRead(store, "GAATACATCTCTAAAGAGCTTGATGCTAATTTGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTT");
    appendRead(store, "TCTAAAGAGCTTGGATGCTAAAATAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAA");
    appendRead(store, "ACATCTCTAAAGAGCTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCT");

    // -----------------------------------------------------------------------
    // Compute Consensus Alignment
    // -----------------------------------------------------------------------

    // Compute consensus alignment.
    seqan::ConsensusAlignmentOptions options;
    options.useGlobalAlignment = true;
    consensusAlignment(store, options);

    // -----------------------------------------------------------------------
    // Check Results
    // -----------------------------------------------------------------------

    SEQAN_ASSERT_EQ(length(store.contigStore), 1u);

    // Print final consensus alignment into buffer.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, /*contigID=*/0, /*beginPos=*/0, /*endPos=*/(int)length(ref), 0, 30);
    
    char const * EXPECTED =
            "GAATACATCTCTAAAGAGCTT-GATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCT---------------------------------------------\n"
            ".....................*...........*...........................................................\n"
            "    ................**.........................................................................\n"
            "         ............G........AA.......................................................\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_BEGIN_TESTSUITE(test_consensus)
{
    SEQAN_CALL_TEST(test_consensus_consensus_alignment_coordinates);
    SEQAN_CALL_TEST(test_consensus_consensus_alignment_contig_ids);
    SEQAN_CALL_TEST(test_consensus_consensus_alignment_no_contig_ids);
    SEQAN_CALL_TEST(test_consensus_consensus_alignment_global_alignment);
}
SEQAN_END_TESTSUITE
