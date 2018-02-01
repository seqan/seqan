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

#include <seqan/basic.h>

#include "test_consensus.h"
#include "test_consensus_realign.h"

SEQAN_BEGIN_TESTSUITE(test_consensus)
{
    SEQAN_CALL_TEST(test_consensus_realign_one_contig_small);

    SEQAN_CALL_TEST(test_consensus_consensus_score_sequence_entry_consensus_score);
    SEQAN_CALL_TEST(test_consensus_consensus_score_sequence_entry_fractional_score);
    SEQAN_CALL_TEST(test_consensus_consensus_score_sequence_entry_weightedconsensus_score);

    SEQAN_CALL_TEST(test_consensus_sequence_entry_for_score_metafunction_consensus_score);
    SEQAN_CALL_TEST(test_consensus_sequence_entry_for_score_metafunction_fractional_score);
    SEQAN_CALL_TEST(test_consensus_sequence_entry_for_score_metafunction_weightedconsensus_score);

    SEQAN_CALL_TEST(test_consensus_sequence_entry_for_score_consensus_score);
    SEQAN_CALL_TEST(test_consensus_sequence_entry_for_score_fractional_score);
    SEQAN_CALL_TEST(test_consensus_sequence_entry_for_score_weightedconsensus_score);

    SEQAN_CALL_TEST(test_consensus_host_consensus_score);
    SEQAN_CALL_TEST(test_consensus_host_fractional_score);
    SEQAN_CALL_TEST(test_consensus_host_weightedconsensus_score);

    SEQAN_CALL_TEST(test_consensus_position_consensus_score);
    SEQAN_CALL_TEST(test_consensus_position_fractional_score);
    SEQAN_CALL_TEST(test_consensus_position_weightedconsensus_score);

    SEQAN_CALL_TEST(test_consensus_value_consensus_score);
    SEQAN_CALL_TEST(test_consensus_value_fractional_score);
    SEQAN_CALL_TEST(test_consensus_value_weightedconsensus_score);

//    SEQAN_CALL_TEST(test_consensus_write_celera_cgb);
//    SEQAN_CALL_TEST(test_consensus_write_celera_frg);
//    SEQAN_CALL_TEST(test_consensus_write_fasta_read_format);
//    SEQAN_CALL_TEST(test_consensus_convert_simple_read_file);
}
SEQAN_END_TESTSUITE
