// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

//#define TEST_DEBUG_OUTPUT


#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_delta_map.h"
#include "test_journaled_string_tree.h"
#include "test_journaled_string_tree_traverse.h"

SEQAN_BEGIN_TESTSUITE(test_journaled_string_tree)
{
    // Tests for delta map.
    SEQAN_CALL_TEST(test_delta_map_value_mf);
    SEQAN_CALL_TEST(test_delta_map_reference_mf);
    SEQAN_CALL_TEST(test_delta_map_delta_value_mf);
    SEQAN_CALL_TEST(test_delta_map_delta_coverage_mf);
    SEQAN_CALL_TEST(test_delta_map_insert);
    SEQAN_CALL_TEST(test_delta_map_length);
    SEQAN_CALL_TEST(test_delta_map_empty);
    SEQAN_CALL_TEST(test_delta_map_coverage_size);
    SEQAN_CALL_TEST(test_delta_map_set_coverage_size);
    SEQAN_CALL_TEST(test_delta_map_iterator);
    SEQAN_CALL_TEST(test_delta_map_iterator_copy_constructor);
    SEQAN_CALL_TEST(test_delta_map_iterator_assign);
    SEQAN_CALL_TEST(test_delta_map_iterator_value);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_type);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_position);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_snp);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_ins);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_del);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_indel);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_coverage);

    // Tests for journaled string tree
    SEQAN_CALL_TEST(test_journaled_string_tree_container_mf);
    SEQAN_CALL_TEST(test_journaled_string_tree_get_string_tree_mf);
    SEQAN_CALL_TEST(test_journaled_string_tree_host_mf);
    SEQAN_CALL_TEST(test_journaled_string_tree_constructor);
    SEQAN_CALL_TEST(test_journaled_string_tree_init);
    SEQAN_CALL_TEST(test_journaled_string_tree_reinit);
    SEQAN_CALL_TEST(test_journaled_string_tree_container);
    SEQAN_CALL_TEST(test_journaled_string_tree_string_set);
    SEQAN_CALL_TEST(test_journaled_string_tree_full_journal_required);
    SEQAN_CALL_TEST(test_journaled_string_tree_set_block_size);
    SEQAN_CALL_TEST(test_journaled_string_tree_block_size);
    SEQAN_CALL_TEST(test_journaled_string_tree_journal_next_block);
    SEQAN_CALL_TEST(test_journaled_string_tree_host);
    SEQAN_CALL_TEST(test_journaled_string_tree_virtual_block_position);

	// ----------------------------------------------------------------------------
    // Test all variants being SNPs.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_journaled_string_tree_jst_traversal_concept);

	// Test all position 0, all snps, different coverages.
	SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_0_journaled_string_tree);
	SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_1_journaled_string_tree);
	SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_2_journaled_string_tree);
	SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_3_journaled_string_tree);
	SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_4_journaled_string_tree);
	SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_5_journaled_string_tree);

	// Test all position 30, all snps, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_5_journaled_string_tree);

    // Test all position 100, all snps, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_5_journaled_string_tree);

	// Test equidistant position including 0, all snps, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_5_journaled_string_tree);

    // Test equidistant position including last position, all snps, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_5_journaled_string_tree);

    // ----------------------------------------------------------------------------
    // Test all variants being deletionss.
    // ----------------------------------------------------------------------------

    // Test all position 0, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_5_journaled_string_tree);

    // Test all position 30, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_5_journaled_string_tree);

    // Test all position 100, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_5_journaled_string_tree);

    // Test equidistant position including 0, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_5_journaled_string_tree);

    // Test equidistant position including last position, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_5_journaled_string_tree);

    // ----------------------------------------------------------------------------
    // Test all variants being insertions.
    // ----------------------------------------------------------------------------

    // Test all position 0, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_5_journaled_string_tree);

    // Test all position 30, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_5_journaled_string_tree);

    // Test all position 100, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_5_journaled_string_tree);

    // Test equidistant position including 0, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_5_journaled_string_tree);

    // Test equidistant position including last position, different coverages.
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_5_journaled_string_tree);

    // ----------------------------------------------------------------------------
    // Test special variant combinations.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_0_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_2_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_5_journaled_string_tree);

    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_6_4_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_7_4_3_journaled_string_tree);

    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_1_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_3_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_4_journaled_string_tree);
    SEQAN_CALL_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_5_journaled_string_tree);
}
SEQAN_END_TESTSUITE
