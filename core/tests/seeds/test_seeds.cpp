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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the seeds module.
// ==========================================================================

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

// #include "test_basic_iter_indirect.h"
#include "test_seeds_combination.h"
#include "test_seeds_extension.h"
#include "test_seeds_global_chaining.h"
#include "test_seeds_seed_base.h"
#include "test_seeds_seed_chained.h"
#include "test_seeds_seed_diagonal.h"
#include "test_seeds_seed_set_base.h"
#include "test_seeds_seed_simple.h"
#include "test_align_banded_chain_impl.h"
#include "test_banded_chain_alignment_interface.h"

SEQAN_BEGIN_TESTSUITE(test_seeds)
{
    // Test indirect iterator.
    // SEQAN_CALL_TEST(test_seeds_basic_iter_indirect_constructors);
    // SEQAN_CALL_TEST(test_seeds_basic_iter_indirect_metafunctions);
    // SEQAN_CALL_TEST(test_seeds_basic_iter_indirect_basic_functions);

    // Tests for seed diagonals.
    SEQAN_CALL_TEST(test_seeds_seed_diagonal_constructors);
    SEQAN_CALL_TEST(test_seeds_seed_diagonal_metafunctions);

    // Tests for seeds.
    SEQAN_CALL_TEST(test_seeds_seed_base_constructors_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_metafunctions_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_getters_setters_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_basic_functions_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_assign_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_constructors_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_metafunctions_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_getters_setters_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_basic_functions_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_assign_chained);

    SEQAN_CALL_TEST(test_seeds_seed_chained_assign);
    SEQAN_CALL_TEST(test_seeds_seed_chained_metafunctions);
    SEQAN_CALL_TEST(test_seeds_seed_chained_append_diagonal);
    SEQAN_CALL_TEST(test_seeds_seed_chained_truncate_diagonals);
    SEQAN_CALL_TEST(test_seeds_seed_chained_iterators);
    SEQAN_CALL_TEST(test_seeds_seed_chained_front_back);

    SEQAN_CALL_TEST(test_seeds_seed_simple_constructors);
    SEQAN_CALL_TEST(test_seeds_seed_simple_setters);

    // Tests for the combination of seeds.
    SEQAN_CALL_TEST(test_seeds_combination_seeds_combineable_merge_chained);
    SEQAN_CALL_TEST(test_seeds_combination_seeds_combineable_simple_chaining_chained);
    SEQAN_CALL_TEST(test_seeds_combination_seeds_combineable_simple_chaos_chaining_chained);
    SEQAN_CALL_TEST(test_seeds_combination_combine_seeds_merge_chained);
    SEQAN_CALL_TEST(test_seeds_combination_combine_seeds_simple_chaining_chained);
    SEQAN_CALL_TEST(test_seeds_combination_combine_seeds_simple_chaos_chaining_chained);

    // Tests for unordered seed sets and simple seeds.
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_functions_simple_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_score_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_score_simple_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_right_merging_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_impossible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_scored_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_scored_simple_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_right_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_impossible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_scored_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_reached_scored_simple_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_right_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_impossible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_scored_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_scored_simple_unordered);

    // Tests for unordered seed sets and chained seeds.

    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_functions_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_score_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_score_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_right_merging_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_impossible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_scored_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_scored_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_right_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_impossible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_not_reached_scored_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_reached_scored_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_right_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_impossible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_scored_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_scored_chained_unordered);


    // Tests for seed extension algorithms
    SEQAN_CALL_TEST(test_seeds_extension_match_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_ungapped_xdrop_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_gapped_xdrop_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_match_extension_chained);
    SEQAN_CALL_TEST(test_seeds_extension_ungapped_xdrop_extension_chained);

    // Test global chaining of seeds.
    SEQAN_CALL_TEST(test_seeds_global_chaining_sparse_length);

    // Disabled the test for now.  Extension function contains a
    // force-failure assertion and instruction show to implement this.
    // See http://trac.mi.fu-berlin.de/seqan/ticket/344 for details.
    // 
    // TODO(holtgrew): Implement this.
    // 
    // SEQAN_CALL_TEST(test_seeds_extension_gapped_xdrop_extension_chained);

    // Tests for the banded chain alignment algorithm.
    SEQAN_CALL_TEST(test_banded_chain_alignment_empty_set_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_empty_set_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_one_seed_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_one_seed_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_two_seeds_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_two_seeds_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_three_seeds_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_three_seeds_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_special_seeds_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_special_seeds_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_band_extensions_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_band_extensions_affine);

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
}
SEQAN_END_TESTSUITE
