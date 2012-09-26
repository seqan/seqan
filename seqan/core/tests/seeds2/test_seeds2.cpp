// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#include <seqan/seeds2.h>  // Include module under test.

#include "test_align_chain_banded.h"
#include "test_align_dynprog_affine.h"
#include "test_align_dynprog_banded_affine.h"
#include "test_align_dynprog_banded_linear.h"
#include "test_align_dynprog_linear.h"
#include "test_basic_iter_indirect.h"
#include "test_seeds_combination.h"
#include "test_seeds_extension.h"
#include "test_seeds_global_chaining.h"
#include "test_seeds_seed_base.h"
#include "test_seeds_seed_chained.h"
#include "test_seeds_seed_diagonal.h"
#include "test_seeds_seed_set_base.h"
#include "test_seeds_seed_simple.h"

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

    // Tests for seed algorithms
    SEQAN_CALL_TEST(test_seeds_extension_match_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_ungapped_xdrop_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_gapped_xdrop_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_match_extension_chained);
    SEQAN_CALL_TEST(test_seeds_extension_ungapped_xdrop_extension_chained);
    // Disabled the test for now.  Extension function contains a
    // force-failure assertion and instruction show to implement this.
    // See http://trac.mi.fu-berlin.de/seqan/ticket/344 for details.
    // 
    // TODO(holtgrew): Implement this.
    // 
    // SEQAN_CALL_TEST(test_seeds_extension_gapped_xdrop_extension_chained);

    // Tests for the banded chain alignment algorithms.
    SEQAN_CALL_TEST(test_align_chain_banded_compute_upper_left_overlap);
    SEQAN_CALL_TEST(test_align_chain_banded_compute_lower_right_overlap);
    //SEQAN_CALL_TEST(test_align_chain_banded_align_linear);
    // TODO(holtgrew): Not testing banded chain alignment with affine gap costs for now.
    // SEQAN_CALL_TEST(test_align_chain_banded_align_affine);

    // Tests for the classic NW dynamic programming.
    SEQAN_CALL_TEST(test_align_dynprog_linear_resize_matrix);
    SEQAN_CALL_TEST(test_align_dynprog_linear_init_gutter_free);
    SEQAN_CALL_TEST(test_align_dynprog_linear_init_gutter_not_free);
    SEQAN_CALL_TEST(test_align_dynprog_linear_fill_matrix);
    //SEQAN_CALL_TEST(test_align_dynprog_linear_traceback);
    // Tests for the banded NW dynamic programming.
    SEQAN_CALL_TEST(test_align_dynprog_banded_linear_resize_matrix);
    SEQAN_CALL_TEST(test_align_dynprog_banded_linear_init_gutter_free);
    SEQAN_CALL_TEST(test_align_dynprog_banded_linear_init_gutter_not_free);
    SEQAN_CALL_TEST(test_align_dynprog_banded_linear_fill_matrix);
    SEQAN_CALL_TEST(test_align_dynprog_banded_linear_traceback);

    // Tests for the classic Gotoh dynamic programming.
    SEQAN_CALL_TEST(test_align_dynprog_affine_resize_matrix);
    SEQAN_CALL_TEST(test_align_dynprog_affine_init_gutter_free);
    SEQAN_CALL_TEST(test_align_dynprog_affine_init_gutter_not_free);
    SEQAN_CALL_TEST(test_align_dynprog_affine_fill_matrix);
    //SEQAN_CALL_TEST(test_align_dynprog_affine_traceback);
    // Tests for the banded Gotoh dynamic programming.
    SEQAN_CALL_TEST(test_align_dynprog_banded_affine_resize_matrix);
    SEQAN_CALL_TEST(test_align_dynprog_banded_affine_init_gutter_free);
    SEQAN_CALL_TEST(test_align_dynprog_banded_affine_init_gutter_not_free);
    SEQAN_CALL_TEST(test_align_dynprog_banded_affine_fill_matrix);
    //SEQAN_CALL_TEST(test_align_dynprog_banded_affine_traceback);

    // Test global chaining of seeds.
    SEQAN_CALL_TEST(test_seeds_global_chaining_sparse_length);
}
SEQAN_END_TESTSUITE
