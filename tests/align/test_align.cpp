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

#include <seqan/basic.h>
#include <seqan/stream.h>

#include "test_align_stream.h"
#include "test_align_fragment.h"
#include "test_align_gaps.h"
#include "test_align_gaps_iterator.h"

// TODO(holtgrew): Test Align<>, AlignCols<>!

#include "test_alignment_dp_band.h"
#include "test_alignment_dp_cell.h"
#include "test_alignment_dp_profile.h"
#include "test_alignment_dp_formula.h"
#include "test_alignment_dp_matrix.h"
#include "test_alignment_dp_matrix_navigator.h"
#include "test_alignment_dp_trace_segment.h"
#include "test_alignment_dp_adapt_tracesegments.h"
#include "test_alignment_dp_traceback.h"

#include "test_alignment_algorithms_band_position.h"
#include "test_alignment_algorithms_global.h"
#include "test_alignment_algorithms_global_banded.h"
#include "test_alignment_algorithms_local.h"
#include "test_alignment_algorithms_local_banded.h"
#include "test_alignment_algorithms_dynamic_gap.h"
#include "test_align_global_alignment_specialized.h"

#include "test_align_alignment_operations.h"
#include "test_evaluate_alignment.h"

SEQAN_BEGIN_TESTSUITE(test_align)
{
    // -----------------------------------------------------------------------
    // Test Gaps Data Structures.
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_gaps_array_gaps_metafunctions);
    SEQAN_CALL_TEST(test_align_gaps_array_constructor_and_source);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_set_source);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_assign_source);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_gap_operations_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_gap_operations_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_gap_operations_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_ungapped);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_iterator_interface_begin);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_iterator_interface_end);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_iterator_interface_iter);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_ungapped);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_ungapped);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clear_clipping);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_copy_gaps);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_copy_clipping);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_is_nothing);
    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clear);

    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_metafunctions);
    SEQAN_CALL_TEST(test_align_gaps_anchor_constructor_and_source);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_set_source);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_assign_source);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_gap_operations_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_gap_operations_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_gap_operations_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_ungapped);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_iterator_interface_begin);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_iterator_interface_end);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_iterator_interface_iter);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_ungapped);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_ungapped);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_gaps_center);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_gaps_leading);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_gaps_trailing);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clear_clipping);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_copy_gaps);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_copy_clipping);
    // TODO(holtgrew): Extend anchor gaps such that this works.
    // SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_is_nothing);
    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clear);

    // -----------------------------------------------------------------------
    // Test for Fragment
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_fragment);

    // -----------------------------------------------------------------------
    // Test Gaps Iterators.
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_gaps_iterator_array_metafunctions);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_trivial_iterator_array_functions);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_rooted_random_access_iterator_array_functions);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_movement);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_relations);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_pointer_arithmetic);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_forward_iteration);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_reverse_iteration);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_count_gaps_count_characters_is_gap);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_clipped_count_gaps_count_characters_is_gap);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_gap_operations_center);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_gap_operations_leading);
    SEQAN_CALL_TEST(test_align_gaps_iterator_array_gap_operations_trailing);

    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_metafunctions);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_trivial_iterator_anchor_functions);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_rooted_random_access_iterator_anchor_functions);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_movement);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_relations);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_pointer_arithmetic);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_forward_iteration);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_reverse_iteration);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_count_gaps_count_characters_is_gap);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_clipped_count_gaps_count_characters_is_gap);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_gap_operations_center);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_gap_operations_leading);
    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_gap_operations_trailing);

    // ----------------------------------------------------------------------------
    // Test DPProfile.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_alignment_dp_profile_is_global_alignment);
    SEQAN_CALL_TEST(test_alignment_dp_profile_is_local_alignment);
    SEQAN_CALL_TEST(test_alignment_dp_profile_is_traceback_enabled);
    SEQAN_CALL_TEST(test_alignment_dp_profile_is_free_end_gaps);

    // ----------------------------------------------------------------------------
    // Test DPBandConfig.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_dp_band_on_constructor);
    SEQAN_CALL_TEST(test_dp_band_on_lower_diagonal);
    SEQAN_CALL_TEST(test_dp_band_on_upper_diagonal);
    SEQAN_CALL_TEST(test_dp_band_on_set_lower_diagonal);
    SEQAN_CALL_TEST(test_dp_band_on_set_upper_diagonal);
    SEQAN_CALL_TEST(test_dp_band_off_band_size);
    SEQAN_CALL_TEST(test_dp_band_on_band_size);

    // ----------------------------------------------------------------------------
    // Test DPCell.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_dp_cell_value);
    SEQAN_CALL_TEST(test_dp_cell_reference);
    SEQAN_CALL_TEST(test_dp_cell_default_infinity);

    SEQAN_CALL_TEST(test_dp_cell_linear_constructor);
    SEQAN_CALL_TEST(test_dp_cell_linear_copy_constructor);
    SEQAN_CALL_TEST(test_dp_cell_linear_assignment);
    SEQAN_CALL_TEST(test_dp_cell_linear_score);

    SEQAN_CALL_TEST(test_dp_cell_affine_constructor);
    SEQAN_CALL_TEST(test_dp_cell_affine_copy_constructor);
    SEQAN_CALL_TEST(test_dp_cell_affine_assignment);
    SEQAN_CALL_TEST(test_dp_cell_affine_score);
    SEQAN_CALL_TEST(test_dp_cell_affine_vertical_score);
    SEQAN_CALL_TEST(test_dp_cell_affine_horizontal_score);

    SEQAN_CALL_TEST(test_dp_cell_dynamic_constructor);
    SEQAN_CALL_TEST(test_dp_cell_dynamic_copy_constructor);
    SEQAN_CALL_TEST(test_dp_cell_dynamic_assignment);
    SEQAN_CALL_TEST(test_dp_cell_dynamic_score);

    // ----------------------------------------------------------------------------
    // Test DPMatrix.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_alignment_dp_matrix_metafunction_data_host);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_metafunction_size_arr);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_data_host);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_data_lengths);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_data_factors);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_check_dimension);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_clear);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_host);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_set_host);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_begin_standard);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_begin_rooted);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_end_standard);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_end_rooted);

    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_constructor);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_copy_constructor);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_assigment);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_reference);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_getvalue);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_position);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_size);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_host);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_iterator_standard);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_metafunction_iterator_rooted);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_resize);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_resize_with_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_value_with_coordinates);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_set_length);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_length_dimension);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_length);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_empty);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_full_coordinate);

    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_constructor);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_copy_constructor);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_assigment);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_reference);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_getvalue);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_position);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_size);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_host);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_iterator_standard);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_metafunction_iterator_rooted);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_resize);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_resize_with_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_value_with_coordinates);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_set_length);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_length_dimension);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_length);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_empty);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_sparse_coordinate);

    // ----------------------------------------------------------------------------
    // Test DPMatrix Navigator.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_init_unbanded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_init_banded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_go_next_cell);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_assign_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_previous_cell_horizontal);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_coordinate);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_container);

    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_init_unbanded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_init_banded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_go_next);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_assign_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_previous_cell_horizontal);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_coordinate);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_container);

    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_init_unbanded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_disabled_init_unbanded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_init_banded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_disabled_init_banded);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_go_next);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_assign_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_value);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_coordinate);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_container);
    SEQAN_CALL_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_to_global_position);

    // ----------------------------------------------------------------------------
    // Test Recursion Formula.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_dp_formula_trace_global_linear_diagonal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_linear_horizontal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_linear_vertical_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_linear_upper_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_linear_lower_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_linear_all_direction);

    SEQAN_CALL_TEST(test_dp_formula_trace_global_affine_diagonal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_affine_horizontal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_affine_vertical_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_affine_upper_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_affine_lower_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_affine_all_direction);

    SEQAN_CALL_TEST(test_dp_formula_trace_global_dynamic_diagonal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_dynamic_horizontal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_dynamic_vertical_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_dynamic_upper_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_dynamic_lower_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_global_dynamic_all_direction);

    SEQAN_CALL_TEST(test_dp_formula_trace_local_linear_diagonal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_linear_horizontal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_linear_vertical_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_linear_upper_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_linear_lower_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_linear_all_direction);

    SEQAN_CALL_TEST(test_dp_formula_trace_local_affine_diagonal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_affine_horizontal_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_affine_vertical_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_affine_upper_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_affine_lower_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_trace_local_affine_all_direction);

    SEQAN_CALL_TEST(test_dp_formula_notrace_diagonal_direction);
    SEQAN_CALL_TEST(test_dp_formula_notrace_horizontal_direction);
    SEQAN_CALL_TEST(test_dp_formula_notrace_vertical_direction);
    SEQAN_CALL_TEST(test_dp_formula_notrace_upper_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_notrace_lower_band_direction);
    SEQAN_CALL_TEST(test_dp_formula_notrace_all_direction);


    // ----------------------------------------------------------------------------
    // Test Trace Segment.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_constructor);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_assignment);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_position);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_size);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_begin_horizontal);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_begin_vertical);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_end_horizontal);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_end_vertical);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_translate_trace_value);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_operator_stream);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_operator_equal);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_operator_unequal);
    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_record_segment);

    // ----------------------------------------------------------------------------
    // Test Adaptor.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align2_trace_adaptor_trace_segment);
    SEQAN_CALL_TEST(test_align2_trace_adaptor_record_trace_segment);
    SEQAN_CALL_TEST(test_align2_trace_adaptor_adapt_file);
    SEQAN_CALL_TEST(test_align2_trace_adaptor_adapt_align);
    SEQAN_CALL_TEST(test_align2_trace_adaptor_adapt_fragments);
    SEQAN_CALL_TEST(test_align2_trace_adaptor_adapt_alignment_graph);

    // ----------------------------------------------------------------------------
    // Test Traceback.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align2_traceback_linear_unbanded_alignment);
    SEQAN_CALL_TEST(test_align2_traceback_linear_normal_banded_alignment);
    SEQAN_CALL_TEST(test_align2_traceback_linear_small_banded_alignment);
    SEQAN_CALL_TEST(test_align2_traceback_linear_wide_banded_alignment);
    SEQAN_CALL_TEST(test_align2_traceback_affine);
    SEQAN_CALL_TEST(test_align2_traceback_gaps_left_linear_gaps);
    SEQAN_CALL_TEST(test_align2_traceback_gaps_right_linear_gaps);
    SEQAN_CALL_TEST(test_align2_traceback_gaps_left_affine_gaps);
    SEQAN_CALL_TEST(test_align2_traceback_gaps_right_affine_gaps);

    // ----------------------------------------------------------------------------
    // Test Band Locations.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case1);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case2);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case3);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case4);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case5);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case6);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case7);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case8);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case9);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case10);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case11);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case12);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case13);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case14);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case15);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case16);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case17);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case18);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case19);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case20);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case21);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case22);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case23);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case24);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case25);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case26);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case27);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case28);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case29);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case30);
    SEQAN_CALL_TEST(test_alignment_algorithms_band_position_case31);

	// Integration tests of dp algorithms.

    // ----------------------------------------------------------------------------
    // Unbanded Alignment Algorithms.
    // ----------------------------------------------------------------------------

    // Global Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_global_linear);
    SEQAN_CALL_TEST(test_align_global_alignment_shorter_interfaces_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_global_linear);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_global_affine);
    SEQAN_CALL_TEST(test_align_global_alignment_shorter_interfaces_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_global_affine);

    // Overlap Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_overlap_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_overlap_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_overlap_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_overlap_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_overlap_linear);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_overlap_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_overlap_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_overlap_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_overlap_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_overlap_affine);

    // Semi-Global Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_semi_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_semi_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_semi_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_semi_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_semi_global_linear);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_semi_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_semi_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_semi_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_semi_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_semi_global_affine);

    // Global Alignment with Differnt Container Types
    SEQAN_CALL_TEST(test_alignment_algorithms_global_different_container);

    // Local Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_local_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_local_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_local_linear);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_local_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_local_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_local_affine);

    // Dynamic Gaps.
    SEQAN_CALL_TEST(test_alignment_algorithms_global_dynamic_cost);
    SEQAN_CALL_TEST(test_alignment_algorithms_local_dynamic_cost);

    // Suboptimal Alignment.

    // This is working on the old  module - should be replaced by the new module at some stage.
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_align);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_gaps);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_fragment);

    // TODO(rmaerker): Here are the tests that should run when the Waterman-Eggert is adapted
//    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_suboptimal_linear);
//    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_suboptimal_linear);
//    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_suboptimal_linear);
//    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_suboptimal_linear);
//
//    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_suboptimal_affine);
//    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_suboptimal_affine);
//    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_suboptimal_affine);
//    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_suboptimal_affine);

    // ----------------------------------------------------------------------------
    // Banded Alignment Algorithms.
    // ----------------------------------------------------------------------------

    // Global Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_global_linear_banded);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_global_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_global_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_global_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_global_affine_banded);

    // Overlap Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_overlap_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_overlap_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_overlap_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_overlap_linear_banded);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_overlap_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_overlap_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_overlap_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_overlap_affine_banded);

    // Semi-Global Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_semi_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_semi_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_semi_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_semi_global_linear_banded);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_semi_global_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_semi_global_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_semi_global_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_semi_global_affine_banded);

    SEQAN_CALL_TEST(test_align_global_alignment_banded_shorter_interfaces_linear);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_shorter_interfaces_affine);

    // Global Alignment with Differnt Container Types
    SEQAN_CALL_TEST(test_alignment_algorithms_global_banded_different_container);

    // Local Alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_local_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_local_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_local_linear_banded);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_local_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_graph_local_affine_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_local_affine_banded);

    // Dynamic Gaps.
    SEQAN_CALL_TEST(test_alignment_algorithms_global_dynamic_cost_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_local_dynamic_cost_banded);

    // Suboptimal Alignment.
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_align);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_gaps);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_fragment);

    // TODO(rmaerker): Here are the tests that should run when the Waterman-Eggert is adapted
//    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_suboptimal_linear_banded);
//    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_suboptimal_linear_banded);
//    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_suboptimal_linear_banded);
//    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_suboptimal_linear_banded);
//
//    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_suboptimal_affine_banded);
//    SEQAN_CALL_TEST(test_alignment_algorithms_gaps_gaps_suboptimal_affine_banded);
//    SEQAN_CALL_TEST(test_alignment_algorithms_graph_gaps_suboptimal_affine_banded);
//    SEQAN_CALL_TEST(test_alignment_algorithms_fragments_gaps_suboptimal_affine_banded);


    // ----------------------------------------------------------------------------
    // Test specialized alignments.
    // ----------------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_align);
    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_gaps);
    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_fragments);
    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_graph);

    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_align);
    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_gaps);
    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_fragments);
    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_graph);

    SEQAN_CALL_TEST(test_align_global_alignment_score_hirschberg);
    SEQAN_CALL_TEST(test_align_global_alignment_score_myers);
    SEQAN_CALL_TEST(test_align_global_alignment_score_myers_hirschberg);
    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_single_character);

    // -----------------------------------------------------------------------
    // Test Operations On Align Objects
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_integrate_align);
    SEQAN_CALL_TEST(test_align_integrate_align_infix_of_infix);

    // -----------------------------------------------------------------------
    // Test Printing.
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_stream_gaps_write);
    SEQAN_CALL_TEST(test_align_stream_gaps_stream);
    SEQAN_CALL_TEST(test_align_stream_align_write);
    SEQAN_CALL_TEST(test_align_stream_align_stream);

    // -----------------------------------------------------------------------
    // Test Alignment Evaluation
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_compute_alignment_stats);
}
SEQAN_END_TESTSUITE
