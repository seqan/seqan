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

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_align_stream.h"
#include "test_align_gaps.h"
#include "test_align_gaps_iterator.h"
// TODO(holtgrew): Test Align<>, AlignCols<>!
#include "test_align_global_alignment.h"
#include "test_align_global_alignment_banded.h"
#include "test_align_global_alignment_specialized.h"
#include "test_align_global_alignment_score.h"
#include "test_align_local_alignment.h"
#include "test_align_alignment_operations.h"

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

    // -----------------------------------------------------------------------
    // Test Unbanded Global Alignment Algorithms
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_notop_left_noright_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_notop_left_noright_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_notop_left_noright_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_notop_left_noright_bottom_nw);

    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_notop_left_noright_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_notop_left_noright_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_notop_left_noright_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_notop_left_noright_bottom_gotoh);

    SEQAN_CALL_TEST(test_align_global_alignment_shorter_interfaces_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_shorter_interfaces_gotoh);

    // -----------------------------------------------------------------------
    // Test Banded Global Alignment Algorithms
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_notop_left_noright_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_notop_left_noright_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_notop_left_noright_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_top_left_right_bottom_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_notop_left_noright_bottom_nw);

    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_notop_left_noright_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_notop_left_noright_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_notop_left_noright_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_top_left_right_bottom_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_notop_left_noright_bottom_gotoh);

    SEQAN_CALL_TEST(test_align_global_alignment_banded_shorter_interfaces_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_banded_shorter_interfaces_gotoh);

    // -----------------------------------------------------------------------
    // Test Specialized Global Alignment Algorithms
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_align);
    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_gaps);
    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_fragments);
    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_graph);

    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_align);
    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_gaps);
    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_fragments);
    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_graph);

    // -----------------------------------------------------------------------
    // Test Global Alignment Score Functions
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_global_alignment_score_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_score_banded_nw);
    SEQAN_CALL_TEST(test_align_global_alignment_score_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_score_banded_gotoh);
    SEQAN_CALL_TEST(test_align_global_alignment_score_hirschberg);
    SEQAN_CALL_TEST(test_align_global_alignment_score_myers);
    SEQAN_CALL_TEST(test_align_global_alignment_score_myers_hirschberg);

    // -----------------------------------------------------------------------
    // Test Unbanded And Banded Local Alignment Algorithms
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_local_alignment_align);
    SEQAN_CALL_TEST(test_align_local_alignment_gaps);
    SEQAN_CALL_TEST(test_align_local_alignment_graph);
    SEQAN_CALL_TEST(test_align_local_alignment_fragment);

    SEQAN_CALL_TEST(test_align_local_alignment_banded_align);
    SEQAN_CALL_TEST(test_align_local_alignment_banded_gaps);
    SEQAN_CALL_TEST(test_align_local_alignment_banded_graph);
    SEQAN_CALL_TEST(test_align_local_alignment_banded_fragment);

    // -----------------------------------------------------------------------
    // Test Unbanded And Banded Local Alignment Enumeration Algorithms
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_align);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_gaps);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_fragment);

    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_align);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_gaps);
    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_fragment);

    // -----------------------------------------------------------------------
    // Test Operations On Align Objects
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_integrate_align);

    // -----------------------------------------------------------------------
    // Test Printing.
    // -----------------------------------------------------------------------
    
    SEQAN_CALL_TEST(test_align_stream_gaps_write);
    SEQAN_CALL_TEST(test_align_stream_gaps_stream);
    SEQAN_CALL_TEST(test_align_stream_align_write);
    SEQAN_CALL_TEST(test_align_stream_align_stream);
}
SEQAN_END_TESTSUITE
