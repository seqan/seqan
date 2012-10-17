// ==========================================================================
//                               fm_index_beta
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_index_fm.h"
#include "test_index_fm_rank_support_bit_string.h"
#include "test_index_fm_rank_support_bit_string_iterator.h"
#include "test_index_fm_prefix_sum_table.h"
#include "test_index_fm_right_array_binary_tree.h"
#include "test_index_fm_right_array_binary_tree_iterator.h"
#include "test_index_fm_wavelet_tree.h"
#include "test_index_fm_sparse_string.h"
#include "test_index_fm_compressed_sa.h"
#include "test_index_fm_compressed_sa_iterator.h"
#include "test_index_fm_stree.h"


SEQAN_BEGIN_TESTSUITE(test_fm_index_beta)
{
    SEQAN_CALL_TEST(test_rsbs_defaultConstructor);
    SEQAN_CALL_TEST(test_rsbs_resize);
    SEQAN_CALL_TEST(test_rsbs_getBuPos);
    SEQAN_CALL_TEST(test_rsbs_getSBuPos);
    SEQAN_CALL_TEST(test_rsbs_getPosInBu);
    SEQAN_CALL_TEST(test_rsbs_getSetBit);
    SEQAN_CALL_TEST(test_rsbs_append_value);
    SEQAN_CALL_TEST(test_rsbs_rank);
    SEQAN_CALL_TEST(test_rsbs_update_ranks_);
    SEQAN_CALL_TEST(test_rsbs_constructor);
    SEQAN_CALL_TEST(test_rsbs_equalOperator);
    SEQAN_CALL_TEST(test_rsbs_assignOperator);
    SEQAN_CALL_TEST(test_rsbs_open_save);
    
    SEQAN_CALL_TEST(test_rsbs_iterator_get_value);

    SEQAN_CALL_TEST(prefix_sum_table_constructor);
    SEQAN_CALL_TEST(prefix_sum_table_get_alphabet_size);
    SEQAN_CALL_TEST(prefix_sum_table_get_character_position);
    SEQAN_CALL_TEST(prefix_sum_table_get_character);
    SEQAN_CALL_TEST(prefix_sum_table_determine_dollar_substitute);
    SEQAN_CALL_TEST(prefix_sum_table_get_pivot_position);
    SEQAN_CALL_TEST(prefix_sum_table_prefix_sum);
    SEQAN_CALL_TEST(prefix_sum_table_get_value);
    SEQAN_CALL_TEST(prefix_sum_table_length);
    SEQAN_CALL_TEST(prefix_sum_table_insert_dollar_);
    SEQAN_CALL_TEST(prefix_sum_table_resize);
    SEQAN_CALL_TEST(prefix_sum_table_set_prefix_sum);
    SEQAN_CALL_TEST(prefix_sum_table_value);
    SEQAN_CALL_TEST(prefix_sum_table_open_save);

    SEQAN_CALL_TEST(wavelet_tree_structure_constructor);
    SEQAN_CALL_TEST(wavelet_tree_structure_clear);
    SEQAN_CALL_TEST(wavelet_tree_structure_empty);
    SEQAN_CALL_TEST(wavelet_tree_structure_get_fibre);
    SEQAN_CALL_TEST(wavelet_tree_structure_length);
    SEQAN_CALL_TEST(wavelet_tree_structure_resize);
    SEQAN_CALL_TEST(wavelet_tree_structure_open_save);

    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_begin);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_container);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_end);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_get_character);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_get_child_pos);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_get_num_child_vertieces);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_get_position);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_go_child);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_go_down);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_go_right);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_go_to_position);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_go_up);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_is_leaf);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_is_root);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_and_go_right);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_character);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_child_vertieces_);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_left_child_pos_);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_position_);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_right_child_pos_);

    SEQAN_CALL_TEST(test_wavelet_tree_clear);
    SEQAN_CALL_TEST(test_wavelet_tree_constructor);
    SEQAN_CALL_TEST(test_wavelet_tree_dollar_position);
    SEQAN_CALL_TEST(test_wavelet_tree_dollar_substitute);
    SEQAN_CALL_TEST(test_wavelet_tree_empty);
    SEQAN_CALL_TEST(test_wavelet_tree_get_character);
    SEQAN_CALL_TEST(test_wavelet_tree_get_fibre);
    SEQAN_CALL_TEST(test_wavelet_tree_get_occ);
    SEQAN_CALL_TEST(test_wavelet_tree_fill_wavelet_tree_);
    SEQAN_CALL_TEST(test_wavelet_tree_open_save);

    SEQAN_CALL_TEST(sparse_string_get_value);
    SEQAN_CALL_TEST(sparse_string_clear_length_resize);
    SEQAN_CALL_TEST(sparse_string_empty);
    SEQAN_CALL_TEST(sparse_string_get_fibre);
    
    SEQAN_CALL_TEST(compressed_sa_clear_length_resize);
    SEQAN_CALL_TEST(compressed_sa_empty);
    SEQAN_CALL_TEST(compressed_sa_compressed_sa_create);
    SEQAN_CALL_TEST(compressed_sa_get_fibre);
    SEQAN_CALL_TEST(compressed_sa_get_next_pos_);
    SEQAN_CALL_TEST(compressed_sa_set_lf_table);
    SEQAN_CALL_TEST(compressed_sa_value_access);
    SEQAN_CALL_TEST(compressed_sa_open_save);
    
    SEQAN_CALL_TEST(compressed_sa_iterator_begin);
    SEQAN_CALL_TEST(compressed_sa_iterator_end);
    
    SEQAN_CALL_TEST(test_lf_table_lf_mapping);

    SEQAN_CALL_TEST(test_fm_index_constructor);
    SEQAN_CALL_TEST(test_fm_index_clear);
    SEQAN_CALL_TEST(test_fm_index_determine_dollar_substitute_);
    SEQAN_CALL_TEST(test_fm_index_empty);
    SEQAN_CALL_TEST(test_fm_index_find_first_index_);
    SEQAN_CALL_TEST(test_fm_index_get_fibre);
    SEQAN_CALL_TEST(test_fm_index_search);
    SEQAN_CALL_TEST(test_fm_index_open_save);

    SEQAN_CALL_TEST(fm_index_iterator_constuctor);
    SEQAN_CALL_TEST(fm_index_iterator_go_down);
    SEQAN_CALL_TEST(fm_index_iterator_is_leaf);
    SEQAN_CALL_TEST(fm_index_iterator_go_right);
    SEQAN_CALL_TEST(fm_index_iterator_go_up);
    SEQAN_CALL_TEST(fm_index_iterator_representative);
    SEQAN_CALL_TEST(fm_index_iterator_path_label);
    SEQAN_CALL_TEST(fm_index_iterator_is_root);
    SEQAN_CALL_TEST(fm_index_iterator_count_occurrences);
    SEQAN_CALL_TEST(fm_index_iterator_range);
}
SEQAN_END_TESTSUITE


