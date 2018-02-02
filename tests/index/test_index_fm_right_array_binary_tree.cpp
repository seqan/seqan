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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#include "test_index_helpers.h"
#include "test_index_fm_right_array_binary_tree.h"
#include "test_index_fm_right_array_binary_tree_iterator.h"

using namespace seqan;

SEQAN_BEGIN_TESTSUITE(test_index_fm_right_array_binary_tree)
{
    SEQAN_CALL_TEST(wavelet_tree_structure_constructor);
    SEQAN_CALL_TEST(wavelet_tree_structure_clear);
    SEQAN_CALL_TEST(wavelet_tree_structure_empty);
    SEQAN_CALL_TEST(wavelet_tree_structure_get_fibre);
    SEQAN_CALL_TEST(wavelet_tree_structure_length);
    SEQAN_CALL_TEST(wavelet_tree_structure_resize);
//    SEQAN_CALL_TEST(wavelet_tree_structure_open_save);

    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_begin);
//    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_container);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_end);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_get_character);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_get_child_pos);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_get_num_child_vertices);
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
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_child_vertices_);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_left_child_pos_);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_position_);
    SEQAN_CALL_TEST(wavelet_tree_structure_iterator_set_right_child_pos_);
}
SEQAN_END_TESTSUITE
