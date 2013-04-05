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
// Tests for the sequence_journaled module.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_sequence_journaled.h"
#include "test_sequence_journaled_iterator.h"


SEQAN_BEGIN_TESTSUITE(test_sequence_journaled) {
    // Call tests of the sequence journal with unbalanced tree journal.
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_assign);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_set);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_host);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_clear);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_erase_position);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_erase_begin_end);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_insert);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_insert_value);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_assign_value);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_subscript_operator);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_assign_infix);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_length);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_virtual_to_host_position);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_host_to_virtual_position);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_copy_constructor);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_begin_end_iterator);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_begin_end_const_iterator);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_subscript_operator_randomized);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_fuzzying);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_segments_read_only);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_segments_read_write);
   
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_iterator_sum);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_iterator_difference);
    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_iterator_relations);
    // TODO(rmaerker): Problems with implementation of end in unbalanced trees, however, unbalanced trees ought be deprecated right now.
//    SEQAN_CALL_TEST(test_sequence_journaled_unbalanced_tree_iterator_decrement);

    // Call tests of the sequence journal with sorted array journals.
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_assign);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_set);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_host);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_clear);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_erase_position);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_erase_begin_end);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_insert);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_insert_value);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_assign_value);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_subscript_operator);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_assign_infix);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_length);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_virtual_to_host_position);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_host_to_virtual_position);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_copy_constructor);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_begin_end_iterator);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_begin_end_const_iterator);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_subscript_operator_randomized);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_fuzzying);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_segments_read_only);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_segments_read_write);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_flatten);

    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_iterator_sum);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_iterator_difference);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_iterator_relations);
    SEQAN_CALL_TEST(test_sequence_journaled_sorted_array_iterator_decrement);
}
SEQAN_END_TESTSUITE
