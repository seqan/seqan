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
// Tests for the sub module basic_fundamental.
// ==========================================================================

#include <seqan/basic/basic_debug.h>
#include <seqan/basic/basic_aggregate.h>
#include <seqan/stream.h>

#include "test_basic_aggregate.h"

SEQAN_BEGIN_TESTSUITE(test_basic_aggregate)
{
    // -----------------------------------------------------------------------
    // Tests for Pairs
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_set);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_pair_base_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_base_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_packed_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_pair_bit_packed_stream_output);

    // -----------------------------------------------------------------------
    // Tests for Triples
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_set);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_triple_base_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_base_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_comparison_same_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_comparison_different_spec);
    SEQAN_CALL_TEST(test_basic_aggregates_triple_packed_stream_output);

    // -----------------------------------------------------------------------
    // Tests for Tuples
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_set);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_shift_left);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_shift_right);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_clear);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_length);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_comparison_same_spec);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_comparison_different_spec);  // TODO(holtgrew): Could be added for completeness case, not supported right now.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_base_stream_output);

    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_metafunctions);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_constructors);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_assign);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_set);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_move);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_value);  // TODO(holtgrew): Need proxy for this.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_get_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_assign_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_set_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_move_value);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_shift_left);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_shift_right);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_clear);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_length);
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_comparison_same_spec);
    // SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_comparison_different_spec);  // TODO(holtgrew): Could be added for completeness case, not supported right now.
    SEQAN_CALL_TEST(test_basic_aggregates_tuple_bit_packed_stream_output);
}
SEQAN_END_TESTSUITE

