// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// This file coordinates the calls of tests to ensure that the sequence
// module fulfills the requirements.
// ==========================================================================
#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_sequence.h"
#include "test_segment_beta.h"
// #include "test_string_set.h"

SEQAN_BEGIN_TESTSUITE(test_sequence_concept)
{

    //unsigned a = 1;
    //int b = -1;

    //std::cerr << (a<b) << std::endl;

    // --------------------------------------------------------------------------
    // Testing Alloc Strings With Simple Types
    // --------------------------------------------------------------------------
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_copy_constructible);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_default_constructible);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_less);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_less_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_greater);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_greater_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_unequal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_append);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_append_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_assignable);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_assign_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_back);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_begin);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_begin_position);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_capacity);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_clear);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_end);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_end_position);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_erase);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_erase_back);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_front);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_get_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_insert);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_insert_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_iter);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_length);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_move_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_replace);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_reserve);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_resize);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_swap);
    SEQAN_CALL_TEST(test_sequence_alloc_string_dna_value);

    // --------------------------------------------------------------------------
    // Testing Alloc Strings With Non Simple Types
    // --------------------------------------------------------------------------
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_copy_constructible);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_default_constructible);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_less);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_less_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_greater);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_greater_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_unequal);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_append);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_append_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_assignable);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_assign_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_back);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_begin);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_begin_position);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_capacity);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_clear);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_end);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_end_position);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_erase);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_erase_back);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_front);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_get_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_insert);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_insert_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_iter);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_length);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_move_value);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_replace);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_reserve);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_swap);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_resize);
    SEQAN_CALL_TEST(test_sequence_alloc_string_counting_char_value);

    // --------------------------------------------------------------------------
    // Testing Alloc Segments With Simple Types
    // (Alloc Segments = Segments of Alloc Strings)
    // --------------------------------------------------------------------------
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_constructible);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_copy_constructible);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_default_constructible);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_less);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_less_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_greater);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_greater_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_equal);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_unequal);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_assignable);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_assign_value);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_back);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_begin);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_begin_position);
//     SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_clear);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_end);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_end_position);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_front);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_get_value);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_iter);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_length);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_move_value);
    SEQAN_CALL_TEST(test_sequence_alloc_segment_dna_value);

/*
 * Not tested: append(), appendValue(), capacity(), erase(), eraseBack(),
 * insert(), insertValue(), reserve().
 */

    // --------------------------------------------------------------------------
    // Testing Alloc Strings With Simple Types
    // --------------------------------------------------------------------------
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_copy_constructible);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_default_constructible);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_less_greater_equal);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_append);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_append_value);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_assign);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_assign_value);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_assign_value_by_id);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_begin);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_begin_position);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_back);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_clear);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_concat);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_end);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_end_position);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_erase);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_erase_back);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_front);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_get_value);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_get_value_by_id);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_infix);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_infix_with_length);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_insert);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_insert_value);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_iter);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_length);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_move_value);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_prefix);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_resize);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_suffix);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_swap);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_value);
    // SEQAN_CALL_TEST(test_sequence_alloc_string_set_dna_value_by_id);
}
SEQAN_END_TESTSUITE
