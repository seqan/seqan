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

#include <seqan/stream.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "test_modifier_alphabet.h"
#include "test_modifier_functors.h"
#include "test_modifier_shortcuts.h"
#include "test_modifier_string.h"
#include "test_modifier_view.h"
#include "test_modifier_string_padding.h"


SEQAN_BEGIN_TESTSUITE(test_modifier)
{
    // Test the modifier shortcuts.
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna_string_reverse);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna5_string_reverse);
    SEQAN_CALL_TEST(test_modifer_shortcuts_rna_string_reverse);
    SEQAN_CALL_TEST(test_modifer_shortcuts_rna5_string_reverse);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna_string_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna5_string_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_rna_string_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_rna5_string_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna_string_reverse_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_dna5_string_reverse_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_rna_string_reverse_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_rna5_string_reverse_complement);
    SEQAN_CALL_TEST(test_modifer_shortcuts_complement_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_complement_in_place_string_set);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_complement_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_complement_in_place_string_set);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_reverse_in_place_string_set);
    SEQAN_CALL_TEST(test_modifier_reverse_iterator_metafunctions);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_lower_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_lower_in_place_string_set);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_upper_in_place_string);
    SEQAN_CALL_TEST(test_modifer_shortcuts_to_upper_in_place_string_set);

    // Test modifier functors.
    SEQAN_CALL_TEST(test_modifier_functors_functor_upcase);
    SEQAN_CALL_TEST(test_modifier_functors_functor_lowcase);
    SEQAN_CALL_TEST(test_modifier_functors_dna_complement);

    // Test modified alphabet.
    SEQAN_CALL_TEST(test_modifier_alphabet_size_metafunctions);
    SEQAN_CALL_TEST(test_modifier_DnaQ);
    SEQAN_CALL_TEST(test_modifier_alphabet_enumerate);
    SEQAN_CALL_TEST(test_modifier_alphabet_convert);
    SEQAN_CALL_TEST(test_modifier_alphabet_ord_value);
    /*
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_eq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_neq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_lt);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_gt);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_leq);
    SEQAN_CALL_TEST(test_modifier_alphabet_operator_geq);
    */

    // Test modifier view string.
    SEQAN_CALL_TEST(test_modifier_view_iterator_metafunctions);
    SEQAN_CALL_TEST(test_modifier_view_iterator);
    // SEQAN_CALL_TEST(test_modifier_view_const_iterator);
    SEQAN_CALL_TEST(test_modifier_view_string_caesar_chiffre);
    SEQAN_CALL_TEST(test_modifier_view_string_upper_case);
    SEQAN_CALL_TEST(test_modifier_view_string_low_case);
    SEQAN_CALL_TEST(test_modifier_view_string_alphabet_conversion);
    SEQAN_CALL_TEST(test_modifier_view_string_nested_modifier);
    SEQAN_CALL_TEST(test_modifier_convert_in_place);

    // Test modified string and iterator in some configurations.
    SEQAN_CALL_TEST(test_modifier_modified_string_metafunctions);
    SEQAN_CALL_TEST(test_modifier_modified_string_construct);
    SEQAN_CALL_TEST(test_modifier_modified_string_assignment);
    SEQAN_CALL_TEST(test_modifier_modified_string_length);
    SEQAN_CALL_TEST(test_modifier_modified_string_cascade);

    // TODO(holtgrew): SEQAN_CALL_TEST(test_modifier_modified_iterator_metafunctions);
    SEQAN_CALL_TEST(test_modifier_modified_iterator_construct);

    SEQAN_CALL_TEST(test_modifier_modified_string_mod_view);
    SEQAN_CALL_TEST(test_modifier_modified_string_mod_view_segment);

    SEQAN_CALL_TEST(test_modifier_modified_string_mod_pos);

//    SEQAN_CALL_TEST(test_modifier_modified_string_literal);
    SEQAN_CALL_TEST(test_modifier_modified_string_const_literal);
    SEQAN_CALL_TEST(test_modifier_modified_string_reverse_segment);

    SEQAN_CALL_TEST(test_modifier_minimal);
    SEQAN_CALL_TEST(test_modifier_reverse_back_front);

    // Test modified string with padding.
    SEQAN_CALL_TEST(test_modified_string_padding_construction);
    SEQAN_CALL_TEST(test_modified_string_padding_expand);
    SEQAN_CALL_TEST(test_modified_string_padding_length);
    SEQAN_CALL_TEST(test_modified_string_padding_begin);
    SEQAN_CALL_TEST(test_modified_string_padding_end);
    SEQAN_CALL_TEST(test_modified_string_padding_difference);
    SEQAN_CALL_TEST(test_modified_string_padding_iterator);
    SEQAN_CALL_TEST(test_modified_string_padding_defect_2190);
}
SEQAN_END_TESTSUITE
