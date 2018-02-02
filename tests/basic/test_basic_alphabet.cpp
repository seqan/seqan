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
#include <seqan/basic/basic_metaprogramming.h>
#include <seqan/basic/basic_alphabet.h>

#include "test_basic_fundamental_helpers.h"

#include "test_basic_alphabet_concepts.h"
#include "test_basic_alphabet_math.h"
#include "test_basic_alphabet_adapt_builtins.h"
#include "test_basic_alphabet_qualities.h"
#include "test_basic_alphabet_bio.h"
#include "test_basic_alphabet_storage.h"
#include "test_basic_alphabet_residue.h"
#include "test_basic_alphabet_profile.h"

SEQAN_BEGIN_TESTSUITE(test_basic_alphabet)
{
    // -----------------------------------------------------------------------
    // Test Math Functions
    // -----------------------------------------------------------------------

    // SEQAN_CALL_TEST(test_basic_alphabet_math_metafunctions); //deprecated
    // SEQAN_CALL_TEST(test_basic_alphabet_math_min_value);     //deprecated
    // SEQAN_CALL_TEST(test_basic_alphabet_math_max_value);     //deprecated

    // -----------------------------------------------------------------------
    // Test Adaptions of Builtin Types
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_metafunction_is_char_type);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_bool);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_char);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_short);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_int);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_long);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_int8);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_uint8);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_int16);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_uint16);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_int32);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_uint32);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_int64);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_uint64);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_float);
    SEQAN_CALL_TEST(test_basic_alphabet_adapt_builtins_concepts_double);

    // -----------------------------------------------------------------------
    // Test Alphabet With Qualities
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_alphabet_qualities_quality_value_size_metafunction);
    SEQAN_CALL_TEST(test_basic_alphabet_qualities_quality_has_qualities_metafunction);
    SEQAN_CALL_TEST(test_basic_alphabet_qualities_convert_quality);

    // -----------------------------------------------------------------------
    // Test Bio Alphabet Features
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_alphabet_bio_gap_value_function);
    SEQAN_CALL_TEST(test_basic_alphabet_bio_unknown_value_function);

    // -----------------------------------------------------------------------
    // Test Alphabet Storage Code
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_alphabet_storage_bits_per_value_metafunction);
    SEQAN_CALL_TEST(test_basic_alphabet_storage_value_size_metafunction);
    SEQAN_CALL_TEST(test_basic_alphabet_storage_value_size_function);
    SEQAN_CALL_TEST(test_basic_alphabet_storage_integral_for_value_metafunction);
    SEQAN_CALL_TEST(test_basic_alphabet_storage_bytes_per_value_metafunction);

    // -----------------------------------------------------------------------
    // Test SimpleType
    // -----------------------------------------------------------------------

    // We test this implicitely through the residue types.

    // -----------------------------------------------------------------------
    // Test Residues (SimpleType Specializations)
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_dna);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_dna5);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_dna_q);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_dna5_q);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_rna);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_rna5);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_iupac);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_metafunctions_amino_acid);

    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_dna);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_dna5);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_dna_q);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_dna5_q);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_rna);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_rna5);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_iupac);
    SEQAN_CALL_TEST(test_basic_alphabet_residue_usage_amino_acid);

    // -----------------------------------------------------------------------
    // Test ProfileChar
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_alphabet_profile_metafunctions);
    SEQAN_CALL_TEST(test_basic_alphabet_profile_constructors);
    SEQAN_CALL_TEST(test_basic_alphabet_profile_relations);
    SEQAN_CALL_TEST(test_basic_alphabet_profile_empty);
}
SEQAN_END_TESTSUITE
