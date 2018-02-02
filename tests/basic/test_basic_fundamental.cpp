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

#include <seqan/basic.h>

#include "test_basic_fundamental_helpers.h"
#include "test_basic_fundamental_tags.h"
#include "test_basic_fundamental_metafunctions.h"
#include "test_basic_fundamental_transport.h"
#include "test_basic_fundamental_comparison.h"
#include "test_basic_fundamental_conversion.h"
#include "test_basic_array_construct_destruct.h"
#include "test_basic_hosted_type_interface.h"

SEQAN_BEGIN_TESTSUITE(test_basic_metaprogramming)
{
    // -----------------------------------------------------------------------
    // Fundamental Tags
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_fundamental_tags_tag);
    SEQAN_CALL_TEST(test_basic_fundamental_tags_tags);
    SEQAN_CALL_TEST(test_basic_fundamental_tags_tag_list);
    SEQAN_CALL_TEST(test_basic_fundamental_tags_tag_selector);
    SEQAN_CALL_TEST(test_basic_fundamental_tags_length_tag_list);

    // -----------------------------------------------------------------------
    // Fundamental Metafunctions
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_value);
    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_get_value);
    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_reference);
    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_size);
    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_difference);
    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_position);
    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_spec);
    SEQAN_CALL_TEST(test_basic_fundamental_metafunctions_deepest_spec);

    // -----------------------------------------------------------------------
    // Fundamental Transport Functions
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_fundamental_transport_assign);
    SEQAN_CALL_TEST(test_basic_fundamental_transport_move);
    SEQAN_CALL_TEST(test_basic_fundamental_transport_set);

    // -----------------------------------------------------------------------
    // Fundamental Comparison Code
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_fundamental_comparison_compare_type);

    // -----------------------------------------------------------------------
    // Fundamental Conversion Code
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_fundamental_convert_metafunction);
    SEQAN_CALL_TEST(test_basic_fundamental_convert_function);

    // -----------------------------------------------------------------------
    // Array Construction Destruction Code
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_array_construct_destruct_metafunction_is_simple);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_construct_value_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_destruct_value_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_construct_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_construct_copy_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_construct_move_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_destruct_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_fill_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_copy_forward_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_copy_backward_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_copy_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_move_forward_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_move_backward_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_move_pointer);
    SEQAN_CALL_TEST(test_basic_array_construct_destruct_array_clear_space_pointer);

    // -----------------------------------------------------------------------
    // Hosted Type Interface
    // -----------------------------------------------------------------------

    // TODO(holtgrew): Move to own sub module, test that then.
}
SEQAN_END_TESTSUITE

