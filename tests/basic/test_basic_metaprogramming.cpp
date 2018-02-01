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

#include <cstdio>
#include <seqan/basic/basic_debug.h>
#include <seqan/basic/basic_metaprogramming.h>

#include "test_basic_metaprogramming_logic.h"
#include "test_basic_metaprogramming_control.h"
#include "test_basic_metaprogramming_math.h"
#include "test_basic_metaprogramming_type.h"
#include "test_basic_metaprogramming_enable_if.h"

SEQAN_BEGIN_TESTSUITE(test_basic_metaprogramming)
{
    // -----------------------------------------------------------------------
    // Metaprogramming Logic
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_bool_type);
    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_eval);
    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_or);
    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_or_c);
    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_and);
    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_and_c);
    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_if);
    SEQAN_CALL_TEST(test_basic_metaprogramming_logic_if_c);

    // -----------------------------------------------------------------------
    // Metaprogramming Control Structures
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_metaprogramming_control_loop_reverse);
    SEQAN_CALL_TEST(test_basic_metaprogramming_control_loop);
    SEQAN_CALL_TEST(test_basic_metaprogramming_control_switch);
    SEQAN_CALL_TEST(test_basic_metaprogramming_control_if);

    // -----------------------------------------------------------------------
    // Metaprogramming Math
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_metaprogramming_math_log2);
    SEQAN_CALL_TEST(test_basic_metaprogramming_math_log2_floor);
    SEQAN_CALL_TEST(test_basic_metaprogramming_math_log2_power);

    // -----------------------------------------------------------------------
    // Metaprogramming Type Queries / Modification
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_metaprogramming_type_same_type);
    SEQAN_CALL_TEST(test_basic_metaprogramming_type_make_signed);
    SEQAN_CALL_TEST(test_basic_metaprogramming_type_make_unsigned);
    SEQAN_CALL_TEST(test_basic_metaprogramming_type_remove_reference);
    SEQAN_CALL_TEST(test_basic_metaprogramming_type_remove_const);
    SEQAN_CALL_TEST(test_basic_metaprogramming_type_is_const);
    SEQAN_CALL_TEST(test_basic_metaprogramming_type_class_identifier);

    // -----------------------------------------------------------------------
    // Metaprogramming Conditional Enabling
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_basic_metaprogramming_enable_if_disable_if);
    SEQAN_CALL_TEST(test_basic_metaprogramming_enable_if2_disable_if2);
}
SEQAN_END_TESTSUITE

