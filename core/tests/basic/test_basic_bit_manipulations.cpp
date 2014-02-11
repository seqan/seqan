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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Test suite for basic bit manipulations.
// ==========================================================================

#include <seqan/basic/basic_debug.h>
#include <seqan/basic/basic_bit_operations.h>

#include "test_basic_bit_manipulations.h"

SEQAN_BEGIN_TESTSUITE(test_basic_bit_manipulations)
{
    SEQAN_CALL_TEST(test_basic_bit_manipulations_bitwise_and);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_bitwise_or);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_bitwise_and_not);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_bitwise_not);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_test_all_zeros);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_set_all_zeros);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_test_all_ones);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_set_all_ones);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_bit_scan_forward);
    SEQAN_CALL_TEST(test_basic_bit_manipulations_bit_scan_reverse);
}
SEQAN_END_TESTSUITE
