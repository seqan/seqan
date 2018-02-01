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
#include <seqan/basic.h>

#include "test_basic_iterator.h"
#include "test_basic_iterator_zip.h"

SEQAN_BEGIN_TESTSUITE(test_basic_iterator)
{
    SEQAN_CALL_TEST(test_basic_iterator_adapt_pointer_metafunctions);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_pointer_transport);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_pointer_transport_value);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_pointer_movement);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_pointer_arithmetics);

    SEQAN_CALL_TEST(test_basic_iterator_adapt_std_iterator_metafunctions);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_std_iterator_constructors);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_std_iterator_transport);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_std_iterator_transport_value);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_std_iterator_movement);
    SEQAN_CALL_TEST(test_basic_iterator_adapt_std_iterator_arithmetics);

    SEQAN_CALL_TEST(test_basic_iterator_adaptor_metafunctions);
    SEQAN_CALL_TEST(test_basic_iterator_adaptor_constructors);
    SEQAN_CALL_TEST(test_basic_iterator_adaptor_transport);
    SEQAN_CALL_TEST(test_basic_iterator_adaptor_transport_value);
    SEQAN_CALL_TEST(test_basic_iterator_adaptor_movement);
    SEQAN_CALL_TEST(test_basic_iterator_adaptor_arithmetics);
    SEQAN_CALL_TEST(test_basic_iterator_adaptor_rooted_functions);

    SEQAN_CALL_TEST(test_basic_iterator_position_metafunctions);
    SEQAN_CALL_TEST(test_basic_iterator_position_constructors);
    SEQAN_CALL_TEST(test_basic_iterator_position_transport);
    SEQAN_CALL_TEST(test_basic_iterator_position_transport_value);
    SEQAN_CALL_TEST(test_basic_iterator_position_movement);
    SEQAN_CALL_TEST(test_basic_iterator_position_arithmetics);
    SEQAN_CALL_TEST(test_basic_iterator_position_rooted_functions);

    SEQAN_CALL_TEST(test_basic_iterator_zip_metafunctions);
    SEQAN_CALL_TEST(test_basic_iterator_zip_constructors);
    SEQAN_CALL_TEST(test_basic_iterator_zip_make_zip_iterator);
    SEQAN_CALL_TEST(test_basic_iterator_zip_transport);
    SEQAN_CALL_TEST(test_basic_iterator_zip_transport_value);
    SEQAN_CALL_TEST(test_basic_iterator_zip_movement);
    SEQAN_CALL_TEST(test_basic_iterator_zip_arithmetics);
}
SEQAN_END_TESTSUITE

