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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Test suite testing delta map.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_delta_store.h"
#include "test_delta_map.h"

SEQAN_BEGIN_TESTSUITE(test_delta_map)
{
    // Tests for delta store.
    SEQAN_CALL_TEST(test_delta_map_delta_store_is_delta_type);
    SEQAN_CALL_TEST(test_delta_map_delta_store_select_delta_type);
    SEQAN_CALL_TEST(test_delta_map_delta_store_get_delta_store);
    SEQAN_CALL_TEST(test_delta_map_delta_store_add_delta_value);
    SEQAN_CALL_TEST(test_delta_map_delta_store_erase_delta_value);
    SEQAN_CALL_TEST(test_delta_map_delta_store_delta_value);
    SEQAN_CALL_TEST(test_delta_map_delta_store_clear);
    SEQAN_CALL_TEST(test_delta_map_delta_store_deletion_size);
    SEQAN_CALL_TEST(test_delta_map_delta_store_insertion_size);
    SEQAN_CALL_TEST(test_delta_map_delta_store_net_size);

    // Tests for delta map.
    SEQAN_CALL_TEST(test_delta_map_insert);
    SEQAN_CALL_TEST(test_delta_map_erase);
    SEQAN_CALL_TEST(test_delta_map_lower_bound);
    SEQAN_CALL_TEST(test_delta_map_upper_bound);
    SEQAN_CALL_TEST(test_delta_map_count);
    SEQAN_CALL_TEST(test_delta_map_equal_range);
    SEQAN_CALL_TEST(test_delta_map_find);
    SEQAN_CALL_TEST(test_delta_map_size);
    SEQAN_CALL_TEST(test_delta_map_empty);

    SEQAN_CALL_TEST(test_delta_map_iterator);
    SEQAN_CALL_TEST(test_delta_map_iterator_copy_constructor);
    SEQAN_CALL_TEST(test_delta_map_iterator_assign);
    SEQAN_CALL_TEST(test_delta_map_iterator_value);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_value);

    // Tests for delta map entry.
    SEQAN_CALL_TEST(test_delta_map_entry_delta_coverage);
    SEQAN_CALL_TEST(test_delta_map_entry_delta_type);
    SEQAN_CALL_TEST(test_delta_map_entry_delta_position);
}
SEQAN_END_TESTSUITE
