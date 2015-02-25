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

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_delta_map.h"
#include "test_journaled_string_tree.h"

SEQAN_BEGIN_TESTSUITE(test_journaled_string_tree)
{
    // Tests for delta map.

    SEQAN_CALL_TEST(test_delta_map_insert);
    SEQAN_CALL_TEST(test_delta_map_length);
    SEQAN_CALL_TEST(test_delta_map_empty);
    SEQAN_CALL_TEST(test_delta_map_coverage_size);
    SEQAN_CALL_TEST(test_delta_map_set_coverage_size);
    SEQAN_CALL_TEST(test_delta_map_iterator);
    SEQAN_CALL_TEST(test_delta_map_iterator_copy_constructor);
    SEQAN_CALL_TEST(test_delta_map_iterator_assign);
    SEQAN_CALL_TEST(test_delta_map_iterator_value);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_type);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_position);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_value);
    SEQAN_CALL_TEST(test_delta_map_iterator_delta_coverage);

    // Tests for journaled string tree
    SEQAN_CALL_TEST(test_journaled_string_tree_container_mf);
    SEQAN_CALL_TEST(test_journaled_string_tree_get_string_tree_mf);
    SEQAN_CALL_TEST(test_journaled_string_tree_host_mf);
    SEQAN_CALL_TEST(test_journaled_string_tree_constructor);
    SEQAN_CALL_TEST(test_journaled_string_tree_init);
    SEQAN_CALL_TEST(test_journaled_string_tree_reinit);
    SEQAN_CALL_TEST(test_journaled_string_tree_container);
    SEQAN_CALL_TEST(test_journaled_string_tree_string_set);
    SEQAN_CALL_TEST(test_journaled_string_tree_full_journal_required);
    SEQAN_CALL_TEST(test_journaled_string_tree_set_block_size);
    SEQAN_CALL_TEST(test_journaled_string_tree_block_size);
    SEQAN_CALL_TEST(test_journaled_string_tree_journal_next_block);
    SEQAN_CALL_TEST(test_journaled_string_tree_host);
    SEQAN_CALL_TEST(test_journaled_string_tree_local_to_global_pos);

//    SEQAN_CALL_TEST(test_journaled_string_tree_save_open);
//    SEQAN_CALL_TEST(test_journaled_string_tree_jst_traversal_concept);
}
SEQAN_END_TESTSUITE
