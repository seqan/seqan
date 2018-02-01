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

//#define DEBUG_JST_TRAVERSAL

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_journaled_string_tree_traverser.h"

SEQAN_BEGIN_TESTSUITE(test_journaled_string_tree_traverser)
{
    // Tests for journaled string tree
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_constructor);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_traverser);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_init);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_context_size);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_branch_size);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_at_end);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_is_base);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_advance);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_context_iterator);

    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_container);
    SEQAN_CALL_TEST(test_journaled_string_tree_traverser_basic_traversal);
}
SEQAN_END_TESTSUITE
