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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for align_extend
// ==========================================================================

// work around obscure FreeBSD issue
#ifndef _GLIBCXX_USE_C99
#define _GLIBCXX_USE_C99 1
#define _GLIBCXX_USE_C99_UNDEF 1
#endif
#include <string>

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_blast.h"

SEQAN_BEGIN_TESTSUITE(test_blast)
{
    SEQAN_CALL_TEST(test_blast_scoring_scheme_conversion);
    SEQAN_CALL_TEST(test_blast_scoring_adapter);
    SEQAN_CALL_TEST(test_blast_blastmatch_stats_and_score);
    SEQAN_CALL_TEST(test_blast_blastmatch_bit_score_e_value);

    SEQAN_CALL_TEST(test_blast_write_tabular);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_header_generationblastplus);
    SEQAN_CALL_TEST(test_blast_write_tabular_customfields);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_header_customfields);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_header_customfields_generationblastplus);
    SEQAN_CALL_TEST(test_blast_write_pairwise);

}
SEQAN_END_TESTSUITE

#ifdef _GLIBCXX_USE_C99_UNDEF
#undefine _GLIBCXX_USE_C99_UNDEF
#undefine _GLIBCXX_USE_C99
#endif