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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Tests for the blast module
// ==========================================================================

#include <seqan/basic.h>

#ifndef COMPILER_MSVC

#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/blast.h>

#include "test_blast_misc.h"
#include "test_blast_statistics.h"
#include "test_blast_output.h"
#include "test_blast_input.h"

SEQAN_BEGIN_TESTSUITE(test_blast)
{
    // MISC
    SEQAN_CALL_TEST(test_blast_program);
    SEQAN_CALL_TEST(test_blast_context_targs);

    // STATISTICS
    SEQAN_CALL_TEST(test_blast_scoring_scheme_conversion);
    SEQAN_CALL_TEST(test_blast_scoring_scheme);
    SEQAN_CALL_TEST(test_blast_blastmatch_stats_and_score);
    SEQAN_CALL_TEST(test_blast_blastmatch_bit_score_e_value);

    // WRITING (lowlevel tag)
    SEQAN_CALL_TEST(test_blast_write_lowlevel);
    // WRITING (tabular tag)
    SEQAN_CALL_TEST(test_blast_write_tabular_without_comments);
    SEQAN_CALL_TEST(test_blast_write_tabular_without_comments_customfields);
    SEQAN_CALL_TEST(test_blast_write_tabular_without_comments_legacy);
    SEQAN_CALL_TEST(test_blast_write_tabular_without_comments_constexpr);
    SEQAN_CALL_TEST(test_blast_write_tabular_without_comments_customfields_constexpr);
    SEQAN_CALL_TEST(test_blast_write_tabular_without_comments_legacy_constexpr);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_comments);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_comments_customfields);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_comments_legacy);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_comments_constexpr);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_comments_customfields_constexpr);
    SEQAN_CALL_TEST(test_blast_write_tabular_with_comments_legacy_constexpr);
    // WRITING (report tag)
    SEQAN_CALL_TEST(test_blast_write_report);
    SEQAN_CALL_TEST(test_blast_write_report_constexpr);
    SEQAN_CALL_TEST(test_blast_write_report_constexpr_dynmatrix);

    // READING (lowlevel tag)
    SEQAN_CALL_TEST(test_blast_read_lowlevel);
    // READING (tabular tag)
    SEQAN_CALL_TEST(test_blast_read_tabular_without_comments);
    SEQAN_CALL_TEST(test_blast_read_tabular_without_comments_customfields);
    SEQAN_CALL_TEST(test_blast_read_tabular_without_comments_legacy);
    SEQAN_CALL_TEST(test_blast_read_tabular_without_comments_constexpr);
    SEQAN_CALL_TEST(test_blast_read_tabular_without_comments_customfields_constexpr);
    SEQAN_CALL_TEST(test_blast_read_tabular_without_comments_legacy_constexpr);
    SEQAN_CALL_TEST(test_blast_read_tabular_with_comments);
    SEQAN_CALL_TEST(test_blast_read_tabular_with_comments_customfields);
    SEQAN_CALL_TEST(test_blast_read_tabular_with_comments_legacy);
    SEQAN_CALL_TEST(test_blast_read_tabular_with_comments_constexpr);
    SEQAN_CALL_TEST(test_blast_read_tabular_with_comments_customfields_constexpr);
    SEQAN_CALL_TEST(test_blast_read_tabular_with_comments_legacy_constexpr);
}
SEQAN_END_TESTSUITE

#else
#pragma message("Due to a bug in Microsoft Visual Studio 2015 the BLAST module is deactivated.")
int main() {}
#endif
