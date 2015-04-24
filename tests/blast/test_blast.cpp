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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
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

#include "test_blast_statistics.h"
#include "test_blast_output.h"
#include "test_blast_input.h"

SEQAN_BEGIN_TESTSUITE(test_blast)
{
    // STATISTICS
    SEQAN_CALL_TEST(test_blast_scoring_scheme_conversion);
    SEQAN_CALL_TEST(test_blast_scoring_adapter);
    SEQAN_CALL_TEST(test_blast_blastmatch_stats_and_score);
    SEQAN_CALL_TEST(test_blast_blastmatch_bit_score_e_value);

    //TODO extra test for different programs

    // WRITING (tests writeMatch() interfaces)
    SEQAN_CALL_TEST(test_blast_write_match_tabular);
    SEQAN_CALL_TEST(test_blast_write_match_tabular_run_time_context_args); // same as first but properties set as vars
    SEQAN_CALL_TEST(test_blast_write_match_tabular_legacy);
    SEQAN_CALL_TEST(test_blast_write_match_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_write_match_tabular_with_header_legacy);
    SEQAN_CALL_TEST(test_blast_write_match_customfields_tabular);
    SEQAN_CALL_TEST(test_blast_write_match_customfields_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_write_match_lowlevel_tabular);

    // WRITING (tests writeMatch() and writeRecordHeader() interfaces)
    SEQAN_CALL_TEST(test_blast_write_header_tabular);
    SEQAN_CALL_TEST(test_blast_write_header_tabular_legacy);
    SEQAN_CALL_TEST(test_blast_write_header_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_write_header_tabular_with_header_legacy);
    SEQAN_CALL_TEST(test_blast_write_header_customfields_tabular);
    SEQAN_CALL_TEST(test_blast_write_header_customfields_tabular_with_header);

    // WRITING (tests writeRecord() interfaces)
    SEQAN_CALL_TEST(test_blast_write_record_tabular);
    SEQAN_CALL_TEST(test_blast_write_record_tabular_legacy);
    SEQAN_CALL_TEST(test_blast_write_record_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_write_record_tabular_with_header_legacy);
    SEQAN_CALL_TEST(test_blast_write_record_customfields_tabular);
    SEQAN_CALL_TEST(test_blast_write_record_customfields_tabular_with_header);

    // WRITING (tests formattedFile's writeRecord() interface)
    SEQAN_CALL_TEST(test_blast_write_formatted_file_tabular);
    SEQAN_CALL_TEST(test_blast_write_formatted_file_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_write_formatted_file_customfields_tabular);
    SEQAN_CALL_TEST(test_blast_write_formatted_file_customfields_tabular_with_header);

    //TODO extra tests for the stuff in the context

    // WRITING (tests writeTop(), writeRecord() and writeBottom() for
    //          PAIRWISE -- this is the only test for pairwise)
    SEQAN_CALL_TEST(test_blast_write_pairwise);
    SEQAN_CALL_TEST(test_blast_write_pairwise_formatted_file);

    // READING (onMatch, readMatch, skipMatch, skipUntilMatch)
    SEQAN_CALL_TEST(test_blast_read_match_tabular);
    SEQAN_CALL_TEST(test_blast_read_match_tabular_legacy);
    SEQAN_CALL_TEST(test_blast_read_match_customfields_tabular);
    SEQAN_CALL_TEST(test_blast_read_match_lowlevel_tabular);
//     SEQAN_CALL_TEST(test_blast_read_match_lowlevel_tabular_legacy);

    SEQAN_CALL_TEST(test_blast_read_match_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_read_match_tabular_with_header_legacy);
    SEQAN_CALL_TEST(test_blast_read_match_customfields_tabular_with_header);
//     SEQAN_CALL_TEST(test_blast_read_match_lowlevel_tabular_with_header);
//     SEQAN_CALL_TEST(test_blast_read_match_lowlevel_tabular_with_header_legacy);

    // READING (readHeader, skipHeader)
    SEQAN_CALL_TEST(test_blast_read_header_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_read_header_tabular_with_header_legacy);
    SEQAN_CALL_TEST(test_blast_read_header_customfields_tabular_with_header);

    // READING (readRecord)
    SEQAN_CALL_TEST(test_blast_read_record_tabular);
    // when there is no header, you can't tell the difference for legacy, so no extra test
    SEQAN_CALL_TEST(test_blast_read_record_customfields_tabular);

    SEQAN_CALL_TEST(test_blast_read_record_tabular_with_header);
    SEQAN_CALL_TEST(test_blast_read_record_tabular_with_header_legacy);
    SEQAN_CALL_TEST(test_blast_read_record_customfields_tabular_with_header);

    // READING (formattedFile)
    SEQAN_CALL_TEST(test_blast_read_formatted_file_tabular);
    SEQAN_CALL_TEST(test_blast_read_formatted_file_tabular_with_header);
}
SEQAN_END_TESTSUITE

#ifdef _GLIBCXX_USE_C99_UNDEF
#undef _GLIBCXX_USE_C99_UNDEF
#undef _GLIBCXX_USE_C99
#endif
