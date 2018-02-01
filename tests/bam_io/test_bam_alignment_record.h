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
// Some basic tests for BamAlignmentRecord and related functions.
// ==========================================================================

#ifndef TESTS_BAM_IO_TEST_BAM_ALIGNMENT_RECORD_H_
#define TESTS_BAM_IO_TEST_BAM_ALIGNMENT_RECORD_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

// Simply instantiate class.
SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_class)
{
    using namespace seqan;

    // The BamAlignmentRecord class is not too interesting, we simply
    // instantiate it and call clear once.

    BamAlignmentRecord record;
    clear(record);
}

// Test hasFlag* functions.

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_multiple)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_MULTIPLE;
    SEQAN_ASSERT(hasFlagMultiple(record));

    record.flag = 0xffff ^ BAM_FLAG_MULTIPLE;
    SEQAN_ASSERT_NOT(hasFlagMultiple(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_unmapped)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_UNMAPPED;
    SEQAN_ASSERT(hasFlagUnmapped(record));

    record.flag = 0xffff ^ BAM_FLAG_UNMAPPED;
    SEQAN_ASSERT_NOT(hasFlagUnmapped(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_next_unmapped)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_NEXT_UNMAPPED;
    SEQAN_ASSERT(hasFlagNextUnmapped(record));

    record.flag = 0xffff ^ BAM_FLAG_NEXT_UNMAPPED;
    SEQAN_ASSERT_NOT(hasFlagNextUnmapped(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_rc)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_RC;
    SEQAN_ASSERT(hasFlagRC(record));

    record.flag = 0xffff ^ BAM_FLAG_RC;
    SEQAN_ASSERT_NOT(hasFlagRC(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_next_rc)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_NEXT_RC;
    SEQAN_ASSERT(hasFlagNextRC(record));

    record.flag = 0xffff ^ BAM_FLAG_NEXT_RC;
    SEQAN_ASSERT_NOT(hasFlagNextRC(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_first)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_FIRST;
    SEQAN_ASSERT(hasFlagFirst(record));

    record.flag = 0xffff ^ BAM_FLAG_FIRST;
    SEQAN_ASSERT_NOT(hasFlagFirst(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_last)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_LAST;
    SEQAN_ASSERT(hasFlagLast(record));

    record.flag = 0xffff ^ BAM_FLAG_LAST;
    SEQAN_ASSERT_NOT(hasFlagLast(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_secondary)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_SECONDARY;
    SEQAN_ASSERT(hasFlagSecondary(record));

    record.flag = 0xffff ^ BAM_FLAG_SECONDARY;
    SEQAN_ASSERT_NOT(hasFlagSecondary(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_qc_no_pass)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_QC_NO_PASS;
    SEQAN_ASSERT(hasFlagQCNoPass(record));

    record.flag = 0xffff ^ BAM_FLAG_QC_NO_PASS;
    SEQAN_ASSERT_NOT(hasFlagQCNoPass(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_duplicate)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_DUPLICATE;
    SEQAN_ASSERT(hasFlagDuplicate(record));

    record.flag = 0xffff ^ BAM_FLAG_DUPLICATE;
    SEQAN_ASSERT_NOT(hasFlagDuplicate(record));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_alignment_record_has_flag_supplementary)
{
    using namespace seqan;

    BamAlignmentRecord record;

    record.flag = BAM_FLAG_SUPPLEMENTARY;
    SEQAN_ASSERT(hasFlagSupplementary(record));

    record.flag = 0xffff ^ BAM_FLAG_SUPPLEMENTARY;
    SEQAN_ASSERT_NOT(hasFlagSupplementary(record));
}

#endif  // TESTS_BAM_IO_TEST_BAM_ALIGNMENT_RECORD_H_
