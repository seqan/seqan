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

#ifndef TESTS_BAM_IO_TEST_BAM_SAM_CONVERSION_H_
#define TESTS_BAM_IO_TEST_BAM_SAM_CONVERSION_H_

// TODO(holtgrew): The code below is my feeble attempt to use SeqAn strings for binary stuff. Maybe pure C-style stuff is more appropriate?

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_two_tags)
{
    using namespace seqan;

    CharString bamTags = "XXAyXYAx";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:A:y\tXY:A:x"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_A)
{
    using namespace seqan;

    CharString bamTags = "XXAy";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:A:y"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_c)
{
    using namespace seqan;

    CharString bamTags = "XXc\xFF";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:i:-1"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_C)
{
    using namespace seqan;

    CharString bamTags = "XXC\xFF";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:i:255"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_s)
{
    using namespace seqan;

    CharString bamTags = "XXs\x2E\xFB";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:i:-1234"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_S)
{
    using namespace seqan;

    CharString bamTags = "XXS\x2E\xFB";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:i:64302"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_i)
{
    using namespace seqan;

    CharString bamTags = "XXi\xFF\xFE\x1D\xC0";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:i:-1071776001"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_I)
{
    using namespace seqan;

    CharString bamTags = "XXI\xFF\xFE\x1D\xC0";
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:i:3223191295"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_f)
{
#if defined (__arm__) && defined(__ARM_PCS_VFP) // NOTE(h-2): armhf CRASHES here for unknown reasons
    return;
#endif
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 7);
    char const * DATA = "XXf\x00\x00\x00\x3f";
    arrayCopy(DATA, DATA + 7, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86) // rounding errors on non-x86
    SEQAN_ASSERT_EQ(CharString("XX:f:0.5"), CharString(samTags));
#endif
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_Z)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 18);
    char const * DATA = "XXZThis is a test\0";
    arrayCopy(DATA, DATA + 18, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:Z:This is a test"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_H)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 12);
    char const * DATA = "XXHFF00FF00\0";
    arrayCopy(DATA, DATA + 12, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:H:FF00FF00"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_Bc)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 10);
    char const * DATA = "XXBc\2\0\0\0\xff\xff";
    arrayCopy(DATA, DATA + 10, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:B:c,-1,-1"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_BC)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 10);
    char const * DATA = "XXBC\2\0\0\0\xff\xff";
    arrayCopy(DATA, DATA + 10, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:B:C,255,255"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_Bs)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 12);
    char const * DATA = "XXBs\2\0\0\0\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 12, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:B:s,-1,-1"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_BS)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 12);
    char const * DATA = "XXBS\2\0\0\0\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 12, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:B:S,65535,65535"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_Bi)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 16);
    char const * DATA = "XXBi\2\0\0\0\xff\xff\xff\xff\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 16, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:B:i,-1,-1"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_BI)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 16);
    char const * DATA = "XXBI\2\0\0\0\xff\xff\xff\xff\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 16, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
    SEQAN_ASSERT_EQ(CharString("XX:B:I,4294967295,4294967295"), CharString(samTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_bam_to_sam_type_Bf)
{
    using namespace seqan;

    CharString bamTags;
    resize(bamTags, 16);
    char const * DATA = "XXBf\2\0\0\0\x00\x00\x00\x3f\x00\x00\x00\x3f";
    arrayCopy(DATA, DATA + 16, &bamTags[0]);
    CharString samTags;
    assignTagsBamToSam(samTags, bamTags);
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86) // rounding errors on non-x86
    SEQAN_ASSERT_EQ(CharString("XX:B:f,0.5,0.5"), CharString(samTags));
#endif
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_two_tags)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:A:y\tXY:A:x";
    assignTagsSamToBam(bamTags, samTags);
    SEQAN_ASSERT_EQ(CharString("XXAyXYAx"), CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_A)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:A:y";
    assignTagsSamToBam(bamTags, samTags);
    SEQAN_ASSERT_EQ(CharString("XXAy"), CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_i)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:i:-123456";
    assignTagsSamToBam(bamTags, samTags);
    SEQAN_ASSERT_EQ(CharString("XXi\xC0\x1D\xFE\xFF"), CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_f)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:f:0.5";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 7);
    char const * DATA = "XXf\x00\x00\x00\x3f";
    arrayCopy(DATA, DATA + 7, begin(expected, Standard()));
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86) // rounding errors on non-x86
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
#endif
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_Z)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:Z:This is a test";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 18);
    char const * DATA = "XXZThis is a test\0";
    arrayCopy(DATA, DATA + 18, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_H)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:H:FF00FF00";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 12);
    char const * DATA = "XXHFF00FF00\0";
    arrayCopy(DATA, DATA + 12, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_Bc)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:B:c,-1,-1";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 10);
    char const * DATA = "XXBc\2\0\0\0\xff\xff";
    arrayCopy(DATA, DATA + 10, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_BC)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:B:C,255,255";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 10);
    char const * DATA = "XXBC\2\0\0\0\xff\xff";
    arrayCopy(DATA, DATA + 10, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_Bs)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:B:s,-1,-1";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 12);
    char const * DATA = "XXBs\2\0\0\0\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 12, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_BS)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:B:S,65535,65535";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 12);
    char const * DATA = "XXBS\2\0\0\0\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 12, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_Bi)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:B:i,-1,-1";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 16);
    char const * DATA = "XXBi\2\0\0\0\xff\xff\xff\xff\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 16, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_BI)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:B:I,4294967295,4294967295";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 16);
    char const * DATA = "XXBI\2\0\0\0\xff\xff\xff\xff\xff\xff\xff\xff";
    arrayCopy(DATA, DATA + 16, begin(expected, Standard()));
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
}

SEQAN_DEFINE_TEST(test_assign_tags_sam_to_bam_type_Bf)
{
    using namespace seqan;

    CharString bamTags;
    CharString samTags = "XX:B:f,0.5,0.5";
    assignTagsSamToBam(bamTags, samTags);
    CharString expected;
    resize(expected, 16);
    char const * DATA = "XXBf\2\0\0\0\x00\x00\x00\x3f\x00\x00\x00\x3f";
    arrayCopy(DATA, DATA + 16, begin(expected, Standard()));
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86) // rounding errors on non-x86
    SEQAN_ASSERT_EQ(expected, CharString(bamTags));
#endif
}

#endif  // TESTS_BAM_IO_TEST_BAM_SAM_CONVERSION_H_
