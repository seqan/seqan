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

#ifndef TESTS_BAM_IO_TEST_BAM_TAGS_DICT_DICT_H_
#define TESTS_BAM_IO_TEST_BAM_TAGS_DICT_DICT_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_type_size)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(1, getBamTypeSize('A'));
    SEQAN_ASSERT_EQ(1, getBamTypeSize('c'));
    SEQAN_ASSERT_EQ(1, getBamTypeSize('C'));
    SEQAN_ASSERT_EQ(2, getBamTypeSize('s'));
    SEQAN_ASSERT_EQ(2, getBamTypeSize('S'));
    SEQAN_ASSERT_EQ(4, getBamTypeSize('i'));
    SEQAN_ASSERT_EQ(4, getBamTypeSize('I'));
    SEQAN_ASSERT_EQ(4, getBamTypeSize('f'));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_type)
{
    using namespace seqan;

    {
        CharString str("XXAa");
        BamTagsDict bamTags(str);
        SEQAN_ASSERT_EQ('A', getTagType(bamTags, 0));
    }
    {
        CharString str("XXI\xff\xff\xff\xff");
        BamTagsDict bamTags(str);
        SEQAN_ASSERT_EQ('I', getTagType(bamTags, 0));
    }
    {
        CharString str("XXAaXYI\xff\xff\xff\xff");
        BamTagsDict bamTags(str);
        SEQAN_ASSERT_EQ('A', getTagType(bamTags, 0));
        SEQAN_ASSERT_EQ(CharString("XX"), getTagKey(bamTags, 0));
        SEQAN_ASSERT_EQ('I', getTagType(bamTags, 1));
        SEQAN_ASSERT_EQ(CharString("XY"), getTagKey(bamTags, 1));
    }
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_length)
{
    using namespace seqan;

    // Empty string.
    {
        CharString str("");
        BamTagsDict bamTags(str);
        SEQAN_ASSERT_EQ(length(bamTags), 0u);
    }
    // One entry.
    {
        CharString str("XXAa");
        BamTagsDict bamTags(str);
        SEQAN_ASSERT_EQ(length(bamTags), 1u);
    }
    // Two entries.
    {
        CharString str("XXAaXXAa");
        BamTagsDict bamTags(str);
        SEQAN_ASSERT_EQ(length(bamTags), 2u);
    }
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_extract_value_type_A)
{
    using namespace seqan;

    CharString str("XXAa");
    BamTagsDict bamTags(str);
    char c;
    SEQAN_ASSERT(extractTagValue(c, bamTags, 0));
    SEQAN_ASSERT_EQ(c, 'a');
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_extract_value_type_c)
{
    using namespace seqan;

    CharString str("XXc\xff");
    BamTagsDict bamTags(str);
    uint8_t x;
    SEQAN_ASSERT(extractTagValue(x, bamTags, 0));
    SEQAN_ASSERT_EQ(x, 255u);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_extract_value_type_C)
{
    using namespace seqan;
    CharString str("XXC\xff");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(length(bamTags), 1u);
    int8_t x;
    SEQAN_ASSERT(extractTagValue(x, bamTags, 0));
    SEQAN_ASSERT_EQ(x, -1);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_extract_value_type_s)
{
    using namespace seqan;
    CharString str("XXs\xff\xff");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(length(bamTags), 1u);
    uint16_t x;
    SEQAN_ASSERT(extractTagValue(x, bamTags, 0));
    SEQAN_ASSERT_EQ(x, 0xffff);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_extract_value_type_S)
{
    using namespace seqan;
    CharString str("XXs\xff\xff");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(length(bamTags), 1u);
    int16_t x;
    SEQAN_ASSERT(extractTagValue(x, bamTags, 0));
    SEQAN_ASSERT_EQ(x, -1);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_extract_value_type_i)
{
    using namespace seqan;
    CharString str("XXi\xff\xff\xff\xff");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(length(bamTags), 1u);
    uint32_t x;
    SEQAN_ASSERT(extractTagValue(x, bamTags, 0));
    SEQAN_ASSERT_EQ(x, 0xffffffff);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_extract_value_type_I)
{
    using namespace seqan;
    CharString str("XXI\xff\xff\xff\xff");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(length(bamTags), 1u);
    int32_t x;
    SEQAN_ASSERT(extractTagValue(x, bamTags, 0));
    SEQAN_ASSERT_EQ(x, -1);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_A)
{
    using namespace seqan;

    CharString str("XXAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("Aa"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_c)
{
    using namespace seqan;

    CharString str("XXc\xffXAAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("c\xff"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_C)
{
    using namespace seqan;

    CharString str("XXC\xffXAAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("C\xff"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_s)
{
    using namespace seqan;

    CharString str("XXs\xff\xffXAAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("s\xff\xff"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_S)
{
    using namespace seqan;

    CharString str("XXS\xff\xffXAAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("S\xff\xff"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_i)
{
    using namespace seqan;

    CharString str("XXi\xff\xff\xff\xffXAAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("i\xff\xff\xff\xff"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_I)
{
    using namespace seqan;

    CharString str("XXI\xff\xff\xff\xffXAAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("I\xff\xff\xff\xff"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_f)
{
    using namespace seqan;

    CharString str("XXf\xff\xff\xff\xffXAAa");
    BamTagsDict bamTags(str);
    SEQAN_ASSERT_EQ(CharString("f\xff\xff\xff\xff"), suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_Z)
{
    using namespace seqan;

    CharString str("XXZthis is a test");
    appendValue(str, '\0');
    append(str, "XAAa");
    BamTagsDict bamTags(str);
    CharString expected("Zthis is a test");
    appendValue(expected, '\0');
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_H)
{
    using namespace seqan;

    CharString str("XXZFFFF");
    appendValue(str, '\0');
    append(str, "XAAa");
    BamTagsDict bamTags(str);
    CharString expected("ZFFFF");
    appendValue(expected, '\0');
    SEQAN_ASSERT_EQ(expected, suffix(bamTags[0], 2));
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_Bc)
{
    using namespace seqan;
    CharString str("XXBc\x02");
    appendValue(str, '\0');
    appendValue(str, '\0');
    appendValue(str, '\0');
    append(str, "\xff\xff");
    BamTagsDict bamTags(str);
    CharString expected("Bc\x02");
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    append(expected, "\xff\xff");
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_BC)
{
    using namespace seqan;
    CharString str("XXBC\x02");
    appendValue(str, '\0');
    appendValue(str, '\0');
    appendValue(str, '\0');
    append(str, "\xff\xff");
    BamTagsDict bamTags(str);
    CharString expected("BC\x02");
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    append(expected, "\xff\xff");
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_Bs)
{
    using namespace seqan;
    CharString str("XXBs\x02");
    appendValue(str, '\0');
    appendValue(str, '\0');
    appendValue(str, '\0');
    append(str, "\xff\xff\xff\xff");
    BamTagsDict bamTags(str);
    CharString expected("Bs\x02");
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    append(expected, "\xff\xff\xff\xff");
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_BS)
{
    using namespace seqan;
    CharString str("XXBS\x02");
    appendValue(str, '\0');
    appendValue(str, '\0');
    appendValue(str, '\0');
    append(str, "\xff\xff\xff\xff");
    BamTagsDict bamTags(str);
    CharString expected("BS\x02");
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    append(expected, "\xff\xff\xff\xff");
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_Bi)
{
    using namespace seqan;
    CharString str("XXBi\x02");
    appendValue(str, '\0');
    appendValue(str, '\0');
    appendValue(str, '\0');
    append(str, "\xff\xff\xff\xff\xff\xff\xff\xff");
    BamTagsDict bamTags(str);
    CharString expected("Bi\x02");
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    append(expected, "\xff\xff\xff\xff\xff\xff\xff\xff");
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_BI)
{
    using namespace seqan;
    CharString str("XXBI\x02");
    appendValue(str, '\0');
    appendValue(str, '\0');
    appendValue(str, '\0');
    append(str, "\xff\xff\xff\xff\xff\xff\xff\xff");
    BamTagsDict bamTags(str);
    CharString expected("BI\x02");
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    append(expected, "\xff\xff\xff\xff\xff\xff\xff\xff");
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_get_value_type_Bf)
{
    using namespace seqan;
    CharString str("XXBf\x02");
    appendValue(str, '\0');
    appendValue(str, '\0');
    appendValue(str, '\0');
    append(str, "\xff\xff\xff\xff\xff\xff\xff\xff");
    BamTagsDict bamTags(str);
    CharString expected("Bf\x02");
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    appendValue(expected, '\0');
    append(expected, "\xff\xff\xff\xff\xff\xff\xff\xff");
    CharString result = suffix(bamTags[0], 2);
    SEQAN_ASSERT_EQ(expected, result);
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_set_tag_value)
{
    using namespace seqan;

    // No tag.
    {
        CharString bamTags;
        CharString samTags = "";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XX", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:o"), CharString(samTags));
    }
    // Single tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XX", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:o"), CharString(samTags));
    }
    // Front tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XX", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:o\tXY:A:y\tXZ:A:z"), CharString(samTags));
    }
    // Last tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XZ", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:A:y\tXZ:A:o"), CharString(samTags));
    }
    // Center tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XY", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:A:o\tXZ:A:z"), CharString(samTags));
    }
    // Append.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XA", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:A:y\tXZ:A:z\tXA:A:o"), CharString(samTags));
    }

    // String value (replace).
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:Z:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XY", "Example string", 'Z');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:Z:Example string\tXZ:A:z"), CharString(samTags));
    }
    // String value (append).
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        setTagValue(tags, "XA", "Example string", 'Z');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:A:y\tXZ:A:z\tXA:Z:Example string"), CharString(samTags));
    }
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_append_tag_value)
{
    using namespace seqan;

    // No tag.
    {
        CharString bamTags;
        CharString samTags = "";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        appendTagValue(tags, "XX", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:o"), CharString(samTags));
    }
    // 2 tags.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        appendTagValue(tags, "XZ", 'o', 'A');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:A:y\tXZ:A:o"), CharString(samTags));
    }
    // String value (append).
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        appendTagValue(tags, "XA", "Example string", 'Z');
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:A:y\tXZ:A:z\tXA:Z:Example string"), CharString(samTags));
    }
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_erase_tag)
{
    using namespace seqan;

    // Single tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        eraseTag(tags, "XX");
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString(""), CharString(samTags));
    }
    // Front tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        eraseTag(tags, "XX");
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XY:A:y\tXZ:A:z"), CharString(samTags));
    }
    // Last tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        eraseTag(tags, "XZ");
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXY:A:y"), CharString(samTags));
    }
    // Center tag.
    {
        CharString bamTags;
        CharString samTags = "XX:A:x\tXY:A:y\tXZ:A:z";
        assignTagsSamToBam(bamTags, samTags);
        BamTagsDict tags(bamTags);
        eraseTag(tags, "XY");
        assignTagsBamToSam(samTags, bamTags);
        SEQAN_ASSERT_EQ(CharString("XX:A:x\tXZ:A:z"), CharString(samTags));
    }
}

SEQAN_DEFINE_TEST(test_bam_tags_dict_const_bam_tags_sequence)
{
    using namespace seqan;
    CharString bamStr, samStr = "AA:Z:value1\tAB:Z:value2\tAC:i:30";
    assignTagsSamToBam(bamStr, samStr);

    CharString const bamStrConst = bamStr;

    BamTagsDict tagsNonConst(bamStr);
    BamTagsDict tagsConst(bamStrConst);

    buildIndex(tagsNonConst);
    buildIndex(tagsConst);

    SEQAN_ASSERT(host(tagsConst) == host(tagsNonConst));
    appendTagValue(tagsConst, "AS", 10, 'i');
    SEQAN_ASSERT(bamStr != host(tagsConst));

    appendTagValue(tagsNonConst, "AS", 10, 'i');
    SEQAN_ASSERT(bamStr == host(tagsConst));

    unsigned tagId;
    SEQAN_ASSERT(findTagKey(tagId, tagsConst, "AS"));
    SEQAN_ASSERT(tagId == 3u);
    SEQAN_ASSERT(findTagKey(tagId, tagsNonConst, "AS"));
    SEQAN_ASSERT(tagId == 3u);
}

#endif  // TESTS_BAM_IO_TEST_BAM_TAGS_DICT_DICT_H_
