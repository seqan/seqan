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
// Basic tests for the BamHeader and BamHeaderRecord classes.
// ==========================================================================

#ifndef TESTS_BAM_IO_TEST_BAM_HEADER_RECORD_H_
#define TESTS_BAM_IO_TEST_BAM_HEADER_RECORD_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

SEQAN_DEFINE_TEST(test_bam_io_bam_header_class)
{
    using namespace seqan;

    // Nothing interesting, just instantiating.

    BamHeader header;
}

SEQAN_DEFINE_TEST(test_bam_io_bam_header_typedefs)
{
    using namespace seqan;

    typedef BamIOContext<>::TLengthStore TLengthStore SEQAN_UNUSED_TYPEDEF;
}

SEQAN_DEFINE_TEST(test_bam_io_bam_header_record_class)
{
    using namespace seqan;

    BamHeaderRecord record;
}

SEQAN_DEFINE_TEST(test_bam_io_bam_header_record_typedefs)
{
    using namespace seqan;

    typedef BamHeaderRecord::TTagName  TTagName  SEQAN_UNUSED_TYPEDEF;
    typedef BamHeaderRecord::TTagValue TTagValue SEQAN_UNUSED_TYPEDEF;
    typedef BamHeaderRecord::TTag      TTag      SEQAN_UNUSED_TYPEDEF;
    typedef BamHeaderRecord::TTags     TTags     SEQAN_UNUSED_TYPEDEF;
}

SEQAN_DEFINE_TEST(test_bam_io_bam_header_record_find_tag_key)
{
    using namespace seqan;

    typedef BamHeaderRecord::TTag      TTag;

    BamHeaderRecord record;
    appendValue(record.tags, TTag("SN", "aa"));
    appendValue(record.tags, TTag("BC", "bb"));
    appendValue(record.tags, TTag("ZX", "cc"));

    unsigned idx = 0;

    SEQAN_ASSERT(findTagKey(idx, "BC", record));
    SEQAN_ASSERT_EQ(idx, 1u);
    SEQAN_ASSERT(findTagKey(idx, "ZX", record));
    SEQAN_ASSERT_EQ(idx, 2u);
    SEQAN_ASSERT(findTagKey(idx, "SN", record));
    SEQAN_ASSERT_EQ(idx, 0u);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_header_record_get_tag_value)
{
    using namespace seqan;

    typedef BamHeaderRecord::TTag      TTag;

    BamHeaderRecord record;
    appendValue(record.tags, TTag("SN", "aa"));
    appendValue(record.tags, TTag("BC", "bb"));
    appendValue(record.tags, TTag("ZX", "cc"));

    CharString res;
    SEQAN_ASSERT(getTagValue(res, "BC", record));
    SEQAN_ASSERT_EQ(res, CharString("bb"));
    clear(res);
    SEQAN_ASSERT(getTagValue(res, 1, record));
    SEQAN_ASSERT_EQ(res, CharString("bb"));

    clear(res);
    SEQAN_ASSERT(getTagValue(res, "ZX", record));
    SEQAN_ASSERT_EQ(res, CharString("cc"));
    clear(res);
    SEQAN_ASSERT(getTagValue(res, 2, record));
    SEQAN_ASSERT_EQ(res, CharString("cc"));

    clear(res);
    SEQAN_ASSERT(getTagValue(res, "SN", record));
    SEQAN_ASSERT_EQ(res, CharString("aa"));
    clear(res);
    SEQAN_ASSERT(getTagValue(res, 0, record));
    SEQAN_ASSERT_EQ(res, CharString("aa"));
}

#endif  // TESTS_BAM_IO_TEST_BAM_HEADER_RECORD_H_
