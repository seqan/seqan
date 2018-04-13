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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for seqan/stream/virtual_stream.h
// ==========================================================================

#ifndef TEST_STREAM_TEST_VIRTUAL_STREAM_H_
#define TEST_STREAM_TEST_VIRTUAL_STREAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <sstream>

using namespace seqan;

template <typename TStream>
class VStreamTest : public Test
{
public:
    typedef TStream Type;
};

SEQAN_TYPED_TEST_CASE(VStreamTest, CompressedFileTypes);

static const char *FASTA_EXAMPLE = ">seq1\n\
CGATCGATAAT\n\
>seq2\n\
CCTCTCTCTCCCT\n\
>seq3\n\
CCCCCCCC\n";

static const char *FASTQ_EXAMPLE = "@seq1\n\
CGATCGATAAT\n\
+\n\
IIIIIIIIIII\n\
@seq2\n\
CCTCTCTCTCCCT\n\
+\n\
IIIIIIIIIIIII\n\
@seq3\n\
CCCCCCCC\n\
+\n\
IIIIIIII\n";

SEQAN_TYPED_TEST(VStreamTest, Construct)
{
    CharString fileName = getAbsolutePath("/tests/seq_io/test_dna.fq");
    append(fileName, FileExtensions<typename TestFixture::Type>::VALUE[0]);
    VirtualStream<char, Input> vstream(toCString(fileName), OPEN_RDONLY);
    SEQAN_ASSERT((bool)vstream);

    append(fileName, "_X_");
    VirtualStream<char, Input> vstream2(toCString(fileName), OPEN_RDONLY);
    SEQAN_ASSERT_NOT((bool)vstream2);
}

SEQAN_TYPED_TEST(VStreamTest, OpenClose)
{
    CharString fileName = getAbsolutePath("/tests/seq_io/test_dna.fq");
    append(fileName, FileExtensions<typename TestFixture::Type>::VALUE[0]);
    VirtualStream<char, Input> vstream;

    SEQAN_ASSERT(open(vstream, toCString(fileName), OPEN_RDONLY));
    SEQAN_ASSERT((bool)vstream);
    SEQAN_ASSERT(close(vstream));

    SEQAN_ASSERT(open(vstream, toCString(fileName), OPEN_RDONLY));
    SEQAN_ASSERT((bool)vstream);
    SEQAN_ASSERT(close(vstream));

    append(fileName, "_X_");
    SEQAN_ASSERT_NOT(open(vstream, toCString(fileName), OPEN_RDONLY));
    SEQAN_ASSERT_NOT((bool)vstream);
    SEQAN_ASSERT(close(vstream));
}

SEQAN_TYPED_TEST(VStreamTest, Decompression)
{
    typedef typename TestFixture::Type TCompressionTag;
    CharString fileName = getAbsolutePath("/tests/seq_io/test_dna.fq");
    append(fileName, FileExtensions<TCompressionTag>::VALUE[0]);
    VirtualStream<char, Input> vstream(toCString(fileName), OPEN_RDONLY);
    SEQAN_ASSERT((bool)vstream);

    std::stringstream sstr;
    sstr << vstream.streamBuf;
    SEQAN_ASSERT_EQ(CharString(sstr.str()), CharString(FASTQ_EXAMPLE));
    close(vstream);
    SEQAN_ASSERT_NOT((bool)vstream);

    fileName = getAbsolutePath("/tests/seq_io/test_dna.fa");
    append(fileName, FileExtensions<TCompressionTag>::VALUE[0]);
    open(vstream, toCString(toCString(fileName)), OPEN_RDONLY);

    sstr.str("");
    sstr << vstream.streamBuf;
    SEQAN_ASSERT_EQ(CharString(sstr.str()), CharString(FASTA_EXAMPLE));
    close(vstream);
}

SEQAN_TYPED_TEST(VStreamTest, Compression)
{
    CharString buffer;
    for (unsigned i = 0; i != 10000; ++i)
    {
        appendNumber(buffer,i);
        append(buffer, FASTQ_EXAMPLE);
    }

    typedef typename TestFixture::Type TCompressionTag;
    CharString fileName = SEQAN_TEMP_FILENAME();
    append(fileName, FileExtensions<TCompressionTag>::VALUE[0]);
    VirtualStream<char, Output> vostream(toCString(fileName), OPEN_WRONLY);
    SEQAN_ASSERT((bool)vostream);

    vostream << buffer;
    close(vostream);
    SEQAN_ASSERT_NOT((bool)vostream);

    VirtualStream<char, Input> vistream(toCString(fileName), OPEN_RDONLY);
    SEQAN_ASSERT((bool)vistream);
    std::stringstream sstr;
    sstr << vistream.streamBuf;
    SEQAN_ASSERT_EQ(CharString(sstr.str()), buffer);
    close(vistream);
}

SEQAN_TYPED_TEST(VStreamTest, AutoDetection)
{
    typedef typename TestFixture::Type TCompressionTag;
    CharString fileName = getAbsolutePath("/tests/seq_io/test_dna.fq");
    append(fileName, FileExtensions<TCompressionTag>::VALUE[0]);
    std::fstream file(toCString(fileName), std::ios::in | std::ios::binary);

    VirtualStream<char, Input> vstream(file);
    SEQAN_ASSERT((bool)vstream);

    std::stringstream sstr;
    sstr << vstream.streamBuf;
    SEQAN_ASSERT_EQ(CharString(sstr.str()), CharString(FASTQ_EXAMPLE));
    close(vstream);
    SEQAN_ASSERT_NOT((bool)vstream);
}

#endif // ndef TEST_STREAM_TEST_VIRTUAL_STREAM_H_
