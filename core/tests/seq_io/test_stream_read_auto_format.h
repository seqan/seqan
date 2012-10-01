// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Tests for reading with automatic file format detection.
// ==========================================================================


#ifndef TEST_STREAM_TEST_STREAM_READ_AUTO_FORMAT_H_
#define TEST_STREAM_TEST_STREAM_READ_AUTO_FORMAT_H_

SEQAN_DEFINE_TEST(test_stream_read_record_auto_format_quals_fasta)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = ">id1\nAAACCC\n>id2\nGGGTTT";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    CharString meta;
    CharString seq;
    CharString qual;
    AutoSeqStreamFormat formatTag;
    int res = readRecord(meta, seq, qual, reader, formatTag);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "id1");
    SEQAN_ASSERT_EQ(seq, "AAACCC");
    SEQAN_ASSERT_EQ(length(qual), 0u);
    SEQAN_ASSERT_EQ(formatTag.tagId, 1);
}

SEQAN_DEFINE_TEST(test_stream_read_record_auto_format_quals_fastq)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = "@id1\nAAACCC\n+\n!!!!!!\n@id2\nGGGTTT+\nIIIIII\n";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    CharString meta;
    CharString seq;
    CharString qual;
    AutoSeqStreamFormat formatTag;
    int res = readRecord(meta, seq, qual, reader, formatTag);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "id1");
    SEQAN_ASSERT_EQ(seq, "AAACCC");
    SEQAN_ASSERT_EQ(qual, "!!!!!!");
    SEQAN_ASSERT_EQ(formatTag.tagId, 2);
}

SEQAN_DEFINE_TEST(test_stream_read_record_auto_format_no_quals_fasta)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = ">id1\nAAACCC\n>id2\nGGGTTT";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    CharString meta;
    CharString seq;
    AutoSeqStreamFormat formatTag;
    int res = readRecord(meta, seq, reader, formatTag);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "id1");
    SEQAN_ASSERT_EQ(seq, "AAACCC");
    SEQAN_ASSERT_EQ(formatTag.tagId, 1);
}

SEQAN_DEFINE_TEST(test_stream_read_record_auto_format_no_quals_fastq)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = "@id1\nAAACCC\n+\n!!!!!!\n@id2\nGGGTTT+\nIIIIII\n";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    CharString meta;
    CharString seq;
    AutoSeqStreamFormat formatTag;
    int res = readRecord(meta, seq, reader, formatTag);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "id1");
    SEQAN_ASSERT_EQ(seq, "AAACCC");
    SEQAN_ASSERT_EQ(formatTag.tagId, 2);
}

SEQAN_DEFINE_TEST(test_stream_read_auto_format_quals_fasta)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = ">id1\nAAACCC\n>id2\nGGGTTT";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    StringSet<CharString> metas;
    StringSet<CharString> seqs;
    StringSet<CharString> quals;
    AutoSeqStreamFormat formatTag;

    int res = read2(metas, seqs, quals, reader, formatTag);

    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(length(metas), 2u);
    SEQAN_ASSERT_EQ(metas[0], "id1");
    SEQAN_ASSERT_EQ(metas[1], "id2");
    SEQAN_ASSERT_EQ(length(seqs), 2u);
    SEQAN_ASSERT_EQ(seqs[0], "AAACCC");
    SEQAN_ASSERT_EQ(seqs[1], "GGGTTT");
    SEQAN_ASSERT_EQ(length(quals), 0u);
}

SEQAN_DEFINE_TEST(test_stream_read_auto_format_quals_fastq)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = "@id1\nAAACCC\n+\n!!!!!!\n@id2\nGGGTTT+\nIIIIII\n";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    StringSet<CharString> metas;
    StringSet<CharString> seqs;
    StringSet<CharString> quals;
    AutoSeqStreamFormat formatTag;

    int res = read2(metas, seqs, quals, reader, formatTag);

    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(length(metas), 2u);
    SEQAN_ASSERT_EQ(metas[0], "id1");
    SEQAN_ASSERT_EQ(metas[1], "id2");
    SEQAN_ASSERT_EQ(length(seqs), 2u);
    SEQAN_ASSERT_EQ(seqs[0], "AAACCC");
    SEQAN_ASSERT_EQ(seqs[1], "GGGTTT");
    SEQAN_ASSERT_EQ(length(quals), 2u);
    SEQAN_ASSERT_EQ(quals[0], "!!!!!!");
    SEQAN_ASSERT_EQ(quals[1], "IIIIII");
}

SEQAN_DEFINE_TEST(test_stream_read_auto_format_no_quals_fasta)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = ">id1\nAAACCC\n>id2\nGGGTTT";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    StringSet<CharString> metas;
    StringSet<CharString> seqs;
    AutoSeqStreamFormat formatTag;

    int res = read2(metas, seqs, reader, formatTag);

    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(length(metas), 2u);
    SEQAN_ASSERT_EQ(metas[0], "id1");
    SEQAN_ASSERT_EQ(metas[1], "id2");
    SEQAN_ASSERT_EQ(length(seqs), 2u);
    SEQAN_ASSERT_EQ(seqs[0], "AAACCC");
    SEQAN_ASSERT_EQ(seqs[1], "GGGTTT");
}

SEQAN_DEFINE_TEST(test_stream_read_auto_format_no_quals_fastq)
{
    using namespace seqan;

    typedef Stream<CharArray<char const *> > TStream;

    char const * fastaString = "@id1\nAAACCC\n+\n!!!!!!\n@id2\nGGGTTT+\nIIIIII\n";
    TStream stream(fastaString, fastaString + strlen(fastaString));
    RecordReader<TStream, SinglePass<> > reader(stream);

    StringSet<CharString> metas;
    StringSet<CharString> seqs;
    AutoSeqStreamFormat formatTag;

    int res = read2(metas, seqs, reader, formatTag);

    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(length(metas), 2u);
    SEQAN_ASSERT_EQ(metas[0], "id1");
    SEQAN_ASSERT_EQ(metas[1], "id2");
    SEQAN_ASSERT_EQ(length(seqs), 2u);
    SEQAN_ASSERT_EQ(seqs[0], "AAACCC");
    SEQAN_ASSERT_EQ(seqs[1], "GGGTTT");
}

#endif  // TEST_STREAM_TEST_STREAM_READ_AUTO_FORMAT_H_
