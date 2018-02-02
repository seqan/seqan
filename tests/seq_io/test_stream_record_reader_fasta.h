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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// FASTA-Tests for seqan/stream/read_fasta_fastq.h
// ==========================================================================


#ifndef TEST_STREAM_TEST_STREAM_READ_FASTA_H_
#define TEST_STREAM_TEST_STREAM_READ_FASTA_H_

template <typename TRecordReader>
void FASTA_TEST(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fasta()));

    CharString meta;
    Dna5String seq;

    int res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seq, "AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGCCCCAGAGTGTCATGCATGTCGAACGTGTTTTTGGGGCGGTTATATATATATATATT");
    SEQAN_ASSERT_EQ(value(reader), '>');

    res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seq, "ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTA_TEST_PROTEIN(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fasta()));

    CharString meta;
    Peptide seq;

    int res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seq, "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG");
    SEQAN_ASSERT_EQ(value(reader), '>');

    res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seq, "GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTA_TEST_ANNOTATED_PROTEIN(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fasta()));

    CharString meta;
    CharString seq;

    int res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seq, "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIP[annotation+-*/]YIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG");
    SEQAN_ASSERT_EQ(value(reader), '>');

    res = readRecord(meta, seq, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seq, "GLMPFLHTSKHRSMMLRPLSQALFW[CTD(5)->gogogo]TLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTA_TEST_BATCH(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fasta()));

    StringSet<CharString> metas;
    StringSet<Dna5String> seqs;

    int res = read2(metas, seqs, reader, Fasta());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(metas[0], " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seqs[0], "AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGCCCCAGAGTGTCATGCATGTCGAACGTGTTTTTGGGGCGGTTATATATATATATATT");

    SEQAN_ASSERT_EQ(metas[1], "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seqs[1], "ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTA_TEST_BATCH_CONCAT(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fasta()));

    StringSet<CharString, Owner<ConcatDirect<> > > metas;
    StringSet<Dna5String, Owner<ConcatDirect<> > > seqs;

    int res = read2(metas, seqs, reader, Fasta());

    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(metas[0], " sequenceID_with special chars an irregular linebreaks");
    SEQAN_ASSERT_EQ(seqs[0], "AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGCCCCAGAGTGTCATGCATGTCGAACGTGTTTTTGGGGCGGTTATATATATATATATT");

    SEQAN_ASSERT_EQ(metas[1], "sequence2... with no linebreaks and no newline at end");
    SEQAN_ASSERT_EQ(seqs[1], "ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG");
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_double_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST_BATCH(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_double_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST_BATCH(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_single_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<StringReader> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_double_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<StringReader> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_single_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<StringReader> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST_BATCH(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_single_concat_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<StringReader> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST_BATCH_CONCAT(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_double_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<StringReader> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST_BATCH(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_batch_double_concat_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<StringReader> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTA_TEST_BATCH_CONCAT(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_protein_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFileProtein(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST_PROTEIN(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fasta_annotated_protein_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastAFileAnnotatedProtein(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTA_TEST_ANNOTATED_PROTEIN(reader);

    file->close();
}

#endif // def TEST_STREAM_TEST_STREAM_READ_FASTA_H_
