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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Stream Benchmark Demo.  This demo shows how to use different stream types
// for reading and writing FASTA files.  The time required for this is
// printed, and thus this demo can be used as a benchmark tool.
// ==========================================================================

// #define SEQAN_NEW_IO

#include <cstdio>
#include <fstream>
#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include <bzlib.h>
#endif  // #if SEQAN_HAS_BZIP2

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

// Setting buffer size to 4MB, such that the overhead for jumping buffers and
// such only occurs every 4M chars.
const unsigned int BUFFER_SIZE = 1024 * 1024 * 4;

typedef Dna5String TSequence;

struct Options
{
    bool doublePass;
    bool cstdio;
    bool fstream;
    bool mmapString;
    bool gzip;
    bool bzip2;
    bool documentMMap;
    bool nonConcat;
    bool multiSeq;

    Options()
            : doublePass(false), cstdio(false), fstream(false),
              mmapString(false), gzip(false), bzip2(false), documentMMap(false),
              nonConcat(false), multiSeq(false)
    {}
};

int readFileMultiSeqFile(char const * filename, Options const & /*options*/)
{
    typedef StringSet<CharString> TSequenceIds;
    typedef StringSet<TSequence> TSequences;
    typedef Iterator<TSequenceIds>::Type TSequenceIdsIter;
    typedef Iterator<TSequences>::Type TSequencesIter;
    TSequenceIds sequenceIds;
    TSequences sequences;

    double before = sysTime();
    std::cerr << "READING\tWHOLE\tOWNER";
    std::cerr << "\tmultiseq" << std::flush;

    MultiSeqFile multiSeqFile;
    if (!open(multiSeqFile.concat, filename, OPEN_RDONLY))
    {
        std::cerr << std::endl << "Could not open mmap file for reading." << std::endl;
        return 1;
    }

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

    unsigned seqCount = length(multiSeqFile);
    StringSet<TSequence> seqs;
    StringSet<CharString> seqIDs;
    reserve(seqs, seqCount, Exact());
    reserve(seqIDs, seqCount, Exact());

    TSequence seq;
    // CharString qual;
    CharString id;
    for (unsigned i = 0; i < seqCount; ++i)
    {
        assignSeq(seq, multiSeqFile[i], format);    // read sequence
        // assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignSeqId(id, multiSeqFile[i], format);   // read sequence id

        // convert ascii to values from 0..62
        // store dna and quality together in Dna5Q
        // for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
        //     assignQualityValue(seq[j], (int)(ordValue(qual[j]) - 33));
        // we use reserve and append, as assign is not supported
        // by StringSet<..., Owner<ConcatDirect<> > >
        appendValue(seqs, seq, Generous());
        appendValue(seqIDs, id, Generous());
    }

    double after = sysTime();
    fprintf(stderr, "\t%f\n", after - before);
    return 0;
}

template <typename TSpec>
int readFileMMapDocument(char const * filename, Options const & /*options*/, TSpec const & /*tag*/)
{
    typedef StringSet<CharString, TSpec> TSequenceIds;
    typedef StringSet<TSequence, TSpec> TSequences;
    typedef typename Iterator<TSequenceIds>::Type TSequenceIdsIter;
    typedef typename Iterator<TSequences>::Type TSequencesIter;
    TSequenceIds sequenceIds;
    TSequences sequences;

    double before = sysTime();
    std::cerr << "READING\tWHOLE\t";
    if (IsSameType<TSpec, Owner<> >::VALUE)
        std::cerr << "OWNER";
    else
        std::cerr << "CONCAT";
    std::cerr << "\tmmap" << std::flush;
    // fprintf(stderr, "\t%f\n", after - before);
    typedef File<Async<> > TFile;
    typedef String<char, MMap< ExternalConfig<TFile> > > TMMapString;
    TMMapString myString;
    if (!open(myString, filename, OPEN_RDONLY)) {
        std::cerr << std::endl << "Could not open mmap file for reading." << std::endl;
        return 1;
    }
    RecordReader<TMMapString, DoublePass<StringReader> > reader(myString, BUFFER_SIZE);
    int res = read2(sequenceIds, sequences, reader, Fasta());
    SEQAN_ASSERT_EQ(length(sequenceIds), length(sequences));

    // TSequenceIdsIter itId = begin(sequenceIds);
    // TSequencesIter itSeq = begin(sequences);
    // for (; !atEnd(itId); ++itId, ++itSeq) {
        //std::cout << value(itId) << "\t" << value(itSeq) << "\n";
    // }

    if (res != 0)
        std::cerr << std::endl << "There was an error reading the FASTA file." << std::endl;
    double after = sysTime();
    fprintf(stderr, "\t%f\n", after - before);
    return 0;
}

int readFileMMapDocument(char const * filename, Options const & options)
{
    if (options.nonConcat)
        return readFileMMapDocument(filename, options, Owner<Default>());
    else
        return readFileMMapDocument(filename, options, Owner<ConcatDirect<> >());
}

template <typename TFile, typename TPass>
int readFastaFile(StringSet<CharString> & sequenceIds,
                  StringSet<TSequence> & sequences,
                  TFile & file,
                  TPass const & /*tag*/)
{
    (void)sequenceIds;
    (void)sequences;

    RecordReader<TFile, TPass> reader(file, BUFFER_SIZE);
//     Pair<CharString, TSequence> record;
    CharString meta;
    TSequence seq;
    // std::cerr << "READING FILE" << std::endl;
    while (!atEnd(reader)) {
        int res = readRecord(meta, seq, reader, Fasta());
        if (res != 0) {
            std::cerr << std::endl << "Error reading a record" << std::endl;
            return res;
        } else {
            //std::cout << record.i1 << "\t" << record.i2 << "\n";
        }
    }
    return 0;
}

template <typename TPass>
int readFileMMap(char const * filename, Options const & /*options*/, TPass const & /*tag*/)
{
    StringSet<CharString> sequenceIds;
    StringSet<TSequence> sequences;
//     StringSet<CharString> sequences;

    double before = sysTime();
    std::cerr << "READING\tRECORD\t" << std::flush;
    if (IsSameType<TPass, SinglePass<StringReader> >::VALUE)
        std::cerr << "SINGLE PASS\tmmap" << std::flush;
    else
        std::cerr << "DOUBLE PASS\tmmap" << std::flush;
    // fprintf(stderr, "\t%f\n", after - before);
    typedef File<Async<> > TFile;
    String<char, MMap<ExternalConfig<TFile> > > myString;
    if (!open(myString, filename, OPEN_RDONLY)) {
        std::cerr << std::endl << "Could not open mmap file for reading." << std::endl;
        return 1;
    }
    int res = readFastaFile(sequenceIds, sequences, myString, TPass());
    if (res != 0)
        std::cerr << std::endl << "There was an error reading the FASTA file." << std::endl;
    double after = sysTime();
    fprintf(stderr, "\t%f\n", after - before);
    return 0;
}

template <typename TPass>
int readFileDefault(char const * filename, Options const & options, TPass const & /*tag*/)
{
    StringSet<CharString> sequenceIds;
    StringSet<TSequence> sequences;

    std::cerr << "READING\tRECORD\t" << std::flush;
    if (IsSameType<TPass, SinglePass<> >::VALUE)
        std::cerr << "SINGLE PASS\t";
    else
        std::cerr << "DOUBLE PASS\t";
    double before = sysTime();
    if (options.cstdio) {
        std::cerr << "cstdio" << std::flush;
        FILE * f = fopen(filename, "rb");
        if (!f) {
            std::cerr << std::endl << "ERROR: Could not open input file!" << std::endl;
            return 1;
        }
        int res = readFastaFile(sequenceIds, sequences, f, TPass());
        if (res != 0)
            std::cerr << "There was an error reading the FASTA file." << std::endl;
        fclose(f);
    } else if (options.fstream) {
        std::cerr << "fstream" << std::flush;
        std::ifstream f(filename, std::ios_base::in | std::ios_base::binary);
        if (!f.is_open()) {
            std::cerr << std::endl << "ERROR: Could not open input file!" << std::endl;
            return 1;
        }
        int res = readFastaFile(sequenceIds, sequences, f, TPass());
        if (res != 0)
            std::cerr << "There was an error reading the FASTA file." << std::endl;
#if SEQAN_HAS_ZLIB
    } else if (options.gzip) {
        std::cerr << "gzip" << std::flush;
        Stream<GZFile> f;
        if (!open(f, filename, "r"))
        {
            std::cerr << "Could not open input file!" << std::endl;
            return 1;
        }
        int res = readFastaFile(sequenceIds, sequences, f, TPass());
        if (res != 0)
            std::cerr << "There was an error reading the FASTA file." << std::endl;
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
    } else if (options.bzip2) {
        std::cerr << "bzip2" << std::flush;
        Stream<BZ2File> f;
        if (!open(f, filename, "r"))
        {
            std::cerr << "Could not open input file!" << std::endl;
            return 1;
        }
        int res = readFastaFile(sequenceIds, sequences, f, TPass());
        if (res != 0)
            std::cerr << "There was an error reading the FASTA file." << std::endl;
#endif  // #if SEQAN_HAS_BZIP2
    } else {
        SEQAN_ASSERT_FAIL("SHOULD NEVER REACH HERE!");
    }
    double after = sysTime();
    fprintf(stderr, "\t%f\n", after - before);
    return 0;
}

int main(int argc, char const ** argv)
{
    Options options;

    // -----------------------------------------------------------------------
    // Setup Command Line Parser
    // -----------------------------------------------------------------------

    CommandLineParser parser;
    addUsageLine(parser, "benchmark_stream [OPTIONS] INPUT OUTPUT");
    addSection(parser, "Read Variant");
    addOption(parser, CommandLineOption("d", "double-pass", "Use double-pass parsing.", OptionType::Bool));
    addSection(parser, "Stream Type Options");
    addOption(parser, CommandLineOption("c", "cstdio", "Use <cstdio> stream.", OptionType::Bool));
    addOption(parser, CommandLineOption("f", "fstream", "Use <fstream> stream.", OptionType::Bool));
#if SEQAN_HAS_ZLIB
    addOption(parser, CommandLineOption("g", "gzip", "Use gzip stream.", OptionType::Bool));
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
    addOption(parser, CommandLineOption("b", "bzip2", "Use bzlib stream.", OptionType::Bool));
#endif  // #if SEQAN_HAS_BZIP2
    addOption(parser, CommandLineOption("m", "memory-mapped", "Use memory mapped I/O.", OptionType::Bool));
    addOption(parser, CommandLineOption("w", "document-mmapped", "Read whole document at once with memory mapped I/O.", OptionType::Bool));
    addOption(parser, CommandLineOption("n", "non-concat", "Do not use concat direct string for document-mmapped version.", OptionType::Bool));
    addOption(parser, CommandLineOption("s", "multi-seq", "Use MultiSeqFile to read input.", OptionType::Bool));
    requiredArguments(parser, 2);

    // -----------------------------------------------------------------------
    // Parse And Check Command Line Parameters
    // -----------------------------------------------------------------------

    bool stop = !parse(parser, argc, argv);
    if (stop) {
        if (isSetLong(parser, "help"))
            return 0;
        return 1;
    }

    if (isSetLong(parser, "double-pass"))
        options.doublePass = true;
    if (isSetLong(parser, "cstdio"))
        options.cstdio = true;
    if (isSetLong(parser, "fstream"))
        options.fstream = true;
#if SEQAN_HAS_ZLIB
    if (isSetLong(parser, "gzip"))
        options.gzip = true;
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
    if (isSetLong(parser, "bzip2"))
        options.bzip2 = true;
#endif  // #if SEQAN_HAS_BZIP2
    if (isSetLong(parser, "memory-mapped"))
        options.mmapString = true;
    if (isSetLong(parser, "document-mmapped"))
        options.documentMMap = true;
    if (isSetLong(parser, "non-concat"))
        options.nonConcat = true;
    if (isSetLong(parser, "multi-seq"))
        options.multiSeq = true;

    if (options.cstdio + options.fstream + options.gzip + options.bzip2 + options.mmapString +
        options.documentMMap + options.multiSeq == 0) {
        std::cerr << "You have to select exactly one stream type!" << std::endl;
        return 1;
    } else if (options.cstdio + options.fstream + options.gzip + options.bzip2 + options.mmapString + options.documentMMap > 1) {
        std::cerr << "Only one stream type can be selected!" << std::endl;
        return 1;
    }
    if (options.documentMMap && options.doublePass) {
        std::cerr << "Double-pass I/O is implicit with document-mapped." << std::endl;
        return 1;
    }

    // -----------------------------------------------------------------------
    // Read And Write FASTA file.
    // -----------------------------------------------------------------------
    if (options.multiSeq) {
        readFileMultiSeqFile(toCString(getArgumentValue(parser, 0)), options);
    } else if (options.documentMMap) {
        readFileMMapDocument(toCString(getArgumentValue(parser, 0)), options);
    } else if (options.doublePass) {
        if (options.mmapString)
            readFileMMap(toCString(getArgumentValue(parser, 0)), options, DoublePass<StringReader>());
        else
            readFileDefault(toCString(getArgumentValue(parser, 0)), options, DoublePass<>());
    } else {
        if (options.mmapString)
            readFileMMap(toCString(getArgumentValue(parser, 0)), options, SinglePass<StringReader>());
        else
            readFileDefault(toCString(getArgumentValue(parser, 0)), options, SinglePass<>());
    }

    // for (unsigned i = 0; i < length(sequenceIds); ++i) {
    //     std::cerr << '>' << sequenceIds[i] << '\n';
    //     std::cerr << sequences[i] << '\n';
    // }
    
    return 0;
}
