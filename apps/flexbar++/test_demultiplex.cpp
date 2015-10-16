// ==========================================================================
//                                SeqAn-Flexbar
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// Author: Benjamin Strauch <b.strauch@fu-berlin.de>
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// This file provides tests for the demultiplexing of seqan-flexbar which is
// based in the implementation of the original flexbar program in [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"
#include "general_processing.h"
#include "read.h"

using namespace seqan;

// Used for loading the sequences. Needed for the io-test.
template<template<typename> class TRead, typename TSeq>
int loadSeqs(char const * path, std::vector<TRead<TSeq>>& reads)
{
    SeqFileIn seqFile;
    unsigned records = 1000;									//number of records to be read

    if (!open(seqFile, path))
    {
        std::cerr << "Error while opening the sequence-file.\n";
        return 1;
    }

    unsigned counter = 0;
    TRead<TSeq> read;
    while (!atEnd(seqFile) && counter < records)
    {
        readRecord(read.id, read.seq, seqFile);
        reads.emplace_back(std::move(read));
        ++counter;
    }
    return 0;
}
int loadSeqs(char const * path, StringSet<String<char> >& ids, StringSet<String<Dna5Q> >& seqs)
{
	SeqFileIn seqFile;
	unsigned records = 1000;									//number of records to be read

	if (!open(seqFile, path))
	{
		std::cerr << "Error while opening the sequence-file.\n";
		return 1;
	}

    String<char> meta;
    String<Dna5Q> seq;
    unsigned counter = 0;
    while (!atEnd(seqFile) && counter < records)
    {
        readRecord(meta, seq, seqFile);
        appendValue(ids, meta);
        appendValue(seqs, seq);
        ++counter;
    }
	return 0;
}
// Used for loading the sequences. Needed for the io-test.
int loadBarcodes(char const * path, std::vector<std::string>& bcids, std::vector<std::string>& bcs)
{
	SeqFileIn bcFile;
	
	if (!open(bcFile, path, OPEN_RDONLY))
	{
		std::cerr << "Error while opening barcode-file.\n";
		return 1;
	}

    while (!atEnd(bcFile))
    {
        std::string bc;
        std::string id;
        readRecord(id, bc, bcFile);
        bcids.emplace_back(id);
        bcs.emplace_back(bc);
    }
	return 0;
}
// Checks the correctness of the check function which checks the size of the barcodes and reads.
SEQAN_DEFINE_TEST(check_test)
{
	GeneralStats generalStats;
    using TRead = Read<seqan::Dna5QString>;
    std::vector<TRead> reads(5);
    
    StringSet<String<Dna5Q> > seqs;
    reads[0].seq = "TACGTAGCTACTGACTGACT";
    reads[1].seq = "G";
    reads[2].seq = "TACGTCTGACT";
    reads[3].seq = "";
    reads[4].seq = "ATCGA";

    reads[0].id = "ErsteSeq";
    reads[1].id = "LoeschenEins";
    reads[2].id = "ZweiteSeq";
    reads[3].id = "LoeschenLeer";
    reads[4].id = "LoeschenFuenf";

	StringSet<String<Dna5Q> > barcodesFalse;
	appendValue(barcodesFalse, "ACGAGT");
	appendValue(barcodesFalse, "TGCATC");
	appendValue(barcodesFalse, "AGCTAAT");
	appendValue(barcodesFalse, "GTGACA");

	StringSet<String<Dna5Q> > barcodesTrue;
	appendValue(barcodesTrue, "ACGAGT");
	appendValue(barcodesTrue, "TGCATC");
	appendValue(barcodesTrue, "AGCTAA");
	appendValue(barcodesTrue, "GTGACA");

    std::vector<TRead> expectedReads(2);
    expectedReads[0].seq = "TACGTAGCTACTGACTGACT";
    expectedReads[1].seq = "TACGTCTGACT";

    expectedReads[0].id = "ErsteSeq";
    expectedReads[1].id = "ZweiteSeq";

	bool resF = check(reads, barcodesFalse, generalStats);
	bool resT = check(reads, barcodesTrue, generalStats);

	SEQAN_ASSERT_EQ(false, resF);
	SEQAN_ASSERT_EQ(true, resT);
	SEQAN_ASSERT_EQ(2u, length(reads));
	SEQAN_ASSERT_EQ(expectedReads[0].seq, reads[0].seq);
    SEQAN_ASSERT_EQ(expectedReads[0].id, reads[0].id);
    SEQAN_ASSERT_EQ(expectedReads[1].seq, reads[1].seq);
    SEQAN_ASSERT_EQ(expectedReads[1].id, reads[1].id);
}
// Checks the correctness of the getPrefix function which extract the prefices from a set of sequences.
SEQAN_DEFINE_TEST(getPrefix_test)
{
    using TRead = Read<seqan::Dna5QString>;
    std::vector<TRead> reads(3);
    reads[0].seq = "GATACAGACTGAGCATGTGATCGAC";
    reads[1].seq = "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    reads[2].seq = "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
	
    std::vector<TRead> expectedReads(3);
    expectedReads[0].seq = "GATACAGACTGAGCATGTGATCGAC";
    expectedReads[1].seq = "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    expectedReads[2].seq = "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
        
	std::vector<std::string> exspected;
	appendValue(exspected, "GATAC");
	appendValue(exspected, "AATTC");
	appendValue(exspected, "GTTGG");
		
	for (unsigned i = 0; i < length(exspected); ++i)
	{
        SEQAN_ASSERT_EQ(exspected[i], getPrefix(reads[i], 5));
	}
}
// Checks the correct production of all possible barcode variations with 1 Error on one single barcode
SEQAN_DEFINE_TEST(buildVariations_test)
{
	Dna5String barcode = "ACTNG";
	StringSet<Dna5String> exspected;
	appendValue(exspected, Dna5String("ACTNG"));
	appendValue(exspected, Dna5String("CCTNG"));
	appendValue(exspected, Dna5String("GCTNG"));
	appendValue(exspected, Dna5String("TCTNG"));
	appendValue(exspected, Dna5String("NCTNG"));
	appendValue(exspected, Dna5String("AATNG"));
	appendValue(exspected, Dna5String("ACTNG"));
	appendValue(exspected, Dna5String("AGTNG"));
	appendValue(exspected, Dna5String("ATTNG"));
	appendValue(exspected, Dna5String("ANTNG"));
	appendValue(exspected, Dna5String("ACANG"));
	appendValue(exspected, Dna5String("ACCNG"));
	appendValue(exspected, Dna5String("ACGNG"));
	appendValue(exspected, Dna5String("ACTNG"));
	appendValue(exspected, Dna5String("ACNNG"));
	appendValue(exspected, Dna5String("ACTAG"));
	appendValue(exspected, Dna5String("ACTCG"));
	appendValue(exspected, Dna5String("ACTGG"));
	appendValue(exspected, Dna5String("ACTTG"));
	appendValue(exspected, Dna5String("ACTNG"));
	appendValue(exspected, Dna5String("ACTNA"));
	appendValue(exspected, Dna5String("ACTNC"));
	appendValue(exspected, Dna5String("ACTNG"));
	appendValue(exspected, Dna5String("ACTNT"));
	appendValue(exspected, Dna5String("ACTNN"));
	StringSet<Dna5String> res;
    buildVariations(res, barcode);
	
	for (unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
}
// Checks the correct execution of buildAllVariations and the updating of the barcodes and barcode-IDs
SEQAN_DEFINE_TEST(buildAllVariations_test)
{
	Dna5String barcode1 = "AN";
	CharString id1 = "barcode1";
	Dna5String barcode2 = "CA";
	CharString id2 = "barcode2";
	StringSet<Dna5String> barcodes;
	StringSet<CharString> bcids;
	appendValue(barcodes, barcode1);
	appendValue(bcids, id1);
	appendValue(barcodes, barcode2);
	appendValue(bcids, id2);
	StringSet<Dna5String> exspected;
	appendValue(exspected, Dna5String("AN"));
	appendValue(exspected, Dna5String("CN"));
	appendValue(exspected, Dna5String("GN"));
	appendValue(exspected, Dna5String("TN"));
	appendValue(exspected, Dna5String("NN"));
	appendValue(exspected, Dna5String("AA"));
	appendValue(exspected, Dna5String("AC"));
	appendValue(exspected, Dna5String("AG"));
	appendValue(exspected, Dna5String("AT"));
	appendValue(exspected, Dna5String("AN"));
	appendValue(exspected, Dna5String("AA"));
	appendValue(exspected, Dna5String("CA"));
	appendValue(exspected, Dna5String("GA"));
	appendValue(exspected, Dna5String("TA"));
	appendValue(exspected, Dna5String("NA"));
	appendValue(exspected, Dna5String("CA"));
	appendValue(exspected, Dna5String("CC"));
	appendValue(exspected, Dna5String("CG"));
	appendValue(exspected, Dna5String("CT"));
	appendValue(exspected, Dna5String("CN"));

	buildAllVariations(barcodes); //, bcids
	for (unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], barcodes[i]); 
	}
}

// Checks the correctness of the findExactIndex function which searches for one piece of sequence in the barcodes. Implicitly checks the construction of the Index.
SEQAN_DEFINE_TEST(findExactIndex_test)
{
    using TRead = Read<seqan::Dna5QString>;
	std::vector<std::string> barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");

    std::vector<TRead> reads(8);
    reads[0].seq = "CCCCCC";
    reads[1].seq = "AAAAAA";
    reads[2].seq = "TTTTTT";
    reads[3].seq = "GGGGGG";
    reads[4].seq = "CCCNCC";
    reads[5].seq = "GATACA";
    reads[6].seq = "ACGTAC";
    reads[7].seq = "ATGACNAANG";  

    BarcodeMatcher BarcodeMatcher(barcodes);

	int exspected[] = {1,0,3,2,-1,-1,4,-1};

	for (unsigned i = 0; i < length(reads); ++i)
	{
		int res = BarcodeMatcher.getMatchIndex(reads[i]);
		SEQAN_ASSERT_EQ(exspected[i], res);
	}
}
// Checks the correctnes of the findAllExactIndex function which searches for many pieces of sequence in the barcodes. Implicitly checks the construction of the Index.
SEQAN_DEFINE_TEST(matchBarcodes_test) 
{
    using TRead = Read<seqan::Dna5QString>;
    std::vector<std::string> barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");
	
    std::vector<TRead> reads(7);
    reads[0].seq = "CCCCCC";
    reads[1].seq = "AAAAAA";
    reads[2].seq = "TTTTTT";
    reads[3].seq = "GGGGGG";
    reads[4].seq = "CCCNCC";
    reads[5].seq = "GATACA";
    reads[6].seq = "ACGTAC";

    BarcodeMatcher BarcodeMatcher(barcodes);

	int exspected[] = {1,0,3,2,-1,-1,4,};
	std::vector<int> res(7);
    const auto numMatched = std::count_if(reads.begin(), reads.end(), [](const auto& read)->auto{return read.demuxResult != 0;});
    MatchBarcodes(reads, BarcodeMatcher);
    for (unsigned i = 0; i < reads.size(); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], reads[i].demuxResult);
	}
	
}
// Checks the correctnes of the clipBarcodes function which erases the first x bases of a sequence.
SEQAN_DEFINE_TEST(clipBarcodes_test)
{
    using TRead = Read<seqan::Dna5QString>;
    std::vector<TRead> reads(2);
    reads[0].seq = "AAAAAAGTGACTGATCGTACGACTG";
    reads[1].seq = "GGGGGGGGGGGGGGGG";

    reads[0].demuxResult = 1;
    reads[1].demuxResult = 0;

	StringSet<String<Dna5Q> > exspected ;
	appendValue(exspected, "GTGACTGATCGTACGACTG");
	appendValue(exspected, "GGGGGGGGGGGGGGGG");

	clipBarcodes(reads, 6, ClipSoft());
	for (unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], reads[i].seq);
	}
	
}
// Checks the correctnes of the clipBarcodes function using the hardClip method and therefore also clipping sequences without matching barcode.
SEQAN_DEFINE_TEST(clipBarcodesStrict_test)
{
    using TRead = Read<seqan::Dna5QString>;
    std::vector<TRead> reads(2);
    reads[0].seq = "GTGACTGATCGTACGACTG";
    reads[1].seq = "GGGGGGGGGGGGGGGG";

    StringSet<String<Dna5Q> > exspected;
    appendValue(exspected, "GATCGTACGACTG");
    appendValue(exspected, "GGGGGGGGGG");

    clipBarcodes(reads, 6, ClipHard());
    for (unsigned i = 0; i < length(exspected); ++i)
    {
        SEQAN_ASSERT_EQ(exspected[i], reads[i].seq);
    }
}

// Checks the correctness of the demultiplex function performing all demultiplexing operations for exact inline barcode matching.
SEQAN_DEFINE_TEST(demultiplex_Exact_test)
{
    std::vector<Read<seqan::Dna5QString>> reads(4);
    reads[0].seq = "AAAAAAGTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG";
    reads[1].seq = "GGGGGGAGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC";
    reads[2].seq = "AAAAAATAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    reads[3].seq = "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC";

    reads[0].id = "Adenin1";
    reads[1].id = "Guanin";
    reads[2].id = "Adenin2";
    reads[3].id = "Unidentifiziert";

    std::vector<std::string> barcodes;
    appendValue(barcodes, "GGGGGG");
    appendValue(barcodes, "CCCCCC");
    appendValue(barcodes, "AAAAAA");

    auto expectedReads = reads;
    expectedReads[0].seq = "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG";
    expectedReads[1].seq = "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC";
    expectedReads[2].seq = "TAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    expectedReads[3].seq = "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC";

    expectedReads[0].demuxResult = 3;
    expectedReads[1].demuxResult = 1;
    expectedReads[2].demuxResult = 3;
    expectedReads[3].demuxResult = 0;

    //Index<StringSet<String<Dna5Q> >, IndexEsa<> > indexSet(barcodes);
    //Finder<Index<StringSet<String<Dna5Q> >, IndexEsa<> > > esaFinder(indexSet);
    //indexRequire(indexSet, FibreSA());
    BarcodeMatcher BarcodeMatcher(barcodes);

    GeneralStats stats(length(expectedReads),0);
    demultiplex(reads, BarcodeMatcher, false, stats, ExactBarcodeMatching(), false);

    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(expectedReads[i].demuxResult, reads[i].demuxResult);
        for (unsigned j = 0; j < length(expectedReads[i].seq); ++j)
        {
            SEQAN_ASSERT_EQ(expectedReads[i].seq[j], reads[i].seq[j]);
        }
    }
}

// Checks the correctness of the demultiplex function performing all demultiplexing operations
// for exact multiplex barcode matching.
SEQAN_DEFINE_TEST(demultiplex_Exact_Multiplex_test)
{
    std::vector<ReadMultiplex<seqan::Dna5QString>> reads(4);
    reads[0].seq = "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG";
    reads[1].seq = "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC";
    reads[2].seq = "TAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    reads[3].seq = "GGGG";

    reads[0].id = "Adenin1";
    reads[1].id = "Guanin";
    reads[2].id = "Adenin2";
    reads[3].id = "Unidentifiziert";

    reads[0].demultiplex = "AAAAAA";
    reads[1].demultiplex = "GGGGGG";
    reads[2].demultiplex = "GGCCGG";
    reads[3].demultiplex = "AAAAAA";

    std::vector<std::string> barcodes;
    appendValue(barcodes, "GGGGGG");
    appendValue(barcodes, "CCCCCC");
    appendValue(barcodes, "AAAAAA");

    auto expectedReads = reads;
    expectedReads[0].seq = "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG";
    expectedReads[1].seq = "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC";
    expectedReads[3].seq = "TAGCTAGCTAGCTAGCTAGCTAGCTAGC";  // matched reads will be partitioned at the top
    expectedReads[2].seq = "GGGG";

    expectedReads[0].demuxResult = 3;
    expectedReads[1].demuxResult = 1;
    expectedReads[3].demuxResult = 0;
    expectedReads[2].demuxResult = 3;

    BarcodeMatcher BarcodeMatcher(barcodes);

    GeneralStats stats(length(expectedReads),0);
    demultiplex(reads, BarcodeMatcher, false, stats, ExactBarcodeMatching(), false);

    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(expectedReads[i].demuxResult, reads[i].demuxResult);
        for (unsigned j = 0; j < length(expectedReads[i].seq); ++j)
        {
            SEQAN_ASSERT_EQ(expectedReads[i].seq[j], reads[i].seq[j]);
        }
    }
}

// Checks the correctness of the functions if they are applied on external data.
SEQAN_DEFINE_TEST(Input_test)
{
    std::vector<std::string> bcids;
    std::vector<std::string> barcodes;
    using TRead = Read<seqan::Dna5QString>;

	CharString seqpath = SEQAN_PATH_TO_ROOT();
    append(seqpath, "/apps/seqan_flexbar/tests/seqs.fa");
	String<char> bcpath = SEQAN_PATH_TO_ROOT();
    append(bcpath, "/apps/seqan_flexbar/tests/barcodes.fa");

    std::vector<TRead> reads;
	SEQAN_ASSERT_EQ(0, loadSeqs(toCString(seqpath), reads));
	SEQAN_ASSERT_EQ(0, loadBarcodes(toCString(bcpath), bcids, barcodes));

	//Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(bcs);
	//Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);
	//indexRequire(indexSet, FibreSA());
    BarcodeMatcher BarcodeMatcher(barcodes);

	StringSet<String<int> > groups;
    GeneralStats stats(barcodes.size(),0);
    demultiplex(reads, BarcodeMatcher, false, stats, ExactBarcodeMatching(), false);
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(check_test);
	SEQAN_CALL_TEST(getPrefix_test);
	SEQAN_CALL_TEST(buildVariations_test);
	SEQAN_CALL_TEST(buildAllVariations_test);
	SEQAN_CALL_TEST(findExactIndex_test); 
	SEQAN_CALL_TEST(matchBarcodes_test); 
	SEQAN_CALL_TEST(clipBarcodes_test);
	SEQAN_CALL_TEST(clipBarcodesStrict_test);
    SEQAN_CALL_TEST(demultiplex_Exact_test);
    SEQAN_CALL_TEST(demultiplex_Exact_Multiplex_test);
    SEQAN_CALL_TEST(Input_test);
}
SEQAN_END_TESTSUITE
