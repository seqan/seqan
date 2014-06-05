#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace seqan;

// Used for loading the sequences. Needed for the io-test.
int loadSeqs(char const * path, StringSet<String<char> >& ids, StringSet<String<Dna5Q> >& seqs)
{
	SequenceStream seqStream(path, SequenceStream::READ);
	unsigned records = 1000;									//number of records to be read
	resize(ids, records);
	resize(seqs, records);
	
	if (!isGood(seqStream))
	{
		std::cerr << "Error while opening the sequence-file.\n";
		return 1;
	}
	
	if (readBatch(ids, seqs, seqStream, records) != 0)
	{
		std::cerr << "Error while reading the sequences.\n";
		return 1;
	}
	return 0;
}
// Used for loading the sequences. Needed for the io-test.
int loadBarcodes(char const * path, StringSet<String<char> >& bcids, StringSet<String<Dna> >& bcs)
{
	SequenceStream bcStream(path, SequenceStream::READ);
	
	if (!isGood (bcStream))
	{
		std::cerr << "Error while opening barcode-file.\n";
		return 1;
	}

	if (readAll(bcids, bcs, bcStream) != 0)
	{
		std::cerr << "Error while reading the barcodes.\n";
		return 1;
	}
	return 0;
}
// Checks the correctness of the check function which checks the size of the barcodes and reads.
SEQAN_DEFINE_TEST(check_test)
{
	GeneralStats generalStats;
    
    StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "TACGTAGCTACTGACTGACT");
	appendValue(seqs, "G");
	appendValue(seqs, "TACGTCTGACT");
	appendValue(seqs, "");
	appendValue(seqs, "ATCGA");

	StringSet<String<char> > ids;
	appendValue(ids, "ErsteSeq");
	appendValue(ids, "LoeschenEins");
	appendValue(ids, "ZweiteSeq");
	appendValue(ids, "LoeschenLeer");
	appendValue(ids, "LoeschenFuenf");

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

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "TACGTAGCTACTGACTGACT");
	appendValue(exspectedSeqs, "TACGTCTGACT");

	StringSet<String<char> > exspectedIds;
	appendValue(exspectedIds, "ErsteSeq");
	appendValue(exspectedIds, "ZweiteSeq");

	bool resF = check(seqs, ids, barcodesFalse, generalStats);
	bool resT = check(seqs, ids, barcodesTrue, generalStats);

	SEQAN_ASSERT_EQ(false, resF);
	SEQAN_ASSERT_EQ(true, resT);
	SEQAN_ASSERT_EQ(2u, length(seqs));
	SEQAN_ASSERT_EQ(2u, length(ids));
	SEQAN_ASSERT_EQ(exspectedSeqs[0], seqs[0]);
	SEQAN_ASSERT_EQ(exspectedIds[0], ids[0]);
	SEQAN_ASSERT_EQ(exspectedSeqs[1], seqs[1]);
	SEQAN_ASSERT_EQ(exspectedIds[1], ids[1]);
}
// Checks the correctness of the getPrefix function which extract the prefices from a set of sequences.
SEQAN_DEFINE_TEST(getPrefix_test)
{
	StringSet<String<Dna5Q> > given;
	appendValue(given, "GATACAGACTGAGCATGTGATCGAC");
	appendValue(given, "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(given, "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	
	StringSet<String<Dna5Q> > exspectedSet;
	appendValue(exspectedSet, "GATACAGACTGAGCATGTGATCGAC");
	appendValue(exspectedSet, "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSet, "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");

	StringSet<String<Dna5Q> > exspected;
	appendValue(exspected, "GATACA");
	appendValue(exspected, "AATTCC");
	appendValue(exspected, "GTTGGA");
		
	StringSet<String<Dna5Q> > res;
    getPrefix(res, given, 6);
	for (unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
		SEQAN_ASSERT_EQ(exspectedSet[i], given[i]);
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
	/*for (unsigned i = 0; i < 10; ++i)
	{
		SEQAN_ASSERT_EQ(id1, bcids[i]);
	}
	for (unsigned i = 10; i < 20; ++i)
	{
		SEQAN_ASSERT_EQ(id2, bcids[i]);
	}*/
}

// Checks the correctness of the findExactIndex function which searches for one piece of sequence in the barcodes. Implicitly checks the construction of the Index.
SEQAN_DEFINE_TEST(findExactIndex_test)
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");

	StringSet<String<Dna5Q> > readPieces;
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "TTTTTT");
	appendValue(readPieces, "GGGGGG");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "GATACA");
	appendValue(readPieces, "ACGTAC");
	appendValue(readPieces, "ATGACNAANG");	//can't happen in the first place...

	Index<StringSet<String<Dna5Q> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna5Q> >, IndexEsa<> > > esaFinder(indexSet);
	indexRequire(indexSet, FibreSA());

	int exspected[] = {1,0,3,2,-1,-1,4,-1};

	for (unsigned i = 0; i < length(readPieces); ++i)
	{
		int res = findExactIndex(readPieces[i], esaFinder);
		SEQAN_ASSERT_EQ(exspected[i], res);
	}
}
// Checks the correctnes of the findAllExactIndex function which searches for many pieces of sequence in the barcodes. Implicitly checks the construction of the Index.
SEQAN_DEFINE_TEST(findAllExactIndex_test) 
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");
	
	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

	StringSet<String<Dna5Q> > readPieces;
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "TTTTTT");
	appendValue(readPieces, "GGGGGG");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "GATACA");
	appendValue(readPieces, "ACGTAC");

	Index<StringSet<String<Dna5Q> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna5Q> >, IndexEsa<> > > esaFinder(indexSet);
	indexRequire(indexSet, FibreSA());

	int exspected[] = {1,0,3,2,-1,-1,4,};
	String<int> res; 
    findAllExactIndex(res, readPieces, esaFinder, demultiplexStats);
	for (unsigned i = 0; i < length(res); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
	
}
// Checks the correctnes of the clipBarcodes function which erases the first x bases of a sequence.
SEQAN_DEFINE_TEST(clipBarcodes_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	String<int> matches;
	appendValue(matches, 0);
	appendValue(matches, -1);

	StringSet<String<Dna5Q> > exspected ;
	appendValue(seqs, "GTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	clipBarcodes(seqs, matches, 6);
	for (unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], seqs[i]);
	}
	
}
// Checks the correctnes of the clipBarcodes function using the hardClip method and therefore also clipping sequences without matching barcode.
SEQAN_DEFINE_TEST(clipBarcodesStrict_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	StringSet<String<Dna5Q> > exspected ;
	appendValue(seqs, "GTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGG");

	clipBarcodes(seqs, 6);
	for (unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], seqs[i]);
	}
	
}
// Checks the correctnes of the group function which builds the groups of barcodes matching the same barcodes. Implicitly checks resizeGroups.
SEQAN_DEFINE_TEST(group_test)
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "GAGAGA"); //barcode without matching sequences 
	appendValue(barcodes, "CCCCCC");
	
	String<int> matches;
	appendValue(matches, 0);
	appendValue(matches, 4);
	appendValue(matches, 1);
	appendValue(matches, -1);
	appendValue(matches, 1);
	appendValue(matches, 2);

	StringSet<String<int> > exspected;
	resize(exspected, 7);
	appendValue(exspected[0], 3);
	appendValue(exspected[1], 0);
	appendValue(exspected[2], 2);
	appendValue(exspected[2], 4);	//expected[4] is left empty since the barcode has no matching sequences
	appendValue(exspected[3], 5);
	appendValue(exspected[5], 1);
	
	StringSet<String<int> > res;
    group(res, matches, barcodes, false);
	
	for (unsigned i = 0; i < length(exspected); ++i)
	{
		for (unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
		}
	}
}
// Checks the correctness of the doAll function performing all demultiplexing operations for exact inline barcode matching.
SEQAN_DEFINE_TEST(doAll_Exact_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(seqs, "GGGGGGAGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(seqs, "AAAAAATAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(seqs, "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");

	StringSet<String<char> > ids;
	appendValue(ids, "Adenin1");
	appendValue(ids, "Guanin");
	appendValue(ids, "Adenin2");
	appendValue(ids, "Unidentifiziert");

	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSeqs, "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");

	StringSet<String<int> > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], 3);
	appendValue(exspected[1], 1);	//exspected[2] is left empty since the associated barcode is not matched
	appendValue(exspected[3], 0);
	appendValue(exspected[3], 2);

	Index<StringSet<String<Dna5Q> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna5Q> >, IndexEsa<> > > esaFinder(indexSet);
	indexRequire(indexSet, FibreSA());

	StringSet<String<int> > res;
    doAll(res, seqs, barcodes, esaFinder, false, demultiplexStats, false);

	for (unsigned i = 0; i < length(exspected); ++i)
	{
		for (unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}
// Checks the correctness of the doAll function performing all demultiplexing operations
// for exact multiplex barcode matching.
SEQAN_DEFINE_TEST(doAll_Exact_Multiplex_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(seqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(seqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(seqs, "GGGG");

	StringSet<String<Dna5Q> > multiplex;
	appendValue(multiplex, "AAAAAA");
	appendValue(multiplex, "GGGGGG");
	appendValue(multiplex, "GGCCGG");
	appendValue(multiplex, "AAAAAA");

	StringSet<String<Dna> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSeqs, "GGGG");

	StringSet<String<int> > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], 2);
	appendValue(exspected[1], 1);	//exspected[2] is left empty since the associated barcode is not matched
	appendValue(exspected[3], 0);
	appendValue(exspected[3], 3);

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);
	indexRequire(indexSet, FibreSA());

	StringSet<String<int> > res;
    doAll(res, multiplex, barcodes, esaFinder, demultiplexStats, false);

	for (unsigned i = 0; i < length(exspected); ++i)
	{
		for (unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}
// Checks the correctness of the functions if they are applied on external data.
SEQAN_DEFINE_TEST(Input_test)
{
	StringSet<String<char> > ids;
	StringSet<String<Dna5Q> > seqs;
	StringSet<String<char> > bcids;
	StringSet<String<Dna> > bcs;

	CharString seqpath = SEQAN_PATH_TO_ROOT();
    append(seqpath, "/extras/apps/seqan_flexbar/test_data/seqs.fasta");
	String<char> bcpath = SEQAN_PATH_TO_ROOT();
    append(bcpath, "/extras/apps/seqan_flexbar/test_data/barcodes.fasta");

	SEQAN_ASSERT_EQ(0, loadSeqs(toCString(seqpath), ids, seqs));
	SEQAN_ASSERT_EQ(0, loadBarcodes(toCString(bcpath), bcids, bcs));

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(bcs)+1);

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(bcs);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);
	indexRequire(indexSet, FibreSA());

	StringSet<String<int> > groups;
    doAll(groups, seqs, bcs, esaFinder, false, demultiplexStats, false);
	String<StringSet<String<Dna5Q> > > gSeqs;
	String<StringSet<String<char> > > gIds;
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	omp_set_num_threads(8);
	int tnum = omp_get_max_threads();
	std::cout<<"\nRunning Tests using " << tnum << " thread(s).\n\n";
	SEQAN_CALL_TEST(check_test);
	SEQAN_CALL_TEST(getPrefix_test);
	SEQAN_CALL_TEST(buildVariations_test);
	SEQAN_CALL_TEST(buildAllVariations_test);
	SEQAN_CALL_TEST(findExactIndex_test); 
	SEQAN_CALL_TEST(findAllExactIndex_test); 
	SEQAN_CALL_TEST(clipBarcodes_test);
	SEQAN_CALL_TEST(clipBarcodesStrict_test);
	SEQAN_CALL_TEST(group_test);
	SEQAN_CALL_TEST(doAll_Exact_test);
	SEQAN_CALL_TEST(doAll_Exact_Multiplex_test);
	SEQAN_CALL_TEST(Input_test);
}
SEQAN_END_TESTSUITE
