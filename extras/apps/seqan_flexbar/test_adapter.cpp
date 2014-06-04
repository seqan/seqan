#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include "adapter_trimming.h"

SEQAN_DEFINE_TEST(get_overlap_test)
{
	typedef seqan::String<seqan::Dna5> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps > TAlign;
	typedef seqan::Row<TAlign>::Type TRow;

    TSequence seq1 = "AAAAAAAAA";
    TSequence seq2 = "TTTTTTT";

    TAlign align;
    resize(rows(align), 2);
    seqan::assignSource(row(align,0),seq1);
    seqan::assignSource(row(align,1),seq2);

    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);

    // Overlap of 4.
    // ---AAAAAAAAA
    // TTTTTTT-----
    seqan::insertGaps(row1, 0, 3);
    seqan::insertGaps(row2, length(row2), 5);

    SEQAN_ASSERT_EQ(getOverlap(align), 4u);

    // Overlap of 0.
    // -------AAAAAAAAA
    // TTTTTTT---------
    seqan::clearGaps(align);
    seqan::insertGaps(row1, 0, 7);
    seqan::insertGaps(row2, length(row2), 9);

    SEQAN_ASSERT_EQ(getOverlap(align), 0u);

    // Overlap of 7.
    // AAAAAAAAA
    // -TTTTTTT-
    seqan::clearGaps(align);
    seqan::insertGap(row2, 0);
    seqan::insertGap(row2, length(row2));

    SEQAN_ASSERT_EQ(getOverlap(align), 7u);
}

SEQAN_DEFINE_TEST(count_gap_test)
{
	seqan::Dna5String seq("AAGTCTATCTA");
	seqan::Gaps<seqan::Dna5String> row(seq);

	// no gaps yet.
	SEQAN_ASSERT_EQ(countTotalGaps(row), 0u);

	// insert a few gaps at the start.
	seqan::insertGaps(row, 0, 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 10u);

	// create multiple continuous gaps

	seqan::insertGaps(row, 15, 10);
	seqan::insertGaps(row, length(row), 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 30u);
}

SEQAN_DEFINE_TEST(insert_size_test)
{
	// The insert method is specified only for actually
	// overlapping sequences so we test those.
	typedef seqan::String<seqan::Dna5> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps > TAlign;
	typedef seqan::Row<TAlign>::Type TRow;

    TSequence seq1 = "AAAAA";
    TSequence seq2 = "TTT";

    TAlign align;
    resize(rows(align), 2);
    seqan::assignSource(row(align,0),seq1);
    seqan::assignSource(row(align,1),seq2);

    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);

    // Insert size is 7
    // AAAAA--
    // ----TTT
    seqan::insertGaps(row1, seqan::length(row1), 2);
    seqan::insertGaps(row2, 0, 4);
    SEQAN_ASSERT_EQ(getInsertSize(align), 7u);

    // Insert size is 5
    // AAAAA
    // --TTT
    seqan::clearGaps(align);
    seqan::insertGaps(row2, 0, 2);
    SEQAN_ASSERT_EQ(getInsertSize(align), 5u);

    // Insert size is 4 (seq1 goes one base into the adapter).
	// AAAA|A
	// -TTT|-
	seqan::clearGaps(align);
	seqan::insertGap(row2, seqan::length(row2));
	seqan::insertGap(row2, 0);
	SEQAN_ASSERT_EQ(getInsertSize(align), 4u);
}

// Returns (num2 - num1 + 1) in a string "xxxx_num1_num2_xxx...".
int insertSize(char* s)
{
	int insert = 0;
	strtok(s, "_");
	insert -= atoi(strtok(NULL, "_"));
	insert += atoi(strtok(NULL, "_")) + 1;
	return insert;
}

SEQAN_DEFINE_TEST(match_test)
{
	Auto a = Auto();
	User u(7, 2);

	// Up to 5 overlap no error is allowed.
	SEQAN_ASSERT(isMatch(5,0,a));
	SEQAN_ASSERT_NOT(isMatch(5,1,a));
	// From 5 to 10 one error is allowed.
	SEQAN_ASSERT(isMatch(10,1,a));
	SEQAN_ASSERT_NOT(isMatch(10,2,a));
	// Otherwise 33% errors are allowed.
	SEQAN_ASSERT(isMatch(100,33,a));
	SEQAN_ASSERT_NOT(isMatch(100,34,a));

	// We need an overlap of length 7 and no more than 2 errors.
	SEQAN_ASSERT_NOT(isMatch(5,3,u));
	SEQAN_ASSERT_NOT(isMatch(7,3,u));
	SEQAN_ASSERT(isMatch(7,2,u));
}

SEQAN_DEFINE_TEST(strip_adapter_test)
{
	typedef seqan::String<seqan::Dna5Q> TSeq;
	typedef seqan::String<seqan::Dna5> TAda;

	TSeq seq("AAAAAAAAAATTTTT");
	TAda ada("TTTTTTTTTTT");

	int len = length(seq);
	int removed = stripAdapter(seq, ada, Auto());
	SEQAN_ASSERT_EQ(removed, 5);
	SEQAN_ASSERT_EQ(len - length(seq), 5u);
}

SEQAN_DEFINE_TEST(align_adapter_test)
{
	typedef seqan::String<seqan::Dna5Q> TSeq;
	typedef seqan::String<seqan::Dna5> TAda;

	TSeq seq("AAAAAAAAAATTTTT");
	TAda ada("TTTTTTTTTTT");

	seqan::Pair<unsigned, seqan::Align<TSeq> > pair;
    alignAdapter(pair, seq, ada);
	SEQAN_ASSERT_EQ(pair.i1, 5u);
}

SEQAN_DEFINE_TEST(strip_pair_test)
{
	typedef seqan::String<seqan::Dna5Q> TSeq;

	TSeq seq1("AAAAAAAAAATTTTTTT");
	TSeq seq2("CCCCCCCCCCAAAAAAA");
	TSeq seq3("AAAAAAAAAAGGGGGGG");

	unsigned insert = stripPair(seq1,seq2);
	SEQAN_ASSERT_EQ(insert, 27u);

	unsigned insert2 = stripPair(seq1,seq3);
	SEQAN_ASSERT_EQ(insert2, 0u);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    int tnum = 1;
#ifdef _OPENMP
	omp_set_num_threads(8);
	tnum = omp_get_max_threads();
#endif
	std::cout<<"\nRunning Tests using " << tnum << " thread(s).\n\n";
    SEQAN_CALL_TEST(get_overlap_test);
	SEQAN_CALL_TEST(count_gap_test);
	SEQAN_CALL_TEST(insert_size_test);
	SEQAN_CALL_TEST(match_test);
	SEQAN_CALL_TEST(strip_adapter_test);
	SEQAN_CALL_TEST(align_adapter_test);
	SEQAN_CALL_TEST(strip_pair_test);
}
SEQAN_END_TESTSUITE
