#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

// Test SWIFT finder with empty pattern.
SEQAN_DEFINE_TEST(test_index_swift_find_empty_pattern)
{
    using namespace seqan2;

    typedef Finder<Dna5String, Swift<SwiftSemiGlobal> > TSwiftFinder;

    typedef Dna5String TReadSeq_;
    typedef StringSet<TReadSeq_> TReadSet;
    typedef typename Value<TReadSet>::Type TReadSeq;
    typedef typename Value<TReadSeq>::Type TAlphabet;
    typedef Shape<TAlphabet, UngappedShape<10> > TShape;
    typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> > TQGramIndex;
    typedef Pattern<TQGramIndex, Swift<SwiftSemiGlobal> > TSwiftPattern;

    Dna5String contigSeq = "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT";

    StringSet<Dna5String> readSeqs;

    TQGramIndex index(readSeqs);
    TSwiftPattern swiftPattern(index);

    TSwiftFinder swiftFinder(contigSeq);
    while (find(swiftFinder, swiftPattern, 0.1))
        continue;
}

// Test SWIFT finder with edge case hash collision
// DNA 32-mers are mapped into 64-bit space
//     enc(AA..AA) = 0
//     enc(TT..TT) = 2^64-1 = uint64_t::max()
// Because an empty bucket (TBucketMap::EMPTY) is also marked by uint64_t::max(), there is no way to know if a bucket
// was filled with the TT..TT k-mer or still empty.
// This test checks that the index constructor throws an exception if the q-gram size is too large.
// It only applies to shapes with variable size (SimpleShape) as the fixed size shapes (UngappedShape, GappedShape)
// are checked via a static assert.
SEQAN_DEFINE_TEST(test_index_swift_edge_case_hash_collision)
{
    using namespace seqan2;

    typedef DnaString TReadSeq_;
    typedef StringSet<TReadSeq_> TReadSet;
    typedef typename Value<TReadSet>::Type TReadSeq;
    typedef typename Value<TReadSeq>::Type TAlphabet;
    typedef Shape<TAlphabet, SimpleShape> TShape;
    typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> > TQGramIndex;

    TReadSet readSeqs;
    TReadSeq lastKmerCollision = "AGGTATTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "CCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

    appendValue(readSeqs, lastKmerCollision);

    try
    {
        TQGramIndex index(readSeqs, TShape{32});
        SEQAN_FAIL("The expected exception was not thrown.");
    }
    catch(const std::runtime_error & e)
    {
        std::string expected = "Incompatible q-gram size. Must be in range [1, 31].";
        SEQAN_ASSERT_EQ(expected, e.what());
        TQGramIndex index(readSeqs, TShape{31});
    }
}

SEQAN_BEGIN_TESTSUITE(test_index_swift)
{
	SEQAN_CALL_TEST(test_index_swift_find_empty_pattern);
    SEQAN_CALL_TEST(test_index_swift_edge_case_hash_collision);
}
SEQAN_END_TESTSUITE
