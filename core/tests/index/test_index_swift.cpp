#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

// Test SWIFT finder with empty pattern.
SEQAN_DEFINE_TEST(test_index_swift_find_empty_pattern)
{
    using namespace seqan;

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

SEQAN_BEGIN_TESTSUITE(test_index_swift)
{
	SEQAN_CALL_TEST(test_index_swift_find_empty_pattern);
}
SEQAN_END_TESTSUITE
