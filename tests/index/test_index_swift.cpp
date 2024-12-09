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
/*
    DNA 32-mers are mapped into 64-bit space
        enc(AA..AA) = 0
        enc(TT..TT) = 2^64-1 = uint64_t::max()

    Because an empty bucket is also marked by uint64_t::max() there is no way to know if a bucket was filled with the TT..TT k-mer or still empty 
*/
SEQAN_DEFINE_TEST(test_index_swift_edge_case_hash_collision)
{
    using namespace seqan2;

    typedef Finder<DnaString, Swift<SwiftSemiGlobal> > TSwiftFinder;

    typedef DnaString TReadSeq_;
    typedef StringSet<TReadSeq_> TReadSet;
    typedef typename Value<TReadSet>::Type TReadSeq;
    typedef typename Value<TReadSeq>::Type TAlphabet;
    typedef Shape<TAlphabet, UngappedShape<32u> > TShape;
    typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> > TQGramIndex;
    typedef Pattern<TQGramIndex, Swift<SwiftSemiGlobal> > TSwiftPattern;

    
    DnaString randomContig = "GTGTGCATTTTTCATTTCCCACGTTTTTCAGTGTTTCCTGCCATATTCCAGACCTAGAGTTTGAGTTTCTCATTTTTCACTTTTTTTCA"
                             "CCTGCAGTGTGTGTGTGTCTCATTTTCTAACTTTTTCATTGTTTCTCCCCATATTTCAGGTCCTACAGTATGTGTGTCTGATTTTCCAC"
                             "GTTTTTCAGTCTTAATCACCATATTCCAGGTCCTAGGTTGTGCATTTCTCATTTTCCACATTTTTCAGTGTTTCT";
    
    StringSet<DnaString> readSeqs;
    DnaString lastKmerCollision = "AGGTATTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "CCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                                  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

    appendValue(readSeqs, lastKmerCollision);

    TQGramIndex index(readSeqs);
    TSwiftPattern swiftPattern(index);

    TSwiftFinder swiftFinder(randomContig);
    while (find(swiftFinder, swiftPattern, 0.0))
        continue;
}
SEQAN_BEGIN_TESTSUITE(test_index_swift)
{
	SEQAN_CALL_TEST(test_index_swift_find_empty_pattern);
    SEQAN_CALL_TEST(test_index_swift_edge_case_hash_collision);
}
SEQAN_END_TESTSUITE
