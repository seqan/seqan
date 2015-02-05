#ifndef TESTS_INDEX_MINIMIZER_TEST_INDEX_MINIMIZER_H_
#define TESTS_INDEX_MINIMIZER_TEST_INDEX_MINIMIZER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
// A test for strings.

using namespace seqan;

SEQAN_DEFINE_TEST(test_index_minimizer_getOccurrences)
{
    typedef CharString                                            TText;
    typedef Infix<TText>::Type                                    TTextInfix;
    typedef Iterator<TText, Standard>::Type                       TTextIter;
    typedef Index<TText, IndexQGram<MinimizerShape<30,10> > >     TIndex;
    typedef Fibre<TIndex, FibreSA>::Type const                    TSA;
    typedef Infix<TSA>::Type                                      TOccurrences;

    TText text = "ACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTA";

    TIndex index(text);

    TTextIter textEnd = end(text, Standard()) - length(indexShape(index)) + 1;
    for (TTextIter textIt = begin(text, Standard()); textIt != textEnd; textIt++)
    {
        hash(indexShape(index), textIt);
        TTextInfix textInfix = infix(text, textIt, textIt + length(indexShape(index)));

        std::cout << textInfix << std::endl;
        std::cout << value(indexShape(index)) << std::endl;

//        Shape<Dna, UngappedShape<10> > tmpShape;
//        tmpShape.

        TOccurrences occ = getOccurrences(index, indexShape(index), textInfix);

        for (unsigned i = 0; i < length(occ); i++)
        {
            std::cout << infix(text, occ[i], occ[i] + length(indexShape(index))) << std::endl;
//            SEQAN_ASSERT_EQ(infix(text, occ[i], occ[i] + length(indexShape(index))), textInfix);
        }

        std::cout << "-----------------" << std::endl;
    }
}

/*SEQAN_DEFINE_TEST(testMinimizerFind)
{
    typedef Index<String<char>, IndexQGram<MinimizerShape<4, 2> > > TQGramIndex;
    TQGramIndex idx("to be or not to be");
    Finder<TQGramIndex> finder(idx);

    SEQAN_ASSERT(find(finder, "be"));
    SEQAN_ASSERT(position(finder) == 3);
    SEQAN_ASSERT(find(finder, "be"));
    SEQAN_ASSERT(position(finder) == 16);
    SEQAN_ASSERT(!find(finder, "be"));

// }

}
*/
#endif  // TESTS_INDEX_MINIMIZER_TEST_INDEX_MINIMIZER_H_
