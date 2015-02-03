#ifndef TESTS_INDEX_MINIMIZER_TEST_INDEX_MINIMIZER_H_
#define TESTS_INDEX_MINIMIZER_TEST_INDEX_MINIMIZER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
// A test for strings.

using namespace seqan;

SEQAN_DEFINE_TEST(test_index_minimizer_getOccurrences)
{
    DnaString text = "ACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTAACGCATTTTGCATTA";
    Index<DnaString, IndexQGram< MinimizerShape<30,10> > > index(text);
    Shape<Dna, MinimizerShape<30, 10> > myShape;
    Shape<Dna, UngappedShape<30> > shape;
//    std::cout << text << std::endl;
//    std::cout << length(myShape) << " " << weight(myShape) << std::endl;
 //   for (int i = 0; i < length(text) - length(myShape) + 1; i++)
  //  {
   //     std::cout << hashNext(myShape, begin(text) + i) << " " << std::endl;
   // } 
//    hashInit(myShape, begin(text)); 
 //   for (int i=0; i < 10; i++) 
  //      std::cout << hash(myShape, begin(text) + i) << " ";
        //std::cout << hashNext(myShape, begin(text) + i) << std::endl;
 /*   resize(indexSA(index), _qgramQGramCount(index), Exact());
    resize(indexDir(index), _fullDirLength(index), Exact());
    _qgramClearDir(indexDir(index), index.bucketMap); 
    _qgramCountQGrams(indexDir(index), indexBucketMap(index), indexText(index), myShape, getStepSize(index));
    for (unsigned i = 0; i < _fullDirLength(index) ; i++)
        std::cout << indexDir(index)[i] << " " ; 
    std::cout << "done" << std::endl;
    _qgramCummulativeSum(index.dir, False());
    for (unsigned i = 0; i < _fullDirLength(index); i++)
       std::cout << indexDir(index)[i] << " "; 
    std::cout << std::endl;
    _qgramFillSuffixArray(indexSA(index), indexText(index), indexShape(index), indexDir(index), indexBucketMap(index), getStepSize(index), False());
*/
    indexRequire(index, FibreSA());
    for (unsigned k = 0; k < length(indexText(index)) - 30 + 1; k++)
    {
        //std::cout << " searched k-mer starts at " << k <<"th position of the text" << std::endl;
        hash(indexShape(index), begin(indexText(index)) + k);
        for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); i++) 
        //std::cout << getOccurrences(index, indexShape(index))[i] <<" ";
        SEQAN_ASSERT( hash(shape, begin(indexText(index)) + getOccurrences(index, indexShape(index))[i]) == indexShape(index).nhValue);
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
