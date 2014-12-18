//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<Dna> TSequence;                 // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;    // align type
//![main]

//![init]
    Align<String<char> > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), "ataagcgtctcg");
    assignSource(row(ali, 1), "tcatagagttgc");
//![init]

//![alignment]
    Score<int> scoring(2, -1, -2, 0);
    LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, 5);
    while (nextLocalAlignment(ali, enumerator))
    {
        std::cout << "Score = " << getScore(enumerator) << std::endl;
        std::cout << ali;
        std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0)) - 1) << "]";
        std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" <<  (clippedEndPosition(row(ali, 1)) - 1) << "]" << std::endl << std::endl;
    }
    return 0;
}
//![alignment]
