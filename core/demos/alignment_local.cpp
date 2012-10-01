///A tutorial about local alignments.
#include <iostream>
#include <seqan/score.h>
#include <seqan/align.h>

using namespace seqan;

int main()
{
///Example 1: This program applies the Smith-Waterman algorithm to compute the best local alignment between two given sequences.
    StringSet<CharString> strings;
    Align<String<char> > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), "aphilologicaltheorem");
    assignSource(row(ali, 1), "bizarreamphibology");
    int score = localAlignment(ali, Score<int>(3,-3,-2, -2));
    std::cout << "Score = " << score << std::endl;
    std::cout << ali;
    unsigned cBeginPos = clippedBeginPosition(row(ali, 0));
    unsigned cEndPos = clippedEndPosition(row(ali, 0)) - 1;
    std::cout << "Aligns Seq1[" << cBeginPos << ":" << cEndPos << "]";
    std::cout << " and Seq2[" << cBeginPos << ":" << cEndPos << "]" << std::endl << std::endl;


///Example 2: This program applies the Waterman-Eggert algorithm to compute all non-overlapping local alignments with score better or equal 2.
    Align<String<Dna> > ali2;
    resize(rows(ali2), 2);
    assignSource(row(ali2, 0), "ataagcgtctcg");
    assignSource(row(ali2, 1), "tcatagagttgc");

    Score<int> scoring(2, -1, -2, 0);
    LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, 2);
    while (nextLocalAlignment(ali2, enumerator))
    {
        std::cout << "Score = " << getScore(enumerator) << std::endl;
        std::cout << ali2;
        unsigned cBeginPos0 = clippedBeginPosition(row(ali2, 0));
        unsigned cEndPos0 = clippedEndPosition(row(ali2, 0)) - 1;
        unsigned cBeginPos1 = clippedBeginPosition(row(ali2, 1));
        unsigned cEndPos1 = clippedBeginPosition(row(ali2, 1)) - 1;
        std::cout << "Aligns Seq1[" << cBeginPos0 << ":" << cEndPos0 << "]";
        std::cout << " and Seq2[" << cBeginPos1 << ":" <<  cEndPos1 << "]";
        std::cout << std::endl << std::endl;
    }

///Example 3
    Align<String<Dna> > ali3;
    resize(rows(ali3), 2);
    assignSource(row(ali3, 0), "cccccc");
    assignSource(row(ali3, 1), "tttttggccccccgg");

    Score<int> scoring3(1, -1, -1, -1);
    LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator3(scoring3, 5);
    while (nextLocalAlignment(ali3, enumerator3))
    {
        std::cout << "Score = " << getScore(enumerator3) << std::endl;
        std::cout << ali3;
        unsigned cBeginPos0 = clippedBeginPosition(row(ali3, 0));
        unsigned cEndPos0 = clippedEndPosition(row(ali3, 0)) - 1;
        unsigned cBeginPos1 = clippedBeginPosition(row(ali3, 1));
        unsigned cEndPos1 = clippedEndPosition(row(ali3, 1)) - 1;
        std::cout << "Aligns Seq1[" << cBeginPos0 << ":" << cEndPos0 << "]";
        std::cout << " and Seq2[" << cBeginPos1 << ":" << cEndPos1 << "]";
        std::cout << std::endl << std::endl;
    }

///Example 4: This program applies the banded Waterman-Eggert algorithm to compute all non-overlapping local alignments with score or equal 5
///           in the band from diagonal -1 to diagonal 8.
    Align<String<Dna5> > ali4;
    resize(rows(ali4), 2);
    assignSource(row(ali4, 0), "AAAAAAANAAAGGGNGGGGGGGGNGGGGGANAA");
    assignSource(row(ali4, 1), "GGGGGGCGGGGGGGA");

    LocalAlignmentFinder<> finder4(ali4);
    Score<int> scoring4(1, -1, -1, -1);
    LocalAlignmentEnumerator<Score<int>, Banded> enumerator4(scoring3, -1, 8, 5);
    while (nextLocalAlignment(ali4, enumerator4))
    {
        std::cout << "Score = " << getScore(enumerator4) << std::endl;
        std::cout << ali4;
        unsigned cBeginPos0 = clippedBeginPosition(row(ali4, 0));
        unsigned cEndPos0 = clippedEndPosition(row(ali4, 0)) - 1;
        unsigned cBeginPos1 = clippedBeginPosition(row(ali4, 1));
        unsigned cEndPos1 = clippedEndPosition(row(ali4, 1)) - 1;
        std::cout << "Aligns Seq1[" << cBeginPos0 << ":" << cEndPos0 << "]";
        std::cout << " and Seq2[" << cBeginPos1 << ":" << cEndPos1 << "]";
        std::cout << std::endl << std::endl;
    }

    return 0;
}
