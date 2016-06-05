//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<AminoAcid> TSequence;             // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;      // align type
//![main]

//![init]
    TSequence seq1 = "TELKDD";
    TSequence seq2 = "LKTEL";
    int gapExtend = -2;
    int gapOpen = -10;


    TAlign alignAffine;
    resize(rows(alignAffine), 2);
    assignSource(row(alignAffine, 0), seq1);
    assignSource(row(alignAffine, 1), seq2);

    TAlign alignDynamic;
    resize(rows(alignDynamic), 2);
    assignSource(row(alignDynamic, 0), seq1);
    assignSource(row(alignDynamic, 1), seq2);
//![init]

//![alignment]
    int scoreAffine = globalAlignment(alignAffine, Blosum62(gapExtend, gapOpen), AffineGaps());
    std::cout << "ScoreAffine: " << scoreAffine << std::endl;
    std::cout << alignAffine << std::endl;

    int scoreDynamic = globalAlignment(alignDynamic, Blosum62(gapExtend, gapOpen), DynamicGaps());
    std::cout << "ScoreDynamic: " << scoreDynamic << std::endl;
    std::cout << alignDynamic << std::endl;

    return 0;
}
//![alignment]
