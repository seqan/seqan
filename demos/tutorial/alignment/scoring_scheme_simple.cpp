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
    int match = -0;
    int mismatch = -1;
    int gap = -1;

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
//![init]

//![alignment]
    int score = globalAlignment(align, Score<int, Simple>(match, mismatch, gap));
    std::cout << "Score: " << score << std::endl;
    std::cout << align << std::endl;

    return 0;
}
//![alignment]
