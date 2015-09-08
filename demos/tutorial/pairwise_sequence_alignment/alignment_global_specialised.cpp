//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;     // align type

    TSequence seq1 = "GARFIELDTHECAT";
    TSequence seq2 = "GARFIELDTHEBIGCAT";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
//![main]

//![alignment]
    int score = globalAlignment(align, Score<int, Simple>(1, -1, -1), Hirschberg());
    std::cout << "Score: " << score << std::endl;
    std::cout << align << std::endl;

    return 0;
}
//![alignment]
