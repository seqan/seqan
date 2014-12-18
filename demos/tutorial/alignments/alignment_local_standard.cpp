//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<Dna> TSequence;                 // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;      // align type
//![main]

//![init]
    TSequence seq1 = "ACGTGACGGGATGTG";
    TSequence seq2 = "ACGGCGGGACTGACTG";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
//![init]

//![alignment]
    int score = localAlignment(align, Score<int, Simple>(1, -1, -1, -1));
    std::cout << "Score: " << score << std::endl;
    std::cout << align << std::endl;

    return 0;
}
//![alignment]
