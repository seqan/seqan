//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<Rna> TSequence;
    typedef Align<TSequence, ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    typedef Iterator<TRow>::Type TRowIterator;
//![main]

//![init]
    TSequence seq1 = "AAGUGACUUAUUG";
    TSequence seq2 = "AGUCGGAUCUACUG";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
//![init]

//![alignment]
    int score = globalAlignment(align, MyersHirschberg());
    std::cout << "Score: " << score << std::endl;
    std::cout << align << std::endl;
//![alignment]

//![view]
    unsigned aliLength = _max(length(row(align, 0)), length(row(align, 1)));
    for (unsigned i = 0; i < length(rows(align)); ++i)
    {
        TRowIterator it = iter(row(align, i), 0);
        TRowIterator itEnd = iter(row(align, i), aliLength);
        unsigned pos = 0;
        std::cout << "Row " << i << " contains gaps at positions: ";
        std::cout << std::endl;
        while (it != itEnd)
        {
            if (isGap(it))
                std::cout << pos << std::endl;
            ++it;
            ++pos;
        }
    }

    return 0;
}
//![view]
