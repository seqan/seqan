//![completeSolution]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<Dna> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef Align<TSequence, ArrayGaps> TAlign;

    TStringSet strings;
    appendValue(strings, "AAATGACATGGATTG");
    appendValue(strings, "AGCGGACTCTACTTG");
    appendValue(strings, "AGTCGATAACTG");
    appendValue(strings, "AGTCGGATCTACTG");
    appendValue(strings, "AGCGGCATTG");

    int bestScore = MinValue<int>::VALUE;
    int bestSeqIdx1 = 0;
    int bestSeqIdx2 = 0;
    for (unsigned i = 0; i < length(strings) - 1; ++i)
    {
        for (unsigned j = i + 1; j < length(strings); ++j)
        {
            int tmpScore = globalAlignmentScore(strings[i], strings[j], MyersBitVector());
            if (tmpScore > bestScore)
            {
                bestScore = tmpScore;
                bestSeqIdx1 = i;
                bestSeqIdx2 = j;
            }

        }
    }

    return 0;
}
//![completeSolution]
