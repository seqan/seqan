//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
//![main]
//![typedef]
    typedef String<Dna> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;
//![typedef]

//![init]
    TSequence seq1 = "TTGT";
    TSequence seq2 = "TTAGT";

    TStringSet strings;
    appendValue(strings, seq1);
    appendValue(strings, seq2);

    TAlignGraph alignG(strings);
//![init]

//![construct]
    std::cout << alignG << std::endl;

    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG), 0), 0, 2),
            addVertex(alignG, positionToId(stringSet(alignG), 1), 0, 2));

    addVertex(alignG, positionToId(stringSet(alignG), 1), 2, 1);

    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG), 0), 2, 2),
            addVertex(alignG, positionToId(stringSet(alignG), 1), 3, 2));

    std::cout << alignG << std::endl;

    return 0;
}
//![construct]
