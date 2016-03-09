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
    typedef typename VertexDescriptor<TAlignGraph>::Type TVertexDescriptor;
//![typedef]

    std::cout << "//![output_init]" << std::endl;
//![init]
    TSequence seq1 = "TTGT";
    TSequence seq2 = "TTAGT";

    TStringSet strings;
    appendValue(strings, seq1);
    appendValue(strings, seq2);

    TAlignGraph alignG(strings);
    std::cout << alignG << std::endl;
//![init]
    std::cout << "//![output_init]" << std::endl;

    std::cout << "//![output_construct]" << std::endl;
//![construct]
    TVertexDescriptor u,v;

    // TT
    u = addVertex(alignG, positionToId(stringSet(alignG), 0), 0, 2);
    v = addVertex(alignG, positionToId(stringSet(alignG), 1), 0, 2);
    addEdge(alignG, u, v);

    // A
    addVertex(alignG, positionToId(stringSet(alignG), 1), 2, 1);

    // GT
    addVertex(alignG, positionToId(stringSet(alignG), 0), 2, 2);
    addVertex(alignG, positionToId(stringSet(alignG), 1), 3, 2);

    u = findVertex(alignG, positionToId(stringSet(alignG), 0), 2);
    v = findVertex(alignG, positionToId(stringSet(alignG), 1), 3);
    addEdge(alignG, u, v);

    std::cout << alignG << std::endl;
//![construct]
    std::cout << "//![output_construct]" << std::endl;
//![construct]

    return 0;
}
//![construct]
