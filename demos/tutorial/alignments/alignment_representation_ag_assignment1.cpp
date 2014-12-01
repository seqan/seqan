// FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    // Define the types we need.
    typedef String<char> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;

    // Initializing the sequences and the string set.
    TSequence seq1 = "GARFIELDTHECAT";
    TSequence seq2 = "GARFIELDTHEBIGCAT";
    TSequence seq3 = "THEBIGCAT";

    TStringSet strings;
    appendValue(strings, seq1);
    appendValue(strings, seq2);
    appendValue(strings, seq3);

    // Load the string set into the Alignment Graph.
    TAlignGraph alignG(strings);

    // Add two vertices covering "GARFIELD" in the first and the second sequence and connect them with an edge.
    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG),0), 0, 8),
                    addVertex(alignG, positionToId(stringSet(alignG),1), 0, 8));

    // Add two vertices covering "THE" in the first and the second sequence and connect them with an edge.
    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG),0), 8, 3),
                    addVertex(alignG, positionToId(stringSet(alignG),1), 8, 3));

    // Find the vertex covering "THE" in the first sequence and add the vertex covering "THE" in the third sequence and connect them with an edge.
    addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG),0), 8),
                    addVertex(alignG, positionToId(stringSet(alignG),2), 0, 3));

    // Find the vertices covering "THE" in the second and the third sequence and connect them with an edge.
    addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG),1), 8),
                    findVertex(alignG, positionToId(stringSet(alignG),2), 0));

    // Add two vertices covering "FAT" in the second and the third sequence and connect them with an edge.
    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG),1), 11, 3),
                    addVertex(alignG, positionToId(stringSet(alignG),2), 3, 3));

    // Add two vertices covering "CAT" in the first and the second sequence and connect them with an edge.
    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG),0), 11, 3),
                    addVertex(alignG, positionToId(stringSet(alignG),1), 14, 3));

    // Find the vertex covering "CAT" in the first sequence and add the vertex covering "CAT" in the third sequence and connect them with an edge.
    addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG),0), 11),
                    addVertex(alignG, positionToId(stringSet(alignG),2), 6, 3));

    // Find the vertices covering "CAT" in the second and the third sequence and connect them with an edge.
    addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG),1), 14),
                    findVertex(alignG, positionToId(stringSet(alignG),2), 6));

    ::std::cout << alignG << ::std::endl;

    return 0;
}
