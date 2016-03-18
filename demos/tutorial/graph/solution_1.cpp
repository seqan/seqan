#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int main()
{
    typedef unsigned int TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TGraph g;

    TVertexDescriptor vertBerlin = addVertex(g);
    TVertexDescriptor vertHamburg = addVertex(g);
    TVertexDescriptor vertHannover = addVertex(g);
    TVertexDescriptor vertMuenchen = addVertex(g);

    addEdge(g, vertBerlin, vertHamburg, 289u);
    addEdge(g, vertHamburg, vertHannover, 289u);
    addEdge(g, vertHannover, vertMuenchen, 572u);

    std::ofstream dotFile("graph.dot");
    writeRecords(dotFile, g, DotDrawing());

    dotFile.close();

    std::cout << g << std::endl;

    return 0;
}
