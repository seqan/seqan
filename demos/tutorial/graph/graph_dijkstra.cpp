//![includes]
#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;
//![includes]

//![main-typedefs]
int main()
{
    typedef unsigned int TCargo;
    typedef Graph<Undirected<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
//![main-typedefs]

//![create-g]
    TGraph g;
//![create-g]

//![create-vertices]
    TVertexDescriptor vertBerlin = addVertex(g);
    TVertexDescriptor vertHamburg = addVertex(g);
    TVertexDescriptor vertHannover = addVertex(g);
    TVertexDescriptor vertMainz = addVertex(g);
    TVertexDescriptor vertMuenchen = addVertex(g);
//![create-vertices]

//![create-edges]
    addEdge(g, vertBerlin, vertHamburg, 289);
    addEdge(g, vertBerlin, vertHannover, 286);
    addEdge(g, vertBerlin, vertMainz, 573);
    addEdge(g, vertBerlin, vertMuenchen, 586);
    addEdge(g, vertHannover, vertMuenchen, 572);
    addEdge(g, vertHamburg, vertMainz, 521);
//![create-edges]

//![main-graph-io]
    std::ofstream dotFile("graph.dot");
    writeRecords(dotFile, g, DotDrawing());
    dotFile.close();
//![main-graph-io]

//![alternatively-graph-io]
    std::cout << g << '\n';
//![alternatively-graph-io]

//![definition-property-map]
    typedef String<char> TCityName;
    typedef String<TCityName> TProperties;
    TProperties cityNames;
    resizeVertexMap(cityNames, g);
//![definition-property-map]

//![enter-properties]
    assignProperty(cityNames, vertBerlin, "Berlin");
    assignProperty(cityNames, vertHamburg, "Hamburg");
    assignProperty(cityNames, vertMuenchen, "Munich");
    assignProperty(cityNames, vertMainz, "Mainz");
    assignProperty(cityNames, vertHannover, "Hannover");
//![enter-properties]

    std::cout << "//![iterate-and-output-properties]\n";
//![iterate-and-output-properties]
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator itV(g);
    for (; !atEnd(itV); goNext(itV))
    {
        std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << std::endl;
    }
//![iterate-and-output-properties]
    std::cout << "//![iterate-and-output-properties]\n";

//![dijkstra-containers]
    typedef Size<TGraph>::Type TSize;
    InternalPropertyMap<TCargo> cargoMap;
    String<TVertexDescriptor> predMap;
    String<TSize> distMap;
//![dijkstra-containers]
//![dijkstra]
    dijkstra(predMap, distMap, g, vertHannover, cargoMap);
//![dijkstra]

//![dijkstra-output]
    TVertexIterator itV2(g);
    while (!atEnd(itV2))
    {
        std::cout << "Shortest path from " << property(cityNames, vertHannover) << " to " << property(cityNames, value(itV2)) << ": ";
        std::cout << property(distMap, value(itV2)) << std::endl;
        goNext(itV2);
    }

    return 0;
}
//![dijkstra-output]
