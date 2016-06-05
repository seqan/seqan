//![includes]
#include <iostream>
#include <seqan/graph_algorithms.h>
using namespace seqan;
//![includes]

//![typedefs]
int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef Size<TGraph>::Type TSize;
//![typedefs]

    std::cout << "//![main-graph-construction]" << std::endl;
//![main-graph-construction]
    TSize numEdges = 14;
    TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
    TGraph g;
    addEdges(g, edges, numEdges);
    std::cout << g << std::endl;
//![main-graph-construction]
    std::cout << "//![main-graph-construction]" << std::endl;

//![vertex-map]
    String<char> nameMap;
    char names[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
    assignVertexMap(nameMap, g, names);
//![vertex-map]

    std::cout << "//![iterate-dfs]" << std::endl;
//![iterate-dfs]
    TVertexDescriptor start = 0;
    typedef Iterator<TGraph, DfsPreorder>::Type TDfsIterator;
    TDfsIterator dfsIt(g, start);

    std::cout << "Iterate from '" << getProperty(nameMap, start) << "' in depth-first-search ordering: ";
    while (!atEnd(dfsIt))
    {
        std::cout << getProperty(nameMap, getValue(dfsIt)) << ", ";
        goNext(dfsIt);
    }
    std::cout << std::endl;
//![iterate-dfs]
    std::cout << "//![iterate-dfs]" << std::endl;

//![connected-components]
    String<unsigned int> component;
    stronglyConnectedComponents(component, g);
//![connected-components]

    std::cout << "//![output-connected-components]" << std::endl;
//![output-connected-components]
    std::cout << "Strongly Connected Components: " << std::endl;
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
        std::cout << "Component = " << getProperty(component, getValue(it)) << std::endl;
        goNext(it);
    }
//![output-connected-components]
    std::cout << "//![output-connected-components]" << std::endl;
//![return]
    return 0;
}
//![return]
