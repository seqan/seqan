///A tutorial about topological sort.
#include <iostream>
#include <seqan/graph_algorithms.h>


using namespace seqan;


int main() {
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
///Graph creation: 9 directed edges (0,3), (0,1), ...
	TSize numEdges = 9;
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
///One external property map: Vertex names	
	String<std::string> nameMap;
	std::string names[] = {"shirt", "tie", "jacket", "belt", "watch", "undershorts", "pants", "shoes", "socks"};
	assignVertexMap(g,nameMap, names);
///Out-parameter: Order of vertices
	String<TVertexDescriptor> order;
///Topological sort
	topologicalSort(g, order);
///Console output
	std::cout << "Topological sort: " << std::endl;
	typedef Iterator<String<TVertexDescriptor> >::Type TStringIterator;
	TStringIterator it = begin(order);
	TStringIterator itEnd = end(order);
	while(it != itEnd) {
		std::cout << getProperty(nameMap, getValue(it)) << ",";
		goNext(it);
	}
	std::cout << std::endl;
	return 0;
}
