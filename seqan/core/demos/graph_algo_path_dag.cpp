///A tutorial about a shortest path search in a directed acyclic graph.
#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() {
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
///Graph creation: 10 directed edges (0,2), (0,1), ...
	TSize numEdges = 10;
	TVertexDescriptor edges[] = {0,2, 0,1, 1,3, 1,2, 2,5, 2,4, 2,3, 3,5, 3,4, 4,5};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
///One external property map: Weight map
	int weights[] =             {3,   5,   6,   2,   2,   4,   7,   1,   -1,  -2};
	String<int> weightMap;
	assignEdgeMap(g,weightMap, weights);
///Out-parameters: Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;
///DAG-Shortest path from vertex 1
	dagShortestPath(g,1,weightMap,predMap,distMap);
///Console Output
	std::cout << "Single-Source Shortest Paths in DAG: " << std::endl;
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 1 to " << getValue(it) << ": ";
		_printPath(g,predMap,(TVertexDescriptor) 1, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << std::endl;
		goNext(it);
	}
	return 0;
}
