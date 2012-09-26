///A tutorial about a shortest path search using Bellman-Ford.
#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() {
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
///Graph creation: 10 directed edges (0,1), (0,3), ...
	TSize numEdges = 10;
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
///One external property map: Weight map
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};
	String<unsigned int> weightMap;
	assignEdgeMap(g,weightMap, weights);
///Out-parameters: Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;
///Bellman-Ford from vertex 0
///Note: Bellman-Ford also detects negative cycles
	bool noNegativeCycle = bellmanFordAlgorithm(g,0,weightMap,predMap,distMap);
///Console Output
	std::cout << "Single-Source Shortest Paths: " << std::endl;
	std::cout << "Graph without negative cycles? " << noNegativeCycle << std::endl;
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_printPath(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << std::endl;
		goNext(it);
	}
	return 0;
}
