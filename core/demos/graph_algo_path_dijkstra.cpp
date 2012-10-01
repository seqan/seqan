///A tutorial about the dijkstra's algorithm, once using an external map and once using an internal map.
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
///One external property map: Weight map	
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};
	String<unsigned int> weightMap;
	assignEdgeMap(g,weightMap, weights);
///Out-parameters: Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;
///Dijkstra from vertex 0
	dijkstra(g,0,weightMap,predMap,distMap);
///Console Output
	std::cout << "Single-Source Shortest Paths: " << std::endl;
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Path from 0 to " << getValue(it) << ": ";
		_printPath(g,predMap,(TVertexDescriptor) 0, getValue(it));
		std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")" << std::endl;
		goNext(it);
	}
///We can achieve the same thing using an internal map that is edge cargos.
	typedef unsigned int TEdgeCargo;
	typedef Directed<TEdgeCargo> TEdges;
	typedef Graph<TEdges> TCargoGraph;
///Graph creation
	TCargoGraph cargo_g;
	addEdges(cargo_g, edges, numEdges);
///One internal property map: Weight map	
	InternalMap<TEdgeCargo> intMap;
	assignEdgeMap(cargo_g, intMap, weights);
///Out parameters of Dijkstra: Predecessor map and distance map
	clear(predMap);
	clear(distMap);
///Dijkstra from vertex 0 using an internal map
	dijkstra(cargo_g,0,intMap,predMap,distMap);
///Console Output
	std::cout << "Single-Source Shortest Paths: " << std::endl;
	typedef Iterator<TCargoGraph, VertexIterator>::Type TCargoVertexIterator;
	TCargoVertexIterator itC(cargo_g);
	while(!atEnd(itC)) {
		std::cout << "Path from 0 to " << getValue(itC) << ": ";
		_printPath(cargo_g,predMap,(TVertexDescriptor) 0, getValue(itC));
		std::cout << " (Distance: " << getProperty(distMap, getValue(itC)) << ")" << std::endl;
		goNext(itC);
	}
	return 0;
}
