// FRAGMENT(includes)
#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;

// FRAGMENT(main-typedefs)
int main ()
{
	typedef unsigned int TCargo;
	typedef Graph<Undirected<TCargo> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

// FRAGMENT(create-g)
	TGraph g;

// FRAGMENT(create-vertices)
	TVertexDescriptor vertBerlin = addVertex(g);
	TVertexDescriptor vertHamburg = addVertex(g);
	TVertexDescriptor vertHannover = addVertex(g);
	TVertexDescriptor vertMainz = addVertex(g);
	TVertexDescriptor vertMuenchen = addVertex(g);

// FRAGMENT(create-edges)
	addEdge(g, vertBerlin, vertHamburg, 289);
	addEdge(g, vertBerlin, vertHannover, 286);
	addEdge(g, vertBerlin, vertMainz, 573);
	addEdge(g, vertBerlin, vertMuenchen, 586);
	addEdge(g, vertHannover, vertMuenchen, 572);
	addEdge(g, vertHamburg, vertMainz, 521);

// FRAGMENT(main-graph-io)
	FILE* strmWrite = fopen("graph.dot", "w");
	write(strmWrite, g, DotDrawing());
	fclose(strmWrite);


// FRAGMENT(definition-property-map)
	typedef String<char> TCityName;
	typedef String<TCityName> TProperties;
	TProperties cityNames;
	resizeVertexMap(g, cityNames);

// FRAGMENT(enter-properties)
	assignProperty(cityNames, vertBerlin, "Berlin");
	assignProperty(cityNames, vertHamburg, "Hamburg");
	assignProperty(cityNames, vertMuenchen, "Munich");
	assignProperty(cityNames, vertMainz, "Mainz");
	assignProperty(cityNames, vertHannover, "Hannover");

// FRAGMENT(iterate-and-output-properties)
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(g);
	for(;!atEnd(itV);goNext(itV)) {
		std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << std::endl;
	}


// FRAGMENT(dijkstra-containers)
	typedef Size<TGraph>::Type TSize;
	InternalMap<TCargo> cargoMap;
	String<TVertexDescriptor> predMap;
	String<TSize> distMap;
// FRAGMENT(dijkstra)
	dijkstra(g,vertHannover,cargoMap,predMap,distMap);

// FRAGMENT(dijkstra-output)
	TVertexIterator itV2(g);
	while(!atEnd(itV2)) {
		std::cout << "Shortest path from " << property(cityNames, vertHannover) << " to " << property(cityNames, value(itV2)) << ": ";
		std::cout << property(distMap, value(itV2)) << std::endl;
		goNext(itV2);
	}

	return 0;
}
