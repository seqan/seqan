///A tutorial about depth-first search.
#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;


int main() 
{
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
///Graph creation: 8 directed edges (0,3), (0,1), ...
	TSize numEdges = 8;
	TVertexDescriptor edges[] = {0,3, 0,1, 1,4, 2,4, 2,5, 3,1, 4,3, 5,5};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
///One external property map: Vertex names
	char names[] = {'u', 'v', 'w', 'x', 'y', 'z'};
	String<char> nameMap;
	assignVertexMap(g,nameMap, names);
///Out-parameters: Predecessor and discovery maps
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;
///Depth-frist search
	depthFirstSearch(g, predMap, discoveryTimeMap, finishingTimeMap);
///Console output
	std::cout << "Depth-First search: " << std::endl;
	typedef Iterator<Graph<>, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		std::cout << "Discovery time = " << getProperty(discoveryTimeMap, getValue(it)) << ",";
		std::cout << "Finishing time = " << getProperty(finishingTimeMap, getValue(it)) << ",";
		typedef Value<String<unsigned int> >::Type TPredVal;
		TPredVal pre = getProperty(predMap, getValue(it));
		if (pre != getNil<TVertexDescriptor>()) {
			std::cout << "Predecessor = " << getProperty(nameMap, pre) << std::endl;
		} else {
			std::cout << "Predecessor = nil" << std::endl;
		}
		goNext(it);
	}
	return 0;
}
