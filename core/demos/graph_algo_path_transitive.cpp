///A tutorial about the transitive closure algorithm.
#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() {
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
///Graph creation: 5 directed edges (3,0), (1,2), ...
	TSize numEdges = 5;
	TVertexDescriptor edges[] = {3,0, 1,2, 2,1, 1,3, 3,2};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
///Out-parameter: Closure matrix
	String<bool> closure;
///Transitive-Closure 
	transitiveClosure(g,closure);
///Console output
	TSize len = (TSize) std::sqrt((double) length(closure));
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			std::cout << getValue(closure, row*len+col) << ",";
		}
		std::cout << std::endl;
	}
	return 0;
}
