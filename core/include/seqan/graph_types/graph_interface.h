// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimerf
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_INTERFACE_H
#define SEQAN_HEADER_GRAPH_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{

// Default directed graph
template<typename TCargo = void, typename TSpec = Default>
struct Directed;

// Default undirected graph
template<typename TCargo = void, typename TSpec = Default>
struct Undirected;

// Default Tree
template<typename TCargo = void, typename TSpec = Default>
struct Tree;

// Default Automaton
template<typename TAlphabet = char, typename TCargo = void, typename TSpec = Default>
struct Automaton;

// Default Hmm
template<typename TAlphabet = Dna, typename TCargo = double, typename TSpec = Default>
struct Hmm;


//////////////////////////////////////////////////////////////////////////////
// Graph
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class Graph
 * @extends ContainerConcept
 * @headerfile <seqan/graph_types.h>
 * @brief Graph class.
 *
 * @signature template <[typename TSpec]>
 *            class Graph;
 *
 * @tparam TSpec The specializing type.  Default: Directed&lt;&gt;.
 *
 * @section Examples
 *
 * This is an example for Dijktra's algorithm on a directed graph with an external property map.  The property map
 * labels the edges with weights.  The xample only outputs distances, not th edetails of the paths.
 *
 * @include demos/graph/graph_algo_dijkstra.cpp
 *
 * The output of the distances is as follows:
 *
 * @include demos/graph/graph_algo_dijkstra.cpp.stdout
 */

/**
.Class.Graph:
..cat:Graph
..summary:Generic graph.
..signature:Graph<TSpec>
..param.TSpec:The specializing type determines the kind of graph, e.g., directed, undirected, tree, or automaton.
...remarks:The default Graph<> corresponds to a directed graph.
...default:Directed<>
..include:seqan/graph_types.h
..example:This is an example for Dijkstra's algorithm on a directed graph with an external property map. The property map adds weights to the edges. The example only outputs distances, not the details of the paths.
...file:demos/graph/graph_algo_dijkstra.cpp
...text:The output of the distances is as follows:
...output:Distance from 0 to 0: 0
Distance from 0 to 1: 8
Distance from 0 to 2: 9
Distance from 0 to 3: 5
Distance from 0 to 4: 7

*/
template<typename TSpec = Directed<> >
class Graph;

//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Graph
///.Metafunction.Spec.class:Class.Graph

template<typename TSpec>
struct Spec<Graph<TSpec> > 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
struct Spec<Graph<TSpec> const>
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeDescriptor.param.T.type:Class.Graph
///.Metafunction.EdgeDescriptor.class:Class.Graph

template<typename TSpec>
struct EdgeDescriptor<Graph<TSpec> > 
{
	typedef typename EdgeType<Graph<TSpec> >::Type* Type;
};

template<typename TSpec>
struct EdgeDescriptor<Graph<TSpec> const>
{
	typedef typename EdgeType<Graph<TSpec> const>::Type* Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.Graph
///.Metafunction.VertexDescriptor.class:Class.Graph

template<typename TSpec>
struct VertexDescriptor<Graph<TSpec> > 
{
	typedef typename Id<Graph<TSpec> >::Type Type;
};

template<typename TSpec>
struct VertexDescriptor<Graph<TSpec> const>
{
	typedef typename Id<Graph<TSpec> >::Type Type;
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeType.param.T.type:Class.Graph
///.Metafunction.EdgeType.class:Class.Graph

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, true, false, true, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, true, false, true, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Directed<TCargo, WithoutEdgeId> > > {
	typedef EdgeStump<TCargo, true, false, false, WithoutEdgeId> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Directed<TCargo, WithoutEdgeId> > const> {
	typedef EdgeStump<TCargo, true, false, false, WithoutEdgeId> const Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Tree<TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, true, false, false, TreeTag> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Tree<TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, true, false, false, TreeTag> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, true, true, true, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, true, true, true, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Undirected<TCargo, WithoutEdgeId> > > {
	typedef EdgeStump<TCargo, true, true, false, WithoutEdgeId> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Undirected<TCargo, WithoutEdgeId> > const> {
	typedef EdgeStump<TCargo, true, true, false, WithoutEdgeId> const Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, false, false, true, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, false, false, true, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, WithoutEdgeId> > > {
	typedef EdgeStump<TCargo, false, false, false, WithoutEdgeId> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, WithoutEdgeId> > const> {
	typedef EdgeStump<TCargo, false, false, false, WithoutEdgeId> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Hmm<TAlphabet, TCargo, TSpec> > const> {
	typedef typename EdgeType<Graph<Directed<TCargo, TSpec> > const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Hmm<TAlphabet, TCargo, TSpec> > > {
	typedef typename EdgeType<Graph<Directed<TCargo, TSpec> > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.Graph
///.Metafunction.Cargo.class:Class.Graph

template<typename TSpec>
struct Cargo<Graph<TSpec> > {
	typedef typename Cargo<typename EdgeType<Graph<TSpec> >::Type>::Type Type;
};


template<typename TSpec>
struct Cargo<Graph<TSpec> const> {
	typedef typename Cargo<typename EdgeType<Graph<TSpec> const>::Type>::Type Type;
};



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeIdHandler.param.T.type:Class.Graph
///.Metafunction.EdgeIdHandler.class:Class.Graph

template<typename TSpec>
struct EdgeIdHandler<Graph<TSpec> const> {
	typedef typename EdgeIdHandler<typename EdgeType<Graph<TSpec> const>::Type>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
struct EdgeIdHandler<Graph<TSpec> > {
	typedef typename EdgeIdHandler<typename EdgeType<Graph<TSpec> >::Type>::Type Type;
};




//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Alphabet.param.T.type:Class.Graph
///.Metafunction.Alphabet.class:Class.Graph

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Automaton<TAlphabet, TCargo, TSpec> > > {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Automaton<TAlphabet, TCargo, TSpec> > const> {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Hmm<TAlphabet, TCargo, TSpec> > > {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Hmm<TAlphabet, TCargo, TSpec> > const> {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////
// Generic Graph Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Graph#transpose
 * @brief Transposes a graph, either in place or from source to dest.
 *
 * @signature void transpose(source, dest);
 * @signature void transpose(g);
 *
 * @param[in]     source The input @link Graph @endlink.
 * @param[in]     dest   The output @link Graph @endlink.
 * @param[in,out] g      The @link Graph @endlink to transpose.
 */

/*!
 * @fn Graph#numEdges
 * @brief Returns the number of edges in a graph.
 *
 * @signature TSize numEdges(g);
 *
 * @param[in] g The @link Graph @endlink to query for its number of edges.
 *
 * @return TSize The number of edges.  TSize is the size type of g.
 */

/*!
 * @fn Graph#numVertices
 * @brief Returns the number of vertices in a graph.
 *
 * @signature TSize numVertices(g);
 *
 * @param[in] g The @link Graph @endlink to query for its number of vertices.
 *
 * @return TSize The number of vertices.  TSize is the size type of g.
 */

/*!
 * @fn Graph#empty
 * @brief Returns whether there are no vertices and edges in the graph or not.
 *
 * @signature bool empty(g);
 *
 * @param[in] g true if g is empty and valse if not.
 */

/*!
 * @fn Graph#clearEdges
 * @brief Removes all edges from a graph.
 *
 * @signature void clearEdges(g);
 *
 * @param[in,out] g The graph to remove the edges from.
 */

/*!
 * @fn Graph#clear
 * @brief Remove all edges and vertices from a graph.
 *
 * @signature void clear(g);
 *
 * @param[in,out] g The graph to remove edges and vertices from.
 */

/*!
 * @fn Graph#outDegree
 * @brief Return out degree of a vertex.
 *
 * @signature TSize outDegree(g, v);
 *
 * @param[in] g The Graph to query.
 * @param[in] v The descriptor of the vertex to query for its out degree.
 *
 * @return TSize The number of out edges from vertex v.
 */

/*!
 * @fn Graph#inDegree
 * @brief Return in degree of a vertex.
 *
 * @signature TSize inDegree(g, v);
 *
 * @param[in] g The Graph to query.
 * @param[in] v The descriptor of the vertex to query for its in degree.
 *
 * @return TSize The number of in edges to vertex v.
 */

/*!
 * @fn Graph#degree
 * @brief Return degree of a vertex.
 *
 * @signature TSize degree(g, v);
 *
 * @param[in] g The Graph to query.
 * @param[in] v The descriptor of the vertex to query for its degree.
 *
 * @return TSize The number of edges adjacent to vertex v.
 */

/*!
 * @fn Graph#addVertex
 * @brief Add a vertex to a graph.
 *
 * @signature TVertextDescriptor addVertex(g);
 *
 * @param[in,out] g The Graph to add a vertex to.
 *
 * @return TVertexDescriptor The descriptor of the added vertex.
 */

/*!
 * @fn Graph#removeVertex
 * @brief Remove a vertex from a Graph.
 *
 * @signature void removeVertex(g, v);
 *
 * @param[in,out] g The Graph to remove the vertex from.
 * @param[in]     v The descriptor of the vertex to remove.
 */

// TODO(holtgrew): Specialize for automatons.x
/*!
 * @fn Graph#addEdge
 * @brief Adds a new edge to a graph, either with or without cargo.
 *
 * @signatures TEdgeDescriptor addEdge(g, source, target, label[, cargo]);
 * @signatures TEdgeDescriptor addEdge(g, source, target, cargo);
 *
 * @param[in,out] g      The Graph to add the edge to.
 * @param[in]     source Descriptor of the source vertex.
 * @param[in]     source Descriptor of the target vertex.
 * @param[in]     label  Label of the edge, of alphabet type of the Automaton.
 * @param[in]     cargo  Cargo object for the edge.
 *
 * @return TEdgeDescriptor Descriptor of the added edge.
 *
 * @section Remarks
 *
 * For Automaton objects, a label is required but the label can only be given for Automatons.
 */

/*!
 * @fn Graph#removeEdge
 * @brief Removes an edge from the graph.  For automatons, a label is required.
 *
 * @signature void removeEdge(g, source, target[, label]);
 * @signature void removeEdge(g, e);
 *
 * @param[in,out] g      The Graph to remove the edge from.
 * @param[in]     source Descriptor of the source vertex.
 * @param[in]     target Descriptor of the target vertex.
 * @param[in]     e      Descriptor of the edge to remove.
 */

/*!
 * @fn Graph#removeInEdges
 * @brief Removes the incoming edges of a given vertex.
 *
 * @signature void removeOutEdges(g, v);
 *
 * @param[in] g The Graph to remove the edges from.
 * @param[in] v The descriptor of the vertex to remove incoming edges from.
 */

/*!
 * @fn Graph#removeOutEdges
 * @brief Removes the outgoing edges of a given vertex.
 *
 * @signature void removeOutEdges(g, v);
 *
 * @param[in] g The Graph to remove the edges from.
 * @param[in] v The descriptor of the vertex to remove outgoing edges from.
 */

/*!
 * @fn Graph#targetVertex
 * @brief Returns the target vertex of an edge.
 *
 * @signature TVertexDescriptor targetVertex(g, e);
 * @signature TVertexDescriptor targetVertex(it);
 *
 * @param[in] g  The Graph the edge is in.
 * @param[in] e  The descriptor of the edge to remove.
 * @param[in] it An edge iterator.
 *
 * @section Remarks
 *
 * In a tree, the target vertex is always the child.  In an undirected graph, the larger vertex descriptor of the two
 * end points i the target.  For an out-edge iterator, the target is always the OutEdgeIterator has <b>not</b> been
 * initialized with.
 */

/*!
 * @fn Graph#sourceVertex
 * @brief Returns the source vertex of an edge.
 *
 * @signature TVertexDescriptor sourceVertex(g, e);
 * @signature TVertexDescriptor sourceVertex(it);
 *
 * @param[in] g  The Graph the edge is in.
 * @param[in] e  The descriptor of the edge to remove.
 * @param[in] it An edge iterator.
 *
 * @section Remarks
 *
 * In a tree, the source vertex is always the child.  In an undirected graph, the larger vertex descriptor of the two
 * end points i the source.  For an out-edge iterator, the source is always the OutEdgeIterator has been initialized
 * with.
 */

/*!
 * @fn Graph#getAdjacencyMatrix
 * @brief Build an adjacency matrix representation of the graph.
 *
 * @signature void getAdjacencyMatrix(g, matrix);
 *
 * @param[in]  g      The Graph to compute ajacency matrix for.
 * @param[out] matrix The Matrix to fill with the adjacency matrix values.
 */

/*!
 * @fn Graph#findEdge
 * @brief Finds an edge.
 *
 * @signature TEdgeDescriptor findEdge(g, v, c);
 * @signature TEdgeDescriptor findEdge(g, v, w);
 *
 * @param[in] g The Graph to query.
 * @param[in] v The descriptor of the source vertex.
 * @param[in] c An edge label.
 * @param[in] w The descriptor of the target descriptor.
 *
 * @return TEdgeDescriptor Edge descriptor of the found edge.  0 if not present.  NB: in automatons, there is
 *                         always a valid edge descriptor but the target may be nil.
 *
 * @section Remarks
 *
 * In an automaton, an edge is uniquely defined by a vertex and a label.  In all other graphs, two adjacent vertices
 * uniquely define an edge.  If tehre are multiple edges between two vertices then the behaviour is undefined.
 */

/*!
 * @fn Graph#clearVertices
 * @brief Removes all vertices from a graph.
 *
 * @signature void clearVertices(g);
 *
 * @param[in,out] g The graph to remove the vertices from.
 */

/*!
 * @fn Graph#getNil
 * @brief Returning a "nil" value for graphs.
 *
 * Usefulf or various graph algorithms, e.g. missing predecessors or vertices that have not been visited.
 *
 * @signature T getNil(ptr);
 * @signature T getNil<T>();
 *
 * @param[in] ptr Pointer to T to select the type.
 *
 * @return T Pseudo nil value for type T.
 */

/**
.Function.getNil
..class:Class.Graph
..cat:Graph
..summary:Utility function returning a value that represents nil.
Useful for various graph algorithms, e.g., missing predecessors, vertices that have not been visited, etc.
..signature:getNil<T>()
..returns:Pseudo nil value for type T.
..include:seqan/graph_types.h
*/
template <typename T>
inline T
getNil(T *)
{
	return ~0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T
getNil()
{
	T * _tag = 0;
	return getNil(_tag);
}


//////////////////////////////////////////////////////////////////////////////
// Purely internal!!! Never compare to _getInfinity()!!!.
// Just returns a very large value.
//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T
_getInfinity()
{
	T * _tag = 0;
	return supremumValueImpl(_tag);
}

//////////////////////////////////////////////////////////////////////////////

template <>
inline double
_getInfinity()
{
	return 1000000000;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TWeightMap>
inline typename Value<TWeightMap>::Type
_getInfinityDistance(TWeightMap const&)
{
	// We need to divide by 2 because of addition in some graph algorithms: infinity + something
	return (_getInfinity<typename Value<TWeightMap>::Type>()/2);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T
_getInfinityDistance()
{
	return (_getInfinity<T>() / 2);
}

//////////////////////////////////////////////////////////////////////////////

// Simple _getId function to get the id for a vertex descriptor which is the id!
template<typename TId>
inline TId
_getId(TId const id)
{
	SEQAN_CHECKPOINT
	return id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor>
inline void
_createVertices(Graph<TSpec>& g,
				TVertexDescriptor const maxId) 
{
		// Create missing vertices
		while (maxId >= getIdUpperBound(g.data_id_managerV)) addVertex(g);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Graph#addEdges
 * @brief Shortcut to add multiple edges at once; vertices are created implicitely.
 *
 * @signature void addEdges(graph, edges, size);
 *
 * @param[in,out] graph The Graph to add the edges to.
 * @param[in]     edges An array of vertex descriptions.
 * @param[in]     size  Size of the array.  Must be a even.
 *
 * It is assumed that the edges in <tt>edges</tt> are stored as an array of vertex ids: <tt>source1, target1, source2,
 * target2, ...</tt>.  For a tree, the root must be the first vertex in this array and the enumeration is <tt>parent,
 * child, parent, child</tt>.
 */

/**
.Function.addEdges
..class:Class.Graph
..cat:Graph
..summary:Shortcut to add multiple edges at once.
Creates vertices implicitly.
..signature:addEdges(g, edges, size)
..param.g:A graph.
...type:Class.Graph
..param.edges:An array of vertex descriptors. It is assumed that the
edges are stored in the following way: Source1, Target1, Source2, Target2, Source3, ...
For a tree the root must be the first vertex in this array and the enumeration is Parent, Child, Parent, Child, ...
...type:Metafunction.VertexDescriptor
..param.size:Size of the array. Must be a multiple of 2.
...type:Metafunction.Size
..returns:void
..see:Function.addEdge
..include:seqan/graph_types.h
*/
template<typename TSpec, typename TEdgeArray, typename TSize>
inline void
addEdges(Graph<TSpec>& dest,
		 TEdgeArray const & edges,
		 TSize const size) 
{
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	for(TSize i=0;i<size;++i) {
		TVertexDescriptor source = edges[2*i];
		TVertexDescriptor target = edges[2*i+1];
		// Create missing vertices
		if (source>target) _createVertices(dest,source);
		else _createVertices(dest,target);
		// Add edge
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, source));
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, target));
		addEdge(dest, source, target);
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TSpec>
inline TStream &
operator << (TStream & target, 
			 Graph<TSpec> const& source)
{
	write(target, source);
	return target;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
