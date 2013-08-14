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
//       notice, this list of conditions and the following disclaimer.
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

#ifndef SEQAN_HEADER_GRAPH_BASE_H
#define SEQAN_HEADER_GRAPH_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#EdgeDescriptor
 * @brief Type of an object that represents an edge descriptor.
 *
 * @signature EdgeDescriptor<T>::Type;
 *
 * @tparam The @link Graph @endlink type to query.
 *
 * @return Type The resulting edge descriptor type.
 *
 * The edge descriptor is a unique handle to a given edge in a graph.  It is used in various graph functions, e.g. to
 * remove edges, to assign cargo to an edge or to get the end points of an edge.  It is also used to attache properties
 * to edges via property maps.
 *
 * @examples
 *
 * @code{.cpp}
 * EdgeDescriptor<Graph<> >::Type ed;  // eD is an edge descriptor
 * @endcode
 */

/**
.Metafunction.EdgeDescriptor
..class:Class.Graph
..cat:Graph
..summary:Type of an object that represents an edge descriptor.
..signature:EdgeDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs use a pointer to an edge stump as an edge descriptor.
..returns.param.Type:EdgeDescriptor type.
..remarks.text:The edge descriptor is a unique handle to a given edge in a graph.
It is used in various graph functions, e.g., to remove edges, to assign a cargo to an edge or to get the endpoints of an edge.
It is also used to attach properties to edges.
..example.code:EdgeDescriptor<Graph<> >::Type eD; //eD is an edge descriptor
..include:seqan/graph_types.h
*/
template<typename T>
struct EdgeDescriptor;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#Cargo
 * @brief Type for the graph's cargo.
 *
 * @signature Cargo<T>::Type;
 *
 * @tparam T The @link Graph @endlink type to query.
 *
 * @return Type The resulting cargo type.
 */

/**
.Metafunction.Cargo
..class:Class.Graph
..cat:Graph
..example.code:Cargo<Graph<Directed<int> > >::Type c; //c has type int
..include:seqan/graph_types.h
*/
template<typename T>
struct Cargo;


//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#EdgeType
 * @brief Edge type of a graph class.
 *
 * @signature EdgeType<T>::Type;
 *
 * @tparam T The @link Graph @endlink type to query.
 *
 * @return Type The resulting edge stump type that is used in the graph.
 *
 * @section Examples
 *
 * @code{.cpp}
 * EdgeType<TGraph>::Type e;  // e is an edge in a TGraph
 * @endcode
 */

/**
.Metafunction.EdgeType:
..class:Class.Graph
..cat:Graph
..summary:Edge type of a graph object.
..signature:EdgeType<T>::Type
..param.T:Type T must be a graph.
..returns.param.Type:Edge type.
..remarks.text:The specific edge stump type that is used in a graph.
..example.code:EdgeType<TGraph>::Type e; //e is an edge in TGraph
..include:seqan/graph_types.h
*/
template<typename T>
struct EdgeType;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @concept GraphOverAlphabetConcept
 * @brief A graph construted over an alphabet.
 *
 * @mfn GraphOverAlphabetConcept#Alphabet
 * @brief Return the Alphabe type of a graph over an alphabet.
 *
 * @signature Alphabet<T>::Type;
 *
 * @tparam T The @link GraphOverAlphabetConcept @endlink type to query.
 *
 * @return Type The alphabe type.
 *
 * @section Examples
 *
 * @code{.cpp}
 * Alphabet<Graph<Automaton<Dna> > >::Type c;  // c is of type Dna
 * @endcode
 */

/**
.Metafunction.Alphabet:
..class:Spec.Word Graph
..class:Spec.Automaton
..cat:Graph
..summary:Access to the Alphabet type.
..signature:Alphabet<T>::Type
..param.T:Type T must be a type that uses some kind of alphabet internally.
..returns.param.Type:Alphabet type.
..remarks.text:Type T can be for example an automaton where the alphabet type describes the domain of the transition labels.
..example.code:Alphabet<Graph<Automaton<Dna> > >::Type alph; //alph is of type Dna
..include:seqan/graph_types.h
*/
template<typename T>
struct Alphabet;


//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#EdgeIdHandler
 * @brief Type of an object that represents an IdManager.
 *
 * @signature EdgeIdHandler<T>::Type;
 *
 * @tparam T The Graph to query.
 *
 * @return Type The IdManager type.
 *
 * The exact IdManager type depends on the edge stump type.  If the edge stump has no ids then the IdManager simply
 * counts edge ids and otherwise it manages a list of free and used ids.
 */

/**
.Metafunction.EdgeIdHandler:
..class:Class.Graph
..cat:Graph
..summary:Type of an object that represents an Id Manager.
..signature:EdgeIdHandler<T>::Type
..param.T:A graph.
...type:Class.Graph
..returns.param.Type:IdManager type.
..remarks.text:The exact IdManager type depends on the edge stump.
If the edge stump is id-free the IdManager simply counts edge ids, 
otherwise it manages a list of free and used ids.
..include:seqan/graph_types.h
*/
template<typename T>
struct EdgeIdHandler;


//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#VertexIdHandler
 * @brief Type of an object that repreestns an IdManager.
 *
 * @signature VertexIdHandler<T>::Type;
 *
 * @tparam T The Graph to query.
 *
 * @return Type The IdManager type.
 */

/**
.Metafunction.VertexIdHandler:
..class:Class.Graph
..cat:Graph
..summary:Type of an object that represents an Id Manager.
..signature:VertexIdHandler<T>::Type
..param.T:A graph.
..returns.param.Type:IdManager type.
..include:seqan/graph_types.h
*/
template<typename T>
struct VertexIdHandler;


//////////////////////////////////////////////////////////////////////////////
// General Graph Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

struct WithoutEdgeId_;
typedef Tag<WithoutEdgeId_> const WithoutEdgeId;

//////////////////////////////////////////////////////////////////////////////

struct TreeTag_;
typedef Tag<TreeTag_> const TreeTag;


//////////////////////////////////////////////////////////////////////////////
// Graph Iterator Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup GraphIteratorTags Graph Iterator Tags
 * @brief Tags that can be used to get iterators on graphs.
 *
 * @tag GraphIteratorTags#VertexIterator
 * @brief Iterate over all vertices of a graph.
 *
 * @signature typedef Tag<VertexIterator_> const VertexIterator;
 *
 * @tag GraphIteratorTags#EdgeIterator
 * @brief Iterate over all edges of a graph.
 *
 * @signature typedef Tag<EdgeIterator_> const EdgeIterator;
 *
 * @tag GraphIteratorTags#OutEdgeIterator
 * @brief Iterate over all out edges of a vertex.
 *
 * @signature typedef Tag<OutEdgeIterator_> const OutEdgeIterator;
 *
 * @tag GraphIteratorTags#AdjacencyIterator
 * @brief Iterate over all adjacent vertices of a given vertex.
 *
 * @signature typedef Tag<AdjacencyIterator_> const AdjacencyIterator;
 *
 * @tag GraphIteratorTags#BfsIterator
 * @brief Iterate over all vertices of a graph in breadth-first fashion starting from a given vertex.
 *
 * @signature typedef Tag<BfsIterator_> const BfsIterator;
 *
 * @tag GraphIteratorTags#DfsPreorder
 * @brief Iterate over all vertices of a graph in depth-first fashion.
 *
 * @signature typedef Tag<DfsPreorder_> const DfsPreorder;
 */

/**
.Tag.Graph Iterator:
..cat:Graph
..summary:A specification of the iterator to traverse a graph.
..include:seqan/graph_types.h
*/

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.VertexIterator:
	Traverses all vertices of a graph.
..include:seqan/graph_types.h
*/
struct VertexIterator_;
typedef Tag<VertexIterator_> const VertexIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.EdgeIterator:
	Traverses all edges of a graph.
..include:seqan/graph_types.h
*/
struct EdgeIterator_;
typedef Tag<EdgeIterator_> const EdgeIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.OutEdgeIterator:
	Traverses all edges of a graph given a vertex.
..include:seqan/graph_types.h
*/
struct OutEdgeIterator_;
typedef Tag<OutEdgeIterator_> const OutEdgeIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.AdjacencyIterator:
	Traverses all neighbors of a graph given a vertex.
..include:seqan/graph_types.h
*/
struct AdjacencyIterator_;
typedef Tag<AdjacencyIterator_> const AdjacencyIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.BfsIterator:
	Traverses all vertices of a graph in Bfs order.
..include:seqan/graph_types.h
*/
struct BfsIterator_;
typedef Tag<BfsIterator_> const BfsIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.DfsPreorder:
	Traverses all vertices of a graph in Dfs order.
..include:seqan/graph_types.h
*/
struct DfsPreorder_;
typedef Tag<DfsPreorder_> const DfsPreorder;





//////////////////////////////////////////////////////////////////////////////
// Graph - Default edge stump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo = void, bool TList = true, bool TSource = false, bool TId = true, typename TSpec = Default>
class EdgeStump;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.EdgeStump
///.Metafunction.VertexDescriptor.class:Class.EdgeStump

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> > 
{
	typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> const> 
{
	typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};




//////////////////////////////////////////////////////////////////////////////
// Graph - Default Id Manager
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TIdType = unsigned int, typename TSpec = Default>
class IdManager;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeIdHandler.param.T.type:Class.EdgeStump
///.Metafunction.EdgeIdHandler.class:Class.EdgeStump

template<typename TCargo, bool TList, bool TSource, typename TSpec>
struct EdgeIdHandler<EdgeStump<TCargo, TList, TSource, false, TSpec> > {
	typedef IdManager<void> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
struct EdgeIdHandler<EdgeStump<TCargo, TList, TSource, true, TSpec> > {
	typedef IdManager<typename Id<EdgeStump<TCargo, TList, TSource, true, TSpec> >::Type> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct VertexIdHandler {
	typedef IdManager<> Type;
};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
