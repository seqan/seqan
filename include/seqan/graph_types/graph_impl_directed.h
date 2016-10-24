// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_GRAPH_IMPL_DIRECTED_H
#define SEQAN_HEADER_GRAPH_IMPL_DIRECTED_H

// TODO(holtgrew): The graph uses linked lists for storing edges. Thus, the graphs are not guaranteed to have good cache locality.

namespace seqan
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Directed
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class DirectedGraph
 * @extends Graph
 * @brief A directed graph that stores the edges in an adjacency list.
 *
 * <img src="directedGraph.png" title="A directed graph." />
 *
 * @signature template <[typename TCargo[, typename TSpec]]>
 *            class Graph<Directed<TCargo, TSpec> >;
 *
 * @tparam TCargo The cargo type that can be attached to the edges.  Default: <tt>void</tt>.
 * @tparam TSpec  The specializing type.  Default: <tt>Default</tt>.  Use <tt>WithoutEdgeId</tt> here to omit edge
 *                ids.  NB: if edges do not store ids then external property maps do not work.
 */

template <typename TCargo, typename TSpec>
class Graph<Directed<TCargo, TSpec> >
{
    public:
        typedef typename VertexIdHandler<Graph>::Type TVertexIdManager_;
        typedef typename EdgeIdHandler<Graph>::Type TEdgeIdManager_;
        typedef typename EdgeType<Graph>::Type TEdgeStump_;
        typedef Allocator<SinglePool<sizeof(TEdgeStump_)> > TAllocator_;

        String<TEdgeStump_*> data_vertex;            // Pointers to EdgeStump lists
        TVertexIdManager_ data_id_managerV;
        TEdgeIdManager_ data_id_managerE;
        TAllocator_ data_allocator;

//____________________________________________________________________________

        Graph() {}

        ~Graph()
        {
            clear(*this);
        }

        Graph(Graph const & _other) :
            data_allocator(_other.data_allocator)
        {
            _copyGraph(*this, _other);
        }

        Graph const & operator = (Graph const & _other)
        {
            if (this == &_other)
                return *this;

            clear(*this);
            data_allocator = _other.data_allocator;
            _copyGraph(*this, _other);
            return *this;
        }
};

template <typename TCargo>
class Graph<Directed<TCargo, WithSourceId> >
{
    public:
        typedef typename VertexIdHandler<Graph>::Type TVertexIdManager_;
        typedef typename EdgeIdHandler<Graph>::Type TEdgeIdManager_;
        typedef typename EdgeType<Graph>::Type TEdgeStump_;
        typedef Allocator<SinglePool<sizeof(TEdgeStump_)> > TAllocator_;

        String<TEdgeStump_*> data_vertex;            // Pointers to EdgeStump lists
        String<TEdgeStump_*> data_vertex_in;         
        TVertexIdManager_ data_id_managerV;
        TEdgeIdManager_ data_id_managerE;
        TAllocator_ data_allocator;

//____________________________________________________________________________

        Graph() {}

        ~Graph()
        {
            clear(*this);
        }

        Graph(Graph const & _other) :
            data_allocator(_other.data_allocator)
        {
            _copyGraph(*this, _other);
        }

        Graph const & operator = (Graph const & _other)
        {
            if (this == &_other)
                return *this;

            clear(*this);
            data_allocator = _other.data_allocator;
            _copyGraph(*this, _other);
            return *this;
        }
};

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Directed<TCargo, TSpec> > >::Type*> &
_getVertexString(Graph<Directed<TCargo, TSpec> > const & g)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    return const_cast<String<TEdgeStump*>&>(g.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline typename VertexIdHandler<Graph<Directed<TCargo, TSpec> > >::Type &
_getVertexIdManager(Graph<Directed<TCargo, TSpec> > const & g)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename VertexIdHandler<TGraph>::Type TVertexIdManager;
    return const_cast<TVertexIdManager&>(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline typename EdgeIdHandler<Graph<Directed<TCargo, TSpec> > >::Type &
_getEdgeIdManager(Graph<Directed<TCargo, TSpec> > const & g)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeIdHandler<TGraph>::Type TEdgeIdManager;
    return const_cast<TEdgeIdManager&>(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo>
inline void
_copyGraph(Graph<Directed<TCargo, WithSourceId> > & dest,
           Graph<Directed<TCargo, WithSourceId> > const & source,
           bool const & transpose)
{
    typedef Graph<Directed<TCargo, WithSourceId> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
    clear(dest);
    resize(dest.data_vertex, length(source.data_vertex));
    resize(dest.data_vertex_in, length(source.data_vertex_in));
    TIter itInit = begin(dest.data_vertex, Standard());
    TIter itInitEnd = end(dest.data_vertex, Standard());
    for (;itInit != itInitEnd; ++itInit)
        *itInit = (TEdgeStump*) 0;

    TIter itInitIn = begin(dest.data_vertex_in, Standard());
    TIter itInitInEnd = end(dest.data_vertex_in, Standard());
    for (;itInitIn != itInitInEnd; ++itInitIn)
        *itInitIn = (TEdgeStump*) 0;

    TIterConst it = begin(source.data_vertex, Standard());
    TIterConst itEnd = end(source.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for (;it != itEnd; ++it, ++pos)
    {
        TEdgeStump* current = *it;
        TVertexDescriptor sourceVertex = pos;
        while (current != (TEdgeStump*) 0)
        {
            TVertexDescriptor targetVertex = current->data_target;
            // Create missing vertices
            if (sourceVertex>targetVertex)
                _createVertices(dest,sourceVertex);
            else
                _createVertices(dest,targetVertex);

            // Add edge
            TEdgeDescriptor e;
            if (!transpose)
                e = addEdge(dest, sourceVertex, targetVertex);
            else
                e = addEdge(dest, targetVertex, sourceVertex);

            _assignId(e, _getId(current));
            assignCargo(e, getCargo(current));
            current = getNextT(current);
        }
    }
    dest.data_id_managerV = source.data_id_managerV;
    dest.data_id_managerE = source.data_id_managerE;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Directed<TCargo, TSpec> > & dest,
           Graph<Directed<TCargo, TSpec> > const & source,
           bool const & transpose)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
    clear(dest);
    resize(dest.data_vertex, length(source.data_vertex));
    TIter itInit = begin(dest.data_vertex, Standard());
    TIter itInitEnd = end(dest.data_vertex, Standard());
    for (;itInit != itInitEnd; ++itInit)
    {
        *itInit = (TEdgeStump*) 0;
    }
    TIterConst it = begin(source.data_vertex, Standard());
    TIterConst itEnd = end(source.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for (;it != itEnd; ++it, ++pos)
    {
        TEdgeStump* current = *it;
        TVertexDescriptor sourceVertex = pos;
        while (current != (TEdgeStump*) 0)
        {
            TVertexDescriptor targetVertex = current->data_target;
            // Create missing vertices
            if (sourceVertex>targetVertex)
                _createVertices(dest,sourceVertex);
            else
                _createVertices(dest,targetVertex);
            // Add edge
            TEdgeDescriptor e;
            if (!transpose)
                e = addEdge(dest, sourceVertex, targetVertex);
            else
                e = addEdge(dest, targetVertex, sourceVertex);

            _assignId(e, _getId(current));
            assignCargo(e, getCargo(current));
            current = getNextT(current);
        }
    }
    dest.data_id_managerV = source.data_id_managerV;
    dest.data_id_managerE = source.data_id_managerE;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Directed<TCargo, TSpec> > & dest,
           Graph<Directed<TCargo, TSpec> > const & source)
{
    _copyGraph(dest, source, false);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline void
transpose(Graph<Directed<TCargo, TSpec> > const & source,
          Graph<Directed<TCargo, TSpec> > & dest)
{
    _copyGraph(dest, source, true);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline void
transpose(Graph<Directed<TCargo, TSpec> > & g)
{
    Graph<Directed<TCargo, TSpec> > dest;
    _copyGraph(dest, g, true);
    g = dest;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type
numEdges(Graph<Directed<TCargo, TSpec> > const & g)
{
    return idCount(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type
numVertices(Graph<Directed<TCargo, TSpec> > const & g)
{
    return idCount(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline bool
empty(Graph<Directed<TCargo, TSpec> > const & g)
{
    return !(idCount(g.data_id_managerV));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Directed<TCargo, TSpec> > & g)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
    TIter it = begin(g.data_vertex, Standard());
    TIter itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for (; it != itEnd; ++it, ++pos)
    {
        if (*it != (TEdgeStump*) 0)
            removeOutEdges(g, pos);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo>
inline void
clearVertices(Graph<Directed<TCargo, WithSourceId> > & g)
{
    clearEdges(g);
    releaseAll(g.data_id_managerV);
    clear(g.data_vertex);
    clear(g.data_vertex_in);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Directed<TCargo, TSpec> > & g)
{
    clearEdges(g);
    releaseAll(g.data_id_managerV);
    clear(g.data_vertex);
}

//////////////////////////////////////////////////////////////////////////////


template <typename TCargo, typename TSpec>
inline void
clear(Graph<Directed<TCargo, TSpec> > & g)
{
    clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type
outDegree(Graph<Directed<TCargo, TSpec> > const & g,
          TVertexDescriptor const & vertex)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TGraph>::Type TSize;
    TSize count = 0;
    TEdgeStump* current = getValue(g.data_vertex, vertex);
    while (current != 0)
    {
        current = getNextT(current);
        ++count;
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TVertexDescriptor>
inline typename Size<Graph<Directed<TCargo, WithSourceId> > >::Type
inDegree(Graph<Directed<TCargo, WithSourceId> > const & g,
         TVertexDescriptor const & vertex)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));

    typedef Graph<Directed<TCargo, WithSourceId> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TGraph>::Type TSize;
    TEdgeStump* current = g.data_vertex_in[vertex];
    TSize count = 0;
    while (current != 0)
    {
        current = getNextS(current);
        ++count;
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type
inDegree(Graph<Directed<TCargo, TSpec> > const & g,
         TVertexDescriptor const & vertex)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());

    TSize count=0;
    for (; it != itEnd; ++it)
    {
        TEdgeStump* current = *it;
        while (current!=0)
        {
            if (static_cast<TVertexDescriptor>(getTarget(current)) == vertex)
                ++count;

            current = getNextT(current);
        }
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type
degree(Graph<Directed<TCargo, TSpec> > const & g,
       TVertexDescriptor const & vertex)
{
    return inDegree(g,vertex) + outDegree(g,vertex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo>
inline typename VertexDescriptor<Graph<Directed<TCargo, WithSourceId> > >::Type
addVertex(Graph<Directed<TCargo, WithSourceId> > & g)
{
    typedef Graph<Directed<TCargo, WithSourceId> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    TVertexDescriptor vd = obtainId(g.data_id_managerV);
    if (vd == length(g.data_vertex))
    {
        appendValue(g.data_vertex, (TEdgeStump*) 0);
        appendValue(g.data_vertex_in, (TEdgeStump*) 0);
    }
    else
    {
        g.data_vertex[vd] = (TEdgeStump*) 0;
        g.data_vertex_in[vd] = (TEdgeStump*) 0;
    }
    return vd;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Directed<TCargo, TSpec> > >::Type
addVertex(Graph<Directed<TCargo, TSpec> > & g)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    TVertexDescriptor vd = obtainId(g.data_id_managerV);
    if (vd == length(g.data_vertex))
        appendValue(g.data_vertex, (TEdgeStump*) 0);
    else
        g.data_vertex[vd] = (TEdgeStump*) 0;

    return vd;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeVertex(Graph<Directed<TCargo, TSpec> > & g,
             TVertexDescriptor const & v)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));

    removeOutEdges(g, v); // Remove all outgoing edges
    removeInEdges(g, v); // Remove all incoming edges
    releaseId(g.data_id_managerV, v); // Release id
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Directed<TCargo, WithSourceId> > >::Type
addEdge(Graph<Directed<TCargo, WithSourceId> > & g,
        TVertexDescriptor const & source,
        TVertexDescriptor const & target)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, source));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, target));

    typedef Graph<Directed<TCargo, WithSourceId> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Id<TGraph>::Type TId;

    TEdgeStump* edge_ptr;
    allocate(g.data_allocator, edge_ptr, 1);
    valueConstruct(edge_ptr);
    assignSource(edge_ptr, source);
    assignTarget(edge_ptr, target);
    assignNextS(edge_ptr, (TEdgeStump*) 0);
    assignNextT(edge_ptr, (TEdgeStump*) 0);
    TId id = obtainId(g.data_id_managerE);
    _assignId(edge_ptr, id);

    if (g.data_vertex[source] != 0)
        assignNextT(edge_ptr, getValue(g.data_vertex, source));

    if (g.data_vertex_in[target] != 0)
        assignNextS(edge_ptr, getValue(g.data_vertex_in, target));

    g.data_vertex[source] = edge_ptr;
    g.data_vertex_in[target] = edge_ptr;
    return edge_ptr;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Directed<TCargo, TSpec> > >::Type
addEdge(Graph<Directed<TCargo, TSpec> > & g,
        TVertexDescriptor const & source,
        TVertexDescriptor const & target)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, source));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, target));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Id<TGraph>::Type TId;

    TEdgeStump* edge_ptr;
    allocate(g.data_allocator, edge_ptr, 1);
    valueConstruct(edge_ptr);
    assignTarget(edge_ptr, target);
    assignNextT(edge_ptr, (TEdgeStump*) 0);
    TId id = obtainId(g.data_id_managerE);
    _assignId(edge_ptr, id);

    if (g.data_vertex[source] != 0)
        assignNextT(edge_ptr, getValue(g.data_vertex, source));

    value(g.data_vertex, source) = edge_ptr;
    return edge_ptr;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Directed<TCargo, TSpec> > >::Type
addEdge(Graph<Directed<TCargo, TSpec> > & g,
        TVertexDescriptor const & source,
        TVertexDescriptor const & target,
        TCargo const cargo)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    TEdgeDescriptor e = addEdge(g, source, target);
    assignCargo(e, cargo);
    return e;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TVertexDescriptor>
inline void
removeEdge(Graph<Directed<TCargo, WithSourceId> > & g,
           TVertexDescriptor const & source,
           TVertexDescriptor const & target)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, source));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, target));

    typedef Graph<Directed<TCargo, WithSourceId> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    // Find edge and predecessor
    TEdgeStump* pred = 0;
    TEdgeStump* current = g.data_vertex[source];
    while (current != (TEdgeStump*) 0)
    {
        if (static_cast<TVertexDescriptor>(getTarget(current)) == target)
            break;

        pred = current;
        current = getNextT(current);
    }

    // Not found?
    if (current == (TEdgeStump*) 0)
        return;

    // Relink the next pointer of predecessor
    if (pred != (TEdgeStump*) 0)
    {
        assignNextT(pred, getNextT(current));
    }
    else
    {
        g.data_vertex[source] = getNextT(current);
        g.data_vertex_in[target] = getNextT(current);
    } 

    // Deallocate
    releaseId(g.data_id_managerE, _getId(current));
    valueDestruct(current);
    deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeEdge(Graph<Directed<TCargo, TSpec> > & g,
           TVertexDescriptor const & source,
           TVertexDescriptor const & target)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, source));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, target));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    // Find edge and predecessor
    TEdgeStump* pred = 0;
    TEdgeStump* current = g.data_vertex[source];
    while (current != (TEdgeStump*) 0)
    {
        if (static_cast<TVertexDescriptor>(getTarget(current)) == target)
            break;

        pred = current;
        current = getNextT(current);
    }

    // Not found?
    if (current == (TEdgeStump*) 0)
        return;

    // Relink the next pointer of predecessor
    if (pred != (TEdgeStump*) 0)
        assignNextT(pred, getNextT(current));
    else
        g.data_vertex[source] = getNextT(current);

    // Deallocate
    releaseId(g.data_id_managerE, _getId(current));
    valueDestruct(current);
    deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Directed<TCargo, TSpec> > & g,
           TEdgeDescriptor const & edge)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, sourceVertex(g,edge)));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, targetVertex(g,edge)));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    // Find edge and predecessor
    TEdgeStump* pred = 0;
    TEdgeStump* current = g.data_vertex[sourceVertex(g,edge)];
    while (current != (TEdgeStump*) 0)
    {
        if (current == edge)
            break;

        pred = current;
        current = getNextT(current);
    }

    // Not found?
    if (current == (TEdgeStump*) 0)
        return;

    // Relink the next pointer of predecessor
    if (pred != (TEdgeStump*) 0)
        assignNextT(pred, getNextT(current));
    else
        g.data_vertex[sourceVertex(g,edge)] = getNextT(current);

    // Deallocate
    releaseId(g.data_id_managerE, _getId(current));
    valueDestruct(current);
    deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Directed<TCargo, TSpec> > & g,
               TVertexDescriptor const & v)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    while (g.data_vertex[v] != (TEdgeStump*) 0)
    {
        TVertexDescriptor target = targetVertex(g, g.data_vertex[v]);
        removeEdge(g, v, target);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Directed<TCargo, TSpec> > & g,
              TVertexDescriptor const & v)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
    TIter it = begin(g.data_vertex, Standard());
    TIter itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for (; it != itEnd; ++it, ++pos)
    {
        TEdgeStump* current = *it;
        TVertexDescriptor const sourceVertex = pos;
        while (current!=0)
        {
            if (static_cast<TVertexDescriptor>(current->data_target) == v)
            {
                removeEdge(g, sourceVertex, v);
                current = g.data_vertex[sourceVertex];
            }
            else
            {
                current = getNextT(current);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Directed<TCargo, TSpec> > >::Type
targetVertex(Graph<Directed<TCargo, TSpec> > const &,
             TEdgeDescriptor const & edge)
{
    return getTarget(edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Directed<TCargo, WithSourceId> > >::Type
sourceVertex(Graph<Directed<TCargo, WithSourceId> > const &,
             TEdgeDescriptor const & edge)
{
    return getSource(edge);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Directed<TCargo, TSpec> > >::Type
sourceVertex(Graph<Directed<TCargo, TSpec> > const & g,
             TEdgeDescriptor const & edge)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for (; it != itEnd; ++it, ++pos)
    {
        TEdgeDescriptor current = *it;
        while (current != (TEdgeDescriptor) 0)
        {
            if (current == edge)
                return pos;

            current = getNextT(current);
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Directed<TCargo, TSpec> > const & g,
                   TMatrix & mat)
{
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TGraphSize;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TMatrix>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    typedef typename Value<TMatrix>::Type TMatValue;
    TSize len = getIdUpperBound(g.data_id_managerV);
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    clear(mat);
    resize(mat, len * len, (TMatValue) 0);
    TVertexDescriptor pos = 0;
    for (; it != itEnd; ++it, ++pos)
    {
        TEdgeStump* current = *it;
        TVertexDescriptor const source = pos;
        while (current != (TEdgeStump*) 0)
        {
            TVertexDescriptor target = targetVertex(g, current);
            mat[source * len + target] = static_cast<TMatValue>(static_cast<TGraphSize>(mat[source * len + target]) + 1);
            current = getNextT(current);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TVector, typename TGraph, typename TVertex>
inline void
_getVertexAdjacencyVector(TVector & vectIn,
                         TVector & vectOut,
                         TGraph const & g,
                         TVertex const & vertex)
{
    typedef typename Size<TGraph>::Type TGraphSize;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TVector>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    typedef typename Value<TVector>::Type TMatValue;

    TSize lenVectIn = inDegree(g, vertex);
    TSize lenVectOut = outDegree(g, vertex);
    clear(vectIn);
    clear(vectOut);
    resize(vectIn, lenVectIn, 0);
    resize(vectOut, lenVectOut, 0);
    TIterConst itIn = begin(g.data_vertex, Standard());
    TIterConst itEndIn = end(g.data_vertex, Standard());
    TSize count = 0;
    for(; itIn != itEndIn; ++itIn)
    {
        TEdgeStump * currentIn = *itIn;
        while(currentIn != 0)
        {
            if ((TVertexDescriptor) getTarget(currentIn) == vertex)
            {
                TVertexDescriptor source = sourceVertex(g, currentIn);
                vectIn[count] = static_cast<TMatValue>(static_cast<TGraphSize>(vectIn[count]) + source);
                ++count;
            }
            currentIn = getNextT(currentIn);
        }
    }
    count = 0;
    TEdgeStump * currentOut = getValue(g.data_vertex, vertex);
    while(currentOut != 0)
    {
        TVertexDescriptor target = targetVertex(g, currentOut);
        vectOut[count] = static_cast<TMatValue>(static_cast<TGraphSize>(vectOut[count]) + target);
        currentOut = getNextT(currentOut);
        ++count;
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TVector, typename TCargo, typename TSpec, typename TVertex>
inline void
getVertexAdjacencyVector(TVector & vectIn,
                         TVector & vectOut,
                         Graph<Directed<TCargo, TSpec> > const & g,
                         TVertex const & vertex)
{
    _getVertexAdjacencyVector(vectIn, vectOut, g, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Directed<TCargo, TSpec> > >::Type
findEdge(Graph<Directed<TCargo, TSpec> > const & g,
         TVertexDescriptor const & v,
         TVertexDescriptor const & w)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, w));

    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    TEdgeStump* current = g.data_vertex[v];
    while (current != (TEdgeStump*) 0)
    {
        if (static_cast<TVertexDescriptor>(getTarget(current)) == w)
            return current;

        current = getNextT(current);
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
write(TFile & target,
      Graph<Directed<TCargo, TSpec> > const & g)
{
//IOREV _nodoc_
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    write(target, "Adjacency list:\n");
    for (; it != itEnd; ++it, ++pos)
    {
        if (!idInUse(_getVertexIdManager(g), pos))
            continue;

        TEdgeStump* current = getValue(it);
        appendNumber(target, static_cast<int>(pos));
        write(target, " -> ");
        while (current != 0)
        {
            appendNumber(target, static_cast<int>(getTarget(current)));
            writeValue(target, ',');
            current = getNextT(current);
        }
        writeValue(target, '\n');
    }
    it = begin(g.data_vertex, Standard());
    pos = 0;
    write(target, "Edge list:\n");
    for (; it != itEnd; ++it, ++pos)
    {
        TEdgeStump* current = getValue(it);
        while (current != 0)
        {
            write(target, "Source: ");
            appendNumber(target, static_cast<int>(pos));
            write(target, ",Target: ");
            appendNumber(target, static_cast<int>(getTarget(current)));
            write(target, " (Id: ");
            appendNumber(target, static_cast<int>(_getId(current)));
            write(target, ")\n");
            current = getNextT(current);
        }
    }
}

}// namespace seqan

#endif //#ifndef SEQAN_HEADER_...
