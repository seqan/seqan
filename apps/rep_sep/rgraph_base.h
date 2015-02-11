/*=========================================================================
  Copyright (C) 2009 by Stephan Aiche

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================
  $Id$
 ==========================================================================*/

#ifndef REPSEP_HEADER_RGRAPH_BASE_H
#define REPSEP_HEADER_RGRAPH_BASE_H

//////////////////////////////////////////////////////////////////////////////

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
struct GraphCargo {
    typedef typename Id < TAlignedReadStoreElement >::Type TId;
    
    // condensed information of the corresponding reads
    TAlignedReadStoreElement alignedRead;
    
    // derived information for the "candidate" columns
    typedef Pair< TPosition , TColumnAlphabet> TColumnInfo;
    typedef String< TColumnInfo > TColumns;
    typedef typename Size< TColumns >::Type TSize;

    TColumns spanned_columns;

    GraphCargo() {}
    GraphCargo(GraphCargo const& other) 
    {
        alignedRead = other.alignedRead;
        spanned_columns = other.spanned_columns;
    }
};

//////////////////////////////////////////////////////////////////////////////
// for compatibility 
template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
typename Id< GraphCargo<TColumnAlphabet,TAlignedReadStoreElement,TPosition> >::Type & 
id(GraphCargo<TColumnAlphabet, TAlignedReadStoreElement,TPosition> & me)
{
    return me.alignedRead.readId;
}

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
typename Id< GraphCargo<TColumnAlphabet, TAlignedReadStoreElement,TPosition> >::Type const & 
id(GraphCargo<TColumnAlphabet, TAlignedReadStoreElement,TPosition> const & me)
{
    return me.alignedRead.readId;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TIndex>
typename Value<TColumnAlphabet,1>::Type const &
_sequenceCharacter(GraphCargo<TColumnAlphabet, TAlignedReadStoreElement,TPosition> & me, TIndex i)
{
    return _sequenceCharacter(me.spanned_columns[i].i2);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
void addColumn(GraphCargo<TColumnAlphabet, TAlignedReadStoreElement,TPosition> & me, TPosition pos, TColumnAlphabet col)
{
	typedef typename GraphCargo<TColumnAlphabet,TAlignedReadStoreElement,TPosition>::TColumnInfo TCargoColumnInfo;
	TCargoColumnInfo col_info(pos,col);
    appendValue(me.spanned_columns, col_info);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
struct ReadGraph {
    typedef typename Id< GraphCargo<TColumnAlphabet, TAlignedReadStoreElement, TPosition> >::Type TId;
    typedef Graph<Undirected<double> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef String< GraphCargo<TColumnAlphabet,TAlignedReadStoreElement,TPosition> > TVertexCargoMap;
    typedef map<TId, TVertexDescriptor> TIdVertexMap;
    typedef typename Size<TGraph>::Type TSize;
    typedef String< GraphScoring::TScoreValue > TEdgeScoreMap;


    TGraph graph;
    TVertexCargoMap vertexCargo;
    TIdVertexMap idVertexMap;
    TEdgeScoreMap edgeScoreMap;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// not the best way
template <typename T>
struct SelectGraph_;

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
struct SelectGraph_< ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> >
{
    typedef Graph<Undirected<double> > Type;
};

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// convenience methods

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
bool hasMultipleComponents(ReadGraph<TColumnAlphabet, TAlignedReadStoreElement, TPosition> & me)
{
    typedef typename SelectGraph_< ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> >::Type TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef String<TSize> TComponentMap;

    TComponentMap components;
    TSize component_count = connectedComponents(components, me.graph);

    return (component_count > 1);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
inline bool
hasRead(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> & me, typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TId readId)
{
    return (me.idVertexMap.count(readId) != 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
inline GraphCargo<TColumnAlphabet,TAlignedReadStoreElement, TPosition> &
getCargo(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> & me, typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TId readId)
{
    return me.vertexCargo[me.idVertexMap[readId]];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
inline typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TVertexDescriptor const
getVertex(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> & me, typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TId readId) 
{
    if( hasRead(me, readId) )
    {
        return me.idVertexMap[readId];
    }
    else
    {
        return 0;
    }    
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
inline GraphCargo<TColumnAlphabet,TAlignedReadStoreElement, TPosition> &
registerRead(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> & me, typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TId readId) 
{
    typedef typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TVertexDescriptor TVertexDescriptor;

    if( hasRead(me, readId) ) 
    {
        // readId is already registered -> do nothing, just return cargo
        return getCargo(me, readId);
    }
    else 
    {
        TVertexDescriptor vd = addVertex(me.graph);
        GraphCargo<TColumnAlphabet,TAlignedReadStoreElement,TPosition> new_cargo;

        if(length(me.vertexCargo) <= vd) resize(me.vertexCargo, vd + 1, Generous());
        assignProperty(me.vertexCargo, vd, new_cargo);

        // remember association between readId and cargo
        insert(me.idVertexMap, readId, vd);

        return getCargo(me, readId);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
inline void
addEdge(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> & me, 
             typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TId readId1, 
             typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TId readId2)
{
    typedef typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition>::TVertexDescriptor TVertexDescriptor;
    
    TVertexDescriptor vd_readId1 = me.idVertexMap[readId1];
    TVertexDescriptor vd_readId2 = me.idVertexMap[readId2];

    // check if this edge already exists
    if(findEdge(me.graph, vd_readId1 , vd_readId2) == 0)
    {    
      addEdge(me.graph, vd_readId1 , vd_readId2);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TEdgeDescriptor>
inline void
setEdgeScore(ReadGraph<TColumnAlphabet, TAlignedReadStoreElement, TPosition> & me,
             TEdgeDescriptor edge, GraphScoring::TScoreValue score)
{
    if (length(me.edgeScoreMap) <= _getId(edge)) resize(me.edgeScoreMap, _getId(edge) +1, Generous());
    assignProperty(me.edgeScoreMap,edge,score);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TEdgeDescriptor>
inline GraphScoring::TScoreValue
getEdgeScore(ReadGraph<TColumnAlphabet, TAlignedReadStoreElement, TPosition> & me,
             TEdgeDescriptor edge)
{
    return getProperty(me.edgeScoreMap, edge);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
inline typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TSize const
numVertices(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> & me) 
{
    return numVertices(me.graph);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition>
inline typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition>::TSize const
numEdges(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement, TPosition> & me) 
{
    return numEdges(me.graph);
}

#endif
