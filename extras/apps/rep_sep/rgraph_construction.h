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

#ifndef REPSEP_HEADER_RGRAPH_CONSTRUCTION_H
#define REPSEP_HEADER_RGRAPH_CONSTRUCTION_H

#define RGRAPH_DEBUG_CONSTRUCTION_SCORING
#define RGRAPH_DEBUG_CONSTRUCTION

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TSpec, typename TConfig, typename TId>
void construct(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> & me,
          String< Pair< TPosition, String<TColumnAlphabet> > > const & column_set,
          FragmentStore<TSpec, TConfig> const& fragStore,
          TId const contigId)
{
    // fragmentstore typedefs
  	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    //typedef typename Size<TFragmentStore>::Type TSize;
    //typedef typename TFragmentStore::TReadPos TReadPos;

    typedef typename Iterator<typename TFragmentStore::TAlignedReadStore const>::Type TAlignIter;
    TAlignIter alignItBegin = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());    

    typedef typename Iterator<String< Pair< TPosition, String<TColumnAlphabet> > > const>::Type TColumnSetIterator;
    typedef typename Iterator<String<TColumnAlphabet> const>::Type TColumnIterator;

    // graph typedefs
    typedef GraphCargo<TColumnAlphabet, TAlignedReadStoreElement,TPosition>  TGraphCargo;

    TColumnSetIterator colSetIter = begin(column_set);
    TColumnSetIterator colSetIterEnd = end(column_set);

    for(;colSetIter != colSetIterEnd;goNext(colSetIter)) 
    {
        TPosition current_column = value(colSetIter).i1;

        // 
        TColumnIterator colIter = begin(value(colSetIter).i2);
        TColumnIterator colIterEnd = end(value(colSetIter).i2);

        String< typename Id<TGraphCargo>::Type > readsInColumn;

        for(;colIter != colIterEnd ; goNext(colIter))
        {
            // check if graph already has this read
            if(!hasRead(me, _readId(value(colIter)))) {
                // add this read to the graph
                TAlignedReadStoreElement alignedRead;
                bool found = false;
                // find AlignedRead in fragStore .. inefficient but necessary
                TAlignIter alignIt = alignItBegin;
                
                for(; alignIt != alignItEnd ; goNext(alignIt)) 
                {
                    if(alignIt->readId == _readId(value(colIter))) {
                        alignedRead = value(alignIt);
                        found = true;
                        break;
                    }
                }
                
                if(found) {
                    TGraphCargo & new_cargo = registerRead(me, _readId(value(colIter)));
                    new_cargo.alignedRead = alignedRead;
                } 
                else
                {
                    cerr << "ERROR: could not find read in FragmentStore (readId = " << _readId(value(colIter)) << " )" << endl;
                    return;
                }
            }

            // add this column to the vertex inside the graph
            TGraphCargo & cargo = getCargo(me, _readId(value(colIter)));          
            addColumn(cargo, current_column, value(colIter));
#ifdef RGRAPH_DEBUG_CONSTRUCTION
            cout << "added column " << current_column << "(" << value(colIter) << ") to vertex " << getVertex(me, _readId(value(colIter))) << endl;
#endif            

            append(readsInColumn, _readId(value(colIter)));
        }

        // associate all reads in this column with each other
        typedef typename Iterator< String< typename Id<TGraphCargo>::Type > >::Type TVertexIterator;
        TVertexIterator v_iter = begin(readsInColumn);
        TVertexIterator v_iterEnd = end(readsInColumn);

        for(; v_iter != v_iterEnd ; goNext(v_iter)) 
        {
            TVertexIterator v_iterInternal = v_iter;
            goNext(v_iterInternal);

            for(; v_iterInternal != v_iterEnd ; goNext(v_iterInternal)) {
                // add edge between nodes
                addEdge(me, value(v_iter), value(v_iterInternal));
            }        
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// add mate pairs
// TODO: need to be implemented

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TSpec, typename TConfig, typename TId>
void add_mate_pairs(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> & me,
          FragmentStore<TSpec, TConfig> const& fragStore,
          const TId contigId)
{
    // Get rid of unused variable warnings.
    (void)me;
    (void)fragStore;
    (void)contigId;
    /*
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename Size<TFragmentStore>::Type TSize;
    typedef typename TFragmentStore::TReadPos TReadPos;

    typedef typename SelectGraph_< ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> >::Type TGraphInternal;
    
    typedef typename Iterator<TGraphInternal, VertexIterator>::Type TVertexIterator;
    */
}

//////////////////////////////////////////////////////////////////////////////

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TSpec, typename TConfig, typename TId>
void scoreGraph_(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> & me,
          FragmentStore<TSpec, TConfig> const& fragStore,
          TId const contigId,
          GraphScoring const& scoring)
{
    // Get rid of unused variable warnings.
    (void)fragStore;
    (void)contigId;

    //typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    //typedef typename Size<TFragmentStore>::Type TSize;
    //typedef typename TFragmentStore::TReadPos TReadPos;

    typedef typename SelectGraph_< ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> >::Type TGraphInternal;
    
    typedef typename Iterator<TGraphInternal, EdgeIterator>::Type TEdgeIterator;
    typedef typename ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition>::TVertexDescriptor TVertexDescriptor;

    TEdgeIterator edgeIter(me.graph);

    for(; !atEnd(edgeIter); goNext(edgeIter))
    {
        TVertexDescriptor source = sourceVertex(me.graph, value(edgeIter) );
        TVertexDescriptor target = targetVertex(me.graph, value(edgeIter) );

        GraphScoring::TScoreValue edge_score;
        _computeScore(me,edge_score,source,target,scoring);
#ifdef RGRAPH_DEBUG_CONSTRUCTION_SCORING
        cout << "assigned score " << edge_score << endl;
#endif        
        setEdgeScore(me,value(edgeIter),edge_score);
    }
}

template<typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TVertexDescriptor>
void _computeScore(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> & me, 
                    GraphScoring::TScoreValue & edge_score,
                    TVertexDescriptor const vd1,
                    TVertexDescriptor const vd2,
                    GraphScoring const& scoring)
{
    typedef GraphCargo<TColumnAlphabet, TAlignedReadStoreElement,TPosition> TGraphCargo;
    typedef typename TGraphCargo::TSize TSize;

    TGraphCargo & cargo1 = me.vertexCargo[vd1];
    TGraphCargo & cargo2 = me.vertexCargo[vd2];
#ifdef RGRAPH_DEBUG_CONSTRUCTION_SCORING
    cout << "Compute score between " << vd1 << " and " << vd2 << endl;
    cout << cargo1.alignedRead.readId << " " << cargo2.alignedRead.readId << endl;
#endif    
    // reset score
    edge_score = GraphScoring::TScoreValue();
#ifdef RGRAPH_DEBUG_CONSTRUCTION_SCORING
    cout << "v1 spans " << length(cargo1.spanned_columns) << " columns" << endl;
    cout << "v2 spans " << length(cargo2.spanned_columns) << " columns" << endl;
#endif
    for (TSize i = 0 ; i < length(cargo1.spanned_columns) ; ++i) 
    {
        for (TSize j = 0 ; j < length(cargo2.spanned_columns) ; ++j) 
        {
            // check if we deal with the same position in the assembly
            if( cargo1.spanned_columns[i].i1 == cargo2.spanned_columns[j].i1 ) 
            {
#ifdef RGRAPH_DEBUG_CONSTRUCTION_SCORING
                cout << "both span column " << cargo1.spanned_columns[i].i1 << " with character " << _sequenceCharacter(cargo1.spanned_columns[i].i2) << " and " << _sequenceCharacter(cargo2.spanned_columns[j].i2) << endl;
#endif            
                if( _sequenceCharacter(cargo1.spanned_columns[i].i2) 
                    == 
                    _sequenceCharacter(cargo2.spanned_columns[j].i2))
                {
                    edge_score += matchScore(scoring);
                }
                else 
                {
                    edge_score += mismatchScore(scoring);
                }
            }   
        }
    }

    // check mate pairs
    
}

//////////////////////////////////////////////////////////////////////////////

#endif
