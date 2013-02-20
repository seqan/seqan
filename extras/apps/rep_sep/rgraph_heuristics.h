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

#ifndef REPSEP_HEADER_RGRAPH_HEURISTICS_H
#define REPSEP_HEADER_RGRAPH_HEURISTICS_H

#include "rgraph.h"

#define DEBUG_RGRAPH_HEURISTICS

//////////////////////////////////////////////////////////////////////////////

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TVertexDescriptor>
double distance(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> & me, String<TVertexDescriptor> & comp1, String<TVertexDescriptor> & comp2)
{
    // graph typedefs
    typedef ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> RGraph;
    typedef typename SelectGraph_< RGraph >::Type  TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type               TEdgeDescriptor;
    typedef typename Size<String<TVertexDescriptor> >::Type     TSize;

    TSize len_comp1 = length(comp1);
    TSize len_comp2 = length(comp2);
    double contribution = 0.0;
    for (TSize i = 0; i < len_comp1; ++i)
    {
        for (TSize j = 0; j < len_comp2; ++j)
        {
            TEdgeDescriptor ed = findEdge(me.graph,value(comp1,i),value(comp2,j));
            if(ed != 0) contribution += getEdgeScore(me,ed);
        }
    }

    return contribution;
}

//////////////////////////////////////////////////////////////////////////////

struct GuidedParallelComponentMerge{};

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TComponentList>
unsigned solve(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> & me,
           TComponentList & components,
           unsigned copy_count,
           GuidedParallelComponentMerge const)
{

    // graph typedefs
    typedef ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> RGraph;
    typedef typename SelectGraph_< RGraph >::Type  TGraph;
    typedef typename VertexDescriptor<TGraph>::Type             TVertexDescriptor;
    //typedef typename EdgeDescriptor<TGraph>::Type               TEdgeDescriptor;
    //typedef typename Iterator<TGraph, OutEdgeIterator>::Type    TOutEdgeIterator;
    typedef typename Iterator<TGraph, VertexIterator>::Type     TVertexIter;
    typedef typename Value<TComponentList>::Type            TComponent;


    // general typedefs
    typedef unsigned TSize;

    // heuristic typedefs
    typedef ::std::set<TVertexDescriptor>        TVisitedSet;    
    typedef ::std::map<TVertexDescriptor, TSize> TVertexComponentMap;

    TVisitedSet visitedNodes;
    TComponentList componentList;

    TSize vertexCount = numVertices(me);
    resize(componentList,vertexCount,Exact());

    TVertexComponentMap vcMap;

    // generate one component for each vertex
    TVertexIter vIt(me.graph);
    TSize currentComponent = 0;
    for(;!atEnd(vIt);goNext(vIt))
    {
        TComponent empty;
        append(empty,*vIt);
        value(componentList,currentComponent) = empty;
        insert(vcMap,*vIt,currentComponent);
        ++currentComponent;
    }
    TSize componentCounter = vertexCount;

    bool merged = true;

    ::std::map< double, Pair<TSize>, ::std::greater<double> > unused_candidates;

    while(merged)
    {
        merged = false;
        clear(unused_candidates);

        // find smallest component
        for(TSize c1 = 0;c1 < vertexCount;++c1)
        {
            // this component was already merged
            if(length(value(componentList,c1)) == 0) continue;
            
            TSize c_merge = c1 + 1;
            // find the first non-empty component
            while(length(value(componentList,c_merge)) == 0 && c_merge < vertexCount)
                ++c_merge;

            if(c_merge == vertexCount) break;

            double best_dist = distance(me,value(componentList,c1),value(componentList,c_merge));

            // all components left where already merged .. so start at c1 + 1
            for(TSize c2 = c_merge + 1; c2 < vertexCount;++c2)
            {
                // check if component c2 was not already merged
                if(length(value(componentList,c2)) > 0)
                {
                    // this is a possible partner for c1   
                    if(distance(me,value(componentList,c1),value(componentList,c2)) < best_dist)
                    {
                        c_merge = c2;
                        best_dist = distance(me,value(componentList,c1),value(componentList,c_merge));
                    }
                }
            }

            // so we found a best match between c1 and anoter component
            if(best_dist < 0.0)
            {
                for (TSize i = 0; i < length(value(componentList,c_merge)); ++i)
                {
                    append(value(componentList,c1),value(value(componentList,c_merge),i));
                    vcMap[value(value(componentList,c_merge),i)] = c1;
                }
                clear(value(componentList,c_merge));
                --componentCounter;
                merged = true; // we merged two components so we can do another run
                if(componentCounter == copy_count) break; // abort if we reached the aimed amout of components
            }
            else
            {
                insert(unused_candidates,best_dist,Pair<TSize>(c1,c_merge));
            }
        }

        // the remaing merges are all not as good as we would need them to be
        if(!merged && componentCounter > copy_count)
        {
            typename Iterator< ::std::map<double, Pair<TSize>, ::std::greater<double> > >::Type unused_candIt;
            unused_candIt = begin(unused_candidates);

            for(;!atEnd(unused_candIt);goNext(unused_candIt))
            {
                if( length(value(componentList,cargo(unused_candIt).i1) ) != 0 // component 1 is not empty
                    && 
                    length(value(componentList,cargo(unused_candIt).i2) ) != 0 // component 2 is not empty
                    )
                {
                    // merge
                    merged = true;
                    TSize c1 = cargo(unused_candIt).i1;
                    TSize c_merge = cargo(unused_candIt).i2;

                    for (TSize i = 0; i < length(value(componentList,c_merge)); ++i)
                    {
                        append(value(componentList,c1),value(value(componentList,c_merge),i));
                        vcMap[value(value(componentList,c_merge),i)] = c1;
                    }
                    clear(value(componentList,c_merge));
                    --componentCounter;
                    break;
                }
            }
        }
        if(componentCounter == copy_count) break; // abort if we reached the aimed amout of components
    }


    // put computed components into final list
    clear(components);
    TSize component = 0;
    for(unsigned c = 0; c < length(componentList);++c)
    {
        if(length(value(componentList,c)) == 0) continue;
        else
        {
            appendValue(components,value(componentList,c));
            ++component;
        }
    }

    return component;
}

//////////////////////////////////////////////////////////////////////////////

struct SingleComponentExpansion{};

template <typename TColumnAlphabet, typename TAlignedReadStoreElement, typename TPosition, typename TComponentList>
unsigned solve(ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> & me,
           TComponentList & components,
           unsigned ,
           SingleComponentExpansion const)
{
    // constants 
    const double seed_score_threshold = 0.0;

    // graph typedefs
    typedef ReadGraph<TColumnAlphabet,TAlignedReadStoreElement,TPosition> RGraph;
    typedef typename SelectGraph_< RGraph >::Type  TGraph;
    typedef typename VertexDescriptor<TGraph>::Type             TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type               TEdgeDescriptor;
    //typedef typename Iterator<TGraph, OutEdgeIterator>::Type    TOutEdgeIterator;
    typedef typename Iterator<TGraph, EdgeIterator>::Type       TEdgeIterator;
    typedef typename Iterator<TGraph, VertexIterator>::Type     TVertexIterator;
    typedef typename Value<TComponentList>::Type TComponent;

    // general typedefs
    typedef unsigned TSize;

    typedef ::std::set<TVertexDescriptor>         TVisitedSet;    
    //typedef ::std::map<double, TVertexDescriptor> TBestMatchMap;

    TVisitedSet visitedVertices;
    TSize vertexCount = numVertices(me);

//    TComponentList components;
    clear(components);
    // as long as we do not have seen all vertices
    // -> build components
    while ( length(visitedVertices) != vertexCount )
    {
        // find best interaction
        TEdgeIterator edgeIt(me.graph);
        TEdgeDescriptor best_edge = 0;
        double best_edge_score = 0.0;
        bool init_best_edge_score = true; // take first possible value as best match

        TVertexDescriptor sV = 0, tV = 0;

        while(!atEnd(edgeIt))
        {
            TEdgeDescriptor ed = *edgeIt;
            if(!hasKey(visitedVertices,sourceVertex(edgeIt)) && !hasKey(visitedVertices,targetVertex(edgeIt)))
            {
                if(init_best_edge_score || getEdgeScore(me,ed) < best_edge_score) 
                {
                    best_edge_score = getEdgeScore(me, ed);
#ifdef DEBUG_RGRAPH_HEURISTICS
                    cout << "best edge score for init set to: " << best_edge_score << endl;
#endif                    
                    best_edge = ed;
                    init_best_edge_score = false;
                    sV = sourceVertex(edgeIt);
                    tV = targetVertex(edgeIt);
                }
            }
            goNext(edgeIt);
        }

        if(best_edge == 0) break; // we can not further proceed since there are no seeds left
        if(getEdgeScore(me,best_edge) > seed_score_threshold) break; // this is not a good seed ..

        // while possible 
        // -> extend current component
        
        // 1st .. init new component
        TComponent comp;
        appendValue(comp,sV);
        appendValue(comp,tV);
        
        // mark the used vertices as used
        insert(visitedVertices,sV);
        insert(visitedVertices,tV);
     
        TComponent temporary_component;
        
        // 2nd .. try to expand as long as possible
        bool was_expanded = true;
        while(was_expanded)
        {
            // reset 
            was_expanded = false;

            TVertexDescriptor best_expansion = 0;
            bool best_expansion_init = true;
            double best_expansion_score = 0;

            // try to find the best partner for the current component
            TVertexIterator vertexIt(me.graph);
            while(!atEnd(vertexIt))
            {
                // this vertex was not used 
                if(!hasKey(visitedVertices,*vertexIt))
                {
                    // candidate for an expansion
                    clear(temporary_component);
                    appendValue(temporary_component,*vertexIt);
                    double current_expansion_score = distance(me,comp,temporary_component);

                    if(best_expansion_init || current_expansion_score < best_expansion_score)
                    {
                        best_expansion = *vertexIt;
                        best_expansion_init = false;
                        best_expansion_score = current_expansion_score;
                    }
                }
                goNext(vertexIt);
            }
            
            // can we expand ???
            if(best_expansion_score < 0.0)
            {
                appendValue(comp,best_expansion);
                insert(visitedVertices,best_expansion);
                was_expanded = true;
            }
            
        }
        appendValue(components,comp);
    }

#ifdef REPEAT_HEURISTIC_DEBUG
    ::std::cout << "computed " << length(components) << " different components" << ::std::endl;
#endif
    return length(components);
}

#endif
