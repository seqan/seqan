// ===========================================================================
//                 SGIP - Solution of Graph Isomorphism Problem
// ===========================================================================
// Copyright (C) 2012 by Jialu Hu
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ===========================================================================
// Author: Jialu Hu <Jialu.Hu@fu-berlin.de>
// ===========================================================================

#ifndef APPS_SGIP_SGIP_BASE_H_
#define APPS_SGIP_SGIP_BASE_H_

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

// ==========================================================================
// Forwards
// ==========================================================================

bool const VISITED = 1;
bool const UNVISITED = 0;
bool const OCCUPYED = 1;
bool const FREE = 0;
size_t const INITIAL_PRIMER = 100;
unsigned short const PRIMER[168] =
{
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
    263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373,
    379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487,
    491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613,
    617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
    743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
    877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997
};

// Declaration of compare_ template function.
template <typename T>
int compare_(T const &, T const &); // < -1; = 0; > 1;

template <typename TValue>
inline int compare_(std::vector<TValue> const & obj1, std::vector<TValue> const & obj2);

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Class SgipDegreeDir_
// --------------------------------------------------------------------------

template <typename TValue>
struct SgipDegreeDir_
{
    TValue inDegree;
    TValue outDegree;
    SgipDegreeDir_()
    {
        inDegree = 0;
        outDegree = 0;
    }

    SgipDegreeDir_(TValue in, TValue out)
    {
        inDegree = in;
        outDegree = out;
    }
};

// --------------------------------------------------------------------------
// Class SgipDegreeUndir_
// --------------------------------------------------------------------------

template <typename TValue>
struct SgipDegreeUndir_
{
    TValue degree;
};

// --------------------------------------------------------------------------
// Class SgipHash
// --------------------------------------------------------------------------

template <typename TString = std::vector<int> >
struct SgipHash;

template <typename TValue, unsigned int LENGTH>
struct SgipHash<seqan::String<TValue, seqan::Array<LENGTH> > >
{
    size_t operator()(seqan::String<TValue, seqan::Array<LENGTH> > const & s) const
    {
        using namespace seqan;
        typedef typename Iterator<String<TValue, Array<LENGTH> > >::Type TIterator;
        size_t i = 0;
        size_t sum = 0;
        TIterator it = begin(s);
        TIterator itEnd = end(s);
        while (i < LENGTH && it != itEnd)
        {
            sum += PRIMER[INITIAL_PRIMER + i] * s[i];
            i++;
            goNext(it);
        }
        return sum;
    }
};

template <typename TValue>
struct SgipHash<std::vector<TValue> >
{
    size_t operator()(std::vector<TValue> const & s) const
    {
        //typedef typename std::vector<TValue>::iterator TIterator;
        size_t i = 0;
        size_t sum = 0;
        size_t len = s.size();
        while (i < len)
        {
            sum += PRIMER[INITIAL_PRIMER + i] * s[i];
            i++;
        }
        return sum;
    }
};

// --------------------------------------------------------------------------
// Class SgipEqualTo
// --------------------------------------------------------------------------

template <typename TString>
struct SgipEqualTo;

template <typename TValue, unsigned int LENGTH>
struct SgipEqualTo<seqan::String<TValue, seqan::Array<LENGTH> > >
{
    bool operator()(seqan::String<TValue, seqan::Array<LENGTH> > const & s1,
                    seqan::String<TValue, seqan::Array<LENGTH> > const & s2) const
    {
        using namespace seqan;
        typedef String<TValue, Array<LENGTH> >     TString;
        typedef typename Iterator<TString>::Type   TIterator;

        TIterator it, ti, itEnd;
        it = begin(s1);
        ti = begin(s2);
        itEnd = end(s1);
        while (it != itEnd)
        {
            if (getValue(it) != getValue(ti))
                return false;
            goNext(it);
            goNext(ti);
        }
        return true;
    }
};

template <typename TValue>
struct SgipEqualTo<std::vector<TValue> >
{
    bool operator()(std::vector<TValue> const & s1, std::vector<TValue> const & s2) const
    {
        //typedef std::vector<TValue> TString;
        if (s1.size() != s2.size())
            return false;

        size_t len = s1.size();
        for (size_t i = 0; i < len; i++)
        {
            if (s1[i] != s2[i])
                return false;
        }
        return true;
    }
};

// --------------------------------------------------------------------------
// Class LessCompare_
// --------------------------------------------------------------------------

template <typename TTag>
struct LessCompare_;

template <typename TValue>
struct LessCompare_<SgipDegreeDir_<TValue> >
{
    bool operator()(SgipDegreeDir_<TValue> const & obj1, SgipDegreeDir_<TValue> const & obj2) const
    {
        return compare_(obj1, obj2) < 0;
    }
};

template <typename TValue, unsigned int LENGTH>
struct LessCompare_<seqan::String<TValue, seqan::Array<LENGTH> > >
{
    bool operator()(
        seqan::String<TValue, seqan::Array<LENGTH> > const & obj1,
        seqan::String<TValue, seqan::Array<LENGTH> > const & obj2
        ) const
    {
        return compare_(obj1, obj2) < 0;
    }
};

template <typename TValue>
struct LessCompare_<std::vector<TValue> >
{
    bool operator()(std::vector<TValue> const & obj1, std::vector<TValue> const & obj2) const
    {
        return compare_(obj1, obj2) < 0;
    }
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function compare_()
// --------------------------------------------------------------------------

// Comparison of various data types.
template <typename T>
int compare_(T const &, T const &); // < -1; = 0; > 1;

template <typename TValue>
int compare_(SgipDegreeDir_<TValue> const &, SgipDegreeDir_<TValue> const &);

template <>
inline int compare_(SgipDegreeDir_<unsigned int> const & obj1, SgipDegreeDir_<unsigned int> const & obj2)
{
    if (obj1.inDegree < obj2.inDegree)
    {
        return -1;
    }
    else if (obj1.inDegree == obj2.inDegree)
    {
        if (obj1.outDegree < obj2.outDegree)
            return -1;
        if (obj1.outDegree == obj2.outDegree)
            return 0;
    }
    return 1;
}

template <typename TValue>
inline int compare_(seqan::String<TValue> const & obj1, seqan::String<TValue> const & obj2)
{
    using namespace seqan;
    typedef String<TValue> const             TString;
    typedef typename Iterator<TString>::Type TIterator;
    TIterator it, itEnd, ti, tiEnd;
    it = begin(obj1);
    ti = begin(obj2);
    itEnd = end(obj1);
    tiEnd = end(obj2);
    while (it < itEnd && ti < tiEnd)
    {
        if (getValue(it) < getValue(ti))
        {
            return -1;
        }
        else if (getValue(it) > getValue(ti))
        {
            return 1;
        }
        else
        {
            it++;
            ti++;
        }
    }
    if (it < itEnd)
        return 1;
    if (ti < tiEnd)
        return -1;
    return 0;
}

template <typename TValue>
inline int compare_(std::vector<TValue> const & obj1, std::vector<TValue> const & obj2)
{
    using namespace seqan;
    typedef std::vector<TValue> const        TString;
    typedef typename Iterator<TString>::Type TIterator;
    TIterator it, itEnd, ti, tiEnd;
    it = obj1.begin();
    ti = obj2.begin();
    itEnd = obj1.end();
    tiEnd = obj2.end();
    while (it < itEnd && ti < tiEnd)
    {
        if (*it < *ti)
        {
            return -1;
        }
        else if (*it > *ti)
        {
            return 1;
        }
        else
        {
            it++;
            ti++;
        }
    }
    if (it < itEnd)
        return 1;
    if (ti < tiEnd)
        return -1;
    return 0;
}

// --------------------------------------------------------------------------
// Function _createUndir()
// --------------------------------------------------------------------------

// Create undirected graphs from directed graphs.
void _createUndir(seqan::Graph<seqan::Directed<> > const & graph, seqan::Graph<seqan::Undirected<> > & ugraph)
{
    using namespace seqan;
    typedef Graph<Directed<> >                    TGraph;
    typedef VertexDescriptor<TGraph>::Type        TVertexDescriptor;
    typedef Iterator<TGraph, EdgeIterator>::Type  TIterator;

    assert(numVertices(ugraph) == 0);
    std::vector<TVertexDescriptor> edges;
    TIterator it(graph);
    goBegin(it);

    while (!atEnd(it))
    {
        edges.push_back(sourceVertex(it));
        edges.push_back(targetVertex(it));
        goNext(it);
    }
    size_t len = length(edges) / 2;
    addEdges(ugraph, edges, len);
}

// --------------------------------------------------------------------------
// Function _getDistanceMatrixBfs()
// --------------------------------------------------------------------------

// Get distance matrix through bfs method.
template <typename TSpec, typename TMatrix>
void _getDistanceMatrixBfs(seqan::Graph<TSpec> const &, TMatrix &);

template <typename TMatrix, typename TCargo, typename TSpec>
void _getDistanceMatrixBfs(seqan::Graph<seqan::Undirected<TCargo, TSpec> > const & graph, TMatrix & mat)
{
    using namespace seqan;
    typedef Graph<Undirected<TCargo, TSpec> >               TGraph;
    typedef typename VertexDescriptor<TGraph>::Type         TVertexDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TIterator;
    typedef typename Size<TMatrix>::Type                    TSize;
    //typedef typename Value<TMatrix>::Type                   TString;

    TSize len = numVertices(graph);
    TMatrix predmap;
    resize(mat, len);
    resize(predmap, len);

    for (TIterator it(graph); !atEnd(it); goNext(it))
    {
        TVertexDescriptor v = getValue(it);
        breadthFirstSearch(predmap[v], mat[v], graph, v);
    }
}

template <typename TMatrix, typename TCargo, typename TSpec>
void _getDistanceMatrixBfs(seqan::Graph<seqan::Directed<TCargo, TSpec> > const & graph, TMatrix & mat)
{
    using namespace seqan;
    Graph<Undirected<TCargo, TSpec> > ugraph;
    _createUndir(graph, ugraph);
    _getDistanceMatrixBfs(ugraph, mat);
}

// --------------------------------------------------------------------------
// Function _caculateDegreeMap()
// --------------------------------------------------------------------------

// Get degree map of an input graph.
template <typename TKey, typename TValue, typename TSpec, typename THKey, typename THValue, typename TGraph>
void _caculateDegreeMap(std::multimap<TKey, TValue, TSpec> &,
                        std::unordered_map<THKey, THValue> &,
                        seqan::Graph<TGraph> const &);

template <typename TKey, typename TValue, typename TSpec>
void _caculateDegreeMap(
    std::multimap<TKey, seqan::VertexDescriptor<seqan::Graph<seqan::Directed<> > >::Type, TSpec> & orderedmap,
    std::unordered_map<seqan::VertexDescriptor<seqan::Graph<seqan::Directed<> > >::Type, TValue> & degreemap,
    seqan::Graph<seqan::Directed<> > const & graph
    )
{
    using namespace seqan;
    typedef Graph<Directed<> >                                                      TGraph;
    typedef typename VertexDescriptor<TGraph>::Type                                 TVertexDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type                         TIterator;
    //typedef typename std::unordered_map<TVertexDescriptor, TValue>::size_type  THElement;
    typedef typename std::multimap<TKey, TVertexDescriptor, TSpec>::value_type      TMElement;

    for (TIterator it(graph); !atEnd(it); goNext(it))
    {
        TValue    data(inDegree(graph, getValue(it)), outDegree(graph, getValue(it)));
        TMElement mele(data, getValue(it));
        orderedmap.insert(mele);
        degreemap[getValue(it)] = data;
    }
}

// --------------------------------------------------------------------------
// Function _getRefinement()
// --------------------------------------------------------------------------

// Get refinement of a new vertexorder.
template <typename TValue, typename TEntry>
void _getRefinement(
    seqan::String<TValue> const & vertexorder,
    seqan::String<TEntry> const & mat,
    seqan::String<TEntry> & tempmat
    )
{
    using namespace seqan;
    //typedef typename Iterator<String<TValue> >::Type TIterator;
    typedef typename Size<String<TEntry> >::Type     TSize;
    TSize     len, i, j, vi, vj;
    len = length(vertexorder);
    for (i = 0; i < len; i++)
    {
        vi = vertexorder[i];
        for (j = 0; j < len; j++)
        {
            vj = vertexorder[j];
            if (mat[vi * len + vj])
                tempmat[i * len + j] = true;
            else
                tempmat[i * len + j] = false;
        }
    }
}

// --------------------------------------------------------------------------
// Function _getHeuristicList()
// --------------------------------------------------------------------------

// Caculate list of candidates through heuristic approach.
template <typename TString, typename TXMap, typename TMatrix, typename TDMap>
bool _getHeuristicList(TString & rset, TXMap & cand, TMatrix & dismat, TString & heuristicString, TDMap & degreemap)
{
    using namespace seqan;
    //typedef typename TString::iterator TIterator;
    typedef typename TXMap::iterator TIter;
    typedef typename TXMap::key_type T;
    typedef std::unordered_map<TString, T, SgipHash<TString>, SgipEqualTo<TString> > THMap;
    //typedef typename THMap::size_type TValue;

    TIter itEnd = cand.end();
    size_t max_capacity(0);
    heuristicString.clear();
    for (TIter it = cand.begin(); it != itEnd; ++it)
    {
        if (it->second)
            continue;
        rset.push_back(it->first);
        THMap hmap;
        size_t len = length(rset);
        for (TIter ti = cand.begin(); ti != itEnd; ++ti)
        {
            TString t_key;
            t_key.push_back(degreemap[ti->first].inDegree);
            t_key.push_back(degreemap[ti->first].outDegree);
            for (size_t i = 0; i < len; i++)
                t_key.push_back(dismat[rset[i]][ti->first]);
            if (hmap.find(t_key) == hmap.end())
                hmap[t_key] = ti->first;
        }

        if (hmap.size() > max_capacity)
        {
            max_capacity = hmap.size();
            heuristicString.clear();
            heuristicString.push_back(it->first);
        }
        else if (hmap.size() == max_capacity)
        {
            if (compare_(degreemap[it->first], degreemap[heuristicString[0]]) == 0)
            {
                heuristicString.push_back(it->first);
            }
            else if (compare_(degreemap[it->first], degreemap[heuristicString[0]]) < 0)
            {
                heuristicString.clear();
                heuristicString.push_back(it->first);
            }
        }
        hmap.clear();
        rset.pop_back();
    }
    if (max_capacity == cand.size())
        return true;
    return false;
}

// --------------------------------------------------------------------------
// Function _getLeastMat()
// --------------------------------------------------------------------------

// Get least mat when a new resolving set is available.
template <typename TValue, typename TCandMap, typename TMatrix, typename THMap, typename TDMap, typename TString,
    typename TMap, typename TMat>
    void _getLeastMat(std::vector<TValue> & rset,
                      TCandMap & cand,
                      TMatrix & dismat,
                      TDMap & degreemap,
                      TMap & paritymap,
                      THMap & checkmap,
                      TString & vertexorder,
                      TString & standardVertexOrder,
                      TMat & mat,
                      TMat & leastmat,
                      TMat & tempmat,
                      bool isSmaller)
{
    using namespace seqan;
    typedef typename TCandMap::iterator              TIterator;
    typedef std::vector<TValue>                      TVector;
    typedef typename std::vector<TValue>::iterator   TVIterator;
    //typedef typename THMap::key_type                 TKey;
    typedef typename THMap::iterator                 TMIterator;
    typedef typename Iterator<String<TValue> >::Type TStrIterator;

    TIterator itEnd = cand.end();
    TVIterator tiEnd = rset.end();
    checkmap.clear();

    for (TIterator it = cand.begin(); it != itEnd; ++it)
    {
        TVector key;
        key.push_back(degreemap[it->first].inDegree);
        key.push_back(degreemap[it->first].outDegree);
        for (TVIterator ti = rset.begin(); ti != tiEnd; ++ti)
            key.push_back(dismat[it->first][*ti]);
        checkmap[key] = it->first;
    }
    size_t i = 0;
    for (TMIterator im = checkmap.begin(); im != checkmap.end(); ++im)
    {
        vertexorder[i++]   = im->second;
        TStrIterator ipEnd = end(paritymap[im->second]);
        for (TStrIterator ip = begin(paritymap[im->second]); ip != ipEnd; ++ip)
            vertexorder[i++] = getValue(ip);
    }
    assert(length(vertexorder) == degreemap.size());
    _getRefinement(vertexorder, mat, tempmat);
    if (isSmaller || compare_(tempmat, leastmat) < 0)
    {
        assign(leastmat, tempmat);
        assign(standardVertexOrder, vertexorder);
    }
}

// --------------------------------------------------------------------------
// Function _greedySearch()
// --------------------------------------------------------------------------

template <typename TString, typename TMap, typename TXMap,
    typename TDMap, typename THMap, typename TMatrix, typename TMat, typename TValue>
    bool _greedySearch(TString & rset,
                       TMat & mat,
                       TMat & leastmat,
                       TMatrix & dismat,
                       TDMap & degreemap,
                       TMap & paritymap,
                       THMap & checkmap,
                       TXMap & cand,
                       TString & heuristicString,
                       seqan::String<TValue> & vertexorder,
                       seqan::String<TValue> & standardVertexOrder,
                       TMat & tempmat,
                       bool & flag,
                       size_t & len)
{
    using namespace seqan;
    unsigned local_len = length(heuristicString);
    for (unsigned i = 0; i < local_len; i++)
    {
        rset.push_back(heuristicString[i]);
        cand[heuristicString[i]] = OCCUPYED;
        if (flag)
        {
            if (len > rset.size())
            {
                len = rset.size();
                _getLeastMat(rset, cand, dismat, degreemap, paritymap, checkmap,
                    vertexorder, standardVertexOrder, mat, leastmat, tempmat, true);
            }
            else
            {
                _getLeastMat(rset, cand, dismat, degreemap, paritymap, checkmap,
                    vertexorder, standardVertexOrder, mat, leastmat, tempmat, false);
            }
            cand[heuristicString[i]] = FREE;
            rset.pop_back();
            continue;
        }
        if (rset.size() == len)
        {
            cand[heuristicString[i]] = FREE;
            rset.pop_back();
            continue;
        }
        TString filter;
        flag = _getHeuristicList(rset, cand, dismat, filter, degreemap);
        _greedySearch(rset, mat, leastmat, dismat, degreemap, paritymap, checkmap, cand,
            filter, vertexorder, standardVertexOrder, tempmat, flag, len);
        cand[heuristicString[i]] = FREE;
        rset.pop_back();
    }
    flag = false;
    return 1;
}

// --------------------------------------------------------------------------
// Function isParity()
// --------------------------------------------------------------------------

// Check whether two vertex are parity nodes.
template <typename TValue, typename TNum, typename TVertex>
bool isParity(seqan::String<TValue> const & mat, TVertex it, TVertex ti, TNum n)
{
    for (unsigned short i = 0; i < n; i++)
    {
        if (i == it || i == ti)
            continue;
        if (mat[i + n * it] != mat[i + n * ti] || mat[it + n * i] != mat[ti + n * i])
            return false;
    }
    return true;
}

// --------------------------------------------------------------------------
// Function _createParityMap()
// --------------------------------------------------------------------------

// Create parity map.
template <typename TEntry, typename TValue, typename TSpec, typename TArrString>
void _createParityMap(
    seqan::String<TEntry> const & mat,
    std::unordered_map<seqan::VertexDescriptor<seqan::Graph<seqan::Directed<> > >::Type,
    seqan::String<seqan::VertexDescriptor<seqan::Graph<seqan::Directed<> > >::Type> > & parityMap,
    std::multimap<TValue, seqan::VertexDescriptor<seqan::Graph<seqan::Directed<> > >::Type, TSpec> & degreeMap,
    seqan::Graph<seqan::Directed<> > const & graph,
    TArrString & cand
    )
{
    using namespace seqan;
    typedef VertexDescriptor<Graph<Directed<> > >::Type                               TVertexDescriptor;
    typedef std::vector<TVertexDescriptor>                                            TString;
    //typedef typename Iterator<TString, Rooted>::Type                                  TStrIterator;
    typedef typename std::multimap<TValue, TVertexDescriptor, TSpec>::iterator        TMapIterator;
    //typedef typename std::unordered_map<TVertexDescriptor, TString>::size_type   TPair;
    typedef String<bool>                                                              TProperties;
    //typedef Iterator<TProperties>::Type                                               TProIterator;
    TProperties visitedRecord;
    unsigned int numVer = numVertices(graph);

    resizeVertexMap(visitedRecord, graph);
    TMapIterator it = degreeMap.begin();
    TMapIterator itEnd = degreeMap.end();

    while (it != itEnd)
    {
        if (!getProperty(visitedRecord, it->second))
        {
            TString str;
            TMapIterator ti = it;
            ti++;
            while(ti != itEnd)
            {
                if(ti->first.inDegree == it->first.inDegree && ti->first.outDegree == it->first.outDegree)
                {
                    if (!getProperty(visitedRecord, ti->second))
                    {
                        if (isParity(mat, it->second, ti->second, numVer))
                        {
                            str.push_back(ti->second);
                            assignProperty(visitedRecord, ti->second, VISITED);
                        }
                    }
                    goNext(ti);
                }
                else
                {
                    break;
                }
            }
            cand.push_back(it->second);
            parityMap[it->second] = str;
        }
        goNext(it);
    }
}

// --------------------------------------------------------------------------
// Function _createCandMap()
// --------------------------------------------------------------------------

// Create candmap which parity nodes are excluded.
template <typename TVal, typename TString>
void _createCandMap(std::unordered_map<TVal, bool> & hmap, TString & cand)
{
    typedef typename TString::iterator TIterator;
    TIterator itEnd = cand.end();
    if (!hmap.empty())
        hmap.clear();
    for (TIterator it = cand.begin(); it != itEnd; ++it)
        hmap[*it] = FREE;
}

// --------------------------------------------------------------------------
// Function getCanonicalLabel()
// --------------------------------------------------------------------------

template <typename TSpec, typename TMat>
bool getCanonicalLabel(TMat & leastmat, seqan::Graph<TSpec> const & graph)
{
    using namespace seqan;
    typedef bool                                                           TFlag;
    typedef Graph<TSpec>                                                   TGraph;
    typedef typename VertexDescriptor<TGraph>::Type                        TVertexDescriptor;
    //typedef typename Iterator<TGraph, VertexIterator>::Type                TIterator;
    typedef String<TVertexDescriptor>                                      TString;
    typedef SgipDegreeDir_<TVertexDescriptor>                              TDegree_dir;
    typedef LessCompare_<TDegree_dir>                                      _Less;
    typedef std::unordered_map<TVertexDescriptor, TDegree_dir>             TMap;
    typedef std::multimap<TDegree_dir, TVertexDescriptor, _Less>           TMulMap;
    typedef std::vector<TVertexDescriptor>                                 T2String;
    typedef LessCompare_<T2String>                                         _LessStr;
    typedef std::map<T2String, TVertexDescriptor, _LessStr>                THMap;
    typedef std::unordered_map<TVertexDescriptor, TString>                 TPMap;
    typedef StringSet<TString>                                             TMatrix;
    typedef std::deque<TVertexDescriptor>                                  TVertexString;
    typedef std::unordered_map<TVertexDescriptor, bool>                    TCandMap;
	typedef typename Value<TMat>::Type                                     TMatValue;

    TFlag    flag = false;
    TMat     mat;
    TMat     tempmat;
    TMap     degreemap;
    TMulMap  orderedmap;
    TPMap    paritymap;
    TMatrix  dismat;
    T2String rset;
    T2String heuristicString;
    TVertexString cand;
    THMap   checkmap;
    TString standardVertexOrder;
    TString vertexorder;
    TCandMap candmap;
    size_t  len = -1;

    getAdjacencyMatrix(graph, mat);
    resize(tempmat, length(mat), TMatValue(0));
    resize(leastmat, length(mat), TMatValue(1));
    resize(vertexorder, numVertices(graph));
    _caculateDegreeMap(orderedmap, degreemap, graph);
    _createParityMap(mat, paritymap, orderedmap, graph, cand);
    _getDistanceMatrixBfs(graph, dismat);
    _createCandMap(candmap, cand);
    flag = _getHeuristicList(rset, candmap, dismat, heuristicString, degreemap);
    if (!_greedySearch(rset,mat,leastmat,dismat,degreemap,paritymap,checkmap,candmap,heuristicString,vertexorder,
        standardVertexOrder,tempmat,flag,len))
        return false;
    return true;
}

// --------------------------------------------------------------------------
// Function checkIsomorphic()
// --------------------------------------------------------------------------

template <typename TSpec>
bool checkIsomorphic(seqan::Graph<TSpec> const & g1, seqan::Graph<TSpec> const & g2)
{
    typedef seqan::String<bool> TMat;
    TMat leastmat1, leastmat2;

    getCanonicalLabel(leastmat1, g1);
    getCanonicalLabel(leastmat2, g2);

    if (compare_(leastmat1, leastmat2) != 0)
        return false;
    return true;
}

#endif  // #ifndef APPS_SGIP_SGIP_BASE_H_
