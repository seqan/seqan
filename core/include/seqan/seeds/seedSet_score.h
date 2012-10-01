// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_SEEDSET_SCORE_H
#define SEQAN_HEADER_SEEDSET_SCORE_H


namespace SEQAN_NAMESPACE_MAIN
{

template<typename TValue, typename TScore>
inline TScore
_calculateScoringValue(TValue qPos1, 
					   TValue dPos1,
					   TValue qPos2,
					   TValue dPos2,
					   Score<TScore, Simple> const &matrix, 
					   Manhattan)
{
	return (qPos2 - qPos1 + dPos2 - dPos1-2)*scoreGapExtend(matrix)+scoreGapOpen(matrix);
}

template<typename TValue, typename TScore>
inline TScore
_calculateScoringValue(TValue qPos1, 
					   TValue&,
					   TValue qPos2, 
					   TValue &, 
					   Score<TScore, Simple> const &matrix, 
					   QueryDistance)
{
	return (qPos2 - qPos1 -1)*scoreGapExtend(matrix)+scoreGapOpen(matrix);
}

template<typename TValue, typename TScore>
inline TScore
_calculateScoringValue(TValue &,
					   TValue dPos1,
					   TValue &, 
					   TValue dPos2, 
					   Score<TScore, Simple> const &matrix, 
					   DatabaseDistance)
{
	return (dPos2 - dPos1-1)*scoreGapExtend(matrix)+scoreGapOpen(matrix);
}

template<typename TValue, typename TSeed, typename TScore>
inline TValue
_calculateScoringValue(TSeed &, 
					   TValue &, 
					   TValue &, 
					   Score<TScore, Simple> &, 
					   NoGapCost)
{
	return 0;
}


template<typename TValue, typename TSeed, typename TScore>
inline TValue
_calculateScoringValue(TSeed const &seed, 
					   TValue qPos, 
					   TValue dPos, 
					   TValue &, 
					   TScore score, 
					   Manhattan)
{
    return score + rightDim0(seed) - qPos + rightDim1(seed) - dPos;
}

template<typename TValue, typename TSeed, typename TScore>
inline TValue
_calculateScoringValue(TSeed const &seed, 
					   TValue qPos, 
					   TValue &,
					   TValue &,
					   TScore score,
					   QueryDistance)
{
    return score + rightDim0(seed) - qPos;
}

template<typename TValue, typename TSeed, typename TScore>
inline TValue
_calculateScoringValue(TSeed const &seed, 
					   TValue &, 
					   TValue dPos, 
					   TValue &,
					   TScore score,
					   DatabaseDistance)
{
    return score + rightDim1(seed) - dPos;
}

template<typename TValue, typename TSeed, typename TScore>
inline TValue
_calculateScoringValue(TSeed &, 
					   TValue &,
					   TValue &, 
					   TValue &, 
					   TScore &, 
					   NoGapCost)
{
    return score;
}


///.Metafunction.Value.param.T.type:Class.SeedSet
///.Metafunction.Value.class:Class.SeedSet
template <typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
struct ScoringScheme<SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> >
{
	typedef TSpecScoring Type;
};


/**
.Spec.Scored SeedSet
..summary:SeedSet that uses scored Seeds.
..cat:Seed Handling
..general:Class.SeedSet
..signature:SeedSet<TPosition, TSeedSpec, TScoringScheme, TSpec>
..param.TPosition: Type that saves the positions and upper/lower bounds.
...remarks: Positive and negative values are needed.
..param.TSeedSpec:The @Class.Seed@ specialization.
..param.TScoringScheme:The scoring sheme to use.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec> 
class SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec>
{

public:
    static const unsigned int BLOCKSIZE = BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, void> >::Value;
    typedef typename Size<String<TValue, Block<BLOCKSIZE> > >::Type TSize;
    typedef PropertyMap<typename ScoreType<TScoringSpec>::Type, Block<BLOCKSIZE> > PropMap;
    MemoryManager<Seed<TValue, TSeedSpec>, Block<BLOCKSIZE>,FreeMemoryInt> manager;
    std::multimap<TValue,  TSize> fragmentMap;
    std::set<TSize> result;
    TValue maxDistance;
    typename ScoreType<TScoringSpec>::Type qualityValue;

    PropMap scoreMap;
    Score<int, Simple> scoreMatrix;
	TValue last;

    SeedSet(){
	SEQAN_CHECKPOINT
		last=0;
    }

	SeedSet(TValue maxDistance, TValue qualityValue):maxDistance(maxDistance),qualityValue(qualityValue){
	//SaEQAN_CHECKPOINT
		last = 0;
    }


    SeedSet(TValue maxDistance, TValue qualityValue, Score<int,Simple> matrix):maxDistance(maxDistance),qualityValue(qualityValue),scoreMatrix(matrix){
	SEQAN_CHECKPOINT
		last=0;
    }

    ~SeedSet(){
	SEQAN_CHECKPOINT
	clear(manager);
    }
};


//Only difference is the destructor.
template<typename TValue, typename TSpec, typename TScoringSpec> 
class SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec>
{

public:
	TValue last;
    static const unsigned int BLOCKSIZE = BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, void> >::Value;
    typedef typename Size<String<TValue, Block<BLOCKSIZE> > >::Type TSize;
    MemoryManager<Seed<TValue, ChainedSeed>, Block<BLOCKSIZE>, FreeMemoryInt > manager;
    std::multimap<TValue, TSize> fragmentMap;
    std::set<TSize> result;
    TValue maxDistance;
    typename ScoreType<TScoringSpec>::Type qualityValue;
    PropertyMap<typename ScoreType<TScoringSpec>::Type, Block<BLOCKSIZE> > scoreMap;
    Score<int, Simple> scoreMatrix;

SeedSet(){
    SEQAN_CHECKPOINT
}

SeedSet(TValue maxDistance, TValue qualityValue, Score<int,Simple> matrix):maxDistance(maxDistance),qualityValue(qualityValue),scoreMatrix(matrix){
    SEQAN_CHECKPOINT
}


~SeedSet()
{
	
    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCKSIZE> > >::Type TSize;
    std::set<TSize> del;

	typedef typename std::multimap<TValue,  TSize>::iterator TIterator;
	TIterator it_end = fragmentMap.end();
	for (TIterator it = fragmentMap.begin(); it != it_end;++it)
	{
		TSize id = it->second;
		manager[id].~Seed<TValue, ChainedSeed>();
		result.erase(id);
	}

	typedef typename std::set<TSize>::iterator TIterator2;
	TIterator2 it_end2 = result.end();
	for (TIterator2 it = result.begin(); it != it_end2;++it)
	{
		manager[*it].~Seed<TValue, ChainedSeed>();
	}
	
	clear(scoreMap);
	clear(manager);
}

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//										   Container Functions	                                                   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//.Function.append.param.object.type:Class.SeedSet
//.Function.append.class:Class.SeedSet
template<typename TValue, typename TSeedSpec, typename TContainer1, typename TContainer2, typename TSpec, typename TScoringSpec>
void
append(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &target, 
	   TContainer1 &seedSource, 
	   TContainer2 &scoreSource)
{
    SEQAN_CHECKPOINT
    typedef typename Iterator<TContainer1>::Type TIterator;
    typename Iterator<TContainer2>::Type itScore = begin(scoreSource);
	
	TIterator it_end = end(seedSource);
    for (TIterator itSeed = begin(seedSource); itSeed!= it_end; ++itSeed)
	{
		addSeed(target, *itSeed, *itScore,Single());
		++itScore;
    }
}

template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec, typename TValue2>
void
appendValue(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &target, 
	    Seed<TValue, TSeedSpec> &seed, 
	    TValue2 score)
{
    SEQAN_CHECKPOINT
    addSeed(target,seed,score,Single());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//									         Standard Methods													   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
inline void
setScoreMatrix(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, Score<int,Simple> matrix)
{
    SEQAN_CHECKPOINT
    set.scoreMatrix = matrix;
}


template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
inline Score<int,Simple>
getScoreMatrix(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set)
{
    SEQAN_CHECKPOINT
    return set.scoreMatrix;
}

template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
void
clear(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set)
{
    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type TSize;
    std::set<TSize> del;
    TSize x = (obtainID(set.manager));
    while (x < length(set.manager)-1) 
	{
		del.insert(x);
		x = obtainID(set.manager);
    }
    del.insert(-1);
    typename std::set<TSize>::iterator it = del.begin();
    for (TSize i = 0; i < length(set.manager)-1; ++i){
		if (i != *it)
			set.manager[i].~Seed<TValue, TSeedSpec>();
		else
			++it;
		}
    set.fragmentMap.clear();
    set.result.clear();
    clear(set.scoreMap);
    clear(set.manager);
}


template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
inline PropertyMap<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> >*
_scoreMap(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set)
{
    return &set.scoreMap;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                    Addition of a Single new Seed                                               //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TScoringSpec, typename TPosition, typename TSize2, typename TScore, typename TSeedSpec>
bool
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
	TPosition qPos, 
	TPosition dPos, 
	TSize2 length, 
	TScore score, 
	Single)
{
    SEQAN_CHECKPOINT
	typedef SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> TSeedSet;
	typedef String<TValue, Block< BLOCK_SIZE<TSeedSet>::Value> > TBlockString;
    typedef typename Size<TBlockString>::Type TSize;
    typedef typename QualityFactor<TScoringSpec>::Type TQuality;

	TSize position = obtainID(set.manager);
    if (set.manager.change)
		raiseMemory(set.scoreMap);

    new (&set.manager[position]) Seed<TValue, TSeedSpec>(qPos,dPos,length);
    set.scoreMap[position] = score;
    if (_qualityReached(set.manager[position],set.scoreMap[position], set.qualityValue, TQuality()))
		set.result.insert(position);

	set.fragmentMap.insert( std::make_pair(static_cast<int>(dPos-qPos), position));
    return true;
}


//Score calculation
template<typename TValue, typename TSpec, typename TScoringSpec, typename TPosition, typename TSize, typename TSeedSpec>
inline bool
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
    TPosition qPos, 
    TPosition dPos, 
    TSize length, 
    Single)
{
	SEQAN_CHECKPOINT
	return addSeed(set, qPos, dPos, length, scoreMatch(set.scoreMatrix)*length, Single());
}


template<typename TValue, typename TSpec, typename TScoringSpec, typename TPosition, typename TScore>
bool
addSeed(SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> &set, 
	TPosition qBeginPos, 
	TPosition dBeginPos, 
	TPosition qEndPos, 
	TPosition dEndPos, 
	TScore score, 
	Single)
{
    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
    typedef typename QualityFactor<TScoringSpec>::Type TQuality;
   
	//Seed<TValue, SimpleSeed>* pSeed;
    TSize position = obtainID(set.manager);
    if (set.manager.change)
		raiseMemory(set.scoreMap);

    new (&set.manager[position]) Seed<TValue, SimpleSeed>(qBeginPos,dBeginPos,qEndPos,dEndPos);
    set.scoreMap[position] = score;
    if (_qualityReached(set.manager[position],set.scoreMap[position], set.qualityValue ,TQuality()))
		set.result.insert(position);

    set.fragmentMap.insert( std::pair<TValue, TSize >( dEndPos-qEndPos, position ) );
    return true;
}


template<typename TValue, typename TSpec, typename TScoringSpec, typename TScore>
bool
addSeed(SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> &set, 
	Seed<TValue, SimpleSeed> const &seed, 
	TScore score, 
	Single)
{
    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
    typedef typename QualityFactor<TScoringSpec>::Type TQuality;
   
	TSize position = obtainID(set.manager);
    if (set.manager.change)
		raiseMemory(set.scoreMap);

    new (&set.manager[position]) Seed<TValue, SimpleSeed>(leftDim0(seed), leftDim1(seed), rightDim0(seed), rightDim1(seed));
    set.scoreMap[position] = score;
    if (_qualityReached(set.manager[position],set.scoreMap[position], set.qualityValue ,TQuality()))
		set.result.insert(position);

    setLeftDiagonal(set.manager[position], leftDiagonal(seed));
    setRightDiagonal(set.manager[position], rightDiagonal(seed));
    set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(seed), position ) );
    return true;
}


template<typename TValue, typename TSpec, typename TScoringSpec, typename TScore>
bool
addSeed(SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> &set, 
	Seed<TValue, ChainedSeed> const &seed, 
	TScore score, 
	Single)
{
    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
    typedef typename QualityFactor<TScoringSpec>::Type TQuality;
    typedef typename std::list<Triple<TValue, TValue, TValue> >::const_iterator TIterator;
	
	TSize position = obtainID(set.manager);
    if (set.manager.change)
		raiseMemory(set.scoreMap);

    new (&set.manager[position]) Seed<TValue, ChainedSeed>();
    set.scoreMap[position] = score;
	TIterator it_end = _getDiagSet(seed).end();
    for (TIterator it = _getDiagSet(seed).begin(); it != it_end ; ++it)
		appendDiag(set.manager[position],*it);

    if (_qualityReached(set.manager[position],set.scoreMap[position], set.qualityValue ,TQuality()))
		set.result.insert(position);

    setLeftDiagonal(set.manager[position], leftDiagonal(seed));
    setRightDiagonal(set.manager[position], rightDiagonal(seed));
    set.fragmentMap.insert( std::pair<TValue, TSize >( endDiagonal(set.manager[position]), position));
    return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//										       Merging												              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSeedSpec, typename TSpec, typename TPosition, typename TSize2, typename TScoringSpec>
bool 
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set,
	TPosition qPos,
	TPosition dPos,
	TSize2 length,
	typename ScoreType<TScoringSpec>::Type score,
	int gapDistance,
	Merge)
{

    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type TSize;
    typedef typename std::multimap<TValue,TSize >::iterator TIterator;
    typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
    typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	TIterator tmpIt = _findSeedsMerge(set, (TValue) qPos, (TValue) dPos, (TValue) length, score, gapDistance);
    if (tmpIt != set.fragmentMap.end())
	{
		TSize position = tmpIt->second;
		set.fragmentMap.erase(tmpIt);

		_mergeTwoSeedsScore(set.manager[position], set.scoreMap[position], (TValue)qPos, (TValue)dPos, (TValue)length, score, getScoreMatrix(set), TGapCosts(), Merge());
		if (_qualityReached(set.manager[position],set.scoreMap[position], set.qualityValue, TQualityFactor()))
			set.result.insert(position);

		set.fragmentMap.insert(std::pair<TValue, TSize>( endDiagonal(set.manager[position]),position));
		return true;
    }
    return false;
}

//Score calculation
template<typename TValue, typename TSeedSpec, typename TSpec, typename TPosition, typename TSize, typename TScoringSpec>
inline bool 
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
	TPosition qPos, 
	TPosition dPos, 
	TSize length, 
	int gapDistance, 
	Merge)
{
	SEQAN_CHECKPOINT
	return addSeed(set, qPos, dPos, length, length*scoreMatch(set.scoreMatrix), gapDistance, Merge());
}

template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
bool
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
	    Seed<TValue, TSeedSpec> const &seed, 
	    typename ScoreType<TScoringSpec>::Type score, 
	    int gapDistance, 
	    Merge)
{
    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type TSize;
    typedef Seed<TValue, TSeedSpec> * pSeed;
    typename std::multimap<TValue,TSize >::iterator tmpIt = _findSeedsMerge(set, leftDim0(seed), leftDim1(seed), length(seed), score, gapDistance);
    typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
    typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;

    if (tmpIt != set.fragmentMap.end()){
		TSize position = tmpIt->second;
		set.fragmentMap.erase(tmpIt);
		_mergeTwoSeedsScore(set.manager[position], set.scoreMap[position], seed, score, getScoreMatrix(set), TGapCosts(),  Merge());
		if (_qualityReached(set.manager[position], set.scoreMap[position], set.qualityValue, TQualityFactor()))
			set.result.insert(position);

		set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[position]),position));
		return true;
    }
    return false;
}


template<typename TValue, typename TSpec, typename TScoringSpec>
bool 
addSeed(SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> &set, 
	TValue qlPos, 
	TValue dlPos, 
	TValue qrPos, 
	TValue drPos, 
	typename ScoreType<TScoringSpec>::Type score, 
	int gapDistance, 
	Merge)
{
    SEQAN_CHECKPOINT
    typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
    typedef Seed<TValue, SimpleSeed> * pSeed;
    typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
    typedef typename std::multimap<TValue,TSize >::iterator TIterator;
    typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
    
	TIterator tmpIt = _findSeedsMerge(set, qlPos, dlPos, qrPos-qlPos+1, score, gapDistance);
    if (tmpIt != set.fragmentMap.end()){
		TSize position = tmpIt->second;
		TValue x = endDiagonal(set.manager[position]);
		_mergeTwoSeedsScore(set.manager[position], set.scoreMap[position], qlPos, dlPos, qrPos, drPos, score, getScoreMatrix(set), TGapCosts(), Merge());
		if (_qualityReached(set.manager[position],set.scoreMap[position], set.qualityValue, TQualityFactor()))
			set.result.insert(position);
		
		if(x != drPos-qrPos)
		{
			set.fragmentMap.erase(tmpIt);
			set.fragmentMap.insert(std::pair<TValue, TSize>( endDiagonal(set.manager[position]),position));
		}
		set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[position]),position));
		return true;
    }
    return false;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Chaining Algorithms                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


								////////////////////////////////////
								//	     Standard Chaining		  //
								////////////////////////////////////



template<typename TValue, typename TSpec, typename TScoringSpec, typename TSeedSpec>
bool 
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
	TValue qPos,
	TValue dPos,
	TValue length, 
	typename ScoreType<TScoringSpec>::Type score, 
	int gapDistance, 
	SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;

	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qPos, dPos, length, score, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		set.scoreMap[id] += score + _calculateScoringValue(rightDim0(set.manager[id]), rightDim1(set.manager[id]),qPos, dPos, getScoreMatrix(set), TGapCosts());
		_mergeTwoSeedsScore(set.manager[id], qPos, dPos, length, SimpleChain());
		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue, TQualityFactor()))
			set.result.insert(id);
		
		//new endDiagonal
		if(x != dPos-qPos)
		{
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>(endDiagonal(set.manager[id]),id));		
		}
		return true;
	}
	return false;
}

//score calculation
template<typename TValue, typename TSpec, typename TScoringSpec, typename TSeedSpec>
inline bool 
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		TValue qPos, 
		TValue dPos, 
		TValue length,
		int gapDistance,
		SimpleChain)
{
	SEQAN_CHECKPOINT
	return addSeed(set, qPos, dPos, length, length*scoreMatch(set.scoreMatrix), gapDistance, SimpleChain());
}


template<typename TValue, typename TSpec, typename TScoringSpec>
bool 
addSeed(SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> &set,
		TValue qlPos,
		TValue dlPos,
		TValue qrPos,
		TValue drPos, 
		typename ScoreType<TScoringSpec>::Type score,
		int gapDistance,
		SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qlPos, dlPos, qrPos-qlPos+1, score, gapDistance);
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		set.scoreMap[id] += score + _calculateScoringValue(rightDim0(set.manager[id]), rightDim1(set.manager[id]),qlPos, dlPos, getScoreMatrix(set), TGapCosts());
		_mergeTwoSeedsScore(set.manager[id], qlPos, dlPos, qrPos, drPos, SimpleChain());
		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue, TQualityFactor()))
			set.result.insert(id);
		if(x != drPos-qrPos){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert(std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	}
	return false;
}

template<typename TValue, typename TSpec, typename TScoringSpec, typename TSeedSpec>
bool 
addSeed(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		Seed<TValue, TSeedSpec> const &seed, 
		typename ScoreType<TScoringSpec>::Type score, 
		int gapDistance, 
		SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, leftDim0(seed), leftDim1(seed), length(seed), score, gapDistance);
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		set.scoreMap[id] += score + _calculateScoringValue(rightDim0(set.manager[id]), rightDim1(set.manager[id]), leftDim0(seed), leftDim1(seed), getScoreMatrix(set), TGapCosts());
		_mergeTwoSeedsScore(set.manager[id], leftDim0(seed), leftDim1(seed), rightDim0(seed), rightDim1(seed), SimpleChain());
		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue, TQualityFactor()))
			set.result.insert(id);
		if(x != endDiagonal(seed)){
			set.fragmentMap.erase(id);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	}
	return false;
}


template<typename TValue, typename TSpec, typename TScoringSpec>
bool 
addSeed(SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> &set, 
		Seed<TValue, ChainedSeed> const &seed, 
		typename ScoreType<TScoringSpec>::Type score, 
		int gapDistance, 
		SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
        TIterator it = _findSeedsChain(set, leftDim0(seed), leftDim1(seed), length(seed), score, gapDistance);
	if (it != set.fragmentMap.end()){
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		set.scoreMap[id] += score + _calculateScoringValue(rightDim0(set.manager[id]), rightDim1(set.manager[id]),leftDim0(seed), leftDim1(seed), getScoreMatrix(set), TGapCosts());
		_mergeTwoSeedsScore(set.manager[id], leftDim0(seed), leftDim1(seed), _getFirstDiag(seed).i3, SimpleChain());
		typedef typename std::list<Triple<TValue,TValue,TValue> >::const_iterator TIterator2;
		TIterator2 it2_end = _getDiagSet(seed).end();
		for (TIterator2 it2 = ++_getDiagSet(seed).begin(); it2 != it2_end; ++it2)
			appendDiag(set.manager[id],*it2);

		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue ,TQualityFactor()))
			set.result.insert(id);

		if(x != endDiagonal(seed)){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	}
	return false;
}


template<typename TValue, typename TSeedSpec>
void
_mergeTwoSeedsScore(Seed<TValue, TSeedSpec> &seed, 
					TValue qPos, 
					TValue dPos, 
					TValue length, 
					SimpleChain)
{
	SEQAN_CHECKPOINT
	setRightDim0(seed, qPos+length-1);
	setRightDim1(seed, dPos+length-1);
	TValue diag = dPos-qPos;
	if (leftDiagonal(seed) < diag)
		setLeftDiagonal(seed, diag);
	if (rightDiagonal(seed) > diag)
		setRightDiagonal(seed, diag);
}


template<typename TValue>
void
_mergeTwoSeedsScore(Seed<TValue, SimpleSeed> &firstSeed, 
					TValue qlPos, 
					TValue dlPos, 
					TValue qrPos, 
					TValue drPos, 
					SimpleChain)
{
	SEQAN_CHECKPOINT
	setRightDim0(firstSeed, qrPos);
	setRightDim1(firstSeed, drPos);
	TValue diag = dlPos-qlPos;
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
	diag = drPos-qrPos;
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
}


								////////////////////////////////////
								//				Chaos			  //
								////////////////////////////////////


template<typename TValue, typename TText, typename TSpec, typename TPosition, typename TSize2, typename TScoringSpec, typename TSpecSeed>
bool 
addSeed(SeedSet<TValue, TSpecSeed, TScoringSpec, TSpec> &set, 
		TPosition qPos, 
		TPosition dPos, 
		TSize2 length, 
		typename ScoreType<TScoringSpec>::Type score, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Chaos)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;

	TIterator it = _findSeedsChain(set, (TValue)qPos, (TValue)dPos, (TValue)length, score, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		typename ScoreType<TScoringSpec>::Type tmpScore = _mergeTwoSeedsScore(set.manager[id], (TValue)qPos, (TValue)dPos, (TValue)length, query, database, getScoreMatrix(set), TGapCosts(), Chaos());
		set.scoreMap[id] += score + tmpScore;
		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue ,TQualityFactor()))
			set.result.insert(id);

		if(x+qPos != dPos){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert(std::pair<TValue, TSize>(endDiagonal(set.manager[id]),id));
		}
		return true;
	}
	return false;
}

template<typename TValue, typename TText, typename TSpec, typename TPosition, typename TSize, typename TScoringSpec, typename TSpecSeed>
inline bool 
addSeed(SeedSet<TValue, TSpecSeed, TScoringSpec, TSpec> &set, 
		TPosition qPos, 
		TPosition dPos, 
		TSize length,
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Chaos)
{
	SEQAN_CHECKPOINT
	return addSeed(set, qPos, dPos, length, length*scoreMatch(set.scoreMatrix), query, database, gapDistance, Chaos());
}


template<typename TValue, typename TText, typename TSpec, typename TScoringSpec>
bool
addSeed(SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> &set, 
		Seed<TValue, ChainedSeed> const &seed, 
		typename ScoreType<TScoringSpec>::Type score, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Chaos)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	
	TIterator it = _findSeedsChain(set, leftDim0(seed), leftDim1(seed), length(seed), score, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		TValue tmpScore = _mergeTwoSeedsScore(set.manager[id], leftDim0(seed), leftDim1(seed), _getFirstDiag(seed).i3, query, database, getScoreMatrix(set), TGapCosts(), Chaos());
		set.scoreMap[id] += score + tmpScore;
		typedef typename std::list<Triple<TValue,TValue,TValue> >::const_iterator TIterator2;
		TIterator2 it2_end = _getDiagSet(seed).end();
		for (TIterator2 it2 = ++_getDiagSet(seed).begin(); it2 != it2_end; ++it2)
			appendDiag(set.manager[id],*it2);

		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue ,TQualityFactor()))
			set.result.insert(id);

		if(x != endDiagonal(seed)){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	}
	return false;
}

template<typename TValue, typename TText, typename TSpec, typename TScoringSpec>
bool
addSeed(SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> &set, 
		Seed<TValue, SimpleSeed> const &seed, 
		typename ScoreType<TScoringSpec>::Type score, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Chaos)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	
	TIterator it = _findSeedsChain(set, leftDim0(seed), leftDim1(seed), length(seed), score, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		TValue tmpScore = _mergeTwoSeedsScore(set.manager[id], leftDim0(seed), leftDim1(seed), 1, query, database, getScoreMatrix(set), TGapCosts(), Chaos());
		set.scoreMap[id] += score + tmpScore;
		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue ,TQualityFactor()))
			set.result.insert(id);

		if(x == endDiagonal(seed)){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		setRightDim0(set.manager[id], rightDim0(seed));
		setRightDim1(set.manager[id], rightDim1(seed));
		return true;
	}
	return false;
}

//SimpleSeed
template<typename TValue, typename TScore, typename TText, typename TGapCost>
TValue
_mergeTwoSeedsScore(Seed<TValue, SimpleSeed>  &firstSeed, 
					TValue qPos, 
					TValue dPos, 
					TValue length,
					String<TText> const &query,
					String<TText> const &database,
					Score<TScore, Simple> const &scoreMatrix,
					TGapCost tag, 
					Chaos)
{
	SEQAN_CHECKPOINT
	typename std::list<Triple <TValue, TValue, TValue> >::iterator begin1, end2, it;
	TValue databaseGap = dPos - rightDim1(firstSeed)-1;
	TValue queryGap = qPos - rightDim0(firstSeed)-1;

	TValue gap = abs((TScore)(databaseGap-queryGap));
	TValue currentScore = 0;
	if (gap == 0)
	{
		TValue de = rightDim1(firstSeed);
		for (int i = rightDim0(firstSeed)+1; i <qPos;++i){
			currentScore += score(scoreMatrix,i,++de,query,database);
		}
	} else {
		TValue lPositionQuery = rightDim0(firstSeed);
		TValue lPositionDatabase = rightDim1(firstSeed);
		TValue rPositionQuery = qPos;
		TValue rPositionDatabase = dPos;

		TValue gap = (databaseGap < queryGap) ? databaseGap : queryGap;
		for (int i = 0; i <gap;++i){
			currentScore += score(scoreMatrix,--rPositionQuery,--rPositionDatabase,query,database);
		}
		TValue tmpScore = currentScore;
		TValue tmpLength = 0;

		for (int i = 0; i < gap ; ++i){
			currentScore += score(scoreMatrix,++lPositionQuery,++lPositionDatabase,query,database) - score(scoreMatrix,rPositionQuery++,rPositionDatabase++,query,database);
			if (currentScore > tmpScore){
				tmpScore = currentScore;

				tmpLength = i+1;
			}
		}
		setRightDim0(firstSeed,rightDim0(firstSeed)+tmpLength);
		currentScore += _calculateScoringValue(rightDim0(firstSeed), rightDim1(firstSeed),qPos-gap+tmpLength, dPos-gap+tmpLength, scoreMatrix, tag);
	}
	setRightDim0(firstSeed,qPos + length-1);
	setRightDim1(firstSeed,dPos + length-1);

	TValue diag = dPos-qPos;
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
	return currentScore;
}


//ChainedSeed
template<typename TValue, typename TScore, typename TText, typename TGapCost>
TValue
_mergeTwoSeedsScore(Seed<TValue, ChainedSeed>  &firstSeed, 
					TValue qPos, 
					TValue dPos, 
					TValue length,
					String<TText> const &query,
					String<TText> const &database,
					Score<TScore, Simple> const &scoreMatrix,
					TGapCost tag, 
					Chaos)
{
	SEQAN_CHECKPOINT
	TValue databaseGap = dPos - rightDim1(firstSeed)-1;
	TValue queryGap = qPos - rightDim0(firstSeed)-1;

	TValue gap = abs(databaseGap-queryGap);
	TValue currentScore = 0;
	if (gap == 0){
		//cout << "interessant" << endl;
		TValue de = rightDim1(firstSeed);
		for (int i = rightDim0(firstSeed)+1; i <qPos;++i){
			currentScore += score(scoreMatrix,i,++de,query,database);
		}
		setRightDim0(firstSeed,qPos + length-1);
	} 
	else 
	{
		TValue lPositionQuery = rightDim0(firstSeed);
		TValue lPositionDatabase = rightDim1(firstSeed);
		TValue rPositionQuery = qPos;
		TValue rPositionDatabase = dPos;

		TValue gap = (databaseGap < queryGap)? databaseGap : queryGap;
		for (int i = 0; i <gap;++i){
			currentScore += score(scoreMatrix,--rPositionQuery,--rPositionDatabase,query,database);
		}
		TValue tmpScore = currentScore;
		TValue tmpLength = 0;

		for (int i = 0; i < gap; ++i){
			currentScore += score(scoreMatrix,++lPositionQuery,++lPositionDatabase,query,database) - score(scoreMatrix,rPositionQuery++,rPositionDatabase++,query,database);
			if (currentScore > tmpScore){
				tmpScore = currentScore;
				tmpLength = i+1;
			}
		}
		setRightDim0(firstSeed,rightDim0(firstSeed)+tmpLength);
		currentScore += _calculateScoringValue(rightDim0(firstSeed), rightDim1(firstSeed),qPos-gap+tmpLength, dPos-gap+tmpLength, scoreMatrix, tag);//abs(databaseGap-queryGap)*scoreGap(scoreMatrix);

		firstSeed.seedSet.push_back(Triple<TValue,TValue,TValue>(qPos-gap+tmpLength, dPos - gap +tmpLength,gap-tmpLength+length));
	}

	TValue diag = dPos-qPos;
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);


	return currentScore;
}





								////////////////////////////////////
								//				Blat			  //
								////////////////////////////////////

template<typename TValue, typename TText, typename TSpec, typename TScoringSpec, typename TSpecSeed>
bool 
addSeed(SeedSet<TValue, TSpecSeed, TScoringSpec, TSpec> &set, 
		TValue qPos, 
		TValue dPos, 
		TValue length, 
		typename ScoreType<TScoringSpec>::Type score, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Blat)
{
	SEQAN_CHECKPOINT
	typedef typename  Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qPos, dPos, score, length, gapDistance);
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		int dLength = dPos - rightDim1(set.manager[id]);
		int qLength = qPos - rightDim0(set.manager[id]);
		int dLog = (int) ceil(log((double)dLength));
		int qLog = (int) ceil(log((double)qLength));
		int k;
		int maxValue = (dLog > qLog) ? dLog : qLog;
		if ((maxValue < dLength) && (maxValue < qLength))
			k = maxValue;
		else 
		    k = (dLog < qLog) ? dLog : qLog;
		typename ScoreType<TScoringSpec>::Type tmpScore = _mergeTwoSeedsScore(set.manager[id], qPos, dPos, length, getScoreMatrix(set), query, database, k, TGapCosts(), Blat());
		set.scoreMap[id] += score + tmpScore;
		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue ,TQualityFactor()))
			set.result.insert(id);

		if(endDiagonal(set.manager[id]) != dPos-qPos){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	}
	return false;
}

//same as above but score is calculated
template<typename TValue, typename TText, typename TSpec, typename TScoringSpec, typename TSpecSeed>
inline bool 
addSeed(SeedSet<TValue, TSpecSeed, TScoringSpec, TSpec> &set,
		TValue qPos, 
		TValue dPos, 
		TValue length, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance,
		Blat)
{
	SEQAN_CHECKPOINT
	return addSeed(set, qPos, dPos, length, length*scoreMatch(set.scoreMatrix), query, database, gapDistance, Blat());
}


template<typename TValue, typename TText, typename TSpec, typename TScoringSpec>
bool 
addSeed(SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> &set,
		Seed<TValue, ChainedSeed> const &seed,
		typename ScoreType<TScoringSpec>::Type score,
		String<TText> const &query, 
		String<TText> const &database,
		int gapDistance, 
		Blat)
{
	SEQAN_CHECKPOINT
	typedef typename  Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	typedef typename ScoreType<TScoringSpec>::Type TScore;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	TValue qPos = leftDim0(seed);
	TValue dPos = leftDim1(seed);
	TValue length_ = _getFirstDiag(seed).i3;
        TIterator it = _findSeedsChain(set, qPos, dPos, score, length(seed), gapDistance);
	if (it != set.fragmentMap.end()){
		TSize id = it->second;
		int dLength = dPos - rightDim1(set.manager[id]);
		int qLength = qPos - rightDim0(set.manager[id]);
		int dLog = (int) ceil(log((double)dLength));
		int qLog = (int) ceil(log((double)qLength));
		int k;
		int maxValue = (dLog > qLog) ? dLog : qLog;
		if ((maxValue < dLength) && (maxValue < qLength))
			k = maxValue;
		else 
			k = (dLog < qLog) ? dLog : qLog;
		TScore tmpScore = _mergeTwoSeedsScore(set.manager[id], qPos, dPos, length_, getScoreMatrix(set), query, database, k, TGapCosts(), Blat());
		set.scoreMap[id] += score + tmpScore;
		typename std::list<Triple<TValue,TValue,TValue> >::const_iterator seedIt = _getDiagSet(seed).begin();
		++seedIt;
		while (seedIt !=_getDiagSet(seed).end()){
			appendDiag(set.manager[id],*seedIt);
			++seedIt;
		}
		if (rightDiagonal(seed) < rightDiagonal(set.manager[id]))
			setRightDiagonal(set.manager[id], rightDiagonal(seed));

		if (leftDiagonal(seed) < leftDiagonal(set.manager[id]))
			setLeftDiagonal(set.manager[id], leftDiagonal(seed));

		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue ,TQualityFactor()))
			set.result.insert(id);

		if(endDiagonal(set.manager[id]) != endDiagonal(seed)){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} 
	return false;
}

template<typename TValue, typename TText, typename TSpec, typename TScoringSpec>
bool 
addSeed(SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> &set,
		Seed<TValue, SimpleSeed> const &seed,
		typename ScoreType<TScoringSpec>::Type score,
		String<TText> const &query, 
		String<TText> const &database,
		int gapDistance, 
		Blat)
{
	SEQAN_CHECKPOINT
	typedef typename  Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	typedef typename ScoreType<TScoringSpec>::Type TScore;
	typedef typename GapCosts<TScoringSpec>::Type TGapCosts;
	typedef typename QualityFactor<TScoringSpec>::Type TQualityFactor;
	TValue qPos = leftDim0(seed);
	TValue dPos = leftDim1(seed);
        TIterator it = _findSeedsChain(set, qPos, dPos, score, length(seed), gapDistance);
	if (it != set.fragmentMap.end()){
		TSize id = it->second;
		int dLength = dPos - rightDim1(set.manager[id]);
		int qLength = qPos - rightDim0(set.manager[id]);
		int dLog = (int) ceil(log((double)dLength));
		int qLog = (int) ceil(log((double)qLength));
		int k;
		int maxValue = (dLog > qLog) ? dLog : qLog;
		if ((maxValue < dLength) && (maxValue < qLength))
			k = maxValue;
		else 
			k = (dLog < qLog) ? dLog : qLog;
		TScore tmpScore = _mergeTwoSeedsScore(set.manager[id], qPos, dPos, 1, getScoreMatrix(set), query, database, k, TGapCosts(), Blat());
		set.scoreMap[id] += score + tmpScore;
	
		setRightDim0(set.manager[id],rightDim0(seed));
		setRightDim1(set.manager[id],rightDim1(seed));

		if (rightDiagonal(seed) < rightDiagonal(set.manager[id]))
			setRightDiagonal(set.manager[id], rightDiagonal(seed));

		if (leftDiagonal(seed) < leftDiagonal(set.manager[id]))
			setLeftDiagonal(set.manager[id], leftDiagonal(seed));

		if (_qualityReached(set.manager[id],set.scoreMap[id], set.qualityValue ,TQualityFactor()))
			set.result.insert(id);

		if(endDiagonal(set.manager[id]) != endDiagonal(seed)){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} 
	return false;
}


template<typename TValue, typename TText, typename TScore, typename TGapCosts>
TScore
_mergeTwoSeedsScore(Seed<TValue, SimpleSeed>  &firstSeed, 
					TValue qPos, 
					TValue dPos, 
					TValue length, 
					Score<TScore, Simple> const &scoreMatrix, 
					String<TText> const &query, 
					String<TText> const &database, 
					TValue q, 
					TGapCosts tag,
					Blat)
{
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple<TValue, TValue, TValue> > TPieceList;
	typedef typename TPieceList::iterator TIterator;
	Triple<TValue, TValue, TValue> r(qPos,dPos,length);
	TValue leftDiag = leftDiagonal(firstSeed);
	TValue rightDiag = rightDiagonal(firstSeed);
	TValue diag;
	TPieceList tmp;
	_gapFill(rightDim0(firstSeed)+1, rightDim1(firstSeed)+1, qPos, dPos, tmp, query, database,q);
	TScore tmpScore = 0;

	if (tmp.size() > 0){
		TIterator it1 = tmp.begin();
		TIterator it2 = ++tmp.begin();
		diag = it1->i2 - it1->i1;
		if (diag >leftDiag) 
			leftDiag = diag;
		if (diag < rightDiag)
			rightDiag = diag;
		tmpScore += _calculateScoringValue(rightDim0(firstSeed), rightDim1(firstSeed), it1->i1, it1->i2, scoreMatrix, tag);
		for (int i = 0; i < it1->i3; ++i){
			tmpScore += score(scoreMatrix,it1->i1+i,it1->i2+i,query,database);
		}
		while (it2 != tmp.end()){
			diag = it2->i2-it2->i1;
			if (diag >leftDiag) 
				leftDiag = diag;
			if (diag < rightDiag)
				rightDiag = diag;
			tmpScore += _calculateScoringValue(it1->i1 + it1->i3 - 1, it1->i2 + it1->i3 - 1, it2->i1, it2->i2, scoreMatrix, tag);
			tmpScore += it2->i3*scoreMatch(scoreMatrix);
			it1++;
			it2++;
		}	
		tmpScore += _calculateScoringValue(it1->i1 + it1->i3 - 1, it1->i2 + it1->i3 - 1, qPos, dPos, scoreMatrix, tag);
	}else {
		tmpScore += _calculateScoringValue(rightDim0(firstSeed), rightDim1(firstSeed), qPos, dPos, scoreMatrix, tag);
	}
	
	setRightDim0(firstSeed, qPos + length - 1);
	setRightDim1(firstSeed, dPos + length - 1);
	setRightDiagonal(firstSeed, rightDiag);
	setLeftDiagonal(firstSeed, leftDiag);
	return tmpScore;
}


template<typename TValue, typename TText, typename TScore, typename TGapCosts>
TScore
_mergeTwoSeedsScore(Seed<TValue, ChainedSeed>  &firstSeed, 
					TValue qPos, 
					TValue dPos, 
					TValue length, 
					Score<TScore, Simple> const &scoreMatrix, 
					String<TText> const &query, 
					String<TText> const &database, 
					TValue q, 
					TGapCosts tag,
					Blat)
{
	SEQAN_CHECKPOINT
	std::list<Triple<TValue, TValue, TValue> > tmp;
	Triple<TValue, TValue, TValue> r(qPos,dPos,length);
	TValue leftDiag = leftDiagonal(firstSeed);
	TValue rightDiag = rightDiagonal(firstSeed);
	TValue diag;
	_gapFill(rightDim0(firstSeed)+1, rightDim1(firstSeed)+1, qPos, dPos, tmp, query, database,q);
	TScore tmpScore = 0;

	typedef typename std::list<Triple<TValue, TValue, TValue> >::iterator TIterator;
	if (tmp.size() > 0){
		TIterator it1 = tmp.begin();
		TIterator it2 = ++tmp.begin();
		tmpScore += _calculateScoringValue(rightDim0(firstSeed), rightDim1(firstSeed), it1->i1, it1->i2, scoreMatrix, tag);
		for (int i = 0; i < it1->i3; ++i){
			tmpScore += score(scoreMatrix,it1->i1+i,it1->i2+i,query,database);
		}
		while (it2 != tmp.end()){
			tmpScore += _calculateScoringValue(it1->i1+it1->i3-1,it1->i2+it1->i3-1,it2->i1,it2->i2, scoreMatrix, tag);
			tmpScore += it2->i3*scoreMatch(scoreMatrix);
			it1++;
			it2++;
		}	
		tmpScore += _calculateScoringValue(it1->i1+it1->i3-1,it1->i2+it1->i3-1,qPos,dPos,scoreMatrix,tag);
	}else {
		tmpScore += _calculateScoringValue(rightDim0(firstSeed), rightDim1(firstSeed), qPos, dPos, scoreMatrix,tag);
	}

	if (tmp.size() > 0){
		Triple<TValue, TValue, TValue> x = *tmp.begin();
		if ((x.i1 == rightDim0(firstSeed)+1)&& (x.i2 == rightDim1(firstSeed)+1)){
			tmp.front().i1 = _getDiagSet(firstSeed).back().i1;
			tmp.front().i2 = _getDiagSet(firstSeed).back().i2;
			tmp.front().i3 += _getDiagSet(firstSeed).back().i3;
			_getDiagSet(firstSeed).pop_back();
		}
		x = tmp.back();
		if ((x.i1+x.i3 == qPos)&& (x.i2+x.i3 == dPos)){
			qPos = x.i1;
			dPos = x.i2;
			length += x.i3;
			tmp.pop_back();
		}
	}
	_getDiagSet(firstSeed).splice(_getDiagSet(firstSeed).end(),tmp); 
	_getDiagSet(firstSeed).push_back(Triple<TValue,TValue,TValue>(qPos,dPos,length));
	typedef typename std::list<Triple<TValue, TValue, TValue> >::iterator TIterator;
	TIterator it_end = _getDiagSet(firstSeed).end();
	for (TIterator it = _getDiagSet(firstSeed).begin(); it != it_end; ++it){
		diag = it->i2-it->i1;
		if (diag >leftDiag) 
			leftDiag = diag;
		if (diag < rightDiag)
			rightDiag = diag;
	}
	setRightDiagonal(firstSeed, rightDiag);
	setLeftDiagonal(firstSeed, leftDiag);
	return tmpScore;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                    Addition of several new Seeds                                               //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TValue, typename TSeedSpec, typename TIterator, typename TIterator2, typename TAlgoSpec, typename TSpec, typename TScoringSpec>
bool
addSeeds(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		 TIterator begin, 
		 TIterator end, 
		 TIterator2 scoreBegin, 
		 int gapDistance, 
		 TAlgoSpec tag)
{
	SEQAN_CHECKPOINT
	std::multimap<TValue, std::pair<TIterator, TIterator2> > tmpMap; //zum sortieren
	for(TIterator it = begin; it !=end; ++it)
	{
		tmpMap.insert(std::pair<TValue, std::pair<TIterator, TIterator2> >(leftDim0(*it),std::pair<TIterator, TIterator2>(it,scoreBegin)));
		++scoreBegin;
	}
	
	typedef typename std::multimap<TValue, std::pair<TIterator, TIterator2> >::iterator TIterator3;
	TIterator3 it_end = tmpMap.end();
	for (TIterator3 it = tmpMap.begin(); it != it_end; ++it)
		if (!addSeed(set, *it->second.first, *it->second.second, gapDistance, tag))
			addSeed(set, *it->second.first, *it->second.second, Single());
	
	tmpMap.clear();
	return true;
}


template<typename TValue, typename TSeedSpec, typename TIterator, typename TIterator2, typename TSpec, typename TScoringSpec>
bool
addSeeds(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		 TIterator begin, 
		 TIterator end,
		 TIterator2 scoreBegin, 
		 Single)
{
	SEQAN_CHECKPOINT
	for(TIterator it = begin; it !=end; ++it){
		addSeed(set, *it, *scoreBegin, Single());
		++scoreBegin;
	}
	return true;
}

template<typename TValue, typename TSeedSpec, typename TIterator, typename TIterator2, typename TText, typename TSpec, typename TScoringSpec>
bool
addSeeds(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		 TIterator begin, 
		 TIterator end, 
		 TIterator2 scoreBegin, 
		 String<TText> const &query, 
		 String<TText> const &database, 
		 int gapDistance, 
		 Chaos)
{
	SEQAN_CHECKPOINT
	std::multimap<TValue, std::pair<TIterator, TIterator2> > tmpMap;
	while(begin!= end){
		tmpMap.insert(std::pair<TValue, std::pair<TIterator, TIterator2> >(leftDim0(*begin),std::pair<TIterator, TIterator2>(begin,scoreBegin)));
		++begin;
	}
	typedef typename std::multimap<TValue, std::pair<TIterator, TIterator2> >::iterator TIterator3;
	TIterator3 it_end = tmpMap.end();
	for (TIterator3 it = tmpMap.begin(); it != it_end; ++it)
		if (!addSeed(set, *it->second.first, *it->second.second, query, database, gapDistance, Chaos()))
			addSeed(set, *it->second.first, *it->second.second, Single());

	tmpMap.clear();
	return true;
}

template<typename TValue, typename TSeedSpec, typename TIterator, typename TIterator2, typename TText, typename TSpec, typename TScoringSpec>
bool
addSeeds(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		 TIterator begin, 
		 TIterator end, 
		 TIterator2 scoreBegin, 
		 String<TText> const &query, 
		 String<TText> const &database, 
		 int gapDistance, 
		 Blat)
{
	SEQAN_CHECKPOINT
	typename std::multimap<TValue, std::pair<TIterator, TIterator2 > > tmpMap; //zum sortieren
	while(begin!=end){
		tmpMap.insert(std::pair<TValue, std::pair<TIterator, TIterator2> >(leftDim0(*begin),std::pair<TIterator, TIterator2>(begin,scoreBegin)));
		++begin;
	}
	typename std::multimap<TValue, std::pair<TIterator, TIterator2> >::iterator tmpIt;
	for (tmpIt = tmpMap.begin(); tmpIt != tmpMap.end(); ++tmpIt){
		if (!addSeed(set, *tmpIt->second.first, *tmpIt->second.second, query, database, gapDistance, Blat()))
			addSeed(set, *tmpIt->second.first, *tmpIt->second.second, Single());
	}
	tmpMap.clear();
	return true;
}

template<typename TValue, typename TSeedSpec, typename TContainer, typename TContainer2, typename TAlgoSpec, typename TSpec, typename TScoringSpec>
bool
addSeeds(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		 TContainer const &source, 
		 TContainer2 const &scoreSource, 
		 int gapDistance, 
		 TAlgoSpec tag)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer2,Standard>::Type TIterator2;
    typedef typename Iterator<TContainer const,Standard>::Type TIterator;
    TIterator it2 = end(source);
	TIterator2 scoreBegin = begin(scoreSource);

	typename std::multimap<TValue, std::pair<TIterator, TIterator2> > tmpMap; //zum sortieren
	for (typename Iterator<TContainer const,Standard>::Type it1 = begin(source); it1 != it2; ++it1)
		tmpMap.insert(std::pair<TValue, std::pair<TIterator, TIterator2> >(leftDim0(*it1),std::pair<TIterator, TIterator2>(it1,scoreBegin)));
	
	typedef typename std::multimap<TValue, std::pair<TIterator, TIterator2> >::iterator TIterator3;
	TIterator3 it_end = tmpMap.end();
	for (TIterator3 it = tmpMap.begin(); it != it_end; ++it)
		if (!addSeed(set, *it->second.first, *it->second.second, gapDistance, tag))
			addSeed(set, *it->second.first, *it->second.second, Single());
	
	tmpMap.clear();
	return true;
}

template<typename TValue, typename TSeedSpec, typename TContainer, typename TContainer2, typename TSpec, typename TScoringSpec>
bool
addSeeds(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
	TContainer const &source, 
	TContainer2 &scoreSource, 
	Single){
    SEQAN_CHECKPOINT
    typedef typename Iterator<TContainer const,Standard>::Type TIterator;
    typename Iterator<TContainer2,Standard>::Type itScore = begin(scoreSource);
    for (TIterator it = begin(source); it != end(source); ++it){
		addSeed(set, *it, *itScore, Single());
		++itScore;
    }
    return true;
}

template<typename TValue, typename TSeedSpec, typename TContainer, typename TContainer2, typename TText, typename TSpec, typename TScoringSpec, typename TAlgoSpec>
bool
addSeeds(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
		 TContainer const &source, 
		 TContainer2 const &scoreSource, 
		 String<TText> const &query, 
		 String<TText> const &database, 
		 int gapDistance,
		 TAlgoSpec algoSpec)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer const,Standard>::Type TIterator;
        typedef typename Iterator<TContainer2 const, Standard>::Type TIterator2;
	TIterator it2 =end(source);
	TIterator2 itScore = begin(scoreSource);
	typename std::multimap<TValue, std::pair<TIterator, TIterator2> > tmpMap; //zum sortieren
	for (TIterator it1 = begin(source); it1 !=it2; ++it1)
		tmpMap.insert(std::pair<TValue, std::pair<TIterator, TIterator2> >(leftDim0(*it1),std::pair<TIterator, TIterator2>(it1,itScore)));
	
	typedef typename std::multimap<TValue, std::pair<TIterator, TIterator2> >::iterator TIterator3;
	TIterator3 it_end = tmpMap.end();
	for (TIterator3 it = tmpMap.begin(); it != it_end; ++it)
		if (!addSeed(set, *it->second.first, *it->second.second, query, database, gapDistance, algoSpec))
			addSeed(set, *it->second.first, *it->second.second, Single());
	
	tmpMap.clear();
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		                                    Addition of a SeedSet	                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
bool
addSeedSet(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &target, 
		   SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> const &source, 
		   Single)
{
	SEQAN_CHECKPOINT
	typedef SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet const,Standard>::Type TIterator;
	for (TIterator it = begin(source); it != end(source); ++it)
		addSeed(target, *it, seedScore(it), Single());
	
	return true;
}


template<typename TValue, typename TSeedSpec, typename TAlgoSpec, typename TSpec, typename TScoringSpec>
bool
addSeedSet(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &target, 
		   SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> const &source, 
		   int gapDistance, 
		   TAlgoSpec tag)
{
	SEQAN_CHECKPOINT
	typedef SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet const,Standard>::Type TIterator;
	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;

	std::multimap<TValue,TIterator> tmpMap;
	TIterator it1_end = end(source);
	for (TIterator it1 = begin(source); it1 != it1_end; ++it1)
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it1),it1));

	TIterator2 it2_end = tmpMap.end();
	for (TIterator2 it2 = tmpMap.begin(); it2 != it2_end; ++it2)
		if (!addSeed(target, *it2->second, seedScore(it2->second), gapDistance, tag))
			addSeed(target, *it2->second, seedScore(it2->second), Single());

	tmpMap.clear();
	return true;
}

template<typename TValue, typename TSeedSpec, typename TText, typename TSpec, typename TScoringSpec, typename TAlgoSpec>
bool
addSeedSet(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &target, 
		   SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> const &source, 
		   String<TText> const &query, 
		   String<TText> const &database, 
		   int gapDistance, 
		   TAlgoSpec tag)
{
	SEQAN_CHECKPOINT
	typedef SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet const,Standard>::Type TIterator;
	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;

	std::multimap<TValue, TIterator > tmpMap; //zum sortieren
	TIterator it1_end = end(source);
	for (TIterator it1 = begin(source); it1 !=it1_end; ++it1)
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it1),it1));
	
	TIterator2 it2_end = tmpMap.end();
	for (TIterator2 it2 = tmpMap.begin(); it2 != it2_end; ++it2)
		if (!addSeed(target, *it2->second, seedScore(it2->second), query, database, gapDistance, tag))
			addSeed(target, *it2->second, seedScore(it2->second), Single());
	
	tmpMap.clear();
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//									Delete to far away														  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename TValue, typename TSeedSpec, typename TQualityFactor, typename TGapCosts, typename TScore>
void
_deleteEverything(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> >, void> &deletionTarget, 
				  TValue currentPos)
{
	typedef typename Size<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> >, void> >::Type TSize;
	typedef std::multimap<TValue, TSize > TMap;
	typename TMap::iterator it_end = deletionTarget.fragmentMap.end();
	TSize pos;
	typedef typename std::set<TSize>::iterator TSetIterator;
	TSetIterator set_end = deletionTarget.result.end();
	for (typename TMap::iterator it = deletionTarget.fragmentMap.begin(); it != it_end;)
	{
		pos = it->second;
		if (currentPos > (deletionTarget.maxDistance + rightDim0(deletionTarget.manager[pos])))
		{
			deletionTarget.fragmentMap.erase(it++);
			if (set_end == deletionTarget.result.find(pos))
			{
				valueDestruct(&deletionTarget.manager[pos]);
				releaseID(deletionTarget.manager,pos);
			}
		}
		else
			++it;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//								Finding seed for merging/chaining											  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
typename std::multimap<TValue, 	typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type >::iterator
_findSeedsChain(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set, 
				TValue qPos, 
				TValue dPos, 
				TValue length, 
				typename ScoreType<TScoringSpec>::Type score, 
				int gapDistance)
{
//	if(set.last != qPos){
//		_deleteEverything(set, qPos);
//		set.last = qPos;
//	}
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	TIterator itUp=set.fragmentMap.upper_bound(dPos-qPos+gapDistance);
	TIterator tmp = set.fragmentMap.end();
	typename ScoreType<TScoringSpec>::Type maxScore = minValue<TValue>();
	typename ScoreType<TScoringSpec>::Type tmpScore;
	for (TIterator it = set.fragmentMap.lower_bound(dPos-qPos-gapDistance); it != itUp; )
	{
		TSize id = it-> second;
		if (qPos > set.maxDistance + rightDim0(set.manager[id]))
		{//delete it
			TIterator itdelete = it;
			++it;
			set.fragmentMap.erase(itdelete);
			if (set.result.end() == set.result.find(id))
			{
				valueDestruct(&set.manager[id]);
				releaseID(set.manager,id);
			}
		}
		else
		{
			if ((qPos > rightDim0(set.manager[id])) && (dPos > rightDim1(set.manager[id])))
			{
				tmpScore = _calculateScoringValue(set.manager[id], qPos, dPos, length, score, typename GapCosts<TScoringSpec>::Type());
				if (tmpScore > maxScore)
				{
					maxScore = tmpScore;
					tmp = it;
				}
			}
			++it;
		}
	}
	return tmp;
}


template<typename TValue, typename TSeedSpec, typename TSpec, typename TScoringSpec>
typename std::multimap<TValue, typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type >::iterator
_findSeedsMerge(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &set,
				TValue qPos, 
				TValue dPos, 
				TValue ,
				typename ScoreType<TScoringSpec>::Type score, 
				int gapDistance)
{
	SEQAN_CHECKPOINT
	if(set.last != qPos){
		_deleteEverything(set, qPos);
		set.last = qPos;
	}
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> >::Value> > >::Type TSize;
	typedef typename ScoreType<TScoringSpec>::Type TScore;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	
	TScore maxScore = minValue<TScore>();
	TScore tmpScore;
	TValue diag = dPos-qPos;

	TIterator tmpIt = set.fragmentMap.end();
	
	//new seed has higher diagonal number
	TIterator itUp=set.fragmentMap.upper_bound(diag); 
	for (TIterator it=set.fragmentMap.lower_bound(diag-gapDistance); it != itUp; ++it)
	{
		if (dPos <= rightDim1(set.manager[it->second]))
		{
			tmpScore = set.scoreMap[it->second] + score - (rightDim1(set.manager[it->second]) - dPos) * scoreMatch(set.scoreMatrix) +abs((TScore)(endDiagonal(set.manager[it->second])-diag))*scoreGap(set.scoreMatrix);
			if (tmpScore > maxScore)// &&(tmpScore > set.scoreMap[it->second]))
			{
				tmpIt = it;
				maxScore = tmpScore;
			}
		}
	}
	
	//new seed has lower diagnoal number
	itUp=set.fragmentMap.upper_bound(diag+ gapDistance);  
	for (TIterator it=set.fragmentMap.lower_bound (diag+1); it != itUp; ++it)
	{
		if (qPos <= rightDim0(set.manager[it->second]))
		{
			tmpScore = set.scoreMap[it->second] + score- (rightDim0(set.manager[it->second]) - qPos) * scoreMatch(set.scoreMatrix) +abs((TScore)(endDiagonal(set.manager[it->second])-diag))*scoreGapExtend(set.scoreMatrix);
			if (tmpScore > maxScore)//&&(tmpScore > set.scoreMap[it->second]))
			{
				tmpIt = it;
				maxScore = tmpScore;
			}
		}
	}

	return tmpIt;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//										Extension Algorithms														//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
.Function.extendSeedScore
..class:Class.Seed
..summary:Extends a seed and increases the score.
..cat:Seed Handling
..signature:extendSeedScore(seed, score, query, database, direction, tag)
..signature:extendSeedScore(seed, score, scoreDropOff, scoreMatrix, query, database, direction, tag)
..param.seed: The seed to extend.
...type:Class.Seed
..param.score: A reference to the score of the seed. This will be increased by the score of the extension.
..param.query: The query sequence.
...type:Class.String
..param.query: The database sequence.
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended. 0 = left, 1 = right, 2 = both
..param.scoreDropOff: The score drop after which the extension should stop. The extension stops if this value is exceeded.
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.scoreMatrix: The scoring scheme.
...type:Spec.Simple Score
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.tag: The algorithm to use.
...type:Tag.Seed Extension.MatchExtend
...type:Tag.Seed Extension.UngappedXDrop
...type:Tag.Seed Extension.GappedXDrop
..include:seqan/seeds.h
*/

/**
.Function.extendSeedsScore
..summary: Extension of seeds with score calculation.
..cat:Seed Handling
..signature:extendSeedsScore(container, scores, query, database, direction, MatchExtend)
..signature:extendSeedsScore(begin, end, itScore query, database, direction, MatchExtend)
..signature:extendSeedsScore(container, scores, scoreDropOff, scoreMatrix, query, database, direction, tag)
..signature:extendSeedsScore(begin, end, itScores, scoreDropOff, scoreMatrix, query, database, direction, tag)
..param.container: The container with the @Class.Seed@ objects to extend.
...type:Concept.Container
..param.scores: A container with the scores of the seeds. The scores will be increased by the score of the extension. It can be omitted if the seeds container is a @Class.SeedSet@ with scores.
...type:Concept.Container
..param.itScore: Iterator pointing to the first seeds score. The scores will be increased by the score of the extension.
..param.begin: Iterator pointing to the first value to add.
..param.end: Iterator pointing just behind the last value to add.
..param.query: The query sequence.
...type:Class.String
..param.query: The database sequence.
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended. 0 = left, 1 = right, 2 = both
..param.scoreDropOff: The score drop after which the extension should stop. The extension stops if this value is exceeded.
..param.scoreMatrix: The scoring sheme.
...type:Spec.Simple Score
..param.tag: The algorithm to use.
...type:Tag.Seed Extension.MatchExtend
...type:Tag.Seed Extension.tag.UngappedXDrop
...type:Tag.Seed Extension.tag.GappedXDrop
..include:seqan/seeds.h
*/

template<typename TContainer, typename TContainer2, typename TText, typename TValue, typename TDirection>
void
extendSeedsScore(TContainer &seedSet, 
				 TContainer2 &scores, 
				 Score<TValue, Simple> const &matrix,
				 String<TText> const &query, 
				 String<TText> const &database, 
				 TDirection direction, MatchExtend)
{
	SEQAN_CHECKPOINT

	typedef typename Iterator<TContainer,Standard>::Type TIterator;
	
	typename Iterator<TContainer2,Standard>::Type itScore = begin(scores);
	TIterator it_end = end(seedSet);
	for (TIterator it = begin(seedSet); it != it_end; ++it)
	{
		extendSeedScore(*it, *itScore, matrix, query, database, direction, MatchExtend());
		++itScore;
	}
}

template<typename TValue, typename TSeedSpec, typename TSpec, typename TText, typename TScoringSpec, typename TDirection>
void
extendSeedsScore(SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> &seedSet,
				 String<TText> const &query,
				 String<TText> const &database, 
				 TDirection direction, 
				 MatchExtend)
{
	SEQAN_CHECKPOINT
	typedef SeedSet<TValue, TSeedSpec, TScoringSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet,Standard>::Type TIterator;
	TIterator it_end = end(seedSet);
	for (TIterator it = begin(seedSet); it != it_end; ++it)
		extendSeedScore(*it, seedScore(it), getScoreMatrix(seedSet), query, database, direction, MatchExtend());
	
}

template<typename TIterator, typename TIterator2, typename TText, typename TValue, typename TDirection>
void
extendSeedsScore(TIterator begin, 
				 TIterator end, 
				 TIterator2 itScore,
				 Score<TValue, Simple> const &matrix,
				 String<TText> const &query, 
				 String<TText> const &database, 
				 TDirection direction,
				 MatchExtend)
{
	SEQAN_CHECKPOINT
	for (TIterator it = begin; it != end; ++it)
		extendSeedScore(*it, *itScore, matrix, query, database, direction, MatchExtend());	

}

template<typename TContainer, typename TText, typename TContainer2, typename TextendSeedSpec, typename TValue, typename TValue2, typename TDirection>
void
extendSeedsScore(TContainer &seedSet, 
				 TContainer2 &scores, 
				 TValue scoreDropOff, 
				 Score<TValue2, Simple> const &scoreMatrix,
				 String<TText> const &query, 
				 String<TText> const &database,
				 TDirection direction, 
				 TextendSeedSpec tag)
{
	SEQAN_CHECKPOINT
	typename Iterator<TContainer2,Standard>::Type itScore = begin(scores);
	typedef typename Iterator<TContainer,Standard>::Type TIterator;
	
	TIterator it_end = end(seedSet);
	for (TIterator it = begin(seedSet); it != it_end; ++it)
	{
		extendSeedScore(*it, *itScore, scoreDropOff, scoreMatrix, query, database, direction, tag);
		++itScore;
	}
}

template<typename TText, typename TExtendSeedSpec, typename TValue, typename TSpec, typename TScoringSpec, typename TDirection>
void
extendSeedsScore(SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> &seedSet,
				 TValue scoreDropOff, 
				 String<TText> const &query,
				 String<TText> const &database,
				 TDirection direction,
				 TExtendSeedSpec tag)
{
	SEQAN_CHECKPOINT
	typedef SeedSet<TValue, SimpleSeed, TScoringSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet,Standard>::Type TIterator;
	
	TIterator it_end = end(seedSet);
	for (TIterator it = begin(seedSet); it != it_end; ++it)
		extendSeedScore(*it, seedScore(it), scoreDropOff, getScoreMatrix(seedSet), query, database, direction, tag);

}

template<typename TText, typename TExtendSeedSpec, typename TValue, typename TSpec, typename TScoringSpec, typename TDirection>
void
extendSeedsScore(SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> &seedSet,
				 TValue scoreDropOff, 
				 String<TText> const &query,
				 String<TText> const &database,
				 TDirection direction, 
				 TExtendSeedSpec tag)
{
	SEQAN_CHECKPOINT
	typedef SeedSet<TValue, ChainedSeed, TScoringSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet,Standard>::Type TIterator;
	
	TIterator it_end = end(seedSet);
	for (TIterator it = begin(seedSet); it != it_end; ++it)
		extendSeedScore(*it, seedScore(it), scoreDropOff, getScoreMatrix(seedSet), query, database, direction, tag);

}

template<typename TIterator, typename TIterator2, typename TExtendSeedSpec, typename TValue, typename TText, typename TDirection>
void
extendSeedsScore(TIterator begin,
				 TIterator end, 
				 TIterator2 itScore,
				 TValue scoreDropOff,
				 Score<TValue, Simple> const &scoreMatrix,
				 String<TText> const &query, 
				 String<TText> const &database,
				 TDirection direction, 
				 TExtendSeedSpec tag)
{
	SEQAN_CHECKPOINT
	for (TIterator it = begin; it != end; ++it)
	{
		extendSeedScore(*it, *itScore, scoreDropOff, scoreMatrix, query, database, direction, tag);		
		++itScore;
	}

}

template<typename TValue, typename TSeedSpec, typename TText, typename TValue2, typename TValue3, typename TDirection>
void
extendSeedScore(Seed<TValue, TSeedSpec> &seed, 
				TValue3 &currentScore, 
				Score<TValue2, Simple> const &matrix, 
				String<TText> const &query,
				String<TText> const &database, 
				TDirection direction, 
				MatchExtend)
{
	SEQAN_CHECKPOINT

	//left extension
	if (direction != 1){
		int tmpLength = 0;
		TValue queryPos =leftDim0(seed) ;
		TValue dataPos = leftDim1(seed);
		while ((queryPos-1>=0) && (dataPos-1>=0) && (query[queryPos-1] == database[dataPos-1])){
			--queryPos;
			--dataPos;
			++tmpLength;
		}
		setLeftDim0(seed,queryPos);
		setLeftDim1(seed,dataPos);
		currentScore += tmpLength * scoreMatch(matrix);
	}

	//right extension
	if (direction != 0){
		int tmpLength = 0;
		int queryLength = length(query);
		int databaseLength = length(database);
		TValue queryPos =rightDim0(seed) ;
		TValue dataPos = rightDim1(seed);
		while ((queryPos+1 < queryLength) && (dataPos+1 < databaseLength) && (query[queryPos+1] == database[dataPos+1])){
			++queryPos;
			++dataPos;
			++tmpLength;
		}
		setRightDim0(seed,queryPos);
		setRightDim1(seed,dataPos);
		currentScore += tmpLength * scoreMatch(matrix);
	}
}


template<typename TValue, typename TSeedSpec, typename TText, typename TScore, typename TDirection>
void  
extendSeedScore(Seed<TValue,TSeedSpec> &seed, 
				TScore &currentScore, 
				TScore scoreDropOff,
				Score<TScore, Simple> const &scoreMatrix, 
				String<TText> const &query, 
				String<TText> const &database,
				TDirection direction, 
				UngappedXDrop)
{
	SEQAN_CHECKPOINT
	scoreDropOff *=-1;
	TValue tmpScore = 0;
	TValue tmp;
	//left extension
	if (direction != 1){
		TValue xPos = leftDim0(seed)-1;
		TValue yPos = leftDim1(seed)-1;
		TValue last = 0;
	
		while ((tmpScore > scoreDropOff) && (xPos >= 0) && (yPos>=0)){
			if (query[xPos] == database[yPos]){
				last = 0;
				tmp = score(scoreMatrix, xPos, yPos, query, database);
				tmpScore += tmp;
				currentScore += tmp;
				if (tmpScore > 0)
					tmpScore = 0;
			} else{
				tmp = score(scoreMatrix, xPos, yPos, query, database);
				tmpScore += tmp;
				currentScore += tmp;
				++last;
			}
			--xPos;
			--yPos;
		}
		currentScore -= last*scoreMismatch(scoreMatrix);
		setLeftDim0(seed,xPos+last+1);
		setLeftDim1(seed,yPos+last+1);
	}

	//right extension
	if (direction != 0){
		TValue xLength = length(query);
		TValue yLength = length(database);
		TValue xPos = rightDim0(seed)+1;
		TValue yPos = rightDim1(seed)+1;

		TValue last = 0;
		tmpScore= 0;
		while ((tmpScore > scoreDropOff) && (xPos < xLength) && (yPos < yLength)){
			if (query[xPos] == database[yPos]){
				last = 0;
				tmp = score(scoreMatrix, xPos, yPos, query, database);
				tmpScore += tmp;
				currentScore += tmp;
				if (tmpScore > 0)
					tmpScore = 0;
			}else{
				tmp = score(scoreMatrix, xPos, yPos, query, database);
				tmpScore += tmp;
				currentScore += tmp;
				++last;
			}
			++xPos;
			++yPos;
		}
		currentScore -= last*scoreMismatch(scoreMatrix);
		setRightDim0(seed,xPos-last-1);
		setRightDim1(seed,yPos-last-1);
	}
}

template<typename TSeedSpec, typename TPosition, typename TQuery, typename TDatabase, typename TScore, typename TSize>
void 
extendSeedScore(Seed<TPosition, TSeedSpec> &seed, 
           TScore &currentScore, 
		   TScore scoreDropOff, 
		   Score<TScore, Simple> const &scoreMatrix, 
		   TQuery const &query, 
		   TDatabase const &database, 
		   TSize direction, 
		   GappedXDrop)
{
	SEQAN_CHECKPOINT
	//left extension
	if ((direction != 1) && (leftDim0(seed)!= (TPosition)beginPosition(query)) && (leftDim1(seed) != (TPosition)beginPosition(database)))
    {
        typename Prefix<TDatabase const>::Type dataSeg = prefix(database, leftDim1(seed));
        typename Prefix<TQuery const>::Type querySeg = prefix(query, leftDim0(seed));
		
        currentScore += _extendSeedOneDirection(seed, dataSeg, querySeg, scoreDropOff, scoreMatrix, 0 /*left*/);
	}

	//right extension
	if ((direction != 0) && (rightDim0(seed)+1 < (TPosition)endPosition(query)) && (rightDim1(seed)+1 < (TPosition)endPosition(database)))
    {
        typename Suffix<TDatabase const>::Type dataSeg = suffix(database, rightDim1(seed)+1);
        typename Suffix<TQuery const>::Type querySeg = suffix(query, rightDim0(seed)+1);

        currentScore += _extendSeedOneDirection(seed, dataSeg, querySeg, scoreDropOff, scoreMatrix, 1 /*right*/);
	}
}

} //namespace Seqan


#endif //#ifndef SEQAN_HEADER_
