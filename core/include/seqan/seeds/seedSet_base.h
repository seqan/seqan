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

#ifndef SEQAN_HEADER_SEEDSET_BASE_H
#define SEQAN_HEADER_SEEDSET_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

// Forward declaration.
struct SeedIterator;

/*
.Function:_calculateScoringValue
..summary: calculates a Value for a gap during chaining.
..param.seed: The first seed.
..param.qPos: The leftPosition in dimension 0.
..param.dPos: The leftPosition in dimension 1.
..param.length: The length of the seed. ?
..param.tag: The tags define the distance to use.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSeed>
inline TValue
_calculateScoringValue(TSeed const &seed, TValue qPos, TValue dPos, TValue &, Manhattan){
	return rightDim0(seed) - qPos + rightDim1(seed) - dPos;
}

template<typename TValue, typename TSeed>
inline TValue
_calculateScoringValue(TSeed const &seed, TValue qPos, TValue &, TValue &, QueryDistance){
	return rightDim0(seed) - qPos;
}

template<typename TValue, typename TSeed>
inline TValue
_calculateScoringValue(TSeed const &seed, TValue &, TValue dPos, TValue &, DatabaseDistance){
	return rightDim1(seed) - dPos;
}

template<typename TValue, typename TSeed>
inline TValue
_calculateScoringValue(TSeed &, TValue &, TValue &, TValue &, NoGapCost){
	return 0;
}


/*
.Function:_qualityReached
..summary: Checks if a seed has reached a certain quality value.
..param.seed: The first seed.
..param.score: The score of the seed. (If no scores are used, the value will be 0.
..param.qualityValue: The value that shall be reached.
..param.tag: The tags define the kind of quality check to use.
..include:seqan/seeds.h
*/

template<typename TSeed, typename TScore>
inline bool
_qualityReached(TSeed &, TScore score, TScore qualityValue, SeedScore){
	return score >= qualityValue;
}

template<typename TSeed, typename TScore>
inline bool
_qualityReached(TSeed const &seed, TScore , TScore qualityValue, SeedLength){
	return length(seed) >= qualityValue;
}



/**
.Class.SeedSet:
..summary:Manages seeds for local chaining and merging algorithms.
..cat:Seed Handling
..signature:SeedSet<TPosition, TSeedSpec, TScoringScheme, TSpec>
..param.TPosition: Type that saves the positions and upper/lower bounds.
...remarks: Positive and negative values are needed.
..param.TSeedSpec:The @Class.Seed@ specialization.
..param.TScoringScheme:The scoring sheme to use.
..include:seqan/seeds.h
*/

template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec = void> 
class SeedSet;


/**
.Spec.NonScored SeedSet
..summary:SeedSet that uses seeds without scores.
..cat:Seed Handling
..general:Class.SeedSet
..signature:SeedSet<TPosition, TSeedSpec, TScoringScheme<TQualityFactor, TGapCosts, void>, TSpec>
..param.TPosition: Type that saves the positions and upper/lower bounds.
...remarks: Positive and negative values are needed.
..param.TSeedSpec:The @Class.Seed@ specialization.
..param.TScoringScheme:The scoring sheme to use.
..param.TQualityFactor: The quality factor to use.
..param.TGapCosts: The kind of GapCosts to use.
..include:seqan/seeds.h
*/

/**
.Memfunc.SeedSet#SeedSet:
..class:Class.SeedSet
..summary:Constructor
..signature: SeedSet<TValue, TSeedSpec> ()
..signature: SeedSet<TValue, TSeedSpec> (maxDistance, qualityValue)
..param.maxDistance:Maximum distance between two seeds.
..param.qualityValue:Minimum length of a seeed.
*/
template<typename TValue, typename TSeedSpec, typename TQualityFactor, typename TGapCosts> 
class SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, void>
{

public:
	static const unsigned int BLOCKSIZE = BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, void> >::Value;
	typedef typename Size<String<TValue, Block<BLOCKSIZE> > >::Type TSize;
	MemoryManager<Seed<TValue, TSeedSpec>, Block<BLOCKSIZE>, FreeMemoryInt > manager;
	std::multimap<TValue,  TSize> fragmentMap;
	std::set<TSize> result;
	TValue maxDistance;
	TValue qualityValue;
	TValue last;

	SeedSet(){
		last=0;
		SEQAN_CHECKPOINT
	}

	SeedSet(TValue maxDistance, TValue qualityValue):maxDistance(maxDistance),qualityValue(qualityValue){
		SEQAN_CHECKPOINT
		last=0;
	}
	
	~SeedSet(){
		SEQAN_CHECKPOINT
		clear(manager);
	}
};


// Only difference is the destructor.
template<typename TValue, typename TQualityFactor, typename TGapCosts> 
class SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, void>
{

public:
	static const unsigned int BLOCKSIZE = BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, void> >::Value;
	typedef typename Size<String<TValue, Block<BLOCKSIZE> > >::Type TSize;
	MemoryManager<Seed<TValue, ChainedSeed>, Block<BLOCKSIZE>, FreeMemoryInt > manager;
	std::multimap<TValue, TSize> fragmentMap;
	std::set<TSize> result;
	TValue maxDistance;
	TValue qualityValue;
	TValue last;

	SeedSet(){
	SEQAN_CHECKPOINT
		last=0;
	}

	SeedSet(TValue maxDistance, TValue qualityValue):maxDistance(maxDistance),qualityValue(qualityValue){
	SEQAN_CHECKPOINT
		last=0;
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
	clear(manager);
	}
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											   Metafunctions	                                                   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
struct BLOCK_SIZE<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> >
{
	enum {Value = 2048};
};

///.Metafunction.Spec.param.T.type:Class.SeedSet
///.Metafunction.Spec.class:Class.SeedSet
template <typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
struct Spec<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec > >
{
	typedef TSeedSpec Type;
};

///.Metafunction.Value.param.T.type:Class.SeedSet
///.Metafunction.Value.class:Class.SeedSet
template <typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
struct Value<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> >
{
	typedef Seed<TValue, TSeedSpec> Type;
};


///.Metafunction.Iterator.param.T.type:Class.SeedSet
///.Metafunction.Iterator.class:Class.SeedSet
template <typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
struct Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec>, Standard >
{
	typedef Iter<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec>, SeedIterator> Type;
};

template <typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
struct Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> const, Standard >
{
	typedef Iter<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> const, SeedIterator> Type;
};

template <typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
struct Size<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> >
{
	typedef typename Size< String<TValue, Block<BLOCK_SIZE<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> >::Value > > >::Type Type;
};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//										   Container Functions	                                                   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//.Function.append.param.object.type:Class.SeedSet
//.Function.append.class:Class.SeedSet
template<typename TValue, typename TSeedSpec, typename TContainer2, typename TSpecScoring, typename TSpec>
void
append(SeedSet<TValue, TSeedSpec, 
	   TSpecScoring, TSpec> &target, 
	   TContainer2 &source)
{
    SEQAN_CHECKPOINT
    typename Iterator<TContainer2, Standard>::Type it;
    for (it = begin(source); it != end(source);++it)
	addSeed(target, *it ,Single());
}

template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline void
appendValue(SeedSet<TValue, 
			TSeedSpec, TSpecScoring, TSpec> &target, 
			Seed<TValue, TSeedSpec> &seed)
{
	SEQAN_CHECKPOINT
	addSeed(target,seed,Single());
}

///.Function.end.param.object.type:Class.SeedSet
///.Function.end.class:Class.SeedSet
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline typename Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec>, Standard >::Type
begin(SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> &set)
{
	SEQAN_CHECKPOINT
	return typename Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec>, Standard >::Type(set, set.result.begin());
}

///.Function.begin.param.object.type:Class.SeedSet
///.Function.begin.class:Class.SeedSet
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline typename Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec>, Standard >::Type
end(SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> &set)
{
    SEQAN_CHECKPOINT
    return typename Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec>, Standard >::Type(set, set.result.end());
}


///.Function.end.param.object.type:Class.SeedSet
///.Function.end.class:Class.SeedSet
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline typename Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> const, Standard>::Type
begin(SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> const &set)
{
	SEQAN_CHECKPOINT
	return typename Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> const, Standard>::Type(set, set.result.begin());
}

///.Function.begin.param.object.type:Class.SeedSet
///.Function.begin.class:Class.SeedSet
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline typename Iterator<SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> const, Standard>::Type
end(SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> const &set)
{
	SEQAN_CHECKPOINT
	return typename Iterator<SeedSet<TValue,TSeedSpec, TSpecScoring, TSpec> const, Standard>::Type(set, set.result.end());
}

template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline int
length(SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> &set)
{
	SEQAN_CHECKPOINT
	return set.result.size();
}

template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline int
length(SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> const &set)
{
	return set.result.size();
}

template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
void
clear(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec,Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> > , TSpec> >::Value> > >::Type TSize;
	std::set<TSize> del;
	TSize x = (obtainID(set.manager));
	while (x < length(set.manager)-1) {
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
	clear(set.manager);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//									         Standard Methods													   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Function.setMaximumDistance:
..class:Class.SeedSet
..summary:Sets the maximal distance between two seed during a chaining process.
..cat:Seed Handling
..signature:setMaximumDistance(set, distance)
..param.set: The set of seeds.
...type:Class.SeedSet
..param.distance: The maximal distance between to seeds.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline void
setMaximumDistance(SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> &set, 
				   TValue distance)
{
	SEQAN_CHECKPOINT
	set.maxDistance = distance;
}


/**
.Function.setQualityValue:
..class:Class.SeedSet
..summary:Sets the minimum length for a seed to be saved permanently.
..cat:Seed Handling
..signature:setQualityValue(set, distance)
..param.set: The set of seeds.
...type:Class.SeedSet
..param.distance: The minimum length of a seed.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline void
setQualityValue(SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> &set, 
				TValue value)
{
	SEQAN_CHECKPOINT
	set.qualityValue = value;
}

/**
.Function.maximumDistance:
..class:Class.SeedSet
..summary:Sets the maximal distance between two seed during a chaining process.
..cat:Seed Handling
..signature:maximumDistance(set)
..param.set: The set of seeds.
...type:Class.SeedSet
..returns:The maximum distance.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline TValue
maximumDistance(SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> &set)
{
	SEQAN_CHECKPOINT
	return set.maxDistance;
}

/**
.Function.qualityValue:
..class:Class.SeedSet
..summary:Sets the minimum length for a seed to be saved permanently.
..cat:Seed Handling
..signature:qualityValue(set)
..param.set: The set of seeds.
...type:Class.SeedSet
..returns: The minimum length.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSeedSpec, typename TSpecScoring, typename TSpec>
inline TValue
qualityValue(SeedSet<TValue, TSeedSpec, TSpecScoring, TSpec> &set)
{
	SEQAN_CHECKPOINT
	return set.qualityValue;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       addition of new Seeds                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Function.addSeed:
..class:Class.SeedSet
..summary:Adds a seed to an existing set.
..cat:Seed Handling
..signature:addSeed(set, qPos, dPos, length, tag)
..signature:addSeed(set, qPos, dPos, qrPos, drPos, tag)
..signature:addSeed(set, seed, tag)
..param.set:The set to which the new seed should be added.
...type:Class.SeedSet
..param.qPos: Start position in sequence1.
..param.dPos: Start position in sequence2.
..param.length: Length of the seed.
..param.tag: The algorithm that should be used to add the new seed.
...type:Tag.Seed Adding
...remark: Note that not every algorithm can be used with each type of @Class.Seed@.
..param.qrPos: End position in sequence1.
..param.drPos: End Position in sequence2.
..param.seed: The new Seed.
...type:Class.Seed
...remarks: The seed is copied and then added.
..returns:Boolean if succesfully added.
...remarks:Always true for Tag Single.
..include:seqan/seeds.h
*/

/**
.Function.addSeeds:
..class:Class.SeedSet
..summary:Adds several seeds to an existing set. If a merging or chaining algorithm is used seeds are added if the merging or chaining fails.
..cat:Seed Handling
..signature:addSeed(set, container, tag)
..signature:addSeed(set, begin, end, tag)
..param.set:The set to which the new seed sould be added.
...type:Class.SeedSet
..param.container: Content is copied to set.
...type:Concept.Container
..param.begin: Iterator pointing to the first value to add.
..param.end: Iterator pointing just behind the last value to add.
..param.tag: The algorithm that should be used to add the new @Class.Seed@.
...type:Tag.Seed Adding
...remarks: Note that not every algorithm can be used with each specialization of @Class.Seed@.
..include:seqan/seeds.h
*/
template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeed(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qPos, 
		TValue dPos,
		TValue length, 
		Single)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	TSize position = obtainID(set.manager);
	new (&set.manager[position]) Seed<TValue, TSeedSpec>(qPos,dPos,length);
	if (_qualityReached(set.manager[position],0,qualityValue(set), TQualityFactor()))           
		set.result.insert(position);

	set.fragmentMap.insert( std::pair<TValue, TSize >(dPos-qPos, position) );
	return true;
}

template<typename TValue, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeed(SeedSet<TValue, SimpleSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qBeginPos, 
		TValue dBeginPos, 
		TValue qEndPos, 
		TValue dEndPos, 
		Single)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	//Seed<TValue, SimpleSeed>* pSeed;
	TSize position = obtainID(set.manager);
	new (&set.manager[position]) Seed<TValue, SimpleSeed>(qBeginPos,dBeginPos,qEndPos,dEndPos);
	if (_qualityReached(set.manager[position],0,qualityValue(set), TQualityFactor()))
		set.result.insert(position);

	set.fragmentMap.insert( std::pair<TValue, TSize >( dEndPos-qEndPos, position ) );
	return true;
}

template<typename TValue, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeed(SeedSet<TValue, SimpleSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		Seed<TValue, SimpleSeed> const &seed, 
		Single)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	TSize position = obtainID(set.manager);
	new (&set.manager[position]) Seed<TValue, SimpleSeed>(leftDim0(seed), leftDim1(seed), rightDim0(seed), rightDim1(seed));
	if (_qualityReached(set.manager[position],0,qualityValue(set), TQualityFactor()))
		set.result.insert(position);

	setLeftDiagonal(set.manager[position], leftDiagonal(seed));
	setRightDiagonal(set.manager[position], rightDiagonal(seed));
	set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(seed), position ) );
	return true;
}

template<typename TValue, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeed(SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		Seed<TValue, ChainedSeed> const &seed, 
		Single)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typedef typename std::list<Triple<TValue, TValue, TValue> >::const_iterator TIterator;
	TSize position = obtainID(set.manager);
	new (&set.manager[position]) Seed<TValue, ChainedSeed>();
	TIterator it_end = _getDiagSet(seed).end();
	for (TIterator it = _getDiagSet(seed).begin(); it != it_end; ++it)
		appendDiag(set.manager[position],*it);

	if (_qualityReached(set.manager[position], 0, qualityValue(set), TQualityFactor()))
		set.result.insert(position);

	setLeftDiagonal(set.manager[position], leftDiagonal(seed));
	setRightDiagonal(set.manager[position], rightDiagonal(seed));
	set.fragmentMap.insert(std::pair<TValue, TSize >( endDiagonal(set.manager[position]), position));
	return true;
}


template<typename TValue, typename TSeedSpec, typename TIterator, typename TAlgoSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeeds(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		 TIterator begin, 
		 TIterator end, 
		 int gapDistance, 
		 TAlgoSpec tag)
{
	SEQAN_CHECKPOINT
	std::multimap<TValue, TIterator> tmpMap; //zum sortieren
	for (TIterator it = begin; it!=end; ++it)
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it),it));

	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;
	TIterator2 it_end = tmpMap.end();
	for (TIterator2 it = tmpMap.begin(); it != it_end; ++it)
		if (!addSeed(set, *it->second, gapDistance, tag))
			addSeed(set, *it->second, Single());
	
	tmpMap.clear();
	return true;
}


template<typename TSeedSet, typename TIterator>//typename TValue, typename TSeedSpec, typename TIterator, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeedsIt(TSeedSet &set,//SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
	TIterator begin, 
	TIterator end,  Single)
{
	SEQAN_CHECKPOINT
	for(TIterator it = begin; it != end; ++it)
		addSeed(set, *it, Single());
	return true;
}


template<typename TValue, typename TSeedSpec, typename TIterator, typename TText, typename TValue2, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeeds(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		 TIterator begin, 
		 TIterator end, 
		 Score<TValue2, Simple> const &scoreMatrix, 
		 String<TText> const &query, 
		 String<TText> const &database, 
		 int gapDistance, 
		 Chaos)
{
	SEQAN_CHECKPOINT
	std::multimap<TValue, TIterator> tmpMap; 
	for(TIterator it = begin; it != end; ++it){
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it), it));
		++begin;
	}
	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;
	TIterator2 it_end = tmpMap.end();
	for (TIterator2 it = tmpMap.begin(); it != it_end; ++it)
		if (!addSeed(set, *it->second, scoreMatrix, query, database, gapDistance, Chaos()))
			addSeed(set, *it->second, Single());

	tmpMap.clear();
	return true;
}

template<typename TValue, typename TSeedSpec, typename TIterator, typename TText, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeeds(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		 TIterator begin, 
		 TIterator end, 
		 String<TText> const &query, 
		 String<TText> const &database, 
		 int gapDistance, 
		 Blat)
{
	SEQAN_CHECKPOINT
	std::multimap<TValue, TIterator> tmpMap; //zum sortieren
	for(TIterator it = begin; it != end; ++it){
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it), it));
		++begin;
	}
	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;
	TIterator2 it_end = tmpMap.end();
	for (TIterator2 tmpIt = tmpMap.begin(); tmpIt != it_end; ++tmpIt){
		if (!addSeed(set, *tmpIt->second, query, database, gapDistance, Blat()))
			addSeed(set, *tmpIt->second, Single());
	}
	tmpMap.clear();
	return true;
}

template<typename TValue, typename TSeedSpec, typename TContainer, typename TAlgoSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeeds(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		 TContainer &source, 
		 int gapDistance, 
		 TAlgoSpec tag)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer,Standard>::Type TIterator;
	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;
    TIterator it_end = end(source);
	std::multimap<TValue, TIterator> tmpMap; //zum sortieren
	for (TIterator it = begin(source); it != it_end; ++it)
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it),it));
	
	TIterator2 it_end2 = tmpMap.end();
	for (TIterator2 it = tmpMap.begin(); it != it_end2; ++it)
		if (!addSeed(set, *it->second, gapDistance, tag))
			addSeed(set, *it->second, Single());
	
	tmpMap.clear();
	return true;
}


template<typename TValue, typename TSeedSpec, typename TContainer, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeeds(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		 TContainer const &source, 
		 Single)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer const,Standard>::Type TIterator;
	TIterator it_end = end(source);
	for (TIterator it = begin(source); it != it_end; ++it)
		addSeed(set, *it, Single());
	
	return true;
}


template<typename TValue, typename TSeedSpec, typename TContainer, typename TText, typename TValue2, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeeds(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		 TContainer const &source, 
		 Score<TValue2, Simple> const &scoreMatrix, 
		 String<TText> const &query, 
		 String<TText> const &database, 
		 int gapDistance, 
		 Chaos)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer const,Standard>::Type TIterator;
    TIterator it1, it2 = end(source);
	std::multimap<TValue, TIterator> tmpMap; //zum sortieren
	for (TIterator it1 = begin(source); it1 != it2; ++it1)
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it1),it1));
	
	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;
	TIterator2 it_end = tmpMap.end();
	for (TIterator2 tmpIt = tmpMap.begin(); tmpIt != it_end; ++tmpIt)
		if (!addSeed(set, *tmpIt->second, scoreMatrix, query, database, gapDistance, Chaos()))
			addSeed(set, *tmpIt->second, Single());

	tmpMap.clear();
	return true;
}

template<typename TValue, typename TSeedSpec, typename TContainer, typename TText, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeeds(SeedSet<TValue, TSeedSpec, 
		 const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		 TContainer const &source, 
		 String<TText> const &query, 
		 String<TText> const &database, 
		 int gapDistance, 
		 Blat)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer const,Standard>::Type TIterator;
	TIterator it2 = end(source);
	std::multimap<TValue, TIterator> tmpMap; //zum sortieren
	for (TIterator it1 = begin(source); it1 !=it2; ++it1)
		tmpMap.insert(std::pair<TValue, TIterator>(leftDim0(*it1),it1));
	
	typedef typename std::multimap<TValue, TIterator>::iterator TIterator2;
	TIterator2 it_end = tmpMap.end();
	for (TIterator2 tmpIt = tmpMap.begin(); tmpIt != it_end; ++tmpIt)
		if (!addSeed(set, *tmpIt->second, query, database, gapDistance, Blat()))
			addSeed(set, *tmpIt->second, Single());
	
	tmpMap.clear();
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//										       Merging												              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qPos, 
		TValue dPos, 
		TValue length_, 
		int gapDistance, 
		Merge)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
    TIterator tmpIt = _findSeedsMerge(set, qPos, dPos, length_, gapDistance);
	bool change = false;
	if (tmpIt != set.fragmentMap.end())
	{
		TSize position = tmpIt->second;
		set.fragmentMap.erase(tmpIt);
		//addSeed(set,set.manager[position],Single());
		
		_mergeTwoSeeds(set.manager[position], qPos, dPos, length_, Merge());

	
		if (_qualityReached(set.manager[position],0,qualityValue(set), TQualityFactor()))
			set.result.insert(position);

		set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[position]),position));
		change = true;
	}
	return change;
}


template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool
addSeed(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		Seed<TValue, TSeedSpec> const &seed, 
		int gapDistance, 
		Merge)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typedef Seed<TValue, TSeedSpec> * pSeed;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	TIterator tmpIt = _findSeedsMerge(set, leftDim0(seed), leftDim1(seed),length(seed), gapDistance);
	bool change = false;
	if (tmpIt != set.fragmentMap.end())
	{
		TSize position = (*tmpIt).second;
		set.fragmentMap.erase(tmpIt);
		//addSeed(set,set.manager[position],Single());
		
		_mergeTwoSeeds(set.manager[position], seed, Merge());
		if (_qualityReached(set.manager[position],0,qualityValue(set), TQualityFactor()))
			set.result.insert(position);
		set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[position]),position));
		change = true;
	}
	return change;
}


template<typename TValue, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, SimpleSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qlPos, 
		TValue dlPos, 
		TValue qrPos, 
		TValue drPos, 
		int gapDistance, 
		Merge)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typedef Seed<TValue, SimpleSeed> * pSeed;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	TIterator tmpIt = _findSeedsMerge(set, qlPos, dlPos, qrPos-qlPos+1, gapDistance);
	bool change = false;
	if (tmpIt != set.fragmentMap.end())
	{
		TSize position = (*tmpIt).second;
		set.fragmentMap.erase(tmpIt);
		//addSeed(set,set.manager[position],Single());
		_mergeTwoSeeds(set.manager[position], qlPos, dlPos, qrPos, drPos, Merge());
		if (_qualityReached(set.manager[position],0, qualityValue(set), TQualityFactor()))
			set.result.insert(position);

		set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[position]),position));
		change = true;
	}
	return change;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Chaining Algorithms                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


								////////////////////////////////////
								//	     Standard Chaining		  //
								////////////////////////////////////

template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qPos, 
		TValue dPos, 
		TValue length, 
		int gapDistance, 
		SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qPos, dPos, length, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		_mergeTwoSeeds(set.manager[id], qPos, dPos, length, SimpleChain());

		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);

		//new endDiagonal
		if(x!= dPos-qPos){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert(std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}

template<typename TValue, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, SimpleSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qlPos, 
		TValue dlPos, 
		TValue qrPos, 
		TValue drPos, 
		int gapDistance, 
		SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, SimpleSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qlPos, dlPos, qrPos-qlPos+1, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		_mergeTwoSeeds(set.manager[id], qlPos, dlPos, qrPos, drPos, SimpleChain());
		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);
		if(x != drPos-qrPos){
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}

template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		Seed<TValue, TSeedSpec> const &seed, 
		int gapDistance, 
		SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, leftDim0(seed), leftDim1(seed), length(seed), gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		_mergeTwoSeeds(set.manager[id], leftDim0(seed), leftDim1(seed), rightDim0(seed), rightDim1(seed), SimpleChain());
		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);
		if(x != endDiagonal(seed))
		{
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}

template<typename TValue, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		Seed<TValue, ChainedSeed> const &seed, 
		int gapDistance, 
		SimpleChain)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, leftDim0(seed), leftDim1(seed), length(seed), gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		_mergeTwoSeeds(set.manager[id], leftDim0(seed), leftDim1(seed), _getFirstDiag(seed).i3, SimpleChain());
		typedef typename std::list<Triple<TValue,TValue,TValue> >::const_iterator TIterator2; 
		
		TIterator2 it2_end = _getDiagSet(seed).end();
		for (TIterator2 it2 = ++_getDiagSet(seed).begin(); it2 != it2_end; ++it2)
			appendDiag(set.manager[id],*it2);

		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);

		if(x != endDiagonal(seed))
		{
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}


template<typename TValue, typename TSeedSpec>
void
_mergeTwoSeeds(Seed<TValue, TSeedSpec> &seed, 
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
_mergeTwoSeeds(Seed<TValue, SimpleSeed> &firstSeed, 
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

template<typename TValue>
void
_mergeTwoSeeds(Seed<TValue, ChainedSeed> &firstSeed, 
			   TValue qPos, 
			   TValue dPos, 
			   TValue length, 
			   SimpleChain)
{
	SEQAN_CHECKPOINT
	(firstSeed).seedSet.push_back(Triple<TValue, TValue, TValue>(qPos, dPos, length));
	TValue diag = dPos-qPos;
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
}


								////////////////////////////////////
								//				Chaos			  //
								////////////////////////////////////

template<typename TValue, typename TSeedSpec, typename TText, typename TValue2, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qPos, 
		TValue dPos, 
		TValue length, 
		Score<TValue2, Simple> const &scoreMatrix, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Chaos)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qPos, dPos, length, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		_mergeTwoSeeds(set.manager[id], qPos, dPos, length, query, database, scoreMatrix, Chaos());

		
		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);

		if(x != dPos-qPos)
		{
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}

template<typename TValue, typename TValue2, typename TText, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		Seed<TValue, ChainedSeed> const &seed, 
		Score<TValue2, Simple> const &scoreMatrix, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Chaos)
{
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, leftDim0(seed), leftDim1(seed),length(seed), gapDistance);
	if (it != set.fragmentMap.end())
	{
		typedef typename std::list<Triple<TValue,TValue,TValue> >::const_iterator TIterator;
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		_mergeTwoSeeds(set.manager[id], leftDim0(seed), leftDim1(seed), _getFirstDiag(seed).i3, query, database, scoreMatrix, Chaos());
 
		TIterator it2_end = _getDiagSet(seed).end(); 
		for (TIterator it2 = ++_getDiagSet(seed).begin(); it2 != it2_end; ++it2)
			appendDiag(set.manager[id],*it2);

		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);

		if(x != endDiagonal(seed))
		{
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}

template<typename TValue, typename TText, typename TValue2>
void
_mergeTwoSeeds(Seed<TValue, ChainedSeed>  &firstSeed, 
			   TValue qPos, 
			   TValue dPos, 
			   TValue length, 
			   String<TText> const &query, 
			   String<TText> const &database, 
			   Score<TValue2, Simple> const &scoreMatrix, 
			   Chaos)
{
	SEQAN_CHECKPOINT
	TValue databaseGap = dPos - rightDim1(firstSeed)-1;
	TValue queryGap = qPos - rightDim0(firstSeed)-1;
	TValue gap = abs(databaseGap-queryGap);

	if (gap == 0)
		setRightDim1(firstSeed,dPos + length-1);
	else 
	{
		TValue tmpLength, tmpScore;
		TValue currenvoid = 0;
		TValue lPositionQuery = rightDim0(firstSeed);
		TValue lPositionDatabase = rightDim1(firstSeed);
		TValue rPositionQuery = qPos;
		TValue rPositionDatabase = dPos;

		TValue gap = (databaseGap < queryGap)? databaseGap : queryGap;
		for (int i = 0; i <gap;++i)
			currenvoid += score(scoreMatrix, sequenceEntryForScore(scoreMatrix, query, --rPositionQuery),
			                    sequenceEntryForScore(scoreMatrix, database, --rPositionDatabase));
		tmpScore = currenvoid;
		tmpLength = 0;
		for (int i = 0; i < gap ; ++i)
		{
			currenvoid += score(scoreMatrix,sequenceEntryForScore(scoreMatrix, query, ++lPositionQuery),
			                    sequenceEntryForScore(scoreMatrix, database, ++lPositionDatabase)) -
                          score(scoreMatrix,sequenceEntryForScore(scoreMatrix, query, rPositionQuery++),
                                sequenceEntryForScore(scoreMatrix, database, rPositionDatabase++));
			if (currenvoid > tmpScore)
			{
				tmpScore = currenvoid;
				tmpLength = i+1;
			}
		}
		setRightDim0(firstSeed,rightDim0(firstSeed)+tmpLength);
		(firstSeed).seedSet.push_back(Triple<TValue,TValue,TValue>(qPos-gap+tmpLength, dPos - gap +tmpLength,gap-tmpLength+length));		
	}
	TValue diag = dPos-qPos;
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
}


								////////////////////////////////////
								//				Blat			  //
								////////////////////////////////////

template<typename TValue, typename TText, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		TValue qPos, 
		TValue dPos, 
		TValue length, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Blat)
{
	SEQAN_CHECKPOINT
	typedef  typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qPos, dPos, length, gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
		int dLength = dPos - rightDim1(set.manager[id]);
		int qLength = qPos - rightDim0(set.manager[id]);
		int dLog = (int) ceil(log((double)dLength));
		int qLog = (int) ceil(log((double)qLength));
		int k;
		int maxValue = (dLog > qLog) ? dLog : qLog;
		if ((maxValue < dLength) && (maxValue < qLength))
			k = maxValue;
		else 
			k = (dLog < qLog)? dLog : qLog;
		_mergeTwoSeeds(set.manager[id], qPos, dPos, length, query, database, k, Blat());
		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);

		if(x != dPos-qPos)
		{
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}

template<typename TValue, typename TText, typename TSpec, typename TQualityFactor, typename TGapCosts>
bool 
addSeed(SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
		Seed<TValue, ChainedSeed> const &seed, 
		String<TText> const &query, 
		String<TText> const &database, 
		int gapDistance, 
		Blat)
{
	SEQAN_CHECKPOINT
	TValue qPos = leftDim0(seed);
	TValue dPos = leftDim1(seed);
	TValue length_ = _getFirstDiag(seed).i3;
	typedef  typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typename std::multimap<TValue,TSize >::iterator it = _findSeedsChain(set, qPos, dPos, length(seed),  gapDistance);
	if (it != set.fragmentMap.end())
	{
		TSize id = it->second;
		TValue x = endDiagonal(set.manager[id]);
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
		_mergeTwoSeeds(set.manager[id], qPos, dPos, length_, query, database, k, Blat());
		typedef typename std::list<Triple<TValue,TValue,TValue> >::const_iterator TIterator; 
		TIterator it2_end = _getDiagSet(seed).end();
		for (TIterator it2 = ++_getDiagSet(seed).begin(); it2 != it2_end; ++it2)
			appendDiag(set.manager[id],*it2);

		if (_qualityReached(set.manager[id],0,qualityValue(set), TQualityFactor()))
			set.result.insert(id);
		
		if(x != endDiagonal(seed))
		{
			set.fragmentMap.erase(it);
			set.fragmentMap.insert( std::pair<TValue, TSize>( endDiagonal(set.manager[id]),id));
		}
		return true;
	} else
		return false;
}

template<typename TValue, typename TText>
void
_mergeTwoSeeds(Seed<TValue, ChainedSeed>  &firstSeed, 
			   TValue qPos,
			   TValue dPos, 
			   TValue length, 
			   String<TText> const &query,
			   String<TText> const &database, 
			   TValue q, 
			   Blat)
{
	SEQAN_CHECKPOINT
	std::list<Triple<TValue, TValue, TValue> > tmp;
	TValue leftDiag = leftDiagonal(firstSeed);
	TValue rightDiag = rightDiagonal(firstSeed);
	TValue diag;
	_gapFill(rightDim0(firstSeed)+1, rightDim1(firstSeed)+1, qPos, dPos, tmp, query, database,q);
	if (tmp.size() > 0)
	{
		Triple<TValue, TValue, TValue> x = *tmp.begin();
		if ((x.i1 == rightDim0(firstSeed)+1)&& (x.i2 == rightDim1(firstSeed)+1))
		{
			tmp.front().i1 = _getDiagSet(firstSeed).back().i1;
			tmp.front().i2 = _getDiagSet(firstSeed).back().i2;
			tmp.front().i3 += _getDiagSet(firstSeed).back().i3;
			_getDiagSet(firstSeed).pop_back();
		}
		x = tmp.back();
		if ((x.i1+x.i3 == qPos)&& (x.i2+x.i3 == dPos))
		{
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
	for (TIterator it = _getDiagSet(firstSeed).begin(); it != it_end ; ++it)
	{
		diag = it->i2 - it->i1;
		if (diag >leftDiag) 
			leftDiag = diag;
		if (diag < rightDiag)
			rightDiag = diag;
	}
	setRightDiagonal(firstSeed, rightDiag);
	setLeftDiagonal(firstSeed, leftDiag);
}


template<typename TValue, typename TText>
void
_extendSeed(Triple<TValue,TValue,TValue> &position, 
			Segment<String<TText> const, InfixSegment> qSeg, 
			Segment<String<TText> const, InfixSegment> dSeg)
{
	SEQAN_CHECKPOINT
	TValue dMax = length(dSeg);
	TValue qMax = length(qSeg);
	TValue qlength = position.i3;
	while ((position.i1+qlength < qMax) && (position.i2+qlength < dMax) && (qSeg[position.i1+qlength] == dSeg[position.i2+qlength]))
		++qlength;

	position.i3 = qlength;
}

template<class TTriple>
struct SortProcess : public std::binary_function<TTriple, TTriple, bool>{

    inline bool 
	operator()(TTriple left, TTriple right) const {//return true if left is logically less then right for given comparison
		return (left.i1+left.i2 < right.i1+right.i2);
	}

};


// Recursive function that fills the gap with smaller alignments
template<typename TValue, typename TText>
void
_gapFill(TValue qlPos, //query sequence left
		 TValue dlPos, //database sequence left
		 TValue qrPos, //query sequence right
		 TValue drPos, //database sequence right
		 std::list<Triple<TValue, TValue, TValue> > &diagList,
		 String<TText> const &query, 
		 String<TText> const &database,
		 TValue q)
{
    SEQAN_CHECKPOINT
    if ( q > 1)
	{
	    bool change =false;
	    Segment<String<TText> const, InfixSegment> qSeg = infix(query ,qlPos, qrPos);
	    Segment<String<TText> const, InfixSegment> dSeg = infix(database ,dlPos, drPos);

	    typedef Index< String<TText>, IndexQGram<SimpleShape > > TQGramIndex;
	  //  TValue qLastLeft = qlPos;
	   // TValue dLastLeft = dlPos;
	    //TValue t;
	    TQGramIndex index_qgram(dSeg);
	    resize(indexShape(index_qgram), q);
	    Finder<TQGramIndex> finder(index_qgram);
	    String<Dna> qgram;
	    std::list<Triple<TValue,TValue,TValue> > listTmp, listTmp2, seeds;
	    Triple<TValue, TValue, TValue> megaTemp(0,0,q);
	    TValue maximum;
	    for (unsigned int i = 0; i <= length(qSeg)-q; ++i)
		{
			qgram= infix(qSeg,i,i+q);
			maximum = 0;
			megaTemp.i1 = i;
			while (find(finder,qgram))
			{
				megaTemp.i3 = q;
				megaTemp.i2 = position(finder); 
				_extendSeed(megaTemp,qSeg,dSeg);
				if (megaTemp.i3 == maximum)
					listTmp2.push_back(megaTemp);
				else
					if (megaTemp.i3 > maximum)
					{
						listTmp2.clear();
						listTmp2.push_back(megaTemp);
						maximum = megaTemp.i3;
					}
			}

			listTmp.splice(listTmp.end(),listTmp2);
			listTmp2.clear();
			clear(finder);
	    }
	    listTmp.sort(SortProcess<Triple<TValue, TValue, TValue > >());

	    TValue qValue = -1;
	    TValue dValue = -1;
		typedef typename std::list<Triple<TValue,TValue,TValue> >::iterator TIterator;
		TIterator it_end = listTmp.end();
	    for (TIterator itInner = listTmp.begin(); itInner != it_end; ++itInner)
		{
		    if ((itInner->i1 > qValue) &&(itInner->i2 >dValue))
			{
			    seeds.push_back(*itInner);
			    qValue = itInner->i1 + itInner->i3 -1;
			    dValue = itInner->i2 + itInner->i3 -1;
		    }
	    }
	    listTmp.clear();

		it_end = seeds.end();
	    for (TIterator itInner = seeds.begin(); itInner != it_end; ++itInner)
		{
		    itInner->i1 += qlPos;
		    itInner->i2 += dlPos;
	    }

	    seeds.push_front(Triple<TValue,TValue,TValue>(qlPos, dlPos,0));
	    seeds.push_back(Triple<TValue,TValue,TValue>(qrPos, drPos,0));
	    TIterator it2 = ++seeds.begin();
	    TIterator it = seeds.begin();
		it_end = seeds.end();
	    while (it2!= it_end)
		{
		    change = true;
			if ((it2->i1 - it->i1 - it->i3 >5)&&(it2->i2 - it->i2 - it->i3 >5))
				_gapFill(it->i1+it->i3,it->i2+it->i3,it2->i1,it2->i2, diagList, query,database,q-1);

			diagList.push_back(Triple<TValue, TValue, TValue>(it2->i1,it2->i2,it2->i3));
		    ++it;
		    ++it2;
	    }
		diagList.pop_back();

		if ((!change) &&(q>2))
			_gapFill(qlPos, dlPos, qlPos, qrPos, diagList, query, database, q-1);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//									Delete to far away														  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename TValue, typename TSeedSpec, typename TQualityFactor, typename TGapCosts>
void
_deleteEverything(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, void> &deletionTarget, 
				  TValue currentPos)
{
	typedef typename Size<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, void> >::Type TSize;
	typedef std::multimap<TValue, TSize > TMap;
	typedef typename std::set<TSize>::iterator TSetIterator;
	TSetIterator set_end = deletionTarget.result.end();
	typename TMap::iterator it_end = deletionTarget.fragmentMap.end();
	TSize pos;
	for (typename TMap::iterator it = deletionTarget.fragmentMap.begin(); it != it_end;)
	{
		pos = it->second;
		if (currentPos - rightDim0(deletionTarget.manager[pos]) > deletionTarget.maxDistance)
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
//				Finding seed for merging/chaining															  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
typename std::multimap<TValue, 	typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, ChainedSeed, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type >::iterator
_findSeedsChain(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
				TValue qPos, 
				TValue dPos, 
				TValue length, 
				int gapDistance)
{

	if(set.last != qPos){
		_deleteEverything(set, qPos);
		set.last = qPos;
	}
	SEQAN_CHECKPOINT
	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	TValue diag = dPos-qPos;
	TIterator itUp=set.fragmentMap.upper_bound(diag+gapDistance);  
	TIterator tmp = set.fragmentMap.end();
	TValue tempDistance1 = minValue<TValue>();
	TValue tmpValue;
	
	for (TIterator it=set.fragmentMap.lower_bound(diag-gapDistance); it != itUp; it++ )
	{
		TSize id = it->second;
		if ((qPos > rightDim0(set.manager[id])) && (dPos > rightDim1(set.manager[id])))
		{
			tmpValue = _calculateScoringValue(set.manager[id],qPos, dPos, length, TGapCosts());
			if (tmpValue > tempDistance1)
			{
				tmp = it;
				tempDistance1 = tmpValue;
			}
		}
	}
	return tmp;
}



template<typename TValue, typename TSeedSpec, typename TSpec, typename TQualityFactor, typename TGapCosts>
typename std::multimap<TValue, typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type >::iterator
_findSeedsMerge(SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> &set, 
				TValue qPos, 
				TValue dPos, 
				TValue length_, 
				int gapDistance)
{	
	SEQAN_CHECKPOINT
	if(set.last != qPos){
		_deleteEverything(set, qPos);
		set.last = qPos;
	}

	typedef typename Size<String<TValue, Block<BLOCK_SIZE<SeedSet<TValue, TSeedSpec, const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, void> >, TSpec> >::Value> > >::Type TSize;
	TValue maxLength = -1;
	TValue tmpLength;
	TValue diag = dPos-qPos;
	typedef typename std::multimap<TValue,TSize >::iterator TIterator;
	TIterator tmpIt = set.fragmentMap.end();
	
	//new seed has higher diagonal number
	TIterator itUp=set.fragmentMap.upper_bound(diag); 
	for (TIterator it=set.fragmentMap.lower_bound(diag-gapDistance); it != itUp; it++ )
	{
		tmpLength = qPos+length_ - leftDim0(set.manager[it->second]) - abs(it->first-diag);
		if ((dPos <= rightDim1(set.manager[it->second])) && (tmpLength > maxLength) &&(tmpLength > length(set.manager[it->second]))){
			tmpIt = it;
			maxLength = tmpLength;
		}
	}
	
	//new seed has lower diagnoal number
	itUp=set.fragmentMap.upper_bound(diag+ gapDistance);  
	for (TIterator it=set.fragmentMap.lower_bound (diag+1); it != itUp; it++ )
	{
		tmpLength = qPos+length_ - 1 - leftDim0(set.manager[it->second])- abs(it->first-diag);
		if ((qPos <= rightDim0(set.manager[it->second])) && (tmpLength > maxLength) &&(tmpLength > length(set.manager[it->second]))){
			tmpIt = it;
			maxLength = tmpLength;
		}
	}

	return tmpIt;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Extension Algorithms                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Function.extendSeeds
..summary: Extension of seeds.
..cat:Seed Handling
..signature:extendSeeds(container, query, database, direction, MatchExtend)
..signature:extendSeeds(begin, end, query, database, direction, MatchExtend)
..signature:extendSeeds(container, scoreDropOff, scoreMatrix, query, database, direction, tag)
..signature:extendSeeds(begin, end, scoreDropOff, scoreMatrix, query, database, direction, tag)
..param.container: The container with the @Class.Seed@ objects to extend.
...type:Concept.Container
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

template<typename TContainer, typename TQuery, typename TDatabase, typename TDirection>
void
extendSeeds(TContainer &seedSet, 
			String<TQuery> &query, 
			String<TDatabase> &database, 
			TDirection direction, 
			MatchExtend)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer,Standard>::Type TIterator;
	TIterator it2 = end(seedSet);
	for (TIterator it1 = begin(seedSet); it1 != it2; ++it1)
		extendSeed(*it1, query, database, direction, MatchExtend());
}

template<typename TIterator, typename TQuery, typename TDatabase, typename TDirection>
void
extendSeeds(TIterator begin, 
			TIterator end, 
			String<TQuery> &query,
			String<TDatabase> &database, 
			TDirection direction, 
			MatchExtend)
{
	SEQAN_CHECKPOINT
	for (TIterator it = begin; it != end; ++it)
		extendSeed(*it, query, database, direction, MatchExtend());		
}

template<typename TContainer, typename TQuery, typename TDatabase, typename TextendSeedSpec, typename TValue, typename TValue2, typename TDirection>
void
extendSeeds(TContainer &seedSet, 
			TValue scoreDropOff, 
			Score<TValue2, Simple> const &scoreMatrix,
			String<TQuery> const &query, 
			String<TDatabase> const &database, 
			TDirection direction, 
			TextendSeedSpec tag)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TContainer,Standard>::Type TIterator;
	TIterator it2 = end(seedSet);
	for (TIterator it1 = begin(seedSet); it1 != it2; ++it1)
		extendSeed(*it1, scoreDropOff, scoreMatrix, query, database, direction, tag);
}


template<typename TIterator, typename TextendSeedSpec, typename TValue, typename TQuery, typename TDatabase, typename TDirection>
void
extendSeeds(TIterator begin, 
			TIterator end, 
			TValue scoreDropOff, 
			Score<TValue, Simple> const &scoreMatrix, 
			String<TQuery> const &query,
			String<TDatabase> const &database, 
			TDirection direction, 
			TextendSeedSpec tag)
{
	SEQAN_CHECKPOINT
	for (TIterator it = begin; it != end; ++it)
		extendSeed(*it, scoreDropOff, scoreMatrix, query, database, direction, tag);		
}


} //namespace Seqan


#endif //#ifndef SEQAN_HEADER_
