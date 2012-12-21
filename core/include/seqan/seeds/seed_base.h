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


#ifndef SEQAN_HEADER_SEED_H
#define SEQAN_HEADER_SEED_H

namespace SEQAN_NAMESPACE_MAIN
{


struct SeedSimple_;
typedef Tag<SeedSimple_> const SimpleSeed;



/**
.Tag.Seed Extension
..cat:Seed Handling
..summary:The algorithms used to extend a seed.
..see:Function.extendSeed
..see:Function.extendSeeds
..see:Function.extendSeedScore
..see:Function.extendSeedsScore
..tag.MatchExtend:
	Extends a seed until a mismatch occurs.
..tag.UngappedXDrop:
	Ungapped extension of a seed until score drops below a Value.
..tag.GappedXDrop:
	Gapped extension of a seed until score drops below a Value. Only @Spec.SimpleSeed@s.
..include:seqan/seeds.h
*/


/**
.Tag.Seed Adding.tag.Merge:
	Merging of Seeds.
..include:seqan/seeds.h
*/
struct ChainMerge_;
typedef Tag<ChainMerge_> const Merge;



struct _extendSeed_Match;
typedef Tag<_extendSeed_Match> const MatchExtend;


struct _extendSeed_UnGappedXDrop;
typedef Tag<_extendSeed_UnGappedXDrop> const UngappedXDrop;

struct ExtendSeedGappedXDrop_;
typedef Tag<ExtendSeedGappedXDrop_> const GappedXDrop;

//template<typename TPosition = int, typename TSpecSeed = SimpleSeed>class Seed;



/**
.Class.Seed:
..summary:Describes a seed.
..cat:Seed Handling
..signature:Seed<TPosition, TSpecSeed>
..param.TPosition:The type number that should be used. Must have negative numbers (e.g. int/long).
..param.TSpec:The seed type used.
..include:seqan/seeds.h
*/

/**
.Spec.SimpleSeed:
..summary:Describes a seed with start and end position and diagonal upper and lower bounds.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, SimpleSeed>
..param.TPosition:The type number that should be used. Must have negative numbers (e.g. int/long).
.Memfunc.SimpleSeed#Seed:
..class:Spec.SimpleSeed
..summary:Constructor
..signature: Seed<TPosition, SimpleSeed> ()
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, length)
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, qEndPos, dEndPos)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.qEndPos: End in query sequence.
..param.dEndPos: End in database sequence.
..param.length: Length of the seed.
..include:seqan/seeds.h
*/
template<typename TPosition = int, typename TSpecSeed = SimpleSeed> 
class Seed{
public:
	TPosition leftDim0;
	TPosition leftDim1;
	TPosition rightDim0;
	TPosition rightDim1;
	TPosition leftDiagonal;
	TPosition rightDiagonal;

	Seed(){
		SEQAN_CHECKPOINT
	}

	Seed(TPosition leftDim0, TPosition leftDim1, TPosition length):leftDim0(leftDim0),leftDim1(leftDim1){
		SEQAN_CHECKPOINT
		rightDim0 = leftDim0 + length-1;
		rightDim1 = leftDim1 + length-1;
		rightDiagonal = leftDiagonal = leftDim1-leftDim0;
	}

	Seed(TPosition leftDim0, TPosition leftDim1, TPosition rightDim0, TPosition rightDim1):leftDim0(leftDim0),leftDim1(leftDim1),rightDim0(rightDim0), rightDim1(rightDim1){
		SEQAN_CHECKPOINT
		leftDiagonal = _max(leftDim1 - leftDim0, rightDim1-rightDim0);
		rightDiagonal = _min(leftDim1 - leftDim0, rightDim1-rightDim0);
	}


	~Seed(){
	}

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											Meta Functions		                                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///.Metafunction.Spec.param.T.type:Class.Seed
///.Metafunction.Spec.class:Class.Seed

template <typename TPosition, typename TSpecSeed>
struct Spec<Seed<TPosition,TSpecSeed> >
{
	typedef TSpecSeed Type;
};

///.Metafunction.Value.param.T.type:Class.Seed
///.Metafunction.Value.class:Class.Seed
template <typename TPosition, typename TSpecSeed>
struct Value<Seed<TPosition,TSpecSeed> >
{
	typedef TPosition Type;
};



template< typename TBorder, typename TSpec >
struct Size< Seed< TBorder, TSpec > >
{
	typedef size_t Type;
};

template< typename TPosition, typename TSpec >
struct Key< Seed< TPosition, TSpec > >
{
	typedef TPosition Type;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Standard Functions                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Function.startDiagonal
..class:Class.Seed
..summary: Returns the diagonal of the start point.
..cat:Seed Handling
..signature:startDiagonal(seed)
..param.seed: The seed whose start diagonal should be returned.
...type:Class.Seed
..returns: The diagonal of the start point.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
startDiagonal(Seed<TPosition, TSpecSeed> const &me)
{
	SEQAN_CHECKPOINT
	return me.leftDim1-me.leftDim0;
}

/**
.Function.endDiagonal
..class:Class.Seed
..summary: Returns the diagonal of the end point.
..cat:Seed Handling
..signature:endDiagonal(seed)
..param.seed: The seed whose end diagonal should be returned.
...type:Class.Seed
..returns: The diagonal of the end point.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
endDiagonal(Seed<TPosition, TSpecSeed> const &me)
{
	SEQAN_CHECKPOINT
	return me.rightDim1-me.rightDim0;
}



/**
.Function.leftPosition
..class:Class.Seed
..summary:The begin position of segment in a seed.
..cat:Seed Handling
..signature:leftPosition(seed, dim)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..returns:Begin position of the $dim$-th segment in $seed$.
*/
template< typename TPosition, typename TSpecSeed, typename TSize> 
inline TPosition 
leftPosition(Seed<TPosition, TSpecSeed>  & me, 
			 TSize dim)
{
	SEQAN_CHECKPOINT
	return (dim)? me.leftDim1 : me.leftDim0;
}

/**
.Function.setLeftPosition
..class:Class.Seed
..summary:Sets begin position of segment in a seed.
..cat:Seed Handling
..signature:setLeftPosition(seed, dim, new_pos)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..param.new_pos:The new begin position of the $dim$-th segment in $seed$.
..see:Function.leftPosition
*/
template< typename TPosition, typename TSpecSeed, typename TSize, typename TPosition2> 
inline void
setLeftPosition(Seed<TPosition, TSpecSeed>  & me, 
				TSize dim,
				TPosition2 new_pos)
{
	SEQAN_CHECKPOINT
	if (dim) me.leftDim1 = new_pos;
	else me.leftDim0 = new_pos;
}

/**
.Function.rightPosition:
..class:Class.Seed
..summary:The end position of segment in a seed.
..cat:Seed Handling
..signature:rightPosition(seed, dim)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..returns:End position of the $dim$-th segment in $seed$.
..see:Function.leftPosition
*/

template< typename TPosition, typename TSpecSeed, typename TSize> 
inline TPosition 
rightPosition(Seed<TPosition, TSpecSeed>  & me, 
			  TSize dim)
{
	SEQAN_CHECKPOINT
	return (dim)? me.rightDim1 : me.rightDim0;
}

/**
.Function.setRightPosition:
..class:Class.Seed
..summary:Sets end position of segment in a seed.
..cat:Seed Handling
..signature:setRightPosition(seed, dim, new_pos)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..param.new_pos:The new end position of the $dim$-th segment in $seed$.
..see:Function.rightPosition
..see:Function.setLeftPosition
*/
template< typename TPosition, typename TSpecSeed, typename TSize, typename TPosition2> 
inline void 
setRightPosition(Seed<TPosition, TSpecSeed>  & me, 
				TSize dim,
				TPosition2 new_pos)
{
	SEQAN_CHECKPOINT
	if (dim) me.rightDim1 = new_pos;
	else me.rightDim0 = new_pos;
}

/**
.Function.dimension:
..class:Class.Seed
..summary:Dimension of a seed.
..cat:Seed Handling
..signature:dimension(seed)
..param.seed:A seed.
...type:Class.Seed
..returns:The number of segments in $seed$.
*/
template< typename TPosition, typename TSpec > inline
typename Size< Seed< TPosition, TSpec > >::Type
dimension( Seed< TPosition, TSpec > & /*me*/ )
{
	return 2;
}

/**
.Function.leftDim0
..class:Class.Seed
..summary: Returns the first position of the seed in the query.
..cat:Seed Handling
..signature:leftDim0(seed)
..param.seed: The seed whose query position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
leftDim0(Seed<TPosition, TSpecSeed> const &seed)
{
	SEQAN_CHECKPOINT
	return seed.leftDim0;
}

/**
.Function.rightDim0
..class:Class.Seed
..summary: Returns the last position of the seed in the query.
..cat:Seed Handling
..signature:rightDim0(seed)
..param.seed: The seed whose last in the query position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
rightDim0(Seed<TPosition,TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDim0;
}

/**
.Function.leftDim1
..class:Class.Seed
..summary: Returns the first position of the seed in the database.
..cat:Seed Handling
..signature:leftDim1(seed)
..param.seed: The seed whose database position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
leftDim1(Seed<TPosition, TSpecSeed> const &seed)
{
	SEQAN_CHECKPOINT
	return seed.leftDim1;
}

/**
.Function.rightDim1
..class:Class.Seed
..summary: Returns the last position of the seed in the database.
..cat:Seed Handling
..signature:rightDim1(seed)
..param.seed: The seed whose last in the database position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
rightDim1(Seed<TPosition,TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDim1;
}

/**
.Function.leftDiagonal
..class:Class.Seed
..summary: Returns the most left diagonal of the seed (maximum diagonal value).
..cat:Seed Handling
..signature:leftDiagonal(seed)
..param.seed: The seed whose database position should be returned.
...type:Class.Seed
..returns: The most left diagonal.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
leftDiagonal(Seed<TPosition, TSpecSeed> const &seed)
{
	SEQAN_CHECKPOINT
	return seed.leftDiagonal;
}

/**
.Function.rightDiagonal
..class:Class.Seed
..summary: Returns the most right diagonal of the seed (minimum diagonal value).
..cat:Seed Handling
..signature:rightDiagonal(seed)
..param.seed: The seed whose database position should be returned.
...type:Class.Seed
..returns: The most right diagonal.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
rightDiagonal(Seed<TPosition,TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDiagonal;
}

///.Function.length.param.object.type:Class.Seed
///.Function.length.class:Class.Seed
template<typename TPosition, typename TSpecSeed>
inline TPosition 
length(Seed<TPosition, TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDim0-seed.leftDim0+1;
}


/**
.Function.setLeftDim0:
..class:Class.Seed
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim0(seed, start)
..param.seed: The seed whose start position should be updated.
...type:Class.Seed
..param.start: The query position where the seed should start.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setLeftDim0(Seed<TPosition, TSpecSeed> &me, 
			  TPosition start)
{
	SEQAN_CHECKPOINT
	me.leftDim0 = start;
}

/**
.Function.setRightDim0:
..class:Class.Seed
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim0(seed, end)
..param.seed: The seed whose end position should be updated.
...type:Class.Seed
..param.end: The query position where the seed should end.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setRightDim0(Seed<TPosition,TSpecSeed> & me, 
			TPosition end)
{
	SEQAN_CHECKPOINT
	me.rightDim0 = end;
}

/**
.Function.setLeftDim1:
..class:Class.Seed
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim1(seed, start)
..param.seed: The seed whose start position should be updated.
...type:Class.Seed
..param.start: The database position where the seed should start.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setLeftDim1(Seed<TPosition, TSpecSeed> &me, 
				 TPosition start)
{
	SEQAN_CHECKPOINT
	me.leftDim1 = start;
}

/**
.Function.setRightDim1:
..class:Class.Seed
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim1(seed, end)
..param.seed: The seed whose end position should be updated.
...type:Class.Seed
..param.end: The database position where the seed should end.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setRightDim1(Seed<TPosition,TSpecSeed> & me, 
			   TPosition end)
{
	SEQAN_CHECKPOINT
	me.rightDim1 = end;
}

/**
.Function.setLeftDiagonal:
..class:Class.Seed
..summary: Sets a new value for the most left diagonal.
..cat:Seed Handling
..signature:setLeftDiagonal(seed, diag)
..param.seed: The seed whose left diagonal value should be updated.
...type:Class.Seed
..param.diag: The new value for the most left diagonal.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setLeftDiagonal(Seed<TPosition, TSpecSeed> &me,
				TPosition diag)
{
	SEQAN_CHECKPOINT
	me.leftDiagonal = diag;
}

/**
.Function.setRightDiagonal:
..class:Class.Seed
..summary: Sets a new value for the most right diagonal.
..cat:Seed Handling
..signature:setRightDiagonal(seed, diag)
..param.seed: The seed whose right diagonal value should be updated.
...type:Class.Seed
..param.diag: The new value for the most right diagonal.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setRightDiagonal(Seed<TPosition,TSpecSeed> & seed, 
				 TPosition diag)
{
	SEQAN_CHECKPOINT
	seed.rightDiagonal = diag;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      Merge Alogrithms                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename TPosition, typename TSpecSeed>
void 
_mergeTwoSeeds(Seed<TPosition, TSpecSeed> &firstSeed, 
			   Seed<TPosition, TSpecSeed> const &secondSeed, 
			   Merge)
{
	SEQAN_CHECKPOINT
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		setRightDim0(firstSeed,rightDim0(secondSeed));
		setRightDim1(firstSeed,rightDim1(secondSeed));
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}


template<typename TPosition, typename TSpecSeed>
void 
_mergeTwoSeeds(Seed<TPosition, TSpecSeed> &firstSeed, 
			   TPosition q,
			   TPosition d, 
			   TPosition l, 
			   Merge)
{
	SEQAN_CHECKPOINT
	setRightDim0(firstSeed,q+l-1);
	setRightDim1(firstSeed,d+l-1);
	if (leftDiagonal(firstSeed) < d-q)
		setLeftDiagonal(firstSeed, d-q);
	if (rightDiagonal(firstSeed) > d-q)
		setRightDiagonal(firstSeed, d-q);
}


template<typename TPosition>
void
_mergeTwoSeeds(Seed<TPosition, SimpleSeed> &firstSeed, 
			   TPosition qlPos, 
			   TPosition dlPos, 
			   TPosition qrPos, 
			   TPosition drPos, 
			   Merge)
{
	SEQAN_CHECKPOINT
	if (qrPos > rightDim0(firstSeed)){
		typename std::list<Triple <TPosition, TPosition, TPosition> >::iterator begin1, end2, it;
	
		setRightDim0(firstSeed,qrPos);
		setRightDim1(firstSeed,drPos);
	
		TPosition diag = dlPos -qlPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
		diag = drPos - qrPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
	}
}


template<typename TPosition, typename TSpecSeed, typename TPosition2, typename TPosition3, typename TGapCost>
void 
_mergeTwoSeedsScore(Seed<TPosition, TSpecSeed> &firstSeed, 
					TPosition3 &score1, 
					Seed<TPosition, TSpecSeed> const &secondSeed, 
					TPosition3 score2, 
					Score<TPosition2,Simple> const &scoreMatrix, 
					TGapCost &, 
					Merge)
{
	SEQAN_CHECKPOINT
	score1 += score2;
	score1 += abs(endDiagonal(firstSeed)-startDiagonal(secondSeed))*scoreGap(scoreMatrix);
	score1 -= (_max(abs(rightDim0(firstSeed)-leftDim0(secondSeed)),abs(rightDim1(firstSeed)-leftDim1(secondSeed)))+1)*scoreMatch(scoreMatrix);
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		setRightDim0(firstSeed,rightDim0(secondSeed));
		setRightDim1(firstSeed,rightDim1(secondSeed));
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}


template<typename TPosition, typename TPosition2, typename TSpecSeed, typename TPosition3, typename TGapCost>
void 
_mergeTwoSeedsScore(Seed<TPosition, TSpecSeed> &firstSeed, 
					TPosition3 &score1, 
					TPosition q, 
					TPosition d, 
					TPosition l, 
					TPosition3 score2, 
					Score<TPosition2, Simple> const &scoreMatrix, 
					TGapCost &, 
					Merge)
{
	SEQAN_CHECKPOINT
	score1 += score2;
	score1 += abs((TPosition2)(endDiagonal(firstSeed) - d + q))*scoreGap(scoreMatrix);
	score1 -= (::std::max<TPosition>(abs((TPosition2)(rightDim0(firstSeed)- q)),abs((TPosition2)(rightDim1(firstSeed)-d)))+1)*scoreMatch(scoreMatrix);
	setRightDim0(firstSeed,q+l-1);
	setRightDim1(firstSeed,d+l-1);
	if (leftDiagonal(firstSeed) < d-q)
		setLeftDiagonal(firstSeed, d-q);
	if (rightDiagonal(firstSeed) > d-q)
		setRightDiagonal(firstSeed, d-q);
}


template<typename TPosition, typename TPosition2, typename TPosition3, typename TGapCost>
void           
_mergeTwoSeedsScore(Seed<TPosition, SimpleSeed> &firstSeed, 
					TPosition3 &score1, 
					TPosition qlPos, 
					TPosition dlPos, 
					TPosition qrPos, 
					TPosition drPos, 
					TPosition3 score2, 
					Score<TPosition2,Simple> const &scoreMatrix, 
					TGapCost &, 
					Merge)
{
	SEQAN_CHECKPOINT
	if (qrPos > rightDim0(firstSeed)){
		
		score1 += score2;
		score1 += abs(endDiagonal(firstSeed) - dlPos + qlPos)*scoreGap(scoreMatrix);
		score1 -= (_max(abs(rightDim0(firstSeed)- qlPos),abs(rightDim1(firstSeed)-dlPos))+1)*scoreMatch(scoreMatrix);

		setRightDim0(firstSeed,qrPos);
		setRightDim1(firstSeed,drPos);
		TPosition diag = dlPos -qlPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
		diag = drPos - qrPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Extension Algorithms                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
.Function.extendSeed
..class:Class.Seed
..summary:Extends a seed.
..cat:Seed Handling
..signature:extendSeed(seed, query, database, direction, tag)
..signature:extendSeed(seed, scoreDropOff, scoreMatrix, query, database, direction, tag)
..param.seed: The seed to extend.
...type:Class.Seed
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


template<typename TPosition, typename TSpecSeed, typename TQuery, typename TDatabase, typename TDirection>
void 
extendSeed(Seed<TPosition, TSpecSeed> &seed, 
		   TQuery const &query, 
		   TDatabase const &database, 
		   TDirection direction, 
		   MatchExtend)
{
	SEQAN_CHECKPOINT
	//left extension
	if (direction != 1){
		TPosition queryPos =leftDim0(seed) ;
		TPosition dataPos = leftDim1(seed);
		while ((queryPos-1>=0) && (dataPos-1>=0) && (query[queryPos-1] == database[dataPos-1])){
			--queryPos;
			--dataPos;
		}
		setLeftDim0(seed,queryPos);
		setLeftDim1(seed,dataPos);

	}

	//right extension
	if (direction != 0){
		int queryLength = length(query);
		int databaseLength = length(database);
		TPosition queryPos =rightDim0(seed) ;
		TPosition dataPos = rightDim1(seed);
		while ((queryPos+1 < queryLength) && (dataPos+1 < databaseLength) && (query[queryPos+1] == database[dataPos+1])){
			++queryPos;
			++dataPos;
		}

		setRightDim0(seed,queryPos);
		setRightDim1(seed,dataPos);
	}
}


template<typename TPosition, typename TSpecSeed, typename TQuery, typename TDatabase, typename TScore, typename TDirection>
void 
extendSeed(Seed<TPosition,TSpecSeed> &seed, 
		   TScore scoreDropOff, 
		   Score<TScore, Simple> const &scoreMatrix,
		   TQuery const &query,
		   TDatabase const &database,
		   TDirection direction, 
		   UngappedXDrop)
{
	SEQAN_CHECKPOINT
	scoreDropOff *=-1;
	TScore tmpScore = 0;

	//left extension
	if (direction != 1){
		TPosition xPos = leftDim0(seed)-1;
		TPosition yPos = leftDim1(seed)-1;
		TPosition last = 0;

		while ((tmpScore > scoreDropOff) && (xPos >= 0) && (yPos>=0)){
			if (query[xPos] == database[yPos]){
				last = 0;
				tmpScore += score(scoreMatrix, sequenceEntryForScore(scoreMatrix, query, xPos),
				                  sequenceEntryForScore(scoreMatrix, database, yPos));
				if (tmpScore > 0)
					tmpScore = 0;
			} else{
				tmpScore += score(scoreMatrix, sequenceEntryForScore(scoreMatrix, query, xPos),
                                  sequenceEntryForScore(scoreMatrix, database, yPos));
				++last;
			}
		--xPos;
		--yPos;
		}

		setLeftDim0(seed,xPos+last+1);
		setLeftDim1(seed,yPos+last+1);
	}

	//right extension
	if (direction != 0){
		TPosition xLength = length(query);
		TPosition yLength = length(database);
		TPosition xPos = rightDim0(seed)+1;
		TPosition yPos = rightDim1(seed)+1;

		TPosition last = 0;
		tmpScore= 0;
		while ((tmpScore > scoreDropOff) && (xPos < xLength) && (yPos < yLength)){
			if (query[xPos] == database[yPos]){
				last = 0;
				tmpScore += score(scoreMatrix, sequenceEntryForScore(scoreMatrix, query, xPos),
				                  sequenceEntryForScore(scoreMatrix, database, yPos));
				if (tmpScore > 0)
					tmpScore = 0;
			}else{
				tmpScore += score(scoreMatrix, sequenceEntryForScore(scoreMatrix, query, xPos),
				                  sequenceEntryForScore(scoreMatrix, database, yPos));
				++last;
			}
			++xPos;
			++yPos;
		}
		setRightDim0(seed,xPos-last-1);
		setRightDim1(seed,yPos-last-1);
	}
}

// Sets the begin and end position as well as the left and right diagonal of a SimpleSeed after seed extension
template<typename TPosition, typename TBound, typename TExtension, typename TSize>
void
_setExtendedSeedDimensions(Seed<TPosition, SimpleSeed> & seed,
                           TBound lowerBound,
                           TBound upperBound,
                           TExtension extLengthQuery,
                           TExtension extLengthDatabase,
                           TSize direction) {
    if(direction == 0) {
	    // set left and right diagonals
        if (leftDiagonal(seed) < startDiagonal(seed)+upperBound)
	        setLeftDiagonal(seed, startDiagonal(seed)+upperBound);
	    if (rightDiagonal(seed) > startDiagonal(seed)-lowerBound)
	       	setRightDiagonal(seed, startDiagonal(seed)-lowerBound);

        // set new start position of seed
        setLeftDim0(seed, leftDim0(seed)-extLengthQuery);
        setLeftDim1(seed, leftDim1(seed)-extLengthDatabase);
    } else {
	    // set left and right diagonals
	    if (rightDiagonal(seed) > endDiagonal(seed)-upperBound)
	        setRightDiagonal(seed, endDiagonal(seed)-upperBound);
        if (leftDiagonal(seed) < endDiagonal(seed)+lowerBound)
            setLeftDiagonal(seed, endDiagonal(seed)+lowerBound);

        // set new end position of seed
        setRightDim0(seed, rightDim0(seed)+extLengthQuery);
        setRightDim1(seed, rightDim1(seed)+extLengthDatabase);
    }
}

template<typename TSeedSpec, typename TPosition, typename DatabaseInfix, typename QueryInfix, typename TScore, typename TSize>
TScore
_extendSeedOneDirection(Seed<TPosition, TSeedSpec/*SimpleSeed*/> & seed,
                       DatabaseInfix const & dataSeg,
                       QueryInfix const & querySeg, 
                       TScore scoreDropOff,
                       Score<TScore, Simple> const & scoreMatrix,
                       TSize direction) {
    TScore gapCost = scoreGap(scoreMatrix);
    TScore infimum = minValue<TScore>()+1-gapCost;

    TPosition upperBound = 0;
    TPosition lowerBound = 0;

    TPosition xLength = length(querySeg);
    TPosition yLength = length(dataSeg);

    // set variables for calculating the sequence positions
    int factor, xSummand, ySummand;
    if(direction == 0) {
        factor = -1;
        xSummand = xLength;
        ySummand = yLength;
    } else {
        factor = 1;
        xSummand = -1;
        ySummand = -1;
    }

    // antidiagonals of DP matrix
    String<TScore> antiDiagonal1;
    String<TScore> antiDiagonal2;
    String<TScore> antiDiagonal3;

    resize(antiDiagonal1, 1, 0);
    resize(antiDiagonal2, 2, infimum);
    resize(antiDiagonal3, _min(3, xLength+1), infimum);

    String<TScore> *antiDiag1 = &antiDiagonal1;	//smallest diagonal
	String<TScore> *antiDiag2 = &antiDiagonal2;
	String<TScore> *antiDiag3 = &antiDiagonal3;	//current diagonal
	String<TScore> *tmpDiag;

	//Matrix initialization
	if (gapCost >= (-1)*scoreDropOff) {
		(*antiDiag2)[0] = gapCost;
		(*antiDiag2)[1] = gapCost;
	}
	if (2*gapCost >= (-1)*scoreDropOff) {
		(*antiDiag3)[0] = 2*gapCost;
		(*antiDiag3)[2] = 2*gapCost;
	}
		
	TPosition b = 1; // lower bound for i
	TPosition u = 0; // upper bound for i
	TPosition k = 1; // current antidiagonal
	TScore tmp;   // for calculating the maximum for a single matrix entry
	TScore tmpMax1 = 0; // maximum score without the current diagonal
	TScore tmpMax2 = 0; // maximum score including the current diagonal

	//Extension as proposed by Zhang et al
	while(true) {
		++k;
		for (TPosition i = b; i <= u+1; ++i) {
            // calculate matrix entry
			tmp = infimum;
			tmp = _max((*antiDiag2)[i-1], (*antiDiag2)[i]) + gapCost;
			tmp = _max(tmp, (*antiDiag1)[i-1] +
			           score(scoreMatrix, sequenceEntryForScore(scoreMatrix, querySeg, xSummand + factor*i),
			                 sequenceEntryForScore(scoreMatrix, dataSeg, ySummand + factor*(k-i))));

			tmpMax2 = _max(tmpMax2, tmp);
			if (tmp < tmpMax1-scoreDropOff)
				(*antiDiag3)[i] = infimum;
			else
				(*antiDiag3)[i] = tmp;
		}
        
        // narrow the relevant matrix region
		while ((b < (TPosition)length(*antiDiag3)-1) && ((*antiDiag3)[b] == infimum) && (*antiDiag2)[b-1] == infimum) {
			++b;
		}
		++u;
		while ((u >= 0) && ((*antiDiag3)[u] == infimum) && (*antiDiag2)[u] == infimum) {
			--u;
        }
			
		//borders for lower triangle of edit matrix
		b = _max(b, k+1-yLength);
		u = _min(u, xLength-1);

        if (b > u+1) break;

	    //calculate upper/lower bound for diagonals
        if (2*((k+1)/2 - b) - (k%2) > lowerBound) {
			lowerBound = 2*((k+1)/2 - b) - (k%2);
        }
        if (2*(u - k/2) - (k%2) > upperBound) {
            upperBound = 2*(u - k/2) - (k%2);
        }

        // swap diagonals
		tmpDiag = antiDiag1;
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag3 = tmpDiag;

        // extend last diagonal diagonal to be the new longest and initialize
		int d = 0;
		while ((d < 3) && ((TPosition)length(*antiDiag3) <= xLength)) {
			appendValue(*antiDiag3, 0);
			++d;
		}
		for (unsigned int eu = 0; eu < length(*antiDiag3); ++eu)
			(*antiDiag3)[eu] = infimum;

		if ((*antiDiag2)[0]+ gapCost >= tmpMax1-scoreDropOff)
			(*antiDiag3)[0] = (*antiDiag2)[0] + gapCost;
		if ((*antiDiag2)[length(*antiDiag2)-1] + gapCost >= tmpMax1-scoreDropOff)
			(*antiDiag3)[length(*antiDiag3)-1] = (*antiDiag2)[length(*antiDiag2)-1] + gapCost;

		tmpMax1 = tmpMax2;
	}

  //  // print anti diagonals
  //  for(int ii = length(*antiDiag1)-1; ii >= 0 ; --ii) {
  //      for(int jj = ii; jj > 0 ; --jj) std::cout << "  ";
  //      std::cout << " ";
  //      if ((*antiDiag1)[ii] == infimum ) std::cout << "i" << " ";
  //      else std::cout << (*antiDiag1)[ii] << " ";
  //      if (length(*antiDiag2) <= ii+1 ) std::cout << " ";
  //      else if ((*antiDiag2)[ii+1] == infimum ) std::cout << "i" << " ";
  //      else std::cout << (*antiDiag2)[ii+1] << " ";
		//if (length(*antiDiag3) <= ii+2 ) std::cout << std::endl;
  //      else if ((*antiDiag3)[ii+2] == infimum ) std::cout << "i" << std::endl;
  //      else std::cout << (*antiDiag3)[ii+2] << std::endl;
  //  }

	//Find seed start/end
	TPosition extLengthQuery = 0; // length of extension in query
    TPosition extLengthDatabase = 0; // length of extension in database
	TScore tmpMax = infimum;
	if ((k >= xLength+yLength) && ((*antiDiag2)[u+1] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of both sequences
		extLengthQuery = xLength;
        extLengthDatabase = yLength;
		tmpMax = (*antiDiag2)[u+1];
    } else if ((b >= xLength) && b < (TPosition)length(*antiDiag2) && ((*antiDiag2)[b] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of query
        tmpMax = (*antiDiag2)[b];
        extLengthQuery = xLength;
        extLengthDatabase = k-(xLength+1);
    } else if ((k-u-1 >= yLength) && u >= 0 && ((*antiDiag2)[u] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of database
        tmpMax = (*antiDiag2)[u];
        extLengthQuery = u;
        extLengthDatabase = yLength;
    } else {
        for (unsigned int eu = 0; eu < length(*antiDiag1); ++eu) {
            if ((*antiDiag1)[eu] > tmpMax) {
                // extension ends with mismatch
		        tmpMax = (*antiDiag1)[eu];
		        extLengthQuery = _min(xLength, (TPosition)eu);
                extLengthDatabase = _min(yLength, (TPosition)(k - (eu+2)));
		    }
        }
	}

    if(tmpMax != infimum) {
        _setExtendedSeedDimensions(seed, lowerBound, upperBound, extLengthQuery, extLengthDatabase, direction);
        return tmpMax;
    } else {
        return 0;
    }
}


template<typename TSeedSpec, typename TPosition, typename TQuery, typename TDatabase, typename TScore, typename TSize>
void 
extendSeed(Seed<TPosition, TSeedSpec> &seed, 
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
		
        _extendSeedOneDirection(seed, dataSeg, querySeg, scoreDropOff, scoreMatrix, 0 /*left*/);
	}

	//right extension
	if ((direction != 0) && (rightDim0(seed)+1 < (TPosition)endPosition(query)) && (rightDim1(seed)+1 < (TPosition)endPosition(database)))
    {
        typename Suffix<TDatabase const>::Type dataSeg = suffix(database, rightDim1(seed)+1);
        typename Suffix<TQuery const>::Type querySeg = suffix(query, rightDim0(seed)+1);

        _extendSeedOneDirection(seed, dataSeg, querySeg, scoreDropOff, scoreMatrix, 1 /*right*/);
	}
}

} //end of Seqan namespace

#endif //#ifndef SEQAN_HEADER_

