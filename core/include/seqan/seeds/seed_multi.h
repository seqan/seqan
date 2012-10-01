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

#ifndef SEQAN_HEADER_MultiSeed_H
#define SEQAN_HEADER_MultiSeed_H

namespace SEQAN_NAMESPACE_MAIN
{

struct SeedMulti_;
typedef Tag<SeedMulti_> const ChainedSeed;

/**
..Spec.ChainedSeed
..summary:Describes a seed with start and end position2 and diagonal upper and lower bounds. Additionaly diagonal segments
between start and end position2 are stored.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, ChainedSeed>
..param.TPosition:The type of number that should be used. Must have negative numbers (e.g. int/long).
.Memfunc.ChainedSeed#Seed:
..class:Spec.ChainedSeed
..summary:Constructor
..signature: Seed<TPosition, SimpleSeed> ()
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, length)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.length: Length of the seed.
*/

template<typename TPosition> 
class Seed<TPosition, ChainedSeed> {
	
public:
	std::list<Triple<TPosition, TPosition, TPosition> > seedSet;
	TPosition leftDiagonal;
	TPosition rightDiagonal;
 
	Seed(){
		SEQAN_CHECKPOINT
	}
	
	Seed(TPosition leftDim0,
		 TPosition leftDim1,
		 TPosition length)
	{
		SEQAN_CHECKPOINT
		seedSet.push_back(Triple<TPosition, TPosition, TPosition> (leftDim0, leftDim1, length));
		rightDiagonal = leftDiagonal = leftDim1 - leftDim0;
	}
	
	~Seed(){
	}

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Standard Functions                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TPosition>
inline TPosition 
startDiagonal(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.front().i2-seed.seedSet.front().i1;
}

template<typename TPosition>
inline TPosition 
endDiagonal(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.back().i2-seed.seedSet.back().i1;
}


template<typename TPosition>
inline TPosition 
leftDim0(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.front().i1;
}


template<typename TPosition>
inline TPosition 
rightDim0(Seed<TPosition, ChainedSeed> const & seed)
{
	return seed.seedSet.back().i1+seed.seedSet.back().i3-1;
}


template<typename TPosition>
inline TPosition 
leftDim1(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.front().i2;
}

template<typename TPosition>
inline TPosition 
rightDim1(Seed<TPosition, ChainedSeed> const & seed)
{
	return seed.seedSet.back().i2+seed.seedSet.back().i3-1;
}


template<typename TPosition>
inline TPosition 
length(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.seedSet.back().i1 + seed.seedSet.back().i3 - seed.seedSet.front().i1;	
}



template<typename TPosition>
inline void 
setLeftDim0(Seed<TPosition, ChainedSeed> &seed, 
			  TPosition start)
{
	SEQAN_CHECKPOINT
	TPosition lengthDiff = seed.seedSet.front().i1 - start;
	seed.seedSet.front().i3 += lengthDiff;
	seed.seedSet.front().i2 -= lengthDiff;
	seed.seedSet.front().i1 = start;
}


template<typename TPosition>
inline void 
setRightDim0(Seed<TPosition,ChainedSeed> & seed, 
			TPosition end)
{
	SEQAN_CHECKPOINT
	seed.seedSet.back().i3 = end - seed.seedSet.back().i1+1;
}


template<typename TPosition>
inline void 
setLeftDim1(Seed<TPosition, ChainedSeed> &seed, 
				 TPosition start)
{
	SEQAN_CHECKPOINT
	TPosition lengthDiff = seed.seedSet.front().i2 - start;
	seed.seedSet.front().i3 += lengthDiff;
	seed.seedSet.front().i1 -= lengthDiff;
	seed.seedSet.front().i2 = start;
}


template<typename TPosition>
inline void 
setRightDim1(Seed<TPosition,ChainedSeed> & seed, 
			   TPosition end)
{
	SEQAN_CHECKPOINT
	seed.seedSet.back().i3 = end - seed.seedSet.back().i2+1;
}

/*
.Function._getDiagSet:
..summary: Returns the set of matching diagonals.
..cat:Seed Handling
..signature:setRightDim1(seed)
..param.seed: The seed whose end position2 should be updated.
...type:Spec.ChainedSeed
..returns: A reference to the list of seed diagonals.
..include:seqan/seeds.h
*/
template<typename TPosition>
inline const std::list<Triple<TPosition, TPosition, TPosition> >&
_getDiagSet(Seed<TPosition,ChainedSeed> const & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet;
}

template<typename TPosition>
inline std::list<Triple<TPosition, TPosition, TPosition> >&
_getDiagSet(Seed<TPosition,ChainedSeed> & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet;
}

/**
.Function.appendDiag
..class:Class.Seed
..summary: Adds diagonal to the seed.
..cat:Seed Handling
..signature:appendDiag(seed, diagonal)
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeed
..param.diag: The diagonal to add.
...type:Class.Triple
...remarks: A diagonal consists of three values: 1: start in 1. sequence, 2: start in 2. sequence, 3: length of match
..include:seqan/seeds.h
*/
template<typename TPosition>
void
appendDiag(Seed<TPosition,ChainedSeed> & seed, 
		   Triple<TPosition, TPosition, TPosition> diag)
{
	SEQAN_CHECKPOINT
	seed.seedSet.push_back(diag);
}

template<typename TPosition>
Triple<TPosition, TPosition, TPosition>&
_getFirstDiag(Seed<TPosition,ChainedSeed> & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.front();
}

template<typename TPosition>
const Triple<TPosition, TPosition, TPosition>&
_getFirstDiag(Seed<TPosition,ChainedSeed> const & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.front();
}

template<typename TPosition>
Triple<TPosition, TPosition, TPosition>&
_getLastDiag(Seed<TPosition,ChainedSeed> & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.back();
}

template<typename TPosition>
const Triple<TPosition, TPosition, TPosition>&
_getLastDiag(Seed<TPosition,ChainedSeed> const & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.back();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      Merge Alogrithms                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TPosition>
void
_mergeTwoSeeds(Seed<TPosition, ChainedSeed> &firstSeed,
			   TPosition qPos,
			   TPosition dPos,
			   TPosition length,
			   Merge)
{
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::iterator TIterator;
    TIterator begin1, end2, it;
	TPosition diag = dPos -qPos;
	//new seed would be longer?
	if (qPos+length-1 > rightDim0(firstSeed)){
		while ((_getLastDiag(firstSeed)).i1 > qPos)// || ((_getLastDiag(firstSeed)).i2 > dPos))
		{
			(firstSeed).seedSet.pop_back();
		}
		if ((rightDim0(firstSeed) < qPos) && (rightDim1(firstSeed) < dPos)) {
			appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
		} else {
			if (diag == endDiagonal(firstSeed)){
				setRightDim1(firstSeed,dPos+length-1);
			} else {
				TPosition tmp = diag - endDiagonal(firstSeed);
				if (tmp < 0){
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- dPos+1);
					if (tmp2 > 0){
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
					}
				} else {
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - qPos+1);
					if (tmp2 > 0){
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
				}
			}
		}
	}

	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
	}
}


template<typename TPosition>
void
_mergeTwoSeeds(Seed<TPosition, ChainedSeed> &firstSeed, 
			   Seed<TPosition, ChainedSeed> const &secondSeed, 
			   Merge)
{
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::const_iterator TIterator;
        TIterator begin1, end2, it;
	begin1 = _getDiagSet(secondSeed).begin();
	end2 = _getDiagSet(secondSeed).end();
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		while (_getLastDiag(firstSeed).i1 > _getFirstDiag(secondSeed).i1) //|| ((_getLastDiag(firstSeed).i2 > _getFirstDiag(secondSeed).i2)))
		{
			(firstSeed).seedSet.pop_back();
		}

		if ((rightDim0(firstSeed) < leftDim0(secondSeed)) && (rightDim1(firstSeed) < leftDim1(secondSeed))) {
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		} else {
		if (startDiagonal(secondSeed) == endDiagonal(firstSeed)){
			setRightDim1(firstSeed,(*begin1).i2+(*begin1).i3-1);
			++begin1;
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		}
		else {
			
			TPosition tmp = startDiagonal(secondSeed) - endDiagonal(firstSeed);
			if (tmp < 0){
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- leftDim1(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			} else {
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - leftDim0(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			}
		}	
		}
		
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}

template<typename TPosition, typename TPosition2, typename TPosition3, typename TGapCost>
void
_mergeTwoSeedsScore(Seed<TPosition, ChainedSeed> &firstSeed,
					TPosition3 &score1,
					TPosition qPos,
					TPosition dPos,
					TPosition length,
					TPosition3 score2,
					Score<TPosition2,Simple> const &scoreMatrix,
					TGapCost &,
				    Merge)
{
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::iterator TIterator;
    TIterator begin1, end2, it;
	TPosition diag = dPos - qPos;
	score1 += score2;

	
	if (qPos+length-1 > rightDim0(firstSeed)){
		while (_getLastDiag(firstSeed).i1 > qPos)// || (_getLastDiag(firstSeed).i2 > dPos))
		{
			score1 -= firstSeed.seedSet.back().i3*scoreMatch(scoreMatrix);
			TPosition x1 = firstSeed.seedSet.back().i1;
			TPosition x2 = firstSeed.seedSet.back().i2;
			firstSeed.seedSet.pop_back();
			score1 -=(abs(rightDim0(firstSeed)-x1)+abs(rightDim1(firstSeed)-x2))*scoreGap(scoreMatrix); 
		}
		score1 += abs(endDiagonal(firstSeed) - dPos + qPos)*scoreGap(scoreMatrix);
		score1 -=(_max(abs(rightDim0(firstSeed)- qPos),abs(rightDim1(firstSeed)-dPos))+1)*scoreMatch(scoreMatrix);

		if ((rightDim0(firstSeed) < qPos) && (rightDim1(firstSeed) < dPos)) {
			appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
		} else {
			if (diag == endDiagonal(firstSeed))
			{
				setRightDim1(firstSeed,dPos+length-1);
			} 
			else 
			{
				TPosition tmp = diag - endDiagonal(firstSeed);
				if (tmp < 0)
				{
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- dPos+1);
					if (tmp2 > 0)
					{
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
					}
				} 
				else 
				{
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - qPos+1);
					if (tmp2 > 0)
					{
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
					}
				}
			}
		}
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
	}
}


template<typename TPosition, typename TPosition2, typename TPosition3, typename TGapCost>
void
_mergeTwoSeedsScore(Seed<TPosition, ChainedSeed> &firstSeed,
					TPosition3 &score1,
					Seed<TPosition, ChainedSeed> const &secondSeed,
					TPosition3 score2,
					Score<TPosition2,Simple> const &scoreMatrix,
					TGapCost &,
					Merge)
{
	SEQAN_CHECKPOINT
	score1 += score2;
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::const_iterator TIterator;
        TIterator begin1, end2, it;
	begin1 = _getDiagSet(secondSeed).begin();
	end2 = _getDiagSet(secondSeed).end();
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		while (_getLastDiag(firstSeed).i1 > _getFirstDiag(secondSeed).i1)// || (_getLastDiag(firstSeed).i2 > _getFirstDiag(secondSeed).i2))
		{
			score1 -= firstSeed.seedSet.back().i3*scoreMatch(scoreMatrix);
			TPosition x1 = firstSeed.seedSet.back().i1;
			TPosition x2 = firstSeed.seedSet.back().i2;
			firstSeed.seedSet.pop_back();
			score1 -= (abs(rightDim0(firstSeed)-x1)+abs(rightDim1(firstSeed)-x2))*scoreGap(scoreMatrix); 
		}
		score1 += abs(endDiagonal(firstSeed) - startDiagonal(secondSeed))*scoreGap(scoreMatrix);
		score1 -= (_max(abs(rightDim0(firstSeed)- leftDim0(secondSeed)),abs(rightDim1(firstSeed)-leftDim1(secondSeed)))+1)*scoreMatch(scoreMatrix);

		if ((rightDim0(firstSeed) < leftDim0(secondSeed)) && (rightDim1(firstSeed) < leftDim1(secondSeed))) {
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		} else {
		if (startDiagonal(secondSeed) == endDiagonal(firstSeed)){
			setRightDim1(firstSeed,(*begin1).i2+(*begin1).i3-1);
			++begin1;
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		}
		else {
			
			TPosition tmp = startDiagonal(secondSeed) - endDiagonal(firstSeed);
			if (tmp < 0){
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- leftDim1(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			} else {
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - leftDim0(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			}
		}	
		}
		
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}

// Sets the begin and end position as well as the left and right diagonal of a ChainedSeed after seed extension
template<typename TPosition, typename TBound, typename TExtension, typename TSize>
void
_setExtendedSeedDimensions(Seed<TPosition, ChainedSeed> & seed,
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
        seed.seedSet.push_front(Triple<TPosition, TPosition, TPosition>
                                       (leftDim0(seed)-extLengthQuery, leftDim1(seed)-extLengthDatabase, 1));
    } else {
	    // set left and right diagonals
	    if (rightDiagonal(seed) > endDiagonal(seed)-upperBound)
	        setRightDiagonal(seed, endDiagonal(seed)-upperBound);
        if (leftDiagonal(seed) < endDiagonal(seed)+lowerBound)
            setLeftDiagonal(seed, endDiagonal(seed)+lowerBound);

        // set new end position of seed
        seed.seedSet.push_back(Triple<TPosition, TPosition, TPosition>
                                      (rightDim0(seed)+extLengthQuery, rightDim1(seed)+extLengthDatabase, 1));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      Alignment Construction Alogrithm										  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
.Function.getAlignment:
..class:Spec.ChainedSeed
..summary: Constructs a alignment from a @Spec.ChainedSeed@.
..cat:Seed Handling
..signature:getAlignment(seed, align, query, database, scoreMatrix)
..param.seed: The alignment foundation.
...type:Spec.ChainedSeed
..param.align: An emtpy alignment object, that stores the constructed alignment.
...type:Class.Align
..param.query:The Query sequence.
...type:Class.String
..param.database:The database sequence.
...type:Class.String
..param.scoreMatrix:The scoring matrix.
...type:Spec.Simple Score
..returns: Score of the alignment.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TText, typename TPosition2>
int
getAlignment(Seed<TPosition,ChainedSeed> &seed,
			 Align<String<TText>, ArrayGaps> &aligned, 
			 String<TText> &query, 
			 String<TText> &database, 
			 Score<TPosition2, Simple> &scoreMatrix)
{
	SEQAN_CHECKPOINT
	int seedScore = 0;
	typename std::list<Triple< TPosition, TPosition, TPosition> > seedList = _getDiagSet(seed);
	typedef typename std::list<Triple< TPosition, TPosition, TPosition> >::iterator TIterator;

	resize(rows(aligned), 2);
	assignSource(row(aligned, 0), query);
	assignSource(row(aligned, 1), database);

	setClippedBeginPosition(row(aligned, 0), leftDim0(seed));
	setClippedBeginPosition(row(aligned, 1), leftDim1(seed));
	setBeginPosition(row(aligned, 0), 0);
	setBeginPosition(row(aligned, 1), 0);
	setClippedEndPosition(row(aligned, 0), rightDim0(seed)+1);
	setClippedEndPosition(row(aligned, 1), rightDim1(seed)+1);

	TIterator it1 = seedList.begin();
	TIterator it2 = ++seedList.begin();
	
	for (int i =0; i<(*it1).i3;++i){
		seedScore += score(scoreMatrix,(*it1).i1+i, (*it1).i2+i, query, database); 
	}
	
	if (seedList.size()>=2){
		TPosition gapLength;
		TPosition position1 =	it1->i1 - leftDim0(seed) + it1->i3;
		TPosition position2 =	it1->i2 - leftDim1(seed) + it1->i3;
		while (it2 != seedList.end()){
			if (it2->i1 == it1->i1 + it1->i3){ //query teile zusammen
				gapLength = it2->i2 - it1->i2 - it1->i3;
				insertGaps(row(aligned,0),position1, gapLength);
				seedScore += gapLength*scoreGap(scoreMatrix);
				position1 += gapLength;
			}
			else
				if ((*it2).i2 == it1->i2+it1->i3)
				{
					gapLength = it2->i1 - it1->i1 - it1->i3;
					insertGaps(row(aligned,1),position2,gapLength);
					position2 += gapLength;
					seedScore += gapLength*scoreGap(scoreMatrix);
				} 
				else 
				{
					Align<String<TText>, ArrayGaps> alignSeg;
					resize(rows(alignSeg), 2);
					assignSource(row(alignSeg, 0), query);
					assignSource(row(alignSeg, 1), database);
					
					setClippedBeginPosition(row(alignSeg, 0), (*it1).i1+(*it1).i3);
					setClippedBeginPosition(row(alignSeg, 1), (*it1).i2+(*it1).i3);
					setBeginPosition(row(alignSeg, 0), 0);
					setBeginPosition(row(alignSeg, 1), 0);
					setClippedEndPosition(row(alignSeg, 0), (*it2).i1);
					setClippedEndPosition(row(alignSeg, 1), (*it2).i2);

					seedScore += globalAlignment(alignSeg,scoreMatrix,NeedlemanWunsch());//needlemanWunsch(alignSeg,scoreMatrix);
			
					unsigned int j;
					bool gap;
					if (row(alignSeg,1).data_arr[0] == 0){
						j=1;
						gap=false;
					} else {
						j=0;
						gap=true;
					}
					while (j < length(row(alignSeg,1).data_arr)){
						if (gap){
							insertGaps(row(aligned,1),position2,row(alignSeg,1).data_arr[j]);
					
							gap = false;
							position2 += row(alignSeg,1).data_arr[j];
						} else {
							gap = true;
							position2 += row(alignSeg,1).data_arr[j];							
						}
						++j;
					}
					if (row(alignSeg,0).data_arr[0] == 0){
						j=1;
						gap=false;
					} else {
						j=0;
						gap=true;
					}
					 while (j < length(row(alignSeg,0).data_arr)){
						if (gap){
							insertGaps(row(aligned,0),position1,row(alignSeg,0).data_arr[j]);
							gap = false;
							position1+= row(alignSeg,0).data_arr[j];
						} else {
							gap = true;
							position1+= row(alignSeg,0).data_arr[j];							
						}
						++j;
					}
					TPosition tmp1, tmp2;
					tmp1 = position1;
					tmp2 = position2;
					if (tmp1 > tmp2){
						insertGaps(row(aligned,1),position2,tmp1-tmp2);
						position2 += tmp1-tmp2;
					} else 
						if (tmp2 > tmp1){
							insertGaps(row(aligned,0),position1,tmp2-tmp1);
							position1 += tmp2-tmp1;
						}
				}
				position1+= it2->i3;
				position2 += it2->i3;
				for (int i =0; i<it2->i3;++i)
					seedScore += score(scoreMatrix, it2->i1+i, it2->i2+i, query, database);

				++it1;
				++it2;
		}
	}
	return seedScore;
}


/**
.Function.scoreSeed:
..class:Spec.ChainedSeed
..summary: Calculates the score of a seed. 
..cat:Seed Handling
..signature:scoreSeed(seed, query, database, scoreMatrix)
..param.seed: A seed.
...type:Spec.ChainedSeed
..param.query:The Query sequence.
...type:Class.String
..param.database:The database sequence.
...type:Class.String
..param.scoreMatrix:The scoring sheme.
...type:Spec.Simple Score
..returns: Score of the seed.
..remarks: Score has not the same value as the resulting alignment. Gaps between diagonals matches are scored as full length gaps.
..include:seqan/seeds.h
*/
template<typename TPosition, typename TText, typename TScore>
TScore
scoreSeed(Seed<TPosition, ChainedSeed> &seed, String<TText> &query, String<TText> &database, Score<TScore, Simple> &matrix){
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple< TPosition, TPosition, TPosition> >::iterator TIterator;
	int tmpScore =0;
	TIterator it = _getDiagSet(seed).begin();
	for (int i = 0; i < it->i3; ++i){
		tmpScore+=score(matrix, it->i1+i, it->i2+i, query, database);
	}

	if (_getDiagSet(seed).size()>=2){
		TIterator it_end = _getDiagSet(seed).end();
		for (TIterator it2 = ++_getDiagSet(seed).begin(); it2!= it_end; it2++){
			for (int i = 0; i < it2->i3; ++i){
				tmpScore+=score(matrix, it2->i1+i, it2->i2+i, query, database);
			}
			tmpScore += scoreGap(matrix)*(it2->i2-(it->i2+it->i3)+it2->i1-(it->i1+it->i3));
			++it;
		}
	}
	return tmpScore;
}


} //namespace Seqan

#endif //#ifndef SEQAN_HEADER_
