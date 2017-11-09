 /*==========================================================================
                RazerS Spliced - Fast Split Read Mapping 
              
 ============================================================================
  Copyright (C) 2008 by Anne-Katrin Emde

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
 ==========================================================================*/

#ifndef SEQAN_HEADER_RAZERS_SPLICED_H
#define SEQAN_HEADER_RAZERS_SPLICED_H

//#define RAZERS_DEBUG
//#define RAZERS_DEBUG_LIGHT

#include <seqan/align.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/graph_align.h>
#include <seqan/seeds.h>

namespace seqan
{



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Match statistics stuff
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//statistics, expected number of matches
double choose(int n, int k)
{
	double cnk = 1.0;
  
	if (k * 2 > n) 
	k = n - k;
    
	for (int i = 1; i <= k; n--, i++)
	{
		cnk /= i;
		cnk *= n;
	}
	return cnk;
}


//#R  - read length
//#T  - total sequence length
//#I  - length of insertion
//#D  - length of deletion
//#M1 - minimal prefix match length
//#M2 - minimal suffix match length


//## first function denotes the random match to an iid sequence
template<typename TSize, typename TValue>
double
_probability(TSize R, TSize M1, TSize M2, TValue d, TValue d_m1, TValue d_m2)
{
	double sum_prob = 0;
	for (unsigned i1 = 0; i1 < d_m1; ++i1)
	{
		for (unsigned i2 = 0; i2 < d_m2; ++i2)
		{
			for (unsigned j = 0; j < d-i1-i2; ++j)
			{
				sum_prob += (double)(choose(M1,i1)* pow((double)0.25,(double)R-i1-i2-j) * choose(M2,i2) * pow((double)0.75,(double)i1+i2+j) * choose(R-M1-M2,j));
			}
		}
	}
	return sum_prob;
}


//#InsCuts describes the possible number of cuts
//#for a read of length R and an insertion of length I
template<typename TSize>
TSize
_insCuts(TSize R, TSize M1, TSize M2, TSize I){
	return (R-I-M1-M2+1);
}



//#function for calculating expected matches with I insertions
template<typename TSize, typename TValue>
double
_expMatchesIns(TSize R, TSize M1, TSize M2, TSize I, TValue d,TValue d_m1,TValue d_m2, TSize S, TSize T)
{
	return((double)T*((S-R-I-1)*_probability((R-I),M1,M2,d,d_m1,d_m2))*_insCuts(R,M1,M2,I));
}


//#DelCuts describes the possible number of cuts
//#for a read of length R in case of a deletion
template<typename TSize>
TSize
_delCuts(TSize R, TSize M1, TSize M2)
{
	return (R-M1-M2+1);
}


//#Del describes the possible number of configurations
//#for a string of length R 
template<typename TSize>
TSize
_del(TSize R, TSize S)
{
	return((S-R)*(S-R+1)/2);
}


//#function for calculating expected matches with deletion of size D
template<typename TSize, typename TValue>
double
_expMatchesDel(TSize R, TSize M1, TSize M2, TValue d, TValue d_m1, TValue d_m2, TSize S, TSize T)
{
	return ((double)T*(_probability(R,M1,M2,d,d_m1,d_m2))*_delCuts(R,M1,M2)*_del(R,S));
}

template<typename TReadSet, typename TSize, typename TOptions>
void
expNumRandomMatches(TReadSet &readSet, TSize genomeLen, TOptions & options)
{
	TSize M1 = options.minMatchLen;
	TSize M2 = options.minMatchLen;
	TSize d_m1 = (int) M1 * options.errorRate;
	TSize d_m2 = (int) M2 * options.errorRate;
	TSize numReads = length(readSet);
	TSize readLen = (numReads > 0) ? length(readSet[0]) : 0;
	TSize numErrors = (int) readLen * options.errorRate;
	
	//expected number of random deletion matches:
	double delMatches = _expMatchesDel(readLen,M1,M2,numErrors,d_m1, d_m2, genomeLen, numReads);
	
	//expected number of random deletion matches:
	double insMatches = 0;
	for(TSize insLen = 1; insLen <=readLen-M1-M2; ++insLen)
		insMatches += _expMatchesIns(readLen,M1,M2,numErrors,insLen,d_m1, d_m2, genomeLen, numReads);
	
	::std::cout << "Expected number of random deletion-indicating matches: " << delMatches << std::endl;
	::std::cout << "Expected number of random insertion-indicating matches: " << insMatches << std::endl;
 }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// SPLIT MAPPING
// We build two q-gram indices, one for the prefix, the other for the suffix of the read

struct LongestPrefix{};
struct LongestSuffix{};

struct OrientationReverse{};
struct OrientationForward{};




//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches # necessary to have an own splicedmatch function? 
// planned specs: SpliceSite, General, ... 
template < typename TMatches, typename TCounts, typename TSpec, typename TSwiftL, typename TSwiftR >
void compactSplicedMatches(TMatches &matches, 
			TCounts & /*cnts*/, 
			RazerSOptions<TSpec> &options, 
			bool compactFinal, 
			TSwiftL &swiftL, 
			TSwiftR &swiftR)
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreDistCutOff = std::numeric_limits<int>::min();
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;
	
	// sort 
	::std::sort(it, itEnd, LessSplicedErrors<TMatch>());
//	::std::sort(it, itEnd, LessSplicedScore<TMatch>());
	int counter = 0;
	for (; it != itEnd; ++it) 
	{
		++counter;
		if ((*it).orientation == '-') { ++it; continue; }
		if (readNo == (*it).rseqNo)
		{ 
			if ((*it).pairScore <= scoreDistCutOff)
			{
				++it;
				continue;
			}
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = - (*it).pairScore - 1;
					if (options.purgeAmbiguous && (options.distanceRange == 0 || maxErrors <= (int) options.distanceRange))
						maxErrors = -1;
					
					setMaxErrors(swiftL, readNo, maxErrors);
					setMaxErrors(swiftR, readNo, maxErrors);
					
					if (maxErrors == -1 && options._debugLevel >= 2)
						::std::cerr << "(read #" << readNo << " disabled)";
					
					if(options.purgeAmbiguous)
	     				{
						if (options.distanceRange == 0 || -(*it).pairScore < (int) options.distanceRange || compactFinal)
							dit = ditBeg;
						else {
							*dit = *it;	++dit; ++it;
							*dit = *it;	++dit;
						}
					}
		
				}
#endif
				++it;
				continue;
			}
		}
		else
		{
			readNo = (*it).rseqNo;// >> 1;
			hitCount = 0;
			if (options.distanceRange > 0)
				scoreDistCutOff = (*it).pairScore - options.distanceRange;
			ditBeg = dit;
		}
		*dit = *it;	++dit; ++it;
		*dit = *it;	++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
}



//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches # necessary to have an own splicedmatch function? 
template < typename TMatches, typename TCounts, typename TSpec, typename TSwiftL, typename TSwiftR >
void compactSplicedMatchesPurgeAmbiguous(TMatches &matches, TCounts & /*cnts*/, RazerSOptions<TSpec> &options, TSwiftL &, TSwiftR &)
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreDistCutOff = std::numeric_limits<int>::min();
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;
	
	// sort
	::std::sort(it, itEnd, LessSplicedErrors<TMatch>());
//	::std::sort(it, itEnd, LessSplicedScore<TMatch>());
	int counter = 0;
	for (; it != itEnd; ++it) 
	{
		++counter;
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo)
		{ 
			if ((*it).pairScore <= scoreDistCutOff)
			{
				++it;
				continue;
			}
			if (++hitCount >= hitCountCutOff)
			{
				if (hitCount == hitCountCutOff)
					dit = ditBeg;
				++it;
				continue;
			}
		}
		else
		{
			readNo = (*it).rseqNo;
			hitCount = 0;
			if (options.distanceRange > 0)
				scoreDistCutOff = (*it).pairScore - options.distanceRange;
			ditBeg = dit;
		}
		*dit = *it;	++dit; ++it;
		*dit = *it;	++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
}


	template <typename TReadMatch>
	struct LessReadNoPairErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// quality
			if (a.pairScore > b.pairScore) return true;
			if (a.pairScore < b.pairScore) return false;
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;
			return a.editDist < b.editDist;
		}
	};



//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template < typename TMatches, typename TCounts >
void countSplitMatches(TMatches &matches, TCounts &cnt)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	typedef typename Value<TCounts>::Type					TRow;
	typedef typename Value<TRow>::Type						TValue;


	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	::std::sort(it, itEnd, LessReadNoPairErrors<TMatch>());
	
	unsigned readNo = -1;
	short editDist = -1;
	int64_t count = 0;
	int64_t maxVal = std::numeric_limits<TValue>::max();

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo &&
			-editDist == (*it).pairScore)
			++count;
		else
		{
			if (readNo != (unsigned)-1 && (unsigned)editDist < length(cnt))
				cnt[editDist][readNo] = (maxVal < count)? maxVal : count;
			readNo = (*it).rseqNo;
			editDist = -(*it).pairScore;
			count = 1;
		}
	}
	if (readNo != (unsigned)-1 && (unsigned)editDist < length(cnt))
		cnt[editDist][readNo] = count;
}


template<typename TAlign, typename TPosition>
int
countErrorsInAlign(TAlign & align, TPosition end_)
{
	
	typedef typename Source<TAlign>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

	TAlignIterator ali_it0 = iter(row(align,0),0);
	TAlignIterator ali_it1 = iter(row(align,1),0);					
	TAlignIterator ali_it0_stop = iter(row(align,0),end_);
	TAlignIterator ali_it1_stop = iter(row(align,1),end_);


	int errors = 0;
	while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
	{
		while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1))
		{
			if(*ali_it1 != *ali_it0)
				++errors;
			++ali_it0;
			++ali_it1;
		}
		while(ali_it0!=ali_it0_stop && isGap(ali_it0))
		{
			++ali_it0;
			++ali_it1;
			++errors;
		}
		while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
		{
			++ali_it0;
			++ali_it1;
			++errors;
		}
	}
	while(ali_it0!=ali_it0_stop)
	{
		++ali_it0;
		++errors;
	}
	while(ali_it1 != ali_it1_stop)
	{
		++ali_it1;
		++errors;
	}
	
	return errors;
}


			
//////////////////////////////////////////////////////////////////////////////
// Edit distance verification for longest suffix/prefix
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TMyersPatterns,
	typename TSpec,
	typename TSufPrefSpec
>
inline bool
matchVerify(
	TMatch &m,							// resulting match
	Segment<TGenome, InfixSegment> inf,	// potential match genome region
	unsigned rseqNo,					// read number
	TReadSet &readSet,	    			// reads
	TMyersPatterns &forwardPatterns,	// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,// RazerS options
	SwiftSemiGlobal,					// Edit distance
	TSufPrefSpec)						// Swift specialization
{
	typedef Segment<TGenome, InfixSegment>				TGenomeInfix;
	typedef typename Value<TReadSet>::Type				TRead;
	
	// find read match end
	typedef Finder<TGenomeInfix>					TMyersFinder;
	typedef typename Value<TMyersPatterns>::Type			TMyersPattern;
	
	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>		TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>			TReadRev;
	typedef Finder<TGenomeInfixRev>					TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>			TMyersPatternRev;

	TMyersFinder myersFinder(inf);
	TMyersPattern &myersPattern = forwardPatterns[rseqNo];  //have to make sure this only contains the prefix
	
#ifdef RAZERS_DEBUG
	::std::cout << "Verify: " << ::std::endl;
	::std::cout << "Genome: " << inf << "\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout << "Read:   " << readSet[rseqNo] << ::std::endl;
#endif
	
	unsigned ndlLength = _min(sequenceLength(rseqNo, readSet),options.minMatchLen);
	int maxScore = std::numeric_limits<int>::min();
	int minScore = -(int)(ndlLength * options.errorRate);
	TMyersFinder maxPos;
	
	// find end of best semi-global alignment
	while (find(myersFinder, myersPattern, minScore))
	{
		if (maxScore <= getScore(myersPattern)) 
		{
			maxScore = getScore(myersPattern);
			maxPos = myersFinder;
		}
	}
	

	if (maxScore < minScore)
		return false;

	m.editDist	= (unsigned)-maxScore;
	TGenomeInfix oriInf = inf;
	setEndPosition(inf, m.gEnd = (beginPosition(inf) + position(maxPos) + 1));
	
	// limit the beginning to needle length plus errors (== -maxScore)
	if (length(inf) > ndlLength - maxScore)
		setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
	
	// find beginning of best semi-global alignment
	TGenomeInfixRev		infRev(inf);
	TMyersFinderRev		myersFinderRev(infRev);
	TReadRev			readRev;
    TRead               readInf;
	if(IsSameType<TSufPrefSpec,LongestSuffix>::VALUE)
		readInf = infix(readSet[rseqNo],length(readSet[rseqNo])-options.minMatchLen,length(readSet[rseqNo]));
	else
		readInf = infix(readSet[rseqNo],0,options.minMatchLen);
    setHost(readRev, readInf);

	TMyersPatternRev	myersPatternRev(readRev);
	
	_patternMatchNOfPattern(myersPatternRev, options.matchN);
	_patternMatchNOfFinder(myersPatternRev, options.matchN);
	while (find(myersFinderRev, myersPatternRev, maxScore))
		m.gBegin = m.gEnd - (position(myersFinderRev) + 1);
	
	m.mScore = ndlLength;
	m.seedEditDist = m.editDist;
	m.gSeedLen = m.gEnd - m.gBegin;

#ifdef RAZERS_DEBUG
	::std::cout << " before extendMatch " << ::std::endl;
	::std::cout << " match: " << ::std::endl;
	::std::cout << " mScore= " <<m.mScore << ::std::endl;
	::std::cout << " gBegin= " <<m.gBegin << ::std::endl;
	::std::cout << " gEnd= " <<m.gEnd << ::std::endl;
	::std::cout << " edit= " <<m.editDist << ::std::endl;
#endif
	
	//TODO: give only part of read to extension!!!

	//now extend the seed
	extendMatch(readSet,rseqNo,oriInf,m,options,TSufPrefSpec());

#ifdef RAZERS_DEBUG
	::std::cout << " match: " << ::std::endl;
	::std::cout << " mScore= " <<m.mScore << ::std::endl;
	::std::cout << " gBegin= " <<m.gBegin << ::std::endl;
	::std::cout << " gEnd= " <<m.gEnd << ::std::endl;
	::std::cout << " edit= " <<m.editDist << ::std::endl;
#endif

	return true;
}

template<typename TReadSet, typename TSize, typename TInf, typename TMatch, typename TOptions>
void
extendMatch(TReadSet &readSet, TSize rseqNo, TInf & inf, TMatch &m, TOptions &options, LongestSuffix)
{
#ifdef RAZERS_DEBUG
	::std::cout << " extending match left" << ::std::endl;
#endif

	unsigned lDim0 = _max(0,(int)length(readSet[rseqNo])-(int)options.minMatchLen);
	unsigned lDim1 = m.gBegin - beginPosition(inf);
	unsigned rDim0 = length(readSet[rseqNo])-1;
	unsigned rDim1 = m.gEnd - beginPosition(inf)-1;
	Seed<int,SimpleSeed> seed(lDim0, lDim1, rDim0, rDim1);
	Score<int> scoreMatrix(0,-1,-1,-1);
	int scoreDropOff = (sequenceLength(rseqNo,readSet) * options.errorRate) - m.editDist;

#ifdef RAZERS_DEBUG
	::std::cout << "beginPos = " << beginPosition(inf) << std::endl;
	::std::cout << "endPos = " << endPosition(inf) << std::endl;
	::std::cout << " lDim0: " << lDim0 << ::std::endl;
	::std::cout << " lDim1: " << lDim1 << ::std::endl;
	::std::cout << " rDim0: " << rDim0 << ::std::endl;
	::std::cout << " rDim1: " << rDim1 << ::std::endl;
	::std::cout << " scoreDropOff: "<<scoreDropOff << ::std::endl;
	::std::cout << " readInf: "<< infix(readSet[rseqNo],lDim0,rDim0+1) << ::std::endl;
	::std::cout << " gInfInf: "<< infix(inf,lDim1,rDim1+1) << ::std::endl;
	::std::cout << " read: "<< readSet[rseqNo] << ::std::endl;
	::std::cout << " gInf: "<< inf << ::std::endl;
#endif

	int extScore = 0;
	extendSeedScore(seed,extScore,scoreDropOff,scoreMatrix, readSet[rseqNo],inf,0,GappedXDrop());
	m.gBegin = leftDim1(seed) + beginPosition(inf);
	m.mScore = rightDim0(seed) - leftDim0(seed) + 1;
	m.editDist -= extScore;

#ifdef RAZERS_DEBUG
	::std::cout << " lDim0: " << leftDim0(seed) << ::std::endl;
	::std::cout << " lDim1: " << leftDim1(seed) << ::std::endl;
	::std::cout << " rDim0: " << rightDim0(seed) << ::std::endl;
	::std::cout << " rDim1: " << rightDim1(seed) << ::std::endl;
	::std::cout << " scoreDropOff: "<<scoreDropOff << ::std::endl;
	::std::cout << " readInf: "<< infix(readSet[rseqNo],leftDim0(seed),rightDim0(seed)+1) << ::std::endl;
	::std::cout << " gInfInf: "<< infix(inf,leftDim1(seed),rightDim1(seed)+1) << ::std::endl;
	::std::cout << " read: "<< readSet[rseqNo] << ::std::endl;
	::std::cout << " gInf: "<< inf << ::std::endl;
#endif
}

template<typename TReadSet, typename TSize, typename TInf, typename TMatch, typename TOptions>
void
extendMatch(TReadSet &readSet, TSize rseqNo, TInf & inf, TMatch &m, TOptions &options, LongestPrefix)
{
#ifdef RAZERS_DEBUG
	::std::cout << " extending match right" << ::std::endl;
#endif

	unsigned lDim0 = 0;
	unsigned lDim1 = m.gBegin - beginPosition(inf);
	unsigned rDim0 = _min(options.minMatchLen,length(readSet[rseqNo])) - 1;
	unsigned rDim1 = m.gEnd - beginPosition(inf) - 1;
	Seed<int,SimpleSeed> seed(lDim0, lDim1, rDim0, rDim1);
	Score<int> scoreMatrix(0,-1,-1,-1);
	int scoreDropOff = (sequenceLength(rseqNo,readSet) * options.errorRate) - m.editDist;

#ifdef RAZERS_DEBUG
	::std::cout << "beginPos = " << beginPosition(inf) << std::endl;
	::std::cout << "endPos = " << endPosition(inf) << std::endl;
	::std::cout << " lDim0: " << lDim0 << ::std::endl;
	::std::cout << " lDim1: " << lDim1 << ::std::endl;
	::std::cout << " rDim0: " << rDim0 << ::std::endl;
	::std::cout << " rDim1: " << rDim1 << ::std::endl;
	::std::cout << " scoreDropOff: "<<scoreDropOff << ::std::endl;
	::std::cout << " readInf: "<< infix(readSet[rseqNo],lDim0,rDim0+1) << ::std::endl;
	::std::cout << " gInfInf: "<< infix(inf,lDim1,rDim1+1) << ::std::endl;
	::std::cout << " read: "<< readSet[rseqNo] << ::std::endl;
	::std::cout << " gInf: "<< inf << ::std::endl;
#endif

	int extScore = 0;
	extendSeedScore(seed,extScore,scoreDropOff,scoreMatrix, readSet[rseqNo],inf,1,GappedXDrop());
	m.gEnd = rightDim1(seed) + 1 + beginPosition(inf);
	m.mScore = rightDim0(seed) - leftDim0(seed) + 1;
	m.editDist -= extScore;

#ifdef RAZERS_DEBUG
	::std::cout << " lDim0: " << leftDim0(seed) << ::std::endl;
	::std::cout << " lDim1: " << leftDim1(seed) << ::std::endl;
	::std::cout << " rDim0: " << rightDim0(seed) << ::std::endl;
	::std::cout << " rDim1: " << rightDim1(seed) << ::std::endl;
	::std::cout << " scoreDropOff: "<<scoreDropOff << ::std::endl;
	::std::cout << " readInf: "<< infix(readSet[rseqNo],leftDim0(seed),rightDim0(seed)+1) << ::std::endl;
	::std::cout << " gInfInf: "<< infix(inf,leftDim1(seed),rightDim1(seed)+1) << ::std::endl;
	::std::cout << " read: "<< readSet[rseqNo] << ::std::endl;
	::std::cout << " gInf: "<< inf << ::std::endl;
#endif
}


//template <
//	typename TMatch, 
//	typename TGenomeSegment, 
//	typename TReadSet, 
//	typename TDummy,
//	typename TSpec >
//inline bool
//matchVerify(
//	TMatch &m,					// resulting match
//	TGenomeSegment inf,				// potential match genome region
//	unsigned rseqNo,				// read number
//	TReadSet &readSet,				// reads
//	TDummy&,
//	RazerSOptions<TSpec> const &options,		// RazerS options
//	LongestSuffix)				// LongestEditPrefix within errorRate
//{
//	::std::cout << "vcejsfhkjdfhksh11111\n";
//	if( matchVerify(m,inf,rseqNo,
//		readSet, false,
//		options,LongestHammingSuffix()))
//	{
//	//	unsigned tmp = m.gBegin;
//	//	m.gBegin = m.gEnd;
//	//	m.gEnd = tmp;
//		return true;
//	}
//	else return false;
//	
//}
//
//
//template <
//	typename TMatch, 
//	typename TGenomeSegment, 
//	typename TReadSet, 
//	typename TDummy,
//	typename TSpec >
//inline bool
//matchVerify(
//	TMatch &m,					// resulting match
//	TGenomeSegment inf,		// potential match genome region
//	unsigned rseqNo,				// read number
//	TReadSet &readSet,				// reads
//	TDummy & ,
//	RazerSOptions<TSpec> const &options,		// RazerS options
//	LongestPrefix)				// LongestEditPrefix within errorRate
//{
//	::std::cout << "vcejsfhkjdfhksh222222\n";
//	return matchVerify(m,inf,rseqNo,
//		readSet, true,
//		options,LongestHammingPrefix());	
//}


template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet,
	typename TPattern,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,									// resulting match
	Segment<TGenome,InfixSegment>  genomeInf,	// potential match genome region
	unsigned rseqNo,							// read number
	TReadSet& readSet,							// original readSet
	TPattern&,					
	RazerSOptions<TSpec> const &options,		// RazerS options
	SwiftSemiGlobalHamming,						// HammingDistance
	LongestPrefix)								// LongestPrefix
{
	
	typedef Segment<TGenome, InfixSegment>                  TGenomeInfix;
	typedef typename Size<TGenomeInfix>::Type               TSize;
	typedef typename Value<TGenomeInfix>::Type              TDna;
	typedef typename Position<TGenomeInfix>::Type           TPosition;
	typedef typename Value<TReadSet>::Type 			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Infix<TRead>::Type 			TReadInf;
	typedef typename Iterator<TReadInf, Standard>::Type	TReadIterator;
	
	if (length(genomeInf) < options.minMatchLen) return false;
	TReadInf read = infix(readSet[rseqNo],0,length(readSet[rseqNo])-options.minMatchLen);
	

	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git	= begin(genomeInf, Standard());
	TGenomeIterator gitEnd	= end(genomeInf, Standard()) - (length(read) - 1);
	
	// this is max number of errors the seed should have
	unsigned maxErrorsSeed = (unsigned)(options.minMatchLen * options.errorRate);	
	unsigned maxTotalErrors = (unsigned)(length(read) * options.errorRate);	
	unsigned minSeedErrors = maxErrorsSeed + 1;
	unsigned minTotalErrors = maxTotalErrors + 1;
	unsigned bestHitLength = 0;

	for (; git < gitEnd; ++git)
	{
		bool hit = true;
		unsigned hitLength = 0;
		unsigned count = 0;
		unsigned seedErrors = 0;
		unsigned totalErrors = 0;
		TGenomeIterator g = git;	
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
		{
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
			{
				if (count < options.minMatchLen)
				{
					++totalErrors;
					if(++seedErrors > maxErrorsSeed)
					{
						hit = false;
						break;
					}
				}
				else
				{
					if(++totalErrors > maxTotalErrors)
					{
						// we are excluding this last error position 
						--totalErrors;
						break;
					}
				}
			}
			++count;
		}
		if (hit) hitLength = count;
		if (hitLength > bestHitLength ) //simply take the longest hit
		{
			minSeedErrors = seedErrors;
			minTotalErrors = totalErrors;
			bestHitLength = hitLength;
			m.gBegin = git - begin(host(genomeInf), Standard());
		}
	}



	if(bestHitLength < options.minMatchLen) 
		return false;
	
	m.gEnd = m.gBegin + bestHitLength;
	m.mScore = bestHitLength;
	m.editDist = minTotalErrors;
		
#ifdef RAZERS_DEBUG
		std::cout << "m.gBeg  =" << m.gBegin << "\n";
		std::cout << "m.gEnd  =" << m.gEnd << "\n";
		std::cout << "m.edit  =" << m.editDist << "\n";
		std::cout << "m.mScore=" << m.mScore << "\n\n";
#endif

	return true;
	    
}


template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TPattern,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,									// resulting match
	Segment<TGenome,InfixSegment>  genomeInf,	// potential match genome region
	unsigned rseqNo,							// read number
	TReadSet & readSet,							// original readSet
	TPattern&,
	RazerSOptions<TSpec> const &options,		// RazerS options
	SwiftSemiGlobalHamming,
	LongestSuffix)								// LongestSuffix
{

	typedef Segment<TGenome, InfixSegment>                  TGenomeInfix;
	typedef typename Size<TGenomeInfix>::Type               TSize;
	typedef typename Value<TGenomeInfix>::Type              TDna;
	typedef typename Position<TGenomeInfix>::Type           TPosition;
	typedef typename Value<TReadSet>::Type 			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Infix<TRead>::Type 			TReadInf;
	typedef typename Iterator<TReadInf, Standard>::Type	TReadIterator;
	
	if (length(genomeInf) < options.minMatchLen) return false;
	bool debug = false;
	TRead read = infix(readSet[rseqNo],options.minMatchLen,length(readSet[rseqNo]));
	
		
	if(debug)
	{
		::std::cout<< "suffixmatching\n";
		::std::cout << "genome=" << genomeInf << "\nread  =" << read <<"\n";
	}

	TReadIterator ritEnd	= end(read, Standard())-1;
	TReadIterator ritBeg	= begin(read, Standard());
	TGenomeIterator git	= end(genomeInf, Standard())-1;
	TGenomeIterator gitBeg	= begin(genomeInf, Standard()) + options.minMatchLen;
	
	// this is max number of errors the seed should have
	unsigned maxErrorsSeed = (unsigned)(options.minMatchLen * options.errorRate);	
	unsigned maxTotalErrors = (unsigned)(length(read) * options.errorRate);	
	unsigned minSeedErrors = maxErrorsSeed + 1;
	unsigned minTotalErrors = maxTotalErrors + 1;
	unsigned bestHitLength = 0;

	for (; git > gitBeg; --git)
	{
		bool hit = true;
		unsigned hitLength = 0;
		unsigned count = 0;
		unsigned seedErrors = 0;
		unsigned totalErrors = 0;
		TGenomeIterator g = git;	
		for(TReadIterator r = ritEnd; r >= ritBeg; --r, --g)
		{
			if(debug)::std::cout << *r << "\t" << *g << "\n";
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
			{
				if (count < options.minMatchLen)
				{
					++totalErrors;
					if(++seedErrors > maxErrorsSeed)
					{
					//	if(debug) ::std::cout << "-->no\n";
						hit = false;
						break;
					}
				}
				else
				{
					if(++totalErrors > maxTotalErrors)
					{
						// we are excluding this last error position 
						--totalErrors;
						break;
					}
				}
			}
			++count;
		}
		if (hit) hitLength = count;
		if (hitLength > bestHitLength ) //simply take the longest hit
		{
			minSeedErrors = seedErrors;
			minTotalErrors = totalErrors;
			bestHitLength = hitLength;
			m.gEnd = git - begin(host(genomeInf), Standard()) + 1;
			if(debug) ::std::cout << "m.gEnd=" << m.gEnd << ::std::endl;
			
		}
	}


	if(bestHitLength < options.minMatchLen) 
		return false;
	
	
	m.gBegin = m.gEnd - bestHitLength;
	m.mScore = bestHitLength;
	m.editDist = minTotalErrors;
		
#ifdef RAZERS_DEBUG
	std::cout << "m.gBeg  =" << m.gBegin << "\n";
	std::cout << "m.gEnd  =" << m.gEnd << "\n";
	std::cout << "m.edit  =" << m.editDist << "\n";
	std::cout << "m.mScore=" << m.mScore << "\n\n";
#endif	

	return true;
	    
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




template <typename TScore>
bool
findBestSplitPosition(String<Pair<TScore,int> > & maxColsL,
					  String<Pair<TScore,int> > & maxColsR,
					  int & rowPosL1,
					  int rowPosL2,
					  int & rowPosR1,
					  int rowPosR2,
					  int seq0Len,
					  OrientationForward,
					  SwiftSemiGlobal)
{

#ifdef RAZERS_DEBUG
	::std::cout << "findBestSplitEditForward\n";
#endif

	TScore maxSum = std::numeric_limits<TScore>::min();
	int bestL = rowPosL1;
	int bestR = rowPosR1;
	while (rowPosL1 <= rowPosL2 && rowPosR1 >= rowPosR2)
	{
		if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 > maxSum 
			&& (maxColsL[rowPosL1].i2 + maxColsR[rowPosR1].i2 <= seq0Len))  // this is to prevent same bases from being used in both prefix and suffix match
		{																	// this works, because we store the FIRST bestScore in each row, i.e. 
			maxSum = maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1;
			bestL = rowPosL1;
			bestR = rowPosR1;
		}
		++rowPosL1;
		--rowPosR1;
	}
	rowPosL1 = bestL;
	rowPosR1 = bestR;

	return true;
}


	

// Edit distance match combination
template <typename TScore>
bool
findBestSplitPosition(String<Pair<TScore,int> > & maxColsL,
					  String<Pair<TScore,int> > & maxColsR,
					  int & rowPosL1,
					  int rowPosL2,
					  int & rowPosR1,
					  int rowPosR2,
					  int seq0Len,
					  OrientationReverse,
					  SwiftSemiGlobal)
{
#ifdef RAZERS_DEBUG
	::std::cout << "findBestSplitEditReverse\n";
#endif

	TScore maxSum = std::numeric_limits<TScore>::min();
	int bestL = rowPosL2;
	int bestR = rowPosR2;

	while (rowPosL1 <= rowPosL2 && rowPosR1 >= rowPosR2)
	{
		if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 >= maxSum
			&& (maxColsL[rowPosL1].i2 + maxColsR[rowPosR1].i2 <= seq0Len))
		{
			maxSum = maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1;
			bestL = rowPosL1;
			bestR = rowPosR1;
		}
		++rowPosL1;
		--rowPosR1;
	}
	rowPosL1 = bestL;
	rowPosR1 = bestR;

	return true;
}






template <typename TMatch, typename TRead, typename TGenome, typename TSpec>
bool
combineLeftRight(TMatch & mR,
		 TMatch & mL,
		 TRead & read,
		 TGenome & genome,
		 RazerSOptions<TSpec> &options,
		 char orientation,
		 SwiftSemiGlobal)
{

#ifdef RAZERS_DEBUG
	::std::cout << "combinLeftRightEdit\n";
#endif	
	Score<int> scoreType(0,-1,-1,-1);
	typedef typename Infix<TGenome>::Type TGenomeInf;
	typedef ModifiedString<TGenomeInf,ModReverse> TGenomeInfRev;

	TGenomeInf genomeInfL = infix(genome, mL.gBegin+mL.gSeedLen, mL.gEnd);
	TGenomeInf readInfL = infix(read,options.minMatchLen,mL.mScore);
	TGenomeInfRev genomeInfR(infix(genome, mR.gBegin, mR.gEnd-mR.gSeedLen));
	TGenomeInfRev readInfR(infix(read,length(read)-mR.mScore,length(read)-options.minMatchLen));
	
#ifdef RAZERS_DEBUG
	bool debug = true;
	std::cout << "incombineleftright\nmL.mScore =" << mL.mScore << "\n";
	std::cout << "mL.gBegin =" << mL.gBegin << "\n";
	std::cout << "mL.gEnd =" << mL.gEnd << "\n";
	std::cout << "mL.editDist =" << mL.editDist << "\n";
				
	std::cout << "mR.mScore =" << mR.mScore << "\n";
	std::cout << "mR.gBegin =" << mR.gBegin << "\n";
	std::cout << "mR.gEnd =" << mR.gEnd << "\n";
	std::cout << "mR.editDist =" << mR.editDist << "\n";
#endif

	int readLength = length(read);
	unsigned maxErrors = readLength * options.errorRate;

	// neither insertion nor deletion
	if(mR.gEnd - mL.gBegin == readLength)
	{
#ifdef RAZERS_DEBUG
		::std::cout << "normal\n";
#endif
		if(true)return false;
		unsigned halfReadLen = readLength >> 1;

		Align<String<Dna5>, ArrayGaps> align;
		assignSource(row(align, 0), infix(read,0,readLength));
		assignSource(row(align, 1), infix(genome, mL.gBegin, mR.gEnd));
		int sc = globalAlignment(align, scoreType, AlignConfig<false,false,false,false>(), NeedlemanWunsch());
		if(-sc > (int)maxErrors) return false;
		
		mL.gEnd = mL.gBegin + toSourcePosition(row(align, 1),toViewPosition(row(align, 0), halfReadLen-1));
		mR.gBegin = mL.gEnd;
		mL.mScore = halfReadLen;
		mR.mScore = readLength - mL.mScore;
		mL.editDist = countErrorsInAlign(align,toViewPosition(row(align, 0), halfReadLen));
		mR.editDist = -sc-mL.editDist;

#ifdef RAZERS_DEBUG
		std::cout << "mL.mScore =" << mL.mScore << "\n";
		std::cout << "mL.gBegin =" << mL.gBegin << "\n";
		std::cout << "mL.gEnd =" << mL.gEnd << "\n";
		std::cout << "mL.editDist =" << mL.editDist << "\n";
	
		std::cout << "mR.mScore =" << mR.mScore << "\n";
		std::cout << "mR.gBegin =" << mR.gBegin << "\n";
		std::cout << "mR.gEnd =" << mR.gEnd << "\n";
		std::cout << "mR.editDist =" << mR.editDist << "\n";
#endif	

	}

	//potential insertion
	if(mR.gEnd - mL.gBegin < readLength)
	{
#ifdef RAZERS_DEBUG
		::std::cout << "insertion\n";
#endif
		if(mR.gEnd - mL.gBegin < 2*options.minMatchLen)  //too close together // actually minus allowed seed errors
			return false; 

		if(mL.gEnd < mR.gBegin)  //prefix and suffix match do not meet
			return false;

		if(mR.mScore + mL.mScore == mR.gEnd - mL.gBegin)
			if(mL.gEnd == mR.gBegin)
				if(mR.editDist + mL.editDist <= maxErrors ) //prefix and suffix match meet and do not overlap --> perfect
					return true;

//		if((mR.gEnd - mL.gBegin <= mL.mScore) || (mR.gEnd - mL.gBegin <= mR.mScore))//too close together // heuristic way to kick out repeat matches early on
//			return false;

		int diag1L = -maxErrors + mL.seedEditDist;
		int diag2L = maxErrors - mL.seedEditDist;
		int diag1R = -maxErrors + mR.seedEditDist;
		int diag2R = maxErrors - mR.seedEditDist;
 		int minColNum = 0;

		// genomeInf is the shorter sequence --> find best split position on genome
		// rows in alignment matrix represent genome position
		StringSet<TGenomeInf,Dependent<> > strSetL;
		appendValue(strSetL,readInfL);
		appendValue(strSetL,genomeInfL);
		Graph<Alignment<StringSet<TGenomeInf,Dependent<> >, void> > alignL(strSetL);
		String<Pair<int,int> > maxColsL;
		_globalAlignment(alignL,strSetL,scoreType,AlignConfig<false,false,false,false>(),diag1L,diag2L,maxColsL,minColNum,BandedNeedlemanWunsch());
	
		StringSet<TGenomeInfRev,Dependent<> > strSetR;
		appendValue(strSetR,readInfR);
		appendValue(strSetR,genomeInfR);
		Graph<Alignment<StringSet<TGenomeInfRev,Dependent<> >, void > > alignR(strSetR);
		String<Pair<int,int> > maxColsR;
		_globalAlignment(alignR,strSetR,scoreType,AlignConfig<false,false,false,false>(),diag1R,diag2R,maxColsR,minColNum,BandedNeedlemanWunsch());

		//::std::cout << alignL;
		//::std::cout << alignR;
	
		// our begin and start positions are defined by the read positions
		// go from read source to view position to genome source position
		int rowPosL1 = 0;
		if (mL.gSeedLen < (int)mR.gBegin - mL.gBegin) rowPosL1 = mR.gBegin - mL.gBegin - mL.gSeedLen;//or from first possible overlap pos
		int rowPosR2 = 0;
		if (mR.gSeedLen < (int)mR.gEnd - mL.gEnd) rowPosR2 = mR.gEnd - mL.gEnd - mR.gSeedLen;
		
		int rowPosL2 = mR.gEnd - mR.gSeedLen - rowPosR2 - mL.gBegin - mL.gSeedLen;		
		int rowPosR1 = mR.gEnd - mR.gSeedLen - rowPosL1 - mL.gBegin - mL.gSeedLen;


#ifdef RAZERS_DEBUG
		::std::cout << "vorher\nrowPosL1=" << rowPosL1 << ::std::endl;		
		::std::cout << "rowPosL2=" << rowPosL2 << ::std::endl;		
		::std::cout << "rowPosR1=" << rowPosR1 << ::std::endl;		
		::std::cout << "rowPosR2=" << rowPosR2 << ::std::endl;		
#endif

		// compute best split position
		int seq0Len = readLength;
		if(orientation == 'R')
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,OrientationReverse(),SwiftSemiGlobal());
		else
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,OrientationForward(),SwiftSemiGlobal());

		mL.editDist = mR.seedEditDist - maxColsL[rowPosL1].i1;	//scores are negative
		mR.editDist = mL.seedEditDist - maxColsR[rowPosR1].i1;
		if(mR.editDist + mL.editDist > maxErrors )
			return false;

#ifdef RAZERS_DEBUG
		::std::cout << "nachher\nrowPosL1=" << rowPosL1 << ::std::endl;		
		::std::cout << "rowPosL2=" << rowPosL2 << ::std::endl;		
		::std::cout << "rowPosR1=" << rowPosR1 << ::std::endl;		
		::std::cout << "rowPosR2=" << rowPosR2 << ::std::endl;		
#endif

		// best split position in genome is saved in rowPosL1 (and rowPosR1)
		mL.mScore = options.minMatchLen + maxColsL[rowPosL1].i2;
		mL.gEnd = mL.gBegin + mL.gSeedLen + rowPosL1; //read position of best genomic split

		mR.mScore = options.minMatchLen + maxColsR[rowPosR1].i2;
		mR.gBegin = mR.gEnd - mR.gSeedLen - rowPosR1; //read position of best genomic split

#ifdef RAZERS_DEBUG
		if(mL.editDist > 50 || mR.mScore < options.minMatchLen || mL.mScore < options.minMatchLen) 
		{
			::std::cout << "-maxColsL[rowPosL1].i1=" << -maxColsL[rowPosL1].i1 << " -maxColsL[rowPosL1].i2=" << -maxColsL[rowPosL1].i2 << " rowPosL1="<<rowPosL1;
			::std::cout << " maxColsLLen="<< length(maxColsL) << ::std::endl;
			::std::cout << "-maxColsR[rowPosR1].i1=" << -maxColsR[rowPosR1].i1 << "-maxColsR[rowPosR1].i2=" << -maxColsR[rowPosR1].i2 << " rowPosR1="<<rowPosR1;
			::std::cout << " maxColsRLen="<< length(maxColsR) << ::std::endl;
		}
#endif
	}

	//potential deletion
	if(mR.gEnd - mL.gBegin > readLength)
	{
#ifdef RAZERS_DEBUG
		::std::cout << "deletion\n";
#endif
		if(mR.mScore + mL.mScore < readLength)
			return false;
		if(mR.mScore + mL.mScore == readLength && mR.editDist + mL.editDist > maxErrors) //the prefix and suffix match meet, but too many errors
			return false;
		
		int diag1L = -maxErrors + mL.seedEditDist;
		int diag2L = maxErrors - mL.seedEditDist;
		int diag1R = -maxErrors + mR.seedEditDist;
		int diag2R = maxErrors - mR.seedEditDist;
		int minColNum = 0;
		
		// readInf is the shorter sequence --> find best split position on read
		// rows in alignment matrix represent read position
		StringSet<TGenomeInf> strL;
		appendValue(strL,genomeInfL);
		appendValue(strL,readInfL);
		String<Pair<int,int> > maxColsL;
		_globalAlignment(strL,scoreType,AlignConfig<false,false,false,false>(),diag1L,diag2L,maxColsL,minColNum,BandedNeedlemanWunsch());
	
		StringSet<TGenomeInfRev> strR;
		appendValue(strR,genomeInfR);
		appendValue(strR,readInfR);
		String<Pair<int,int> > maxColsR;
		_globalAlignment(strR,scoreType,AlignConfig<false,false,false,false>(),diag1R,diag2R,maxColsR,minColNum,BandedNeedlemanWunsch());
	
		int rowPosL1 = ((int)options.minMatchLen > readLength-mR.mScore) ? (int)0 : readLength-mR.mScore-options.minMatchLen;
		int rowPosL2 = (readLength-(int)options.minMatchLen < mL.mScore) ? readLength-(int)2*options.minMatchLen : mL.mScore - options.minMatchLen;
		int rowPosR1 = (int)readLength - 2*options.minMatchLen - rowPosL1;
		int rowPosR2 = (int)readLength - 2*options.minMatchLen - rowPosL2;
		
#ifdef RAZERS_DEBUG
		::std::cout << "rowPosL1=" << rowPosL1 << ::std::endl;		
		::std::cout << "rowPosL2=" << rowPosL2 << ::std::endl;		
		::std::cout << "rowPosR1=" << rowPosR1 << ::std::endl;		
		::std::cout << "rowPosR2=" << rowPosR2 << ::std::endl;	

		std::cout << "before split:\nmL.mScore =" << mL.mScore << "\n";
		std::cout << "mL.gBegin =" << length(genome)-mL.gBegin << "\n";
		std::cout << "mL.gEnd =" << length(genome)-mL.gEnd << "\n";
		std::cout << "mL.editDist =" << mL.editDist << "\n";

		std::cout << "mR.mScore =" << mR.mScore << "\n";
		std::cout << "mR.gBegin =" << length(genome)-mR.gBegin << "\n";
		std::cout << "mR.gEnd =" << length(genome)-mR.gEnd << "\n";
		std::cout << "mR.editDist =" << mR.editDist << "\n";
#endif

		// compute best split position
		int seq0Len = mR.gEnd - mL.gBegin;
		if(orientation == 'R')
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,OrientationReverse(),SwiftSemiGlobal());
		else
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,OrientationForward(),SwiftSemiGlobal());


		// best split position in read is saved in rowPosL1 (and rowPosR1)
		mL.editDist = mL.seedEditDist - maxColsL[rowPosL1].i1;
		mR.editDist = mR.seedEditDist - maxColsR[rowPosR1].i1;


		if(mR.editDist + mL.editDist > maxErrors)
			return false;

		mR.mScore = options.minMatchLen + rowPosR1;
		mR.gBegin = mR.gEnd - maxColsR[rowPosR1].i2 - mR.gSeedLen; //genomic position of best read split
		mL.mScore = options.minMatchLen + rowPosL1;
		mL.gEnd = mL.gBegin + maxColsL[rowPosL1].i2 + mL.gSeedLen; //genomic position of best read split

	}

	return true;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hamming distance match combination

template <typename TMatch, typename TOriRead, typename TGenome, typename TSpec>
bool
combineLeftRight(TMatch & mR,
		 TMatch & mL,
		 TOriRead const& oriRead,
		 TGenome & genome,
		 RazerSOptions<TSpec> &options,
		 char orientation,
		 SwiftSemiGlobalHamming)
{

#ifdef RAZERS_DEBUG
	::std::cout << "combineLeftRightHamming\n";
#endif
	

	typedef typename Infix<TOriRead const>::Type TRead;
	TRead read = infix(oriRead,0,length(oriRead));
	
	typedef typename Infix<TGenome>::Type TGenomeInf;
	TGenomeInf genomeInf = infix(genome, mL.gBegin, mR.gEnd);
	int readLength = length(read);
	unsigned maxErrors = readLength * options.errorRate;


#ifdef RAZERS_DEBUG
	std::cout << "readLength=" << readLength << "\n";
	std::cout << "sumLen=" << mL.mScore + mR.mScore << "\n";
	std::cout << "gInfLength=" << length(genomeInf) << "\n";
	std::cout << "gInf=" << genomeInf << "\n";
	
	std::cout << "incombineleftright\nmL.mScore =" << mL.mScore << "\n";
	std::cout << "mL.gBegin =" << mL.gBegin << "\n";
	std::cout << "mL.gEnd =" << mL.gEnd << "\n";
	std::cout << "mL.editDist =" << mL.editDist << "\n";
	std::cout << "readPref=" << prefix(oriRead,mL.mScore) << "\n";

	std::cout << "mR.mScore =" << mR.mScore << "\n";
	std::cout << "mR.gBegin =" << mR.gBegin << "\n";
	std::cout << "mR.gEnd =" << mR.gEnd << "\n";
	std::cout << "mR.editDist =" << mR.editDist << "\n";
	std::cout << "readSuff=" << suffix(oriRead,length(oriRead)-mR.mScore) << "\n";
#endif


	// neither insertion nor deletion
	if(mR.gEnd - mL.gBegin == readLength)
	{
#ifdef RAZERS_DEBUG
		::std::cout << "normal\n";
#endif
		unsigned halfReadLen = readLength >> 1;
		unsigned prefixErrors = 0;
		unsigned suffixErrors = 0;
		for (unsigned i = 0 ; i < length(genomeInf); ++i)
		{
			if ((options.compMask[ordValue(read[i])] & options.compMask[ordValue(genomeInf[i])]) == 0)
			{
				if(i >= halfReadLen)
					++suffixErrors;
				else
					++prefixErrors;
			}
		}
	
		if (suffixErrors+prefixErrors <= maxErrors)
		{
			mL.mScore = halfReadLen;
			mR.mScore = length(read)- halfReadLen;
			mR.gBegin = mR.gEnd - mR.mScore;
			mL.gEnd = mL.gBegin + mL.mScore;
			mL.editDist = prefixErrors;
			mR.editDist = suffixErrors;
		}
		else return false;	
		
	}
	//potential insertion
	if(mR.gEnd - mL.gBegin < readLength)
	{
#ifdef RAZERS_DEBUG
		::std::cout << "insertion\n";
#endif

		if(mR.gEnd - mL.gBegin < 2*options.minMatchLen)//too close together
			 return false; 

		if(mR.mScore + mL.mScore < mR.gEnd - mL.gBegin) //prefix and suffix match do not meet
			return false;

		if((mR.mScore + mL.mScore == mR.gEnd - mL.gBegin) && (mR.editDist + mL.editDist > maxErrors)) //prefix and suffix match meet but too many errors
			return false;
		
//		if((mR.gEnd - mL.gBegin <= mL.mScore) || (mR.gEnd - mL.gBegin <= mR.mScore))//too close together 
//			 return false; 
			 
		bool result = findBestSplitPosition(read,genomeInf,mL.mScore,mR.mScore,mL.editDist,mR.editDist, options, orientation, SwiftSemiGlobalHamming());
		if(!result || mR.editDist + mL.editDist > maxErrors) return false;
		mR.gBegin = mR.gEnd - mR.mScore;
		mL.gEnd = mL.gBegin + mL.mScore;

	}
	//potential deletion 	
	if(mR.gEnd - mL.gBegin > readLength)
	{
#ifdef RAZERS_DEBUG
		::std::cout << "deletion\n";
#endif

		if(mR.mScore + mL.mScore < readLength) 
			return false;

		if((mR.mScore + mL.mScore == readLength) && (mR.editDist + mL.editDist > maxErrors)) //the prefix and suffix match meet and do not overlap --> perfect
			return false;

		bool result = findBestSplitPosition(genomeInf,read,mL.mScore,mR.mScore,mL.editDist,mR.editDist, options, orientation,SwiftSemiGlobalHamming());

		if(!result || mR.editDist + mL.editDist > maxErrors) return false;

		mR.gBegin = mR.gEnd - mR.mScore;
		mL.gEnd = mL.gBegin + mL.mScore;
		
	}
#ifdef RAZERS_DEBUG
	std::cout << "incombineleftright\nmL.mScore =" << mL.mScore << "\n";
	std::cout << "mL.gBegin =" << mL.gBegin << "\n";
	std::cout << "mL.gEnd =" << mL.gEnd << "\n";
	std::cout << "mL.editDist =" << mL.editDist << "\n";
				
	std::cout << "mR.mScore =" << mR.mScore << "\n";
	std::cout << "mR.gBegin =" << mR.gBegin << "\n";
	std::cout << "mR.gEnd =" << mR.gEnd << "\n";
	std::cout << "mR.editDist =" << mR.editDist << "\n";
#endif


	return true;
}


// find the best split position for a split match
// positions are relative to shorter sequence 
// (deletion --> read is shorter, insertion --> genomeInf is shorter)
template <typename TSize, typename TLongerSegment, typename TShorterSegment, typename TErrors, typename TOptions>
bool
findBestSplitPosition(TLongerSegment &longSeg, 
			TShorterSegment &shortSeg, 
			TSize & mLmScore, 
			TSize & mRmScore, 
			TErrors & errorsLeft, 
			TErrors & errorsRight, 
			TOptions & options, 
			char orientation, 
			SwiftSemiGlobalHamming)
{
	
#ifdef RAZERS_DEBUG
	::std::cout << "findBestSplitHamming"<<orientation<<"\n";
#endif

	// usually, both types should be the same, but you never know...
	typedef typename Iterator<TLongerSegment const>::Type TLongIterator;
	typedef typename Iterator<TShorterSegment const>::Type TShortIterator;
	typedef typename Size<TShorterSegment const>::Type TShortSize;
	typedef typename Size<TLongerSegment const>::Type TLongSize;


	TShortSize shortLen = length(shortSeg);
	TLongSize  longLen  = length(longSeg);
	TLongSize  lenDiff  = longLen - shortLen;

	TLongSize leftLongBegin = _max((int)options.minMatchLen,(int)shortLen-mRmScore);
	TLongSize leftLongEnd = _min((int)mLmScore,(int)shortLen-(int)options.minMatchLen);
	TLongSize leftLongPos = leftLongBegin;
	
	TLongSize rightLongPos = leftLongBegin + lenDiff;
	TShortSize shortPos = leftLongBegin;

	int bestSumErrors = 0;
	int bestLErrors = 0;
	int bestRErrors = 0;
	int bestPos = shortPos;
	int errorsL = 0;
	int errorsR = 0;
	int errorsPosL = 0;
	int errorsPosR = 0;
	
#ifdef RAZERS_DEBUG
	std::cout << "before\nerrorsLeft= " << errorsLeft << std::endl;
	std::cout << "errorsRight= " << errorsRight << std::endl;
	std::cout << "mLmScore= " << mLmScore << std::endl;
	std::cout << "mRmScore= " << mRmScore << std::endl;

	std::cout << "leftLongPos " << leftLongPos << std::endl;
	std::cout << "leftLongEnd " << leftLongEnd << std::endl;
	std::cout << "rightLongPos " << rightLongPos << std::endl;
	std::cout << "shortPos " << shortPos << std::endl;

	std::cout << longSeg << "=longSeg\n";
	std::cout << shortSeg << "=shortSeg\n";
#endif

	//find best split position --> min. sum errors
	while(leftLongPos < leftLongEnd)
	{
		if((options.compMask[ordValue(shortSeg[shortPos])] & options.compMask[ordValue(longSeg[leftLongPos])]) == 0)
		{
			++errorsL; 
			++errorsPosL;
		}
		if((options.compMask[ordValue(shortSeg[shortPos])] & options.compMask[ordValue(longSeg[rightLongPos])]) == 0)
		{
			--errorsPosR;
			++errorsR;
		}
		if(errorsPosL+errorsPosR < bestSumErrors 
			|| (orientation == 'R' && errorsPosL+errorsPosR == bestSumErrors))
		{
			bestSumErrors = errorsPosL+errorsPosR;
			bestLErrors = errorsPosL;
			bestRErrors = errorsPosR;
			bestPos = shortPos + 1;
		}
		++leftLongPos;
		++rightLongPos;
		++shortPos;
	}
	
	//update to new match lengths
	mLmScore = bestPos;
	mRmScore = shortLen - bestPos; 
	
	//count numErrors for left and right match   
	//(have to count completely new, because mScore may be longer than shortLen --> no able to track errors outside segment)
	errorsRight = errorsLeft = 0;
	for(leftLongPos = 0, shortPos = 0; leftLongPos < (unsigned)mLmScore; ++leftLongPos, ++shortPos)
	{
		if((options.compMask[ordValue(shortSeg[shortPos])] & options.compMask[ordValue(longSeg[leftLongPos])]) == (unsigned) 0)
			++errorsLeft; 
		
	}
	for(rightLongPos = 0, shortPos = 0; rightLongPos < (unsigned)mRmScore; ++rightLongPos, ++shortPos)
	{
		if((options.compMask[ordValue(shortSeg[shortLen-1-shortPos])] & options.compMask[ordValue(longSeg[longLen-1-rightLongPos])]) ==  (unsigned) 0)
			++errorsRight; 
		
	}
	
#ifdef RAZERS_DEBUG
	std::cout << "bestSumErrors= " << bestSumErrors << std::endl;
	std::cout << "errorsPosR= " << errorsPosR << std::endl;
	std::cout << "errorsPosL= " << errorsPosL << std::endl;
	std::cout << "after\nerrorsLeft= " << errorsLeft << std::endl;
	std::cout << "errorsRight= " << errorsRight << std::endl;
	std::cout << "mLmScore= " << mLmScore << std::endl;
	std::cout << "mRmScore= " << mRmScore << std::endl;
#endif

	return true;
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TMatches, 
	typename TReadSet_, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSplicedReads(
	TMatches &		matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet_ & 		readSet,
	TCounts &		cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &		shape,
	Swift<TSwiftSpec> const)
{


	typedef typename Value<TReadSet_>::Type								TRead;
	typedef StringSet<typename Infix<TRead>::Type>						TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape, TQGramIndexSpec> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >							TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>								TMyersPattern;	// verifier
	
	
	// split reads over two indices, one for prefixes, the other for suffixes
	TReadSet readSetL, readSetR;
	unsigned readCount = length(readSet);
	resize(readSetL, readCount);
	resize(readSetR, readCount);

	if(options._debugLevel > 0)
	{
		int64_t genomeLen = 3000000000 * 2;					// ufff make that an option
		expNumRandomMatches(readSet, genomeLen, options);
	}
	
	if(options._debugLevel > 0 ) 
		std::cout << "Performing spliced mapping.\n";
	for (unsigned i = 0; i < readCount; ++i)
	{
		if(length(readSet[i])>=2*options.minMatchLen)
		{
			assign(readSetL[i], infix(readSet[i],0,options.minMatchLen), Exact());
			assign(readSetR[i], infix(readSet[i],length(readSet[i])-options.minMatchLen,length(readSet[i])), Exact());
		}
		else
		{
			assign(readSetL[i], infix(readSet[i],0,0), Exact());
			assign(readSetR[i], infix(readSet[i],0,0), Exact());
		}
	}
	
	
	if(options._debugLevel > 1)::std::cout << "Make index left right\n";
	// configure q-gram index
	TIndex swiftIndexL(readSetL, shape);
	TIndex swiftIndexR(readSetR, shape);
	
#ifdef RAZERS_OPENADDRESSING
	swiftIndexL.alpha = 2;
	swiftIndexR.alpha = 2;
#endif
	
	cargo(swiftIndexL).abundanceCut = options.abundanceCut;
	cargo(swiftIndexR).abundanceCut = options.abundanceCut;
	cargo(swiftIndexL)._debugLevel = 0;
	cargo(swiftIndexR)._debugLevel = options._debugLevel;
	
	// configure Swift
	TSwiftPattern swiftPatternL(swiftIndexL);
	TSwiftPattern swiftPatternR(swiftIndexR);
	swiftPatternL.params.minThreshold = options.threshold;
	swiftPatternR.params.minThreshold = options.threshold;
	swiftPatternL.params.tabooLength = options.tabooLength;
	swiftPatternR.params.tabooLength = options.tabooLength;
	
	// init edit distance verifiers
	String<TMyersPattern> forwardPatternsL;
	String<TMyersPattern> forwardPatternsR;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		resize(forwardPatternsL, readCount, Exact());
		resize(forwardPatternsR, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatternsL[i], indexText(swiftIndexL)[i]);
			setHost(forwardPatternsR[i], indexText(swiftIndexR)[i]);
			_patternMatchNOfPattern(forwardPatternsL[i], options.matchN);
			_patternMatchNOfPattern(forwardPatternsR[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsL[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsR[i], options.matchN);
		}
	}
	
	if(options._debugLevel > 1)::std::cout << "Patterns created\n";

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
	
	unsigned filecount = 0;
	unsigned numFiles = length(genomeFileNameList);
	unsigned gseqNo = 0;
	
	// open genome files, one by one	
	while (filecount < numFiles)
	{
		// open genome file	
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;
	
		// remove the directory prefix of current genome file
		::std::string genomeFile(toCString(genomeFileNameList[filecount]));
		size_t lastPos = genomeFile.find_last_of('/') + 1;
		if (lastPos == genomeFile.npos) lastPos = genomeFile.find_last_of('\\') + 1;
		if (lastPos == genomeFile.npos) lastPos = 0;
		::std::string genomeName = genomeFile.substr(lastPos);
		
		CharString	id;
		//Dna5String	genome;
		String<Dna5Q> genome;
		unsigned gseqNoWithinFile = 0;
		SEQAN_PROTIMESTART(find_time);

		// iterate over genome sequences
		for(; !_streamEOF(file); ++gseqNo)
		{
			if (options.genomeNaming == 0)
			{
				readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
				appendValue(genomeNames, id, Generous());
			}
			read(file, genome, Fasta());			// read Fasta sequence
			
			gnoToFileMap.insert(::std::make_pair(gseqNo,::std::make_pair(genomeName,gseqNoWithinFile)));
			
			if (options.forward)
				mapSplicedReads(matches, genome, gseqNo, readSet, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'F', options);
	
			if (options.reverse)
			{
				reverseComplement(genome);
				mapSplicedReads(matches, genome, gseqNo, readSet, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
		file.close();
		++filecount;
	}

	compactSplicedMatches(matches, cnts, options, false, swiftPatternL, swiftPatternR);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Swift FP: " << options.FP << ::std::endl;
		::std::cerr << "Swift TP: " << options.TP << ::std::endl;
	}
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TOriReadSet,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TCounts,
	typename TSpec >
void mapSplicedReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	TOriReadSet &readSet,			// reads
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternL,	// left index
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternR,	// right index
	TVerifier &forwardPatternsL,
	TVerifier &forwardPatternsR,
	TCounts & cnts,
	char orientation,
	RazerSOptions<TSpec> &options)
{
	typedef typename Value<TOriReadSet>::Type TRead;
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	typedef typename Size<TGenome>::Type			TSize;
	typedef typename Position<TGenome>::Type		TGPos;
	typedef typename MakeSigned_<TGPos>::Type		TSignedGPos;
	typedef typename Value<TMatches>::Type			TMatch;
	typedef typename Infix<TGenome>::Type			TGenomeInf;
	
	// Prefix-Suffix filtration
	typedef Finder<TGenome, Swift<TSwiftSpec> >		TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >	TSwiftFinderR;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >	TSwiftPattern;
	
	typedef Pair<int64_t, TMatch>				TDequeueValue;
	typedef Dequeue<TDequeueValue>				TDequeue;
	typedef typename TDequeue::TIter			TDequeueIterator;
	
	const unsigned NOT_VERIFIED = 1u << (8*sizeof(unsigned)-1);
	
	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}
	
	TReadSet &readSetL = host(host(swiftPatternL));
	TReadSet &readSetR = host(host(swiftPatternR));
	
	if (empty(readSetL) || empty(readSetR))
		return;
	
	
	// Check?
	TSignedGPos maxDistance = options.maxDistance;// + 2 * (int)options.minMatchLength - (int)length(indexShape(host(swiftPatternL)));
	TSignedGPos minDistance = options.minDistance + 2*options.minMatchLen;

	// exit if contig is shorter than minDistance
	if (length(genome) <= (unsigned) minDistance)
		return;
	
	TGenomeInf genomeInf = infix(genome, 0, length(genome));
	TSwiftFinderL swiftFinderL(genome, options.repeatLength, 1);
	TSwiftFinderR swiftFinderR(genomeInf, options.repeatLength, 1);
	
	TDequeue fifo;				// stores left-mate potential matches
	String<int64_t> lastPotMatchNo;		// last number of a left-mate potential
	int64_t lastNo = 0;			// last number over all left-mate pot. matches in the queue
	int64_t firstNo = 0;			// first number over all left-mate pot. match in the queue
	Pair<TGPos> gPair;
	
	resize(lastPotMatchNo, length(host(swiftPatternL)), (int64_t)-1, Exact());
	
	String<Pair<TGPos> > lastRightMatch;		// begin and end of last verified right match
	resize(lastRightMatch, length(host(swiftPatternL)), Pair<TGPos>(0,0), Exact());
	
	TSize gLength = length(genome);
	TMatch mR = {	// to supress uninitialized warnings
		0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0};
	TDequeueValue fL(-1, mR);	// to supress uninitialized warnings
	fL.i2.gseqNo = gseqNo;
	mR.gseqNo = gseqNo;
	fL.i2.orientation = orientation;
	mR.orientation = orientation;
	
	double maxErrorRate = options.errorRate;
	double prefixErrorRate = (double)floor(options.minMatchLen * options.errorRate)/options.minMatchLen;
	if(prefixErrorRate > maxErrorRate) maxErrorRate = prefixErrorRate;
	if(!empty(readSet))
	{
		double extPrefixErrorRate = (double) floor(length(readSet[0]) * options.errorRate)/(length(readSet[0]) - options.minMatchLen);
		if(extPrefixErrorRate > maxErrorRate) maxErrorRate = extPrefixErrorRate;
	}
	Pair<TGPos,TGPos> lastLeftMatch(0,0);
	
	// iterate all verification regions returned by SWIFT
	while (find(swiftFinderR, swiftPatternR, maxErrorRate)) 
	{
		unsigned rseqNo = swiftPatternR.curSeqNo;
		TGPos rEndPos = endPosition(swiftFinderR);	
		TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;
		
		// remove out-of-window left mates from fifo
		while (!empty(fifo) && front(fifo).i2.gBegin + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
		{
			popFront(fifo);
			++firstNo;
		}
		
		// add within-window left mates to fifo
		while (empty(fifo) || back(fifo).i2.gBegin + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
		{
			if (find(swiftFinderL, swiftPatternL, maxErrorRate)) 
			{
				gPair = positionRange(swiftFinderL);
				if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
				{
					// link in
					fL.i1 = lastPotMatchNo[swiftPatternL.curSeqNo]; //link to last previous potential match
					lastPotMatchNo[swiftPatternL.curSeqNo] = lastNo++; //++ general counter and remember last pot match of curSeqNo
					
					fL.i2.rseqNo = swiftPatternL.curSeqNo | NOT_VERIFIED; //das erste bit wird gestetzt
					fL.i2.gBegin = gPair.i1;	//begin und end auf die potential match region
					fL.i2.gEnd = gPair.i2;
					pushBack(fifo, fL);
				}
			} 
			else
				break;
		}
		
		TDequeueIterator it;
		int64_t lastPositive = (int64_t)-1;

		TRead const &read = readSet[rseqNo];
		TSize counter = 0;
		bool noMatchRight = false;
		bool notYetVerifiedRight = true;
		lastLeftMatch.i1 = 0;
		lastLeftMatch.i2 = 0;
		
		// walk through all potential left matches, verify if not verfied already, mark as positive or negative
		for (int64_t i = lastPotMatchNo[rseqNo]; firstNo <= i; i = (*it).i1)
		{
			it = &value(fifo, i - firstNo);
		
			// verify left mate (equal seqNo), if not done already
			if ((*it).i2.rseqNo & NOT_VERIFIED)
			{
				TGPos maxEndPos = (*it).i2.gEnd + length(read)-2*options.minMatchLen + floor(options.errorRate*length(read));
				if(maxEndPos > length(genome)) maxEndPos = length(genome);
				if (matchVerify( (*it).i2,
						infix(genome, (*it).i2.gBegin, maxEndPos),
						rseqNo, 
						readSet, //readSetL 
						forwardPatternsL, 
						options, 
						TSwiftSpec(),
						LongestPrefix()) && 
						!(lastLeftMatch.i1 == (TGPos)(*it).i2.gBegin && lastLeftMatch.i2 == (TGPos)(*it).i2.gEnd ))
				{
					(*it).i2.rseqNo &= ~NOT_VERIFIED; // has been verified positively // go back to regular rseqNo
					lastLeftMatch.i1 = (*it).i2.gBegin;
					lastLeftMatch.i2 = (*it).i2.gEnd;
					// short-cut negative matches
					if (lastPositive == (int64_t)-1)
						lastPotMatchNo[rseqNo] = i;
					else
						value(fifo, lastPositive - firstNo).i1 = i;
					lastPositive = i;
				} 
				else
				{
					(*it).i2.rseqNo = ~NOT_VERIFIED;		// has been verified negatively // 01111111....111
				}
			}
			
			// verify right mate, if left match was found
			if ((*it).i2.rseqNo == rseqNo)
			{
				lastLeftMatch.i1 = (*it).i2.gBegin;
				lastLeftMatch.i2 = (*it).i2.gEnd;

				// dont want to shortcut too much
				if (lastPositive == (int64_t)-1 || i < lastPositive)
					lastPositive = i;

				// first left match --> do right verification
				if(notYetVerifiedRight)
				{
					notYetVerifiedRight = false;
					TGPos minBeginPos = beginPosition(infix(swiftFinderR, genomeInf));
					if((int)minBeginPos - length(read)+2*options.minMatchLen - floor(options.errorRate*length(read)) > 0)
						minBeginPos = minBeginPos - length(read)+2*options.minMatchLen - floor(options.errorRate*length(read));
					else minBeginPos = 0;
					if (!matchVerify(mR, 
						infix(genome,minBeginPos,endPosition(infix(swiftFinderR, genomeInf))),
						rseqNo, 
						readSet,//readSetR, 
						forwardPatternsR,
						options, 
						TSwiftSpec(),
						LongestSuffix()))
					{
						noMatchRight = true;
						continue;
					}
					else 
					{
						if (lastRightMatch[rseqNo].i1 == (TGPos)mR.gBegin && lastRightMatch[rseqNo].i2 == (TGPos)mR.gEnd)
							noMatchRight = true;
						lastRightMatch[rseqNo].i1 = mR.gBegin;
						lastRightMatch[rseqNo].i2 = mR.gEnd;
					}
				}

				//else check if left and right match fit together
				if(!noMatchRight)
				{ 
#ifdef RAZERS_DEBUG
					std::cout << "before\nmL.mScore =" << (*it).i2.mScore << "\n";
					std::cout << "mL.gBegin =" << (*it).i2.gBegin << "\n";
					std::cout << "mL.gEnd =" << (*it).i2.gEnd << "\n";
					std::cout << "mL.editDist =" << (*it).i2.editDist << "\n";
				
					std::cout << "mR.mScore =" << mR.mScore << "\n";
					std::cout << "mR.gBegin =" << mR.gBegin << "\n";
					std::cout << "mR.gEnd =" << mR.gEnd << "\n";
					std::cout << "mR.editDist =" << mR.editDist << "\n";
#endif
					int outerDistance = mR.gEnd - (*it).i2.gBegin;
					if (outerDistance < (int)(2 * options.minMatchLen))
						continue;
//					::std::cout << options.minMatchLen << "<-minMatchLen  outerDistance->" << outerDistance << std::endl;
					int outerDistanceError = length(readSet[rseqNo]) - (int)(outerDistance);
					if (outerDistanceError < 0) outerDistanceError = -outerDistanceError;
					if ((outerDistanceError > (int) options.maxDistance) ||
						(outerDistanceError < (int) options.minDistance))
						continue;

//					::std::cout <<"hier!"<<options.minDistance <<" "<<options.maxDistance<<"\n";

					TMatch mRtmp = mR;
					TMatch mLtmp = (*it).i2;
					if(!combineLeftRight(mRtmp,mLtmp,read,genome,options,orientation,TSwiftSpec()))
					{
						++options.FP;
						continue;
					}
					else ++options.TP;

					++counter;

					//assign an alignment score
					//mRtmp.alignmentScore = (options.matchScore * (sumLength - sumDist)) + (options.mismatchScore * sumDist);
					//if(outerDistanceError != 0) mRtmp.alignmentScore += ((outerDistanceError - 1) * options.gapExtend) + options.gapOpen;
					//mLtmp.alignmentScore = mRtmp.alignmentScore;

					//mRtmp.alignmentScore = (options.matchScore * (sumLength - sumDist)) + (options.mismatchScore * sumDist);
					//if(outerDistanceError != 0)
					//{
					//	int maxErrors = (int)(options.errorRate*length(readSet[rseqNo]));
					//	double gapScoreConvergeTo = ((options.minMatchLen*2 - maxErrors) * options.matchScore) + (maxErrors * options.mismatchScore) ;
					//	double gapPenalty =	gapScoreConvergeTo*(1-exp((double)-options.indelLambda*outerDistanceError));
					//	mRtmp.alignmentScore -= (int)gapPenalty;
					//}
					//mLtmp.alignmentScore = mRtmp.alignmentScore;

					if (orientation == 'R') 
					{
						TSize temp = mLtmp.gBegin;
						mLtmp.gBegin = gLength - mLtmp.gEnd;
						mLtmp.gEnd = gLength - temp;
						temp = mRtmp.gBegin;
						mRtmp.gBegin = gLength - mRtmp.gEnd;
						mRtmp.gEnd = gLength - temp;
						outerDistance = -outerDistance;
					}
					// set a unique pair id
					mLtmp.pairId = mRtmp.pairId = options.nextMatePairId;
					if (++options.nextMatePairId == 0)
						options.nextMatePairId = 1;
						
					// score the whole match pair
					mLtmp.pairScore = mRtmp.pairScore = 0 - mLtmp.editDist - mRtmp.editDist;

#ifdef RAZERS_DEBUG_LIGHT
					bool falsch = false;
					if((unsigned)(mRtmp.mScore + mLtmp.mScore) > length(readSet[rseqNo])) falsch = true;
					if(mRtmp.mScore> length(readSet[rseqNo])-options.minMatchLen || mRtmp.mScore < options.minMatchLen) falsch = true;
					if((int)mRtmp.gEnd - (int)mRtmp.gBegin < 0 || (int)mRtmp.gEnd - (int)mRtmp.gBegin > 500) falsch = true;
					if(mRtmp.editDist > 200) falsch = true;

					if(mLtmp.mScore> length(readSet[rseqNo])-options.minMatchLen ||mLtmp.mScore < options.minMatchLen) falsch = true;
					if((int)mLtmp.gEnd - (int)mLtmp.gBegin < 0 || (int)mLtmp.gEnd - (int)mLtmp.gBegin > 500) falsch = true;
					if(mLtmp.editDist > 200) falsch = true;

					if(falsch)
					{
						::std::cout << "rseqNo="<<rseqNo <<::std::endl;
							std::cout << "mL.mScore =" << mLtmp.mScore << "\n";
							std::cout << "mL.gBegin =" << mLtmp.gBegin << "\n";
							std::cout << "mL.gEnd =" << mLtmp.gEnd << "\n";
							std::cout << "mL.editDist =" << mLtmp.editDist << "\n";
					
							std::cout << "mR.mScore =" << mRtmp.mScore << "\n";
							std::cout << "mR.gBegin =" << mRtmp.gBegin << "\n";
							std::cout << "mR.gEnd =" << mRtmp.gEnd << "\n";
							std::cout << "mR.editDist =" << mRtmp.editDist << "\n";
					
					
					
					if(-mLtmp.pairScore > (int) (options.errorRate * length(readSet[rseqNo])))
					{
						::std::cout << "mLtmp.pairScore = " << mLtmp.pairScore;
						::std::cout << "\treadLen = " << length(readSet[rseqNo]);
						::std::cout << "\trseqNo = " << rseqNo << ::std::endl;
					}
					std::cout << "mScoreSum =" << mLtmp.mScore + mRtmp.mScore << "\t";
					std::cout << "readLength =" << length(readSet[rseqNo]) << "\t";
					std::cout << "mL.mScore =" << mLtmp.mScore << "\t";
					std::cout << "mL.gBegin =" << mLtmp.gBegin << "\t";
					std::cout << "mL.gEnd =" << mLtmp.gEnd << "\t";
					std::cout << "mL.editDist =" << mLtmp.editDist << "\n";
					
					std::cout << "mR.mScore =" << mRtmp.mScore << "\t";
					std::cout << "mR.gBegin =" << mRtmp.gBegin << "\t";
					std::cout << "mR.gEnd =" << mRtmp.gEnd << "\t";
					std::cout << "mR.editDist =" << mRtmp.editDist << "\n";

					Align<String<Dna5Q>, ArrayGaps> alignL;
					Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
					if (options.hammingOnly)
						scoreType.data_mismatch = -1;
					resize(rows(alignL), 2);
					assignSource(row(alignL, 0), infix(readSet[rseqNo],0,mLtmp.mScore));
					assignSource(row(alignL, 1), infix(genome,mLtmp.gBegin, mLtmp.gEnd));
					if (orientation == 'R')
						reverseComplement(source(row(alignL, 1)));

					globalAlignment(alignL, scoreType);
					std::cout << alignL;

					Align<String<Dna5Q>, ArrayGaps> alignR;
					resize(rows(alignR), 2);
					assignSource(row(alignR, 0), infix(readSet[rseqNo],length(readSet[rseqNo])-mRtmp.mScore,length(readSet[rseqNo])));
					assignSource(row(alignR, 1), infix(genome,mRtmp.gBegin, mRtmp.gEnd));
					if (orientation == 'R')
						reverseComplement(source(row(alignR, 1)));

					globalAlignment(alignR, scoreType);
					std::cout << alignR;

					std::cout << "Komplettali:\n";
					Align<String<Dna5Q>, ArrayGaps> align;
					resize(rows(align), 2);
					assignSource(row(align, 0), infix(readSet[rseqNo],0,length(readSet[rseqNo])));
					assignSource(row(align, 1), infix(genome,_min(mRtmp.gBegin,mLtmp.gBegin), _max(mRtmp.gEnd,mLtmp.gEnd)));
					if (orientation == 'R')
						reverseComplement(source(row(align, 1)));

					globalAlignment(align, scoreType);
					std::cout << align;
					}
#endif

				

					// relative positions
					mLtmp.mateDelta = outerDistance;
					mRtmp.mateDelta = -outerDistance;
					
					mLtmp.rseqNo = mRtmp.rseqNo = rseqNo;
					
					if (!options.spec.DONT_DUMP_RESULTS)
					{
						appendValue(matches, mLtmp, Generous());
						appendValue(matches, mRtmp, Generous());
						if (length(matches) > options.compactThresh)
						{
							typename Size<TMatches>::Type oldSize = length(matches);
//							maskDuplicates(matches);	// overlapping parallelograms cause duplicates //TODO: implement!
							compactSplicedMatches(matches, cnts, options, false, swiftPatternL, swiftPatternR);
							options.compactThresh += (options.compactThresh >> 1);
							if (options._debugLevel >= 2)
								::std::cerr << '(' << oldSize - length(matches) << " matches removed)";
						}
					}
				}
			}
		}
			
		// short-cut negative matches
		if (lastPositive == (int64_t)-1)
			lastPotMatchNo[rseqNo] = (int64_t)-1;
		else
			value(fifo, lastPositive - firstNo).i1 = (int64_t)-1; // the first positive's link to previous is removed

		
	}//swiftFinderR



}//function





}//namespace

#endif
