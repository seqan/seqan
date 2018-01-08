 /*==========================================================================
                SplitRazerS - Fast Split Read Mapping 

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

#include <seqan/store.h>
#include <seqan/align.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/align.h>
#include <seqan/seeds.h>
#include <seqan/bam_io.h>

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
	for (unsigned i1 = 0; i1 <= d_m1; ++i1)
	{
		for (unsigned i2 = 0; i2 <= d_m2; ++i2)
		{
			for (unsigned j = 0; j <= d-i1-i2; ++j)
			{
				if(R-M1-M2>=j)
					sum_prob += 
					(double)(choose(M1,i1)* pow((double)0.25,(double)R-i1-i2-j) *
						 choose(M2,i2) * pow((double)0.75,(double)i1+i2+j) * 
						 choose(R-M1-M2,j));
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
//#for a string of length R in a sequence of length S
template<typename TSize>
TSize
_del(TSize R, TSize S, TSize maxD, TSize minD)
{
	//return((S-R)*(S-R+1)/2);
	return (maxD-minD+1)*(S-R+1) - ((maxD+1)*maxD)/2 + ((minD-1)*minD)/2;
}


//#function for calculating expected matches with deletion of size D
template<typename TSize, typename TValue>
double
_expMatchesDel(TSize R, TSize M1, TSize M2, TValue d, TValue d_m1, TValue d_m2, TSize S, TSize T, TSize maxD, TSize minD )
{
	return ((double)T*(_probability(R,M1,M2,d,d_m1,d_m2))*_delCuts(R,M1,M2)*_del(R,S,maxD,minD));
}

template<typename TReadSet, typename TSize, typename TOptions>
void
expNumRandomMatches(TReadSet &readSet, TSize genomeLen, TOptions & options)
{
	TSize M1 = options.minMatchLen;
	TSize M2 = options.minMatchLen;
	TSize d_m1 = (int) options.maxPrefixErrors;
	TSize d_m2 = (int) options.maxSuffixErrors;
	TSize numReads = length(readSet);
	genomeLen = 1000000 * options.specifiedGenomeLen; // go from Mb to bp
	if(options.anchored) genomeLen = options.maxGap + 2*options.libraryError; // upper bound
	TSize readLen = (numReads > 0) ? length(readSet[0]) : 0;
	TSize numErrors = static_cast<int>(readLen * options.errorRate);
	TSize maxD = options.maxGap;
	TSize minD = options.minGap;
	//expected number of random deletion matches:
	double delMatches = _expMatchesDel(readLen,M1,M2,numErrors,d_m1, d_m2, genomeLen, numReads, maxD, minD);
	if (options.reverse) delMatches *= 2;
	
	//expected number of random deletion matches:
	double insMatches = 0;
	for(TSize insLen = 1; insLen <=readLen-M1-M2; ++insLen)
	{
		numErrors = static_cast<int>((readLen-insLen) * options.errorRate);
		insMatches += _expMatchesIns(readLen,M1,M2,insLen,numErrors,d_m1, d_m2, genomeLen, numReads);
	}
	if (options.reverse) insMatches *= 2;

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
// Load anchored reads from SAM file
template <typename TReadSet, typename TNameSet,typename TReadRegions, typename TRazerSOptions>
bool loadReadsSam(
	TReadSet &reads, 
	TNameSet &fastaIDs, 
	TReadRegions &readRegions,
	const char *fileName, 
	TRazerSOptions &options)
{
	typedef typename Value<TReadRegions>::Type	TRegion;
	//typedef typename Value<TReadSet>::Type		TRead;
	typedef typename Value<TNameSet>::Type		TReadName;
	//typedef typename Value<TRegion,2>::Type		TContigPos;
	typedef typename Value<TRegion,2>::Type		TFlagPos;
	typedef typename Value<TFlagPos,2>::Type	TContigPos;

	bool countN = !(options.matchN || options.outputFormat == 1 );
	if (!empty(CharString(options.outputUnmapped))) countN = false;

	BamFileIn file;
	if (!open(file, fileName)) return false;

	options.maxReadRegionsEnd = 0;
	options.minReadRegionsStart = std::numeric_limits<int>::max();
	TContigPos regionBegin = options.minReadRegionsStart;
	TContigPos regionEnd = options.maxReadRegionsEnd;

		//typename std::ifstream::pos_type lineStart = (*file).tellg();
		//lineStart = lineStart - (std::ifstream::pos_type) 1;
		//	(*file).seekp(lineStart);
 

	int i = 0;
    //int lastFlag = -1;
    //TReadName lastQname;
    //bool lastWasAnchored = false;
	int kickoutcount = 0;

    // Read header.
    BamHeader header;
    readHeader(header, file);

    // Read records.
    BamAlignmentRecord record;
    while (!atEnd(file))
    {
        readRecord(record, file);
        if (record.rID == -1)
            continue;  // Skip if orphan.

        // Get the query name, remove everything after the first space.
        TReadName qname = record.qName;
        cropAfterFirst(qname, IsWhitespace());

        // Evaluate flag.
        if (!hasFlagMultiple(record) || !hasFlagUnmapped(record))
            continue;  // Skip if not unmapped mate.
        bool reverse = hasFlagRC(record); // if reverse this read is expected to match downstream of its mate
        int bitFlag = 0;    // has two bits
        if (reverse) bitFlag += 2;  // first bits says whether reversed (or not ->0)
        if ((record.flag & (1 << 7)) == (1 << 7)) bitFlag +=1; // second bit says whether second mate (or first ->0)

        // Read reference name.  Same behaviour as for query name:  Read up to
        // the first whitespace character and skip to next tab char.
        String<char> chrname = contigNames(context(file))[record.rID];
        //need gnameToIdMap !!

        // Get read begin position.
        TContigPos beginPos = record.beginPos;

        // Get CIGAR string.
        String<CigarElement<> > cigar = record.cigar;

		// read in sequence
        String<Dna5Q> readSeq = record.seq;
        SEQAN_ASSERT_GT(length(readSeq), 0u);
		int readLength = (int)length(readSeq);
		if (countN)
		{
			int count = 0;
			int cutoffCount = (int)(options.errorRate * readLength);
			int j = 0;
			for (; j < readLength; ++j)
				if (getValue(readSeq, j) == 'N')
					if (++count > cutoffCount)
					{
						++kickoutcount;
						break;
					}
			if(j<readLength) continue; // read is skipped because of low quality
		}

		// and associated qualities
        assignQualities(readSeq, record.qual);

	
		// now store the information
		if (options.readNaming == 0
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			|| options.fastaIdQual
#endif
			)  //15578976
			appendValue(fastaIDs, qname);	// append read name Fasta id
		
		
//		if (options.trimLength > 0 && readLength > (unsigned)options.trimLength)
//			resize(readSeq, options.trimLength);
		if((int)length(readSeq) > options.maxReadLength) 
 			options.maxReadLength = length(readSeq);

      //reverseComplement(readSeq);
		appendValue(reads, readSeq, Generous());

		appendValue(readRegions,TRegion(0,TFlagPos(0,0)));
		readRegions[i].i2.i1 = bitFlag; 
		readRegions[i].i2.i2 = (signed)beginPos; 
        if (reverse)
		{
			// libraryLength is outer fragment distance 
			regionBegin = _max((int)beginPos + readLength,(int)beginPos+options.libraryLength-readLength-options.libraryError);
			regionEnd = (signed)beginPos+options.libraryLength+options.libraryError+options.maxGap;
			// expected begin position of read
	//		readRegions[i].i2 = (signed)beginPos + options.libraryLength - readLength; 

		}
		else  // read should map upstream of mapped mate --> set the vorzeichen bit
		{
			regionBegin = _max((int)0,(int)beginPos-options.libraryLength+readLength-options.libraryError-options.maxGap);
			regionEnd = _min((int)beginPos,(int)beginPos-options.libraryLength+2*readLength+options.libraryError);
			// expected end position of read
			//readRegions[i].i2 = - _min((int)beginPos,(int)beginPos-options.libraryLength+2*readLength);//+options.maxGap);
		}
//		readRegions[i].i1 = 0; //currently just ment for one chr at a time... need to add genomeId map
		if(regionEnd > (TContigPos) options.maxReadRegionsEnd) options.maxReadRegionsEnd = (unsigned) regionEnd;
		if(regionBegin < (TContigPos)options.minReadRegionsStart) options.minReadRegionsStart = (unsigned) regionBegin;
        //lastFlag = flag;  
        //lastQname = qname;
		++i;
    }
	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality reads.\n";
	return (i > 0);

}



//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches # necessary to have an own splicedmatch function? 
// planned specs: SpliceSite, General, ... 
template < typename TMatches, typename TCounts, typename TSpec, typename TSwiftL, typename TSwiftR >
void compactSplicedMatches(TMatches &matches, 
			TCounts & /*cnts*/, 
			RazerSOptions<TSpec> &options, 
			bool , 
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
//	::std::sort(it, itEnd, LessSplicedErrors<TMatch>());
	::std::sort(it, itEnd, LessSplicedScore<TMatch>());
	int counter = 0;
	for (; it != itEnd; ++it) 
	{
		++counter;
		if ((*it).orientation == '-') { ++it; continue; }
		if (readNo == (*it).rseqNo)
		{ 
			if ((int)(*it).pairScore <= scoreDistCutOff)
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
					int minScore = (*it).pairScore - 1;
					if (options.purgeAmbiguous && 
						(options.distanceRange == 0 || minScore >= options.maxReadLength - (int)options.distanceRange))
					{
						setMaxErrors(swiftL, readNo, -1);
						setMaxErrors(swiftR, readNo, -1);
						if (options._debugLevel >= 2) ::std::cerr << "(read #" << readNo << " disabled)";
						dit = ditBeg;
					}
					else
					{
						// enough optimal read matches found, disable further search but keep maxHits+1 matches
						if(minScore == options.maxReadLength - 1)
						{
							setMaxErrors(swiftL, readNo, -1);
							setMaxErrors(swiftR, readNo, -1);
							if (options._debugLevel >= 2) ::std::cerr << "(read #" << readNo << " disabled)";
						}
					
						*dit = *it;	++dit; ++it;
						*dit = *it;	++dit;
						continue;
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
// final compacting of spliced matches
template < typename TMatches, typename TCounts, typename TSpec>
void compactAndCountSplicedMatches(TMatches &matches, 
			TCounts & states, 
			RazerSOptions<TSpec> &options, 
			bool )
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreDistCutOff = std::numeric_limits<int>::min();
	
	clear(states);
	
	int bestScore = 0;
	int state = -1; // 0 = unique, 1 = multi, 2 = suboptimal
	int numSuccesful = 0;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;
	// int lastPairId = 0;
	// sort 
	::std::sort(it, itEnd, LessSplicedScoreGPos<TMatch>());
//	::std::sort(it, itEnd, LessSplicedErrorsGPos<TMatch>());
	int counter = 0;
	for (; it != itEnd; ++it) 
	{
		++counter;
		if ((*it).orientation == '-') { ++it; continue; }
		if (readNo == (*it).rseqNo) // current match is either multi or suboptimal
		{ 
			if ((int)(*it).pairScore <= scoreDistCutOff)
			{
				++it;
				continue;
			}
			if (++hitCount >= hitCountCutOff)
			{
				if (hitCount == hitCountCutOff)
				{
					if ((int)(*it).pairScore == bestScore)
					{
						if(state == 0)
							state = 1; // the match before is multi
					}
	
					if(options.purgeAmbiguous)
	     				{
						dit = ditBeg;
						resize(states,numSuccesful);
						state  = -1;
					}
		
				}
				++it;
				continue;
			}
			if ((int)(*it).pairScore == bestScore)
			{
				if(state == 0)
				{
					appendValue(states, 1); // the match before is multi
					state = 1;		// current one is multi too
				}
				else // state > 0
					appendValue(states, state); // either current match is multi or suboptimal, same state
				//std::cout << "state = "<< state << " for PairId = " << lastPairId << std::endl;
			}
			if ((int)(*it).pairScore < bestScore)
			{
				appendValue(states, state); // either current match is suboptimal, state before is what it was
				//std::cout << "state = "<< state << " for PairId = " << lastPairId << std::endl;
				state = 2;		// the current one is suboptimal
			}
			// lastPairId = (*it).pairId;
		}
		else
		{
			if(state != -1)
			{
				appendValue(states, state); // append state of the match before
				//std::cout << "state = "<< state << " for PairId = " << lastPairId << std::endl;
			}
			readNo = (*it).rseqNo;// >> 1;
			hitCount = 0;
			if (options.distanceRange > 0)
				scoreDistCutOff = (*it).pairScore - options.distanceRange;
			ditBeg = dit;
			numSuccesful = (dit - begin(matches, Standard()))/2;
			//std::cout << "numSuccesful=" << numSuccesful << std::endl;
			bestScore = (*it).pairScore;
			// lastPairId = (*it).pairId;
			state = 0;
		}
		*dit = *it;	++dit; ++it;
		*dit = *it;	++dit;
		
	}
	if(state != -1)
	{
		appendValue(states, state);
		//std::cout << "state = "<< state << " for PairId = " << lastPairId << std::endl;
	}
	
	resize(matches, dit - begin(matches, Standard()));
//	std::cout << "lengthmatches = " << length(matches) << " numStates = " << length(states) << std::endl;
}


template<typename TAlign, typename TPosition>
int
countErrorsInAlign(TAlign & align, TPosition end_)
{
	
	//typedef typename Source<TAlign>::Type TSource;
	//typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

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
	typedef Segment<TGenome, InfixSegment>              TGenomeInfix;
	typedef typename Value<TReadSet>::Type              TRead;
	
	// find read match end
	typedef Finder<TGenomeInfix>                        TMyersFinder;
	typedef typename Value<TMyersPatterns>::Type        TMyersPattern;
	
	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>    TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>           TReadRev;
	typedef Finder<TGenomeInfixRev>                     TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>       TMyersPatternRev;

	TMyersFinder myersFinder(inf);
	TMyersPattern &myersPattern = forwardPatterns[rseqNo];  //have to make sure this only contains the prefix
	
#ifdef RAZERS_DEBUG
	::std::cout << "Verify: " << ::std::endl;
	::std::cout << "Genome: " << inf << "\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout << "Read:   " << readSet[rseqNo] << ::std::endl;
#endif
	
	unsigned ndlLength = _min(sequenceLength(rseqNo, readSet),options.minMatchLen);
	int maxScore = std::numeric_limits<int>::min();
	int minScore = - maxNumSeedErrors(options,TSufPrefSpec());

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
    TReadRev            readRev;
    TRead               readInf;  // Needs to be a global variable, since ModifiedString cannot hold a pointer to a temporary.
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

template<typename TOptions>
int
maxNumSeedErrors(TOptions &options, LongestSuffix)
{
	return (int)options.maxSuffixErrors;
}

template<typename TOptions>
int
maxNumSeedErrors(TOptions &options, LongestPrefix)
{
	return (int)options.maxPrefixErrors;
}



template<typename TReadSet, typename TSize, typename TInf, typename TMatch, typename TOptions>
void
extendMatch(TReadSet &readSet, TSize rseqNo, TInf & inf, TMatch &m, TOptions &options, LongestSuffix)
{
#ifdef RAZERS_DEBUG
	::std::cout << " extending match left" << ::std::endl;
#endif

//  XXXIMPROV
//	unsigned lDim0 = _max(0,(int)length(readSet[rseqNo])- 2 * (int)options.minMatchLen);
	unsigned lDim0 = _max(0,(int)length(readSet[rseqNo])-(int)options.minMatchLen);
	unsigned lDim1 = m.gBegin - beginPosition(inf);
//  XXXIMPROV
//	unsigned rDim0 = length(readSet[rseqNo])-1-options.minMatchLen;
	unsigned rDim0 = length(readSet[rseqNo])-1;
	unsigned rDim1 = m.gEnd - beginPosition(inf)-1;

//    Seed<int,SimpleSeed> seed(lDim0, lDim1, rDim0, rDim1);
    Seed<Simple> seed(lDim0, lDim1, rDim0+1, rDim1+1);
	Score<int> scoreMatrix(0,-1,-1,-1);
	int scoreDropOff = static_cast<int>((sequenceLength(rseqNo,readSet) * options.errorRate) - m.editDist);

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
	//extendSeedScore(seed,extScore,scoreDropOff,scoreMatrix, readSet[rseqNo],inf,0,GappedXDrop());
	    
    typedef typename Prefix<typename Value<TReadSet>::Type >::Type TQueryPrefix;
    typedef typename Prefix<TInf>::Type TDatabasePrefix;

    TQueryPrefix queryPrefix = prefix(readSet[rseqNo], beginPositionH(seed));
    TDatabasePrefix databasePrefix = prefix(inf, beginPositionV(seed));
    extScore = _extendSeedGappedXDropOneDirection(seed, databasePrefix, queryPrefix, EXTEND_LEFT, scoreMatrix, scoreDropOff);
    
//	m.gBegin = leftDim1(seed) + beginPosition(inf);
//	m.mScore = rightDim0(seed) - leftDim0(seed) + 1;

	m.gBegin = beginPositionV(seed) + beginPosition(inf);
	m.mScore = endPositionH(seed) - beginPositionH(seed);
	m.editDist -= extScore;

#ifdef RAZERS_DEBUG
	::std::cout << " lDim0: " << beginPositionH(seed) << ::std::endl;
	::std::cout << " lDim1: " << beginPositionV(seed) << ::std::endl;
	::std::cout << " rDim0: " << endPositionH(seed) << ::std::endl;
	::std::cout << " rDim1: " << endPositionV(seed) << ::std::endl;
	::std::cout << " scoreDropOff: "<<scoreDropOff << ::std::endl;
	::std::cout << " readInf: "<< infix(readSet[rseqNo],beginPositionH(seed),endPositionH(seed)+1) << ::std::endl;
	::std::cout << " gInfInf: "<< infix(inf,beginPositionV(seed),endPositionV(seed)+1) << ::std::endl;
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
//	Seed<int,SimpleSeed> seed(lDim0, lDim1, rDim0, rDim1);
	Seed<Simple> seed(lDim0, lDim1, rDim0+1, rDim1+1);
	Score<int> scoreMatrix(0,-1,-1,-1);
	int scoreDropOff = static_cast<int>((sequenceLength(rseqNo,readSet) * options.errorRate) - m.editDist);

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
//  XXXIMPROV
//	extendSeedScore(seed,extScore,scoreDropOff,scoreMatrix, prefix(readSet[rseqNo],length(readSet[rseqNo])-options.minMatchLen),inf,1,GappedXDrop());
    
    typedef typename Suffix<typename Value<TReadSet>::Type >::Type TQuerySuffix;
    typedef typename Suffix<TInf>::Type TDatabaseSuffix;

    TQuerySuffix querySuffix = suffix(readSet[rseqNo], endPositionH(seed));
    TDatabaseSuffix databaseSuffix = suffix(inf, endPositionV(seed));
    extScore = _extendSeedGappedXDropOneDirection(seed, databaseSuffix, querySuffix, EXTEND_RIGHT, scoreMatrix, scoreDropOff);


	//extendSeedScore(seed,extScore,scoreDropOff,scoreMatrix, readSet[rseqNo],inf,1,GappedXDrop());
	//m.gEnd = rightDim1(seed) + 1 + beginPosition(inf);
	//m.mScore = rightDim0(seed) - leftDim0(seed) + 1;
	m.gEnd = endPositionV(seed)  + beginPosition(inf);
	m.mScore = endPositionH(seed) - beginPositionH(seed);
	m.editDist -= extScore;

#ifdef RAZERS_DEBUG
	::std::cout << " lDim0: " << beginPositionH(seed) << ::std::endl;
	::std::cout << " lDim1: " << beginPositionV(seed) << ::std::endl;
	::std::cout << " rDim0: " << endPositionH(seed) << ::std::endl;
	::std::cout << " rDim1: " << endPositionV(seed) << ::std::endl;
	::std::cout << " scoreDropOff: "<<scoreDropOff << ::std::endl;
	::std::cout << " readInf: "<< infix(readSet[rseqNo],beginPositionH(seed),endPositionH(seed)) << ::std::endl;
	::std::cout << " gInfInf: "<< infix(inf,beginPositionV(seed),endPositionV(seed)) << ::std::endl;
	::std::cout << " read: "<< readSet[rseqNo] << ::std::endl;
	::std::cout << " gInf: "<< inf << ::std::endl;
#endif
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
	TReadSet& readSet,							// original readSet
	TPattern&,					
	RazerSOptions<TSpec> const &options,		// RazerS options
	SwiftSemiGlobalHamming,						// HammingDistance
	LongestPrefix)								// LongestPrefix
{
	
	typedef Segment<TGenome, InfixSegment>                  TGenomeInfix;
	//typedef typename Size<TGenomeInfix>::Type               TSize;
	//typedef typename Value<TGenomeInfix>::Type              TDna;
	//typedef typename Position<TGenomeInfix>::Type           TPosition;
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
//	unsigned maxErrorsSeed = (unsigned)(options.minMatchLen * options.errorRate);	
	unsigned maxTotalErrors = (unsigned)(length(readSet[rseqNo]) * options.errorRate);	
	unsigned maxErrorsSeed = options.maxPrefixErrors;	
	if(maxErrorsSeed > maxTotalErrors) maxErrorsSeed = maxTotalErrors;
	// unsigned minSeedErrors = maxErrorsSeed + 1;
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
						// we exclude this last error position 
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
			// minSeedErrors = seedErrors;
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
	//typedef typename Size<TGenomeInfix>::Type               TSize;
	//typedef typename Value<TGenomeInfix>::Type              TDna;
	//typedef typename Position<TGenomeInfix>::Type           TPosition;
	typedef typename Value<TReadSet>::Type 			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Infix<TRead>::Type 			TReadInf;
	typedef typename Iterator<TReadInf, Standard>::Type	TReadIterator;
	
	if (length(genomeInf) < options.minMatchLen) return false;
	TRead read = infix(readSet[rseqNo],options.minMatchLen,length(readSet[rseqNo]));
	
		
#ifdef RAZERS_DEBUG
	bool debug = true;
	if(debug)
	{
		::std::cout<< "suffixmatching\n";
		::std::cout << "genome=" << genomeInf << "\nread  =" << read <<"\n";
	}
#endif

	TReadIterator ritEnd	= end(read, Standard())-1;
	TReadIterator ritBeg	= begin(read, Standard());
	TGenomeIterator git	= end(genomeInf, Standard())-1;
	TGenomeIterator gitBeg	= begin(genomeInf, Standard()) + options.minMatchLen;
	
	// this is max number of errors the seed should have
//	unsigned maxErrorsSeed = (unsigned)(options.minMatchLen * options.errorRate);	
	unsigned maxTotalErrors = (unsigned)(length(readSet[rseqNo]) * options.errorRate);	
	unsigned maxErrorsSeed = options.maxSuffixErrors;	
	if(maxErrorsSeed > maxTotalErrors) maxErrorsSeed = maxTotalErrors;
	// unsigned minSeedErrors = maxErrorsSeed + 1;
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
			//if(debug)::std::cout << *r << "\t" << *g << "\n";
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
						// we exclude this last error position 
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
			// minSeedErrors = seedErrors;
			minTotalErrors = totalErrors;
			bestHitLength = hitLength;
			m.gEnd = git - begin(host(genomeInf), Standard()) + 1;
			
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



// Edit distance match combination forward
template <typename TScore>
bool
findBestSplitPosition(String<Pair<TScore,int> > & maxColsL,
					  String<Pair<TScore,int> > & maxColsR,
					  int & rowPosL1,
					  int rowPosL2,
					  int & rowPosR1,
					  int rowPosR2,
					  int seq0Len,
					  int & traceExt,
					  OrientationForward,
					  SwiftSemiGlobal)
{

#ifdef RAZERS_DEBUG
	::std::cout << "findBestSplitEditForward\n";
#endif

	TScore maxSum = std::numeric_limits<TScore>::min();
	int bestL = rowPosL1;
	int bestR = rowPosR1;
	int bestTraceExtL = rowPosL1;
	int bestTraceExtR = rowPosL1;
	while (rowPosL1 <= rowPosL2 && rowPosR1 >= rowPosR2)
	{
		// this is to prevent same bases from being used in both prefix and suffix match 
		// this works, because we store the FIRST bestScore in each row
		if (!(maxColsL[rowPosL1].i2 + maxColsR[rowPosR1].i2 <= seq0Len))
		{
			++rowPosL1;
			--rowPosR1;
			continue;
		}
		if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 > maxSum) 
		{
			maxSum = maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1;
			bestL = rowPosL1;
			bestR = rowPosR1;
			bestTraceExtL = rowPosL1;
		}
		else if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 == maxSum)
			bestTraceExtR = rowPosL1;
			
		++rowPosL1;
		--rowPosR1;
	}
	traceExt = bestTraceExtR - bestTraceExtL;
	rowPosL1 = bestL;
	rowPosR1 = bestR;

	return true;
}

// ---------------------------------------------------------------------------
// BEGIN BREAKPOINT COMPUTATION CODE FROM ALIGN MODULE
// ---------------------------------------------------------------------------

// TODO(holtgrew): These functions also has to be copied to msplazers.
// TODO(holtgrew): This breakpoint computation code is basically reusable, but needs some love (i.e. documentation, demos, interfaces)

template <typename TAlign, typename TStringSet, typename TTrace, typename TValPair, typename TIndexPair, typename TDiagonal>
inline void
_alignBandedNeedlemanWunschTrace(TAlign& align,
					   TStringSet const& str,
					   TTrace const& trace,
					   TValPair const& overallMaxValue,
					   TIndexPair const& overallMaxIndex,
					   TDiagonal const diagL,
					   TDiagonal const diagU)
{
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;

	// Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize lo_row = (diagU <= 0) ? -1 * diagU : 0;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	
	//// Debug stuff
	//TColumn originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cout << count << ',';
	//		++count;
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	// Start the trace from the cell with the max value
	TSize row = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[0] : overallMaxIndex[2];
	TSize col = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[1] : overallMaxIndex[3];

	// Handle the skipped sequence parts
	TSize actualRow = row + lo_row;
	TSize actualCol = col + diagL + actualRow;
	if (actualCol + 1 < len1) _alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, (len1 - (actualCol + 1)),  Horizontal);
	if (actualRow + 1 < len2) _alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, (len2 - (actualRow + 1)),  Vertical);

	if ((actualRow != 0) && (actualCol != 0)) {
		// Find initial direction
		TTraceValue tv = trace[row * diagonalWidth + col];
		if (tv == Horizontal) --col;
		else if (tv == Vertical) {--row; ++col;} 
		else --row;
	
		// Walk until we hit a border
		TSize seqLen = 1;
		TTraceValue newTv = tv;
		while(true) {
			actualRow = row + lo_row;
			actualCol = col + diagL + actualRow;
			newTv = trace[row * diagonalWidth + col];

			// Check if we hit a border
			if ((actualRow == 0) || (actualCol == 0)) break;
			else {
				//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
				if (tv == Diagonal) {
					if (newTv == Horizontal) {
						_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
						--col; seqLen = 1;
					} else if (newTv == Vertical) {
						_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
						--row; ++col; seqLen = 1;
					} else {
						--row; ++seqLen;
					}
				} else {
					if (tv == Horizontal) { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Vertical) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
							--row; ++col; seqLen = 1;
						} else {
							--col; ++seqLen;
						}
					} else { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Horizontal) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
							--col; seqLen = 1;
						} else {
							--row; ++col; ++seqLen;
						}
					}
				}
				tv = newTv;
			}
		}
	
		// Align left overs
		if (seqLen) _alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
	}

	// Handle the remaining sequence
	if (actualCol != 0) _alignTracePrint(align, str[0], str[1], (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _alignTracePrint(align, str[0], str[1], (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);
}

template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TDiagonal, typename TAlignConfig>
inline typename Value<TScore>::Type
_alignBandedNeedlemanWunsch(TTrace& trace,
			TStringSet const& str,
			TScore const & sc,
			TValPair& overallMaxValue,
			TIndexPair& overallMaxIndex,
			TDiagonal diagL,
			TDiagonal diagU,
			TAlignConfig const,
			String<Pair<typename Value<TScore>::Type,int> > & maxCols,
			unsigned minColNum)
{
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TTrace>::Type TSize;

//      ::std::cout << "hierrrrrrrrrr banded ali\n";
	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	TString const& str1 = str[0];
	TString const& str2 = str[1];
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
//      ::std::cout << "len1="<<len1 << " len2=" << len2 << ::std::endl;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	TSize hi_diag = diagonalWidth;
	TSize lo_diag = 0;
	if (diagL > 0) lo_diag = 0;
	else lo_diag = (diagU < 0) ? hi_diag : (TSize) (1 - diagL);
	TSize lo_row = (diagU <= 0) ? -diagU : 0;
	TSize hi_row = len2;
	if (len1 - diagL < hi_row) hi_row = len1 - diagL;
	TSize height = hi_row - lo_row;

	clear(maxCols);

	typedef String<TScoreValue> TRow;
	TRow mat;
	resize(mat, diagonalWidth);
	resize(trace, height * diagonalWidth);
//      ::std::cout <<height << "<-hieght\n";
	overallMaxValue[0] = std::numeric_limits<TScoreValue>::min();
	overallMaxValue[1] = std::numeric_limits<TScoreValue>::min();
	overallMaxIndex[0] = diagonalWidth;     overallMaxIndex[1] = height;
	overallMaxIndex[2] = diagonalWidth;     overallMaxIndex[3] = height;

	//// Debug stuff
	//String<TScoreValue> originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//      for(TSize j=0; j<len1; ++j) {
	//              value(originalMat, i * len1 + j) = count;
	//              std::cerr << count << ',';
	//              ++count;
	//      }
	//      std::cerr << std::endl;
	//}
//std::cerr << std::endl;

	// Classical DP with affine gap costs
	typedef typename Iterator<TRow, Standard>::Type TRowIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TSize actualCol = 0;
	TSize actualRow = 0;
	TScoreValue verti_val = 0;
	TScoreValue hori_val = 0;
	for(TSize row = 0; row < height; ++row) {
		TScoreValue maxRowVal = std::numeric_limits<TScoreValue>::min();
		unsigned maxRowCol = 0;
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if (row + lo_row >= len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TRowIter matIt = begin(mat, Standard()) + lo_diag;
		hori_val = std::numeric_limits<TScoreValue>::min();
		for(TSize col = lo_diag; col<hi_diag; ++col, ++matIt, ++traceIt) {
			actualCol = col + diagL + actualRow;
			//std::cerr << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for mat
				*matIt += score(const_cast<TScore&>(sc), sequenceEntryForScore(const_cast<TScore&>(sc), str1, ((int) actualCol - 1)),
				                sequenceEntryForScore(const_cast<TScore&>(sc), str2, ((int) actualRow - 1)));
				*traceIt = Diagonal;
				if ((verti_val = (col < diagonalWidth - 1) ? *(matIt+1) +
				    scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, ((int) actualCol - 1)),
				                           sequenceEntryForScore(sc, str2, ((int) actualRow - 1))) : std::numeric_limits<TScoreValue>::min()) > *matIt)
				{
					*matIt = verti_val;
					*traceIt = Vertical;
				}
				if ((hori_val = (col > 0) ? hori_val +
				    scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, ((int) actualCol - 1)),
				                             sequenceEntryForScore(sc, str2, ((int) actualRow - 1))) : std::numeric_limits<TScoreValue>::min()) > *matIt)
				{
					*matIt = hori_val;
					*traceIt = Horizontal;
				}
				hori_val = *matIt;
			} else {
				// Usual initialization for first row and column
				if (actualRow == 0)
				    _initFirstRow(TAlignConfig(), *matIt, (TScoreValue) actualCol *
				                  scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, std::max(0,((int) actualCol - 1))),
				                                           sequenceEntryForScore(sc, str2, 0)));
				else {
					_initFirstColumn(TAlignConfig(), *matIt, (TScoreValue) actualRow *
					                 scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, 0),
					                                        sequenceEntryForScore(sc, str2, std::max(0,((int) actualRow - 1)))));
					hori_val = *matIt;
				}
			}
			if(*matIt > maxRowVal && actualCol >= minColNum)
			{
				maxRowVal = *matIt;
				maxRowCol = actualCol;
			}
			// Store the maximum
			if (actualCol == len1 - 1) _lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			if (actualRow == len2 - 1) _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			//std::cerr << row << ',' << col << ':' << *matIt << std::endl;
		}
		appendValue(maxCols,Pair<typename Value<TScore>::Type,int>(maxRowVal,maxRowCol));
	//      std::cout << maxRowVal << ","<<maxRowCol << std::endl;
	}
	return (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxValue[0] : overallMaxValue[1];
}

template<typename TAlign, typename TStringSet, typename TScore, typename TAlignConfig, typename TDiagonal>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
			TStringSet const& str,
			TScore const& sc,
			TAlignConfig const,
			TDiagonal diag1,
			TDiagonal diag2,
			String<Pair<typename Value<TScore>::Type,int> > &maxCols,
			unsigned minColNum,
			NeedlemanWunsch)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	
	// Maximum value
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];
	
	// Create the trace
	String<TraceBack> trace;
	TScoreValue maxScore = _alignBandedNeedlemanWunsch(trace, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig(),maxCols, minColNum);
	
	// Follow the trace and create the graph
	_alignBandedNeedlemanWunschTrace(align, str, trace, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2);
	
	return maxScore;
}

template<typename TStringSet, typename TScore, typename TAlignConfig, typename TDiagonal>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
			TScore const& sc,
			TAlignConfig const,
			TDiagonal diag1,
			TDiagonal diag2,
			String<Pair<typename Value<TScore>::Type,int> > &maxCols,
			unsigned minColNum,
			NeedlemanWunsch)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	
	// Maximum value
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];
	
	// Calculate the score
	String<TraceBack> trace;
	return _alignBandedNeedlemanWunsch(trace, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig(),maxCols,minColNum);
}

// ---------------------------------------------------------------------------
// END BREAKPOINT COMPUTATION CODE FROM ALIGN MODULE
// ---------------------------------------------------------------------------

// Edit distance match combination reverse
template <typename TScore>
bool
findBestSplitPosition(String<Pair<TScore,int> > & maxColsL,
					  String<Pair<TScore,int> > & maxColsR,
					  int & rowPosL1,
					  int rowPosL2,
					  int & rowPosR1,
					  int rowPosR2,
					  int seq0Len,
					  int & traceExt,
					  OrientationReverse,
					  SwiftSemiGlobal)
{
#ifdef RAZERS_DEBUG
	::std::cout << "findBestSplitEditReverse\n";
#endif

	TScore maxSum = std::numeric_limits<TScore>::min();
	int bestL = rowPosL2;
	int bestR = rowPosR2;
	int bestTraceExtR = rowPosL1;
	int bestTraceExtL = rowPosL1;
	
	while (rowPosL1 <= rowPosL2 && rowPosR1 >= rowPosR2)
	{
		if(!(maxColsL[rowPosL1].i2 + maxColsR[rowPosR1].i2 <= seq0Len))
		{
			++rowPosL1;
			--rowPosR1;
			continue;
		}
		if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 >= maxSum)
		{
			maxSum = maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1;
			bestL = rowPosL1;
			bestR = rowPosR1;
			bestTraceExtR = rowPosL1;
			if(maxColsL[rowPosL1].i1 + maxColsR[rowPosR1].i1 > maxSum) bestTraceExtL = rowPosL1;
		}
		++rowPosL1;
		--rowPosR1;
	}
	rowPosL1 = bestL;
	rowPosR1 = bestR;
	traceExt = bestTraceExtR - bestTraceExtL;
	return true;
}





// Edit distance match combination wrapper
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
	unsigned maxErrors = static_cast<unsigned>(readLength * options.errorRate);

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

	//potential insertion // for edit distance this is a heuristic
	if(mR.gEnd - mL.gBegin < readLength)
	{
#ifdef RAZERS_DEBUG
		::std::cout << "insertion\n";
#endif
		if(mR.gEnd - mL.gBegin < static_cast<long int>(2*options.minMatchLen))  //too close together // actually minus allowed seed errors
			return false; 

		if(mL.gEnd < mR.gBegin)  //prefix and suffix match do not meet
			return false;

		if(mR.mScore + mL.mScore == mR.gEnd - mL.gBegin)
			if(mL.gEnd == mR.gBegin)
				if(mR.editDist + mL.editDist <= maxErrors ) //prefix and suffix match meet and do not overlap --> perfect
					return true;

		int diag1L = -static_cast<int>(maxErrors) + mL.seedEditDist;
		int diag2L = maxErrors - mL.seedEditDist;
		int diag1R = -static_cast<int>(maxErrors) + mR.seedEditDist;
		int diag2R = maxErrors - mR.seedEditDist;
 		int minColNum = 0;

		// genomeInf is the shorter sequence --> find best split position on genome
		// rows in alignment matrix represent genome position
		StringSet<TGenomeInf,Dependent<> > strSetL;
		appendValue(strSetL,readInfL);
		appendValue(strSetL,genomeInfL);
		Graph<Alignment<StringSet<TGenomeInf,Dependent<> >, void> > alignL(strSetL);
		String<Pair<int,int> > maxColsL;
		_globalAlignment(alignL,strSetL,scoreType,AlignConfig<false,false,false,false>(),diag1L,diag2L,maxColsL,minColNum,NeedlemanWunsch());
	
		StringSet<TGenomeInfRev,Dependent<> > strSetR;
		appendValue(strSetR,readInfR);
		appendValue(strSetR,genomeInfR);
		Graph<Alignment<StringSet<TGenomeInfRev,Dependent<> >, void > > alignR(strSetR);
		String<Pair<int,int> > maxColsR;
		_globalAlignment(alignR,strSetR,scoreType,AlignConfig<false,false,false,false>(),diag1R,diag2R,maxColsR,minColNum,NeedlemanWunsch());

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
		int seq0Len = readLength - 2*options.minMatchLen;
		int traceExt = 0;
		if(orientation == 'R')
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,traceExt,OrientationReverse(),SwiftSemiGlobal());
		else
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,traceExt,OrientationForward(),SwiftSemiGlobal());

#ifdef RAZERS_DEBUG
		::std::cout << "nachher\nrowPosL1=" << rowPosL1 << ::std::endl;		
		::std::cout << "rowPosL2=" << rowPosL2 << ::std::endl;		
		::std::cout << "rowPosR1=" << rowPosR1 << ::std::endl;		
		::std::cout << "rowPosR2=" << rowPosR2 << ::std::endl;		
		::std::cout << "mR.editDist=" << mR.editDist << ::std::endl;
		::std::cout << "mL.editDist=" << mL.editDist << ::std::endl;
		::std::cout << "maxErros=" << maxErrors << ::std::endl;
#endif
		mL.editDist = mR.seedEditDist - maxColsL[rowPosL1].i1;	//scores are negative
		mR.editDist = mL.seedEditDist - maxColsR[rowPosR1].i1;
		if(mR.editDist + mL.editDist > maxErrors )
			return false;

		mR.traceExtension = mL.traceExtension = traceExt;

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
		
		int diag1L = -static_cast<int>(maxErrors) + mL.seedEditDist;
		int diag2L = maxErrors - mL.seedEditDist;
		int diag1R = -static_cast<int>(maxErrors) + mR.seedEditDist;
		int diag2R = maxErrors - mR.seedEditDist;
		int minColNum = 0;
		
		// readInf is the shorter sequence --> find best split position on read
		// rows in alignment matrix represent read position
		StringSet<TGenomeInf> strL;
		appendValue(strL,genomeInfL);
		appendValue(strL,readInfL);
		String<Pair<int,int> > maxColsL;
		_globalAlignment(strL,scoreType,AlignConfig<false,false,false,false>(),diag1L,diag2L,maxColsL,minColNum,NeedlemanWunsch());
	
		StringSet<TGenomeInfRev> strR;
		appendValue(strR,genomeInfR);
		appendValue(strR,readInfR);
		String<Pair<int,int> > maxColsR;

		_globalAlignment(strR,scoreType,AlignConfig<false,false,false,false>(),diag1R,diag2R,maxColsR,minColNum,NeedlemanWunsch());
	
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
		int traceExt = 0;
		if(orientation == 'R')
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,traceExt,OrientationReverse(),SwiftSemiGlobal());
		else
			findBestSplitPosition(maxColsL,maxColsR,rowPosL1,rowPosL2,rowPosR1,rowPosR2,seq0Len,traceExt,OrientationForward(),SwiftSemiGlobal());


		// best split position in read is saved in rowPosL1 (and rowPosR1)
		mL.editDist = mL.seedEditDist - maxColsL[rowPosL1].i1;
		mR.editDist = mR.seedEditDist - maxColsR[rowPosR1].i1;


		if(mR.editDist + mL.editDist > maxErrors)
			return false;

		mR.mScore = options.minMatchLen + rowPosR1;
		mR.gBegin = mR.gEnd - maxColsR[rowPosR1].i2 - mR.gSeedLen; //genomic position of best read split
		mL.mScore = options.minMatchLen + rowPosL1;
		mR.traceExtension = mL.traceExtension = traceExt;

		mL.gEnd = mL.gBegin + maxColsL[rowPosL1].i2 + mL.gSeedLen; //genomic position of best read split
#ifdef RAZERS_DEBUG
		::std::cout << "rowPosL1=" << rowPosL1 << ::std::endl;		
		::std::cout << "rowPosL2=" << rowPosL2 << ::std::endl;		
		::std::cout << "rowPosR1=" << rowPosR1 << ::std::endl;		
		::std::cout << "rowPosR2=" << rowPosR2 << ::std::endl;	

		std::cout << "after split:\nmL.mScore =" << mL.mScore << "\n";
		std::cout << "mL.gBegin =" << length(genome)-mL.gBegin << "\n";
		std::cout << "mL.gEnd =" << length(genome)-mL.gEnd << "\n";
		std::cout << "mL.editDist =" << mL.editDist << "\n";

		std::cout << "mR.mScore =" << mR.mScore << "\n";
		std::cout << "mR.gBegin =" << length(genome)-mR.gBegin << "\n";
		std::cout << "mR.gEnd =" << length(genome)-mR.gEnd << "\n";
		std::cout << "mR.editDist =" << mR.editDist << "\n";
#endif

	}
	maxErrors = static_cast<unsigned>((mR.mScore + mL.mScore) * options.errorRate);
	if(mR.editDist + mL.editDist > maxErrors) return false; //make sure percent identity critrium is fulfilled on matched part of the read

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
	unsigned maxErrors = static_cast<unsigned>(readLength * options.errorRate);


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

		if(mR.gEnd - mL.gBegin < static_cast<long int>(2*options.minMatchLen))//too close together
			 return false; 

		if(mR.mScore + mL.mScore < mR.gEnd - mL.gBegin) //prefix and suffix match do not meet
			return false;

		if((mR.mScore + mL.mScore == mR.gEnd - mL.gBegin) && (mR.editDist + mL.editDist > maxErrors)) //prefix and suffix match meet but too many errors
			return false;
		
//		if((mR.gEnd - mL.gBegin <= mL.mScore) || (mR.gEnd - mL.gBegin <= mR.mScore))//too close together 
//			 return false; 
		int traceExt = 0;	 
		bool result = findBestSplitPosition(read,genomeInf,mL.mScore,mR.mScore,mL.editDist,mR.editDist, traceExt, options, orientation, SwiftSemiGlobalHamming());
		if(!result || mR.editDist + mL.editDist > maxErrors) return false;
		mR.gBegin = mR.gEnd - mR.mScore;
		mL.gEnd = mL.gBegin + mL.mScore;
		mR.traceExtension = mL.traceExtension = traceExt;


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

		int traceExt = 0;
		bool result = findBestSplitPosition(genomeInf,read,mL.mScore,mR.mScore,mL.editDist,mR.editDist, traceExt, options, orientation,SwiftSemiGlobalHamming());
		
		if(!result || mR.editDist + mL.editDist > maxErrors) return false;

		mR.traceExtension = mL.traceExtension = traceExt;
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

	maxErrors = static_cast<unsigned>((mR.mScore + mL.mScore) * options.errorRate);
	if(mR.editDist + mL.editDist > maxErrors) return false; //make sure percent identity critrium is fulfilled on matched part of the read

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
			int & traceExt, 
			TOptions & options,
			char orientation, 
			SwiftSemiGlobalHamming)
{
	
#ifdef RAZERS_DEBUG
	::std::cout << "findBestSplitHamming"<<orientation<<"\n";
#endif

	// usually, both types should be the same, but you never know...
	//typedef typename Iterator<TLongerSegment const>::Type TLongIterator;
	//typedef typename Iterator<TShorterSegment const>::Type TShortIterator;
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
	// int bestLErrors = 0;
	// int bestRErrors = 0;
	int bestPos = shortPos;
	int errorsL = 0;
	int errorsR = 0;
	int errorsPosL = 0;
	int errorsPosR = 0;
	
	// determine trace extensions
	int bestTraceExtL = shortPos;
	int bestTraceExtR = shortPos;

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
			// bestLErrors = errorsPosL;
			// bestRErrors = errorsPosR;
			bestPos = shortPos + 1;
			if(errorsPosL+errorsPosR < bestSumErrors) bestTraceExtL = bestPos;
		}
		if(errorsPosL+errorsPosR == bestSumErrors)
		{
			bestTraceExtR = shortPos + 1;
		}
		++leftLongPos;
		++rightLongPos;
		++shortPos;
	}
	
	// trace extension:
	traceExt = bestTraceExtR - bestTraceExtL;
	
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
	std::cout << "traceExt= " << traceExt << std::endl;
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
	typename TShapeL,
	typename TShapeR,
	typename TSwiftSpec >
int mapSplicedReads(
	TMatches &		matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet_ & 		readSet,
	TReadRegions &		readRegions,
	TCounts &		cnts,
	RazerSOptions<TSpec> &	options,
	TShapeL const &		shapeL,
	TShapeR const &		shapeR,
	Swift<TSwiftSpec> const)
{


	typedef typename Value<TReadSet_>::Type								TRead;
	typedef StringSet<typename Infix<TRead>::Type>						TReadSet;
	typedef Index<TReadSet, IndexQGram<TShapeL, TQGramIndexSpec> >		TIndexL;			// q-gram index left
	typedef Index<TReadSet, IndexQGram<TShapeR, TQGramIndexSpec> >		TIndexR;			// q-gram index right
	typedef Pattern<TIndexL, Swift<TSwiftSpec> >							TSwiftPatternL;	// filter	//should be the same thing for left and right
	typedef Pattern<TIndexR, Swift<TSwiftSpec> >							TSwiftPatternR;	// filter
	typedef Pattern<TRead, MyersUkkonen>								TMyersPattern;	// verifier
	
	
	// split reads over two indices, one for prefixes, the other for suffixes
	TReadSet readSetL, readSetR;
	unsigned readCount = length(readSet);
	resize(readSetL, readCount);
	resize(readSetR, readCount);

	if(options._debugLevel > 0)
	{
		int64_t genomeLen = static_cast<int64_t>(3000000000lu) * 2;					// ufff make that an option
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
	TIndexL swiftIndexL(readSetL, shapeL);
	TIndexR swiftIndexR(readSetR, shapeR);
	
#ifdef RAZERS_OPENADDRESSING
	swiftIndexL.alpha = 2;
	swiftIndexR.alpha = 2;
#endif
	
	cargo(swiftIndexL).abundanceCut = options.abundanceCut;
	cargo(swiftIndexR).abundanceCut = options.abundanceCut;
	cargo(swiftIndexL)._debugLevel = 0;
	cargo(swiftIndexR)._debugLevel = options._debugLevel;
	
	// configure Swift
	TSwiftPatternL swiftPatternL(swiftIndexL);
	TSwiftPatternR swiftPatternR(swiftIndexR);
	swiftPatternL.params.minThreshold = options.threshold;
	swiftPatternR.params.minThreshold = options.thresholdR;
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
	
	unsigned numFiles = length(genomeFileNameList);
	unsigned gseqNo = 0;
	
	// open genome files, one by one	
	for (unsigned filecount = 0; filecount < numFiles; ++filecount)
	{
		// open genome file	
		SeqFileIn file;
		if (!open(file, toCString(genomeFileNameList[filecount])))
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
		for(; !atEnd(file); ++gseqNo)
		{
            readRecord(id, genome, file);               // read Fasta id and sequence
			if (options.genomeNaming == 0)
            {
                cropAfterFirst(id, IsWhitespace());     // crop id after the first whitespace
				appendValue(genomeNames, id, Generous());
            }
			gnoToFileMap.insert(::std::make_pair(gseqNo,::std::make_pair(genomeName,gseqNoWithinFile)));
			
			if (options.forward)
				mapSplicedReads(matches, genome, gseqNo, readSet, readRegions, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'F', options);
	
			if (options.reverse)
			{
				reverseComplement(genome);
				mapSplicedReads(matches, genome, gseqNo, readSet, readRegions, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
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


// returns start pos furthest to the left
template<typename TValue, typename TRegion, typename TOptions>
inline typename Value<typename Value<TRegion,2>::Type,2>::Type
regionStartPos(TRegion & readRegion,
			   TValue readLength,
			   TOptions & options)
{
	typedef typename Value<typename Value<TRegion,2>::Type,2>::Type TSignedPos;

	if(readRegion.i2.i1 < 2)	// i2 stores expected end pos of match
		return _max((TSignedPos)0,(TSignedPos)(readRegion.i2.i2 - options.libraryLength + readLength - options.libraryError - options.maxGap));
	else				 	    // i2 stores expected start pos
		return (readRegion.i2.i2  + options.libraryLength - readLength - options.libraryError);
	
}

// returns end pos furthest to the right
template<typename TValue, typename TRegion, typename TOptions>
inline typename Value<typename Value<TRegion,2>::Type,2>::Type
regionEndPos(TRegion & readRegion,
			   TValue readLength,
			   TOptions & options)
{
	typedef typename Value<typename Value<TRegion,2>::Type,2>::Type TSignedPos;

	if(readRegion.i2.i1 < 2)	// i2 stores expected end pos of match
		return (TSignedPos)readRegion.i2.i2-options.libraryLength+2*readLength + options.libraryError;
	else					// i2 stores expected start pos
		return readRegion.i2.i2+ options.libraryLength + options.libraryError + options.maxGap;
	
}



template<typename TMatch, typename TValue, typename TRegion, typename TOptions>
inline bool
isValidRegion(TMatch & mL,
			   TMatch & mR,
			   TValue readLength,
			   TRegion & readRegion,
			   TOptions & options)
{ 
	typedef typename Value<typename Value<TRegion,2>::Type,2>::Type TSignedPos;

	if(readRegion.i2.i1 < 2)	// mR.gEnd needs to lie within a specific region
	{
		TSignedPos expEndPos = readRegion.i2.i2 - options.libraryLength + 2*readLength;
		if (mR.gEnd < expEndPos - options.libraryError ||
			mR.gEnd > expEndPos + options.libraryError )
			return false;
	}
	else					// mL.gBegin needs to lie within a specific region
	{
		TSignedPos expStartPos = readRegion.i2.i2 + options.libraryLength - readLength;
		if (mL.gBegin < expStartPos - options.libraryError ||
			mL.gBegin > expStartPos + options.libraryError )
			return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TOriReadSet,
	typename TReadIndexL, 
	typename TReadIndexR, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TCounts,
	typename TSpec >
void mapSplicedReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	TOriReadSet &readSet,			// reads
	TReadRegions & readRegions,
	Pattern<TReadIndexL, Swift<TSwiftSpec> > &swiftPatternL,	// left index
	Pattern<TReadIndexR, Swift<TSwiftSpec> > &swiftPatternR,	// right index
	TVerifier &forwardPatternsL,
	TVerifier &forwardPatternsR,
	TCounts & cnts,
	char orientation,
	RazerSOptions<TSpec> &options)
{
	typedef typename Value<TOriReadSet>::Type TRead;
	typedef typename Fibre<TReadIndexL, FibreText>::Type	TReadSetL;
	typedef typename Fibre<TReadIndexR, FibreText>::Type	TReadSetR;
	typedef typename Size<TGenome>::Type			TSize;
	typedef typename Position<TGenome>::Type		TGPos;
	typedef typename MakeSigned_<TGPos>::Type		TSignedGPos;
	typedef typename Value<TMatches>::Type			TMatch;
	typedef typename Infix<TGenome>::Type			TGenomeInf;
	
	// Prefix-Suffix filtration
	//typedef Finder<TGenome, Swift<TSwiftSpec> >		TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >		TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >	TSwiftFinderR;
	//typedef Pattern<TReadIndexL, Swift<TSwiftSpec> >	TSwiftPatternL;
	//typedef Pattern<TReadIndexR, Swift<TSwiftSpec> >	TSwiftPatternR;
	
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
	
	TReadSetL &readSetL = host(host(swiftPatternL));
	TReadSetR &readSetR = host(host(swiftPatternR));
	
	if (empty(readSetL) || empty(readSetR))
		return;
	
	
	TSignedGPos maxDistance = options.maxGap;
//	raus:TSignedGPos maxDistance = options.maxGap + (int)options.maxReadLength;
	//TSignedGPos minDistance = options.minGap;// + 2*options.minMatchLen ;
	TSignedGPos minDistance = 2*options.minMatchLen ;
	if(!options.hammingOnly) 
		minDistance -= (options.maxPrefixErrors + options.maxSuffixErrors);

	// exit if contig is shorter than minDistance
	if (length(genome) <= (unsigned)2*options.minMatchLen)
		return;
	
	TGPos scanBegin = 0;
	TGPos scanEnd = length(genome);
	if(!empty(readRegions))
	{
		if(options._debugLevel > 1) 
		{
			std::cout << "MaxRegionEndPos=" << options.maxReadRegionsEnd << std::endl;
			std::cout << "MinRegionStartPos=" << options.minReadRegionsStart << std::endl;
		}
		scanBegin = options.minReadRegionsStart;
		if (scanBegin > length(genome))
			scanBegin = length(genome);
		scanEnd = options.maxReadRegionsEnd;
		if (scanEnd > length(genome))
			scanEnd = length(genome);
	}
	TGenomeInf genomeInf = infix(genome, scanBegin, scanEnd);
	//TSwiftFinderL swiftFinderL(genome, options.repeatLength, 1);
	TSwiftFinderL swiftFinderL(genomeInf, options.repeatLength, 1);
	TSwiftFinderR swiftFinderR(genomeInf, options.repeatLength, 1);
	
	TDequeue fifo;						// stores potential prefix matches
	String<int64_t> lastPotMatchNo;		// last number of a potential prefix match
	int64_t lastNo = 0;					// last number over all potential prefix matches in the queue
	int64_t firstNo = 0;				// first number over all potential prefix matches in the queue
	Pair<TGPos> gPair;
	
	resize(lastPotMatchNo, length(host(swiftPatternL)), (int64_t)-1, Exact());
	
	String<Pair<TGPos> > lastRightMatch;		// begin and end of last verified suffix match
	resize(lastRightMatch, length(host(swiftPatternL)), Pair<TGPos>(0,0), Exact());
	
	TSize gLength = length(genome);
	if(orientation == 'R') scanBegin = gLength - scanEnd;
	TMatch mR = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	TMatch temp = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	TDequeueValue fL(-1, mR);	
	fL.i2.gseqNo = gseqNo;
	mR.gseqNo = gseqNo;
	fL.i2.orientation = orientation;
	mR.orientation = orientation;
	
	double maxErrorRateL = (double)options.maxPrefixErrors/options.minMatchLen;
	double maxErrorRateR = (double)options.maxSuffixErrors/options.minMatchLen;

	double totalTime = sysTime();
	double totalTimeCombine = 0.0;
	double totalTimeLeftVerify = 0.0;
	double totalTimeLeftFilter = 0.0;
	double totalTimeRightVerify = 0.0;

	Pair<TGPos,TGPos> lastLeftMatch(0,0);
	// iterate all verification regions returned by SWIFT
	while (find(swiftFinderR, swiftPatternR, maxErrorRateR)) 
	{
		unsigned rseqNo = swiftPatternR.curSeqNo;
		TGPos rEndPos = endPosition(swiftFinderR) + scanBegin;	
		TGPos rBeginPos = beginPosition(swiftFinderR) + scanBegin;
	//	std::cout << "rEnd=" << rEndPos << "\t";	
		//TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;
		TRead const &read = readSet[rseqNo];
		unsigned readLength = length(readSet[rseqNo]);
		TSignedGPos extensionOffset = length(read)-2*options.minMatchLen;
		if(!options.hammingOnly)
			extensionOffset += static_cast<TGPos>(floor(options.errorRate*length(read)));
		// passed all valid regions
		if(!empty(readRegions) && (TSignedGPos)rBeginPos > (TSignedGPos)options.maxReadRegionsEnd)
			break;
		if(!empty(readRegions) &&
			((TSignedGPos)rEndPos < (TSignedGPos)regionStartPos(readRegions[rseqNo],readLength,options) || (TSignedGPos)rBeginPos > (TSignedGPos)regionEndPos(readRegions[rseqNo],readLength,options))) //parallelogram must lie within possible mapping region
			continue;
				
		// Check this again... 
		// remove out-of-window prefixes from fifo
		while(!empty(fifo) &&  front(fifo).i2.gEnd + extensionOffset + maxDistance < (TSignedGPos) rBeginPos)
		//raus:while (!empty(fifo) && front(fifo).i2.gBegin + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
		{
			popFront(fifo);
			++firstNo;
		}
		
		// Check this again... 
		// add within-window prefixes to fifo
		double vTimeBegin1 = sysTime();
		while (empty(fifo) || back(fifo).i2.gBegin + minDistance <= (TSignedGPos)rEndPos )
		{
			if (find(swiftFinderL, swiftPatternL, maxErrorRateL)) 
			{
				gPair = positionRange(swiftFinderL);
				gPair.i1 += scanBegin;
				gPair.i2 += scanBegin;
			//	std::cout << "lBegin=" << gPair.i1 << "\t";	
				if(!empty(readRegions) && (TSignedGPos)gPair.i1 > (TSignedGPos) options.maxReadRegionsEnd)
					break;
				if ((TSignedGPos)gPair.i2 + maxDistance + extensionOffset >= (TSignedGPos)rEndPos
					&& (empty(readRegions) ||
					((TSignedGPos)gPair.i2 > (TSignedGPos)regionStartPos(readRegions[swiftPatternL.curSeqNo],readLength,options)
					&& (TSignedGPos)gPair.i1 < (TSignedGPos)regionEndPos(readRegions[swiftPatternL.curSeqNo],readLength,options)))) //parallelogram must lie within possible mapping region
				//raus:if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
				{
					// link in
					fL.i1 = lastPotMatchNo[swiftPatternL.curSeqNo]; //link to last previous potential match
					lastPotMatchNo[swiftPatternL.curSeqNo] = lastNo++; //++ general counter and remember last pot match of curSeqNo
					
					fL.i2.rseqNo = swiftPatternL.curSeqNo | NOT_VERIFIED; // set first bit
					fL.i2.gBegin = gPair.i1;			
					fL.i2.gEnd = gPair.i2;
					pushBack(fifo, fL);
				}
#ifdef RAZERS_DEBUG
				::std::cout << "Discard potential match, out of window\n";
#endif				
			} 
			else
				break;
		}
		double vTimeEnd1 = sysTime();
		totalTimeLeftFilter += (vTimeEnd1 - vTimeBegin1);
		
		
		TDequeueIterator it;
		int64_t lastPositive = (int64_t)-1;

		TSize counter = 0;
		bool noMatchRight = false;
		bool notYetVerifiedRight = true;
		lastLeftMatch.i1 = 0;
		lastLeftMatch.i2 = 0;
		
		// walk through all potential prefix matches
		// if suffix is positive, verify prefixes (if not verfied already), mark as positive or negative
		for (int64_t i = lastPotMatchNo[rseqNo]; firstNo <= i; i = (*it).i1)
		{
			//CHECK HIER raus do suffix match verification only once 
			if(notYetVerifiedRight)
			{
				notYetVerifiedRight = false;
				TGPos minBeginPos = beginPosition(infix(swiftFinderR, genomeInf));
				if((int)minBeginPos - extensionOffset > 0)
					minBeginPos = minBeginPos - extensionOffset;
				else minBeginPos = 0;
				double vTimeBegin = sysTime();
				if (!matchVerify(mR, 
					infix(genome,minBeginPos,endPosition(infix(swiftFinderR, genomeInf))),
					rseqNo, 
					readSet,//readSetR, 
					forwardPatternsR,
					options, 
					TSwiftSpec(),
					LongestSuffix()))
				{	
					double vTimeEnd = sysTime();
					totalTimeRightVerify += (vTimeEnd - vTimeBegin);
					noMatchRight = true;
//					continue;
				}
				else 
				{
					double vTimeEnd = sysTime();
					totalTimeRightVerify += (vTimeEnd - vTimeBegin);
					if (lastRightMatch[rseqNo].i1 == (TGPos)mR.gBegin && lastRightMatch[rseqNo].i2 == (TGPos)mR.gEnd)
						noMatchRight = true;
					lastRightMatch[rseqNo].i1 = mR.gBegin;
					lastRightMatch[rseqNo].i2 = mR.gEnd;

					// check here if the match lies within the candidate region!!!
				}
			}



		//	std::cout << "i=" <<i << "\t";
			it = &value(fifo, i - firstNo);
			//CHECK HIER raus noMatchRight --> \FCberspringen korrekt?
			if (noMatchRight || (*it).i2.gBegin + minDistance > (TSignedGPos)rEndPos) 
			{ 	
				if (lastPositive == (int64_t)-1)
					lastPotMatchNo[rseqNo] = i;
				else
					value(fifo, lastPositive - firstNo).i1 = i;
				lastPositive = i;
				continue;
			}

			// verify left mate (equal seqNo), if not done already
			if ((*it).i2.rseqNo & NOT_VERIFIED)
			{
				TGPos maxEndPos = (*it).i2.gEnd + extensionOffset;
				if(maxEndPos > length(genome)) maxEndPos = length(genome);
				double vTimeBegin = sysTime();
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
					double vTimeEnd = sysTime();
					totalTimeLeftVerify += (vTimeEnd - vTimeBegin);
					(*it).i2.rseqNo &= ~NOT_VERIFIED; // has been verified positively // unset first bit
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
					double vTimeEnd = sysTime();
					totalTimeLeftVerify += (vTimeEnd - vTimeBegin);
					(*it).i2.rseqNo = ~NOT_VERIFIED;		// has been verified negatively
				}
				
			}
			
			// prefix match was found
			if ((*it).i2.rseqNo == rseqNo)
			{
				lastLeftMatch.i1 = (*it).i2.gBegin;
				lastLeftMatch.i2 = (*it).i2.gEnd;

				// dont shortcut too much
				if (lastPositive == (int64_t)-1 || i < lastPositive)
					lastPositive = i;

/*				// CHECK HIER rein
				// do suffix match verification 
				if(notYetVerifiedRight)
				{
					notYetVerifiedRight = false;
					TGPos minBeginPos = beginPosition(infix(swiftFinderR, genomeInf));
					if((int)minBeginPos - extensionOffset > 0)
						minBeginPos = minBeginPos - extensionOffset;
					else minBeginPos = 0;
					double vTimeBegin = sysTime();
					if (!matchVerify(mR, 
						infix(genome,minBeginPos,endPosition(infix(swiftFinderR, genomeInf))),
						rseqNo, 
						readSet,//readSetR, 
						forwardPatternsR,
						options, 
						TSwiftSpec(),
						LongestSuffix()))
					{	
						double vTimeEnd = sysTime();
						totalTimeRightVerify += (vTimeEnd - vTimeBegin);
						noMatchRight = true;
						continue;
					}
					else 
					{
						double vTimeEnd = sysTime();
						totalTimeRightVerify += (vTimeEnd - vTimeBegin);
						if (lastRightMatch[rseqNo].i1 == (TGPos)mR.gBegin && lastRightMatch[rseqNo].i2 == (TGPos)mR.gEnd)
							noMatchRight = true;
						lastRightMatch[rseqNo].i1 = mR.gBegin;
						lastRightMatch[rseqNo].i2 = mR.gEnd;
					}
				}
*/
				//else check if prefix and suffix match fit together
				if(!noMatchRight)
				{
					int outerDistance = mR.gEnd - (*it).i2.gBegin;
#ifdef RAZERS_DEBUG
					std::cout << "before\nmL.mScore =" << (*it).i2.mScore << "\n";
					std::cout << "mL.gBegin =" << (*it).i2.gBegin << "\n";
					std::cout << "mL.gEnd =" << (*it).i2.gEnd << "\n";
					std::cout << "mL.editDist =" << (*it).i2.editDist << "\n";
				
					std::cout << "mR.mScore =" << mR.mScore << "\n";
					std::cout << "mR.gBegin =" << mR.gBegin << "\n";
					std::cout << "mR.gEnd =" << mR.gEnd << "\n";
					std::cout << "mR.editDist =" << mR.editDist << "\n";
					std::cout << "outerDist =" << outerDistance << "\n";
#endif
					if (outerDistance < (int)(2 * options.minMatchLen)) //subtract suffix/prefix-errors in case of edit distance
						continue;
//					::std::cout << "outerDistance->" << outerDistance << std::endl;
					int outerDistanceError = length(readSet[rseqNo]) -(int)(outerDistance);
					if (outerDistanceError < 0) outerDistanceError = -outerDistanceError;
					if ((outerDistanceError > (int) options.maxGap) || // f\FCr edit + static_cast<TGPos>(floor(options.errorRate*length(read))
						(outerDistanceError < (int) options.minGap))   // f\FCr edit - static_cast<TGPos>(floor(options.errorRate*length(read))
						continue;

					TMatch mRtmp = mR;
					TMatch mLtmp = (*it).i2;
					if(!empty(readRegions) &&
						!isValidRegion(mLtmp,mRtmp,readLength,readRegions[rseqNo],options))
//						((TSignedGPos)mLtmp.gBegin < (TSignedGPos)readRegions[rseqNo].i2 ||
//					 	(TSignedGPos)mRtmp.gEnd > (TSignedGPos)readRegions[rseqNo].i3)) //match must lie within possible mapping region
						continue; 
					double t1 = sysTime();
					if(options.spec.DONT_VERIFY || 
					  !combineLeftRight(mRtmp,mLtmp,read,genome,options,orientation,TSwiftSpec()))
					{
						++options.FP;
						double t2 = sysTime();
						totalTimeCombine += (t2 - t1) ;
						continue;
					}
					else {
						double t2 = sysTime();
						totalTimeCombine += (t2 - t1) ;
						++options.TP;
					}
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
					//mLtmp.pairScore = mRtmp.pairScore = 0 - mLtmp.editDist - mRtmp.editDist;
					// score the whole match pair by # of matches bases - # mismatched bases (not entirely correct for edit distance)
#ifdef TRY_SCORES
					double identityScore = ( (100.00 - (100.00* (double)(mLtmp.editDist + mRtmp.editDist)/(mRtmp.mScore + mLtmp.mScore))) - 80.0 ) * 5.0;
					if(options._debugLevel > 1 )std::cout << "identityScore: " << identityScore << std::endl;
					
					double aliScore = mRtmp.mScore + mLtmp.mScore - 2* mLtmp.editDist - 2* mRtmp.editDist;
					aliScore = (double)(aliScore*100)/readLength;
					mLtmp.pairScore = mRtmp.pairScore = (int)(identityScore+aliScore)/2.0;
#else
					mLtmp.pairScore = mRtmp.pairScore = mRtmp.mScore + mLtmp.mScore - 2* mLtmp.editDist - 2* mRtmp.editDist;
#endif


					if(outerDistanceError != 0) 
					{ //subtract one if there is an indel in the middle (such that perfect matches are better than indel matches..)
						mLtmp.pairScore -= (int)((double)options.penaltyC * length(readSet[rseqNo])/100) ; 
						mRtmp.pairScore -= (int)((double)options.penaltyC * length(readSet[rseqNo])/100) ; 
					}
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

					if (!empty(readRegions) && options.anchored && options.outputFormat != 4)
					{
						if(readRegions[rseqNo].i2.i1 > 1) // match actually is on reverse strand
						{
							temp = mLtmp;
							mLtmp = mRtmp;
							mLtmp = temp;
							mLtmp.orientation = mRtmp.orientation = 'R';
						}
					}
					
					if (!options.spec.DONT_DUMP_RESULTS)
					{
						appendValue(matches, mLtmp, Generous());
						appendValue(matches, mRtmp, Generous());
					}
				}
			}
		}
		if (length(matches) > options.compactThresh)
		{
			typename Size<TMatches>::Type oldSize = length(matches);
//			maskDuplicates(matches);	// overlapping parallelograms cause duplicates for edit distance //TODO: implement!
			compactSplicedMatches(matches, cnts, options, false, swiftPatternL, swiftPatternR);
			options.compactThresh += (options.compactThresh >> 1);
			if (options._debugLevel >= 2)
				::std::cerr << '(' << oldSize - length(matches) << " matches removed)";
		}
			
		// short-cut negative matches
		if (lastPositive == (int64_t)-1)
			lastPotMatchNo[rseqNo] = (int64_t)-1;
		else
			value(fifo, lastPositive - firstNo).i1 = (int64_t)-1; // the first positive's link to previous is removed

		
	}//swiftFinderR
	double totalTimeEnd = sysTime();
	if(options._debugLevel > 1)
	{
		std::cout << "TotalTime:"<< totalTimeEnd - totalTime << std::endl;
		std::cout << "TotalTimeLeftFilter:" <<totalTimeLeftFilter << std::endl;
		std::cout << "TotalTimeLeftVerify:"<< totalTimeLeftVerify << std::endl;
		std::cout << "TotalTimeRightVerify:"<< totalTimeRightVerify << std::endl;
		std::cout << "TotalTimeCombine:"<< totalTimeCombine << std::endl;
	}


}//function




//////////////////////////////////////////////////////////////////////////////
// Find split read matches in many genome sequences (given as StringSet)
template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet,
	typename TCounts, 
	typename TSpec, 
	typename TShapeL,
	typename TShapeR,
	typename TSwiftSpec >
int mapSplicedReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	TReadRegions &			readRegions,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShapeL const &			shapeL,
	TShapeR const &			shapeR,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type							TRead;
	typedef Index<TReadSet, IndexQGram<TShapeL, TQGramIndexSpec> >	TIndexL;			// q-gram index
	typedef Index<TReadSet, IndexQGram<TShapeR, TQGramIndexSpec> >	TIndexR;			// q-gram index
	typedef Pattern<TIndexL, Swift<TSwiftSpec> >						TSwiftPatternL;	// filter
	typedef Pattern<TIndexR, Swift<TSwiftSpec> >						TSwiftPatternR;	// filter
	typedef Pattern<TRead, MyersUkkonen>							TMyersPattern;	// verifier

	unsigned readCount = length(readSet);

	// configure q-gram index
	TIndexL swiftIndexL(readSet, shapeL);
	cargo(swiftIndexL).abundanceCut = options.abundanceCut;
	cargo(swiftIndexL)._debugLevel = options._debugLevel;

	TIndexR swiftIndexR(readSet, shapeR);
	cargo(swiftIndexR).abundanceCut = options.abundanceCut;
	cargo(swiftIndexR)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPatternL swiftPatternL(swiftIndexL);
	swiftPatternL.params.minThreshold = options.threshold;
	swiftPatternL.params.tabooLength = options.tabooLength;

	TSwiftPatternR swiftPatternR(swiftIndexR);
	swiftPatternL.params.minThreshold = options.thresholdR;
	swiftPatternL.params.tabooLength = options.tabooLength;

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

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	CharString	id;
	
	// iterate over genome sequences
	SEQAN_PROTIMESTART(find_time);
	for(unsigned gseqNo = 0; gseqNo < length(genomeSet); ++gseqNo)
	{
		if (options.forward)
			mapSplicedReads(matches, genomeSet[gseqNo], gseqNo, readSet, readRegions, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'F', options);

		if (options.reverse)
		{
			reverseComplement(genomeSet[gseqNo]);
			mapSplicedReads(matches, genomeSet[gseqNo], gseqNo, readSet, readRegions, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'R', options);
			reverseComplement(genomeSet[gseqNo]);
		}

	}
	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);

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


//////////////////////////////////////////////////////////////////////////////
// Wrapper for different template specializations
template <typename TMatches, typename TReadSet, typename TCounts, typename TSpec>
int mapSplicedReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet &		readSet, 
	TReadRegions &			readRegions,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{

	Shape<Dna, SimpleShape>		ungappedL;
	Shape<Dna, OneGappedShape>	onegappedL;
	Shape<Dna, GenericShape>	gappedL;

	Shape<Dna, SimpleShape>		ungappedR;
	Shape<Dna, OneGappedShape>	onegappedR;
	Shape<Dna, GenericShape>	gappedR;

	// 2x3x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting combination of shape
		if (stringToShape(ungappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, ungappedL, ungappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, ungappedL, onegappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, ungappedL, gappedR, Swift<SwiftSemiGlobalHamming>());
		}

		if (stringToShape(onegappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, onegappedL, ungappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, onegappedL, onegappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, onegappedL, gappedR, Swift<SwiftSemiGlobalHamming>());
		}

		if (stringToShape(gappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, gappedL, ungappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, gappedL, onegappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, gappedL, gappedR, Swift<SwiftSemiGlobalHamming>());
		}
	} 
	else 
	{
		// select best-fitting combination of shape
		if (stringToShape(ungappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, ungappedL, ungappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, ungappedL, onegappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, ungappedL, gappedR, Swift<SwiftSemiGlobal>());
		}

		if (stringToShape(onegappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, onegappedL, ungappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, onegappedL, onegappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, onegappedL, gappedR, Swift<SwiftSemiGlobal>());
		}

		if (stringToShape(gappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, gappedL, ungappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, gappedL, onegappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, readRegions, cnts, options, gappedL, gappedR, Swift<SwiftSemiGlobal>());
		}
	}

	return RAZERS_INVALID_SHAPE;
}


template <typename TMatches, typename TGenomeSet, typename TReadSet, typename TCounts, typename TSpec>
int mapSplicedReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet &				readSet, 
	TReadRegions &			readRegions,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{


	Shape<Dna, SimpleShape>		ungappedL;
	Shape<Dna, OneGappedShape>	onegappedL;
	Shape<Dna, GenericShape>	gappedL;


	Shape<Dna, SimpleShape>		ungappedR;
	Shape<Dna, OneGappedShape>	onegappedR;
	Shape<Dna, GenericShape>	gappedR;

	// 2x3x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting combination of shape
		if (stringToShape(ungappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, ungappedL, ungappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, ungappedL, onegappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, ungappedL, gappedR, Swift<SwiftSemiGlobalHamming>());
		}

		if (stringToShape(onegappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, onegappedL, ungappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, onegappedL, onegappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, onegappedL, gappedR, Swift<SwiftSemiGlobalHamming>());
		}

		if (stringToShape(gappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, gappedL, ungappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, gappedL, onegappedR, Swift<SwiftSemiGlobalHamming>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, gappedL, gappedR, Swift<SwiftSemiGlobalHamming>());
		}
	} 
	else 
	{
		// select best-fitting combination of shape
		if (stringToShape(ungappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, ungappedL, ungappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, ungappedL, onegappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, ungappedL, gappedR, Swift<SwiftSemiGlobal>());
		}

		if (stringToShape(onegappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, onegappedL, ungappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, onegappedL, onegappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, onegappedL, gappedR, Swift<SwiftSemiGlobal>());
		}

		if (stringToShape(gappedL, options.shape))
		{
			if (stringToShape(ungappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, gappedL, ungappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(onegappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, gappedL, onegappedR, Swift<SwiftSemiGlobal>());
			if (stringToShape(gappedR, options.shapeR))
				return mapSplicedReads(matches, genomeSet, readSet, readRegions, cnts, options, gappedL, gappedR, Swift<SwiftSemiGlobal>());
		}
	}

	return RAZERS_INVALID_SHAPE;
}



}//namespace

#endif
