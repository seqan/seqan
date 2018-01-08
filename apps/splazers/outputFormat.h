 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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

#ifndef SEQAN_HEADER_OUTPUT_FORMAT_H
#define SEQAN_HEADER_OUTPUT_FORMAT_H

#include <iostream>
#include <fstream>
#include <sstream>

#include "razers.h"
#include <seqan/align.h>

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Quality-based score

	template <typename TQualityString = CharString>
	struct Quality;

	template <typename TValue, typename TQualityString>
	class Score<TValue, Quality<TQualityString> >
	{
	public:
		TValue data_match;
		TValue data_mismatch;
		TValue data_gap_extend;
		TValue data_gap_open;

		TQualityString const *data_qual;

	public:
		Score():
			data_match(0),
			data_mismatch(-1),
			data_gap_extend(-1),
			data_gap_open(-1),
			data_qual(NULL)
		{
		}
		Score(TValue _match, TValue _mismatch, TValue _gap):
			data_match(_match),
			data_mismatch(_mismatch),
			data_gap_extend(_gap),
			data_gap_open(_gap),
			data_qual(NULL)
		{
		}
		Score(TValue _match, TValue _mismatch, TValue _gap_extend, TValue _gap_open, TQualityString const &_qual):
			data_match(_match),
			data_mismatch(_mismatch),
			data_gap_extend(_gap_extend),
			data_gap_open(_gap_open),
			data_qual(&_qual)
		{
		}

		Score(Score const & other):
			data_match(other.data_match),
			data_mismatch(other.data_mismatch),
			data_gap_extend(other.data_gap_extend),
			data_gap_open(other.data_gap_open),
			data_qual(other.data_qual)
		{
		}
		~Score()
		{
		}

		Score & operator = (Score const & other)
		{
			data_match = other.data_match;
			data_mismatch = other.data_mismatch;
			data_gap_extend = other.data_gap_extend;
			data_gap_open = other.data_gap_open;
			data_qual = other.data_qual;
			return *this;
		}
	};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TQualityString, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, Quality<TQualityString> > const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &seq1,
	  TSeq2 const &seq2)
{
	if (seq1[pos1] != seq2[pos2])
		if (me.data_qual)
			return (*me.data_qual)[pos2];
		else
			return scoreMismatch(me);
	else
		return scoreMatch(me);
}


//////////////////////////////////////////////////////////////////////////////
// Less-operators ...

	// ... to sort matches and remove duplicates with equal gBegin
	template <typename TReadMatch>
	struct LessGPosRNo : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

//////////////////////////////////////////////////////////////////////////////
// Determine error distribution
template <typename TErrDistr, typename TMatches, typename TReads, typename TGenomes, typename TOptions>
inline unsigned
getErrorDistribution(
	TErrDistr &posError, 
	TMatches &matches, 
	TReads &reads, 
	TGenomes &genomes, 
	TOptions &options)
{
	typename Iterator<TMatches, Standard>::Type	it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type	itEnd = end(matches, Standard());

	Dna5String genome;
	unsigned unique = 0;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;

		Dna5String const &read = reads[(*it).rseqNo];
		genome = infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd);
		if ((*it).orientation == 'R')
			reverseComplement(genome);
		for (unsigned i = 0; i < length(posError) && i < length(read); ++i)
			if ((options.compMask[ordValue(genome[i])] & options.compMask[ordValue(read[i])]) == 0)
				++posError[i]; 
		++unique;
	}
	return unique;
}

template <typename TErrDistr, typename TCount1, typename TCount2, typename TMatches, typename TReads, typename TGenomes, typename TSpec>
inline unsigned
getErrorDistribution(
	TErrDistr &posError,
	TCount1 &insertions,
	TCount2 &deletions,
	TMatches &matches, 
	TReads &reads, 
	TGenomes &genomes, 
	RazerSOptions<TSpec> &options)
{
	typedef Align<String<Dna5>, ArrayGaps> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TIter;

	//typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	typename Iterator<TMatches, Standard>::Type	it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type	itEnd = end(matches, Standard());

	Align<Dna5String, ArrayGaps> align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	unsigned unique = 0;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;

		assignSource(row(align, 0), reads[(*it).rseqNo]);
		assignSource(row(align, 1), infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd));
		if ((*it).orientation == 'R')
			reverseComplement(source(row(align, 1)));
		globalAlignment(align, scoreType);
		
		TRow& row0 = row(align, 0);
		TRow& row1 = row(align, 1);
		
		TPosition begin = beginPosition(cols(align));
		TPosition end = endPosition(cols(align));
		
		TIter it0 = iter(row0, begin);
		TIter it1 = iter(row1, begin);
		TIter end0 = iter(row0, end);
		
		unsigned pos = 0;
		for (; it0 != end0 && pos < length(posError); ++it0, ++it1)
		{
			if (isGap(it0))
				++insertions;
			else
			{
				if (isGap(it1))
					++deletions;
				else
					if ((options.compMask[ordValue(getValue(it0))] & options.compMask[ordValue(getValue(it1))]) == 0)
						++posError[pos];
				++pos;
			}
		}
		++unique;
	}
	return unique;
}


//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TSource, typename TSpec, typename TPosition>
inline void
dumpAlignment(TFile & target, Align<TSource, TSpec> const & source, TPosition begin_, TPosition end_)
{
	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;

	TRowsPosition row_count = length(rows(source));

	// Print sequences
	for (TRowsPosition i = 0; i < row_count; ++i)
	{
		if (i == 0)
			target << "#Read:   ";
		else
			target << "#Genome: ";
		TRow& row_ = row(source, i);
		typedef typename Iterator<typename Row<TAlign>::Type const>::Type TIter;
		TIter begin1_ = iter(row_, begin_);
		TIter end1_ = iter(row_, end_);
		for (; begin1_ != end1_; ++begin1_) {
			if (isGap(begin1_))
			    target << gapValue<char>();
			else
			    target << *begin1_;
		}
		target << '\n';
	}
}

template<typename TMatches, typename TCounts, typename TOptions>
void
countCoocurrences(TMatches & matches, TCounts & cooc, TOptions & options)
{
	clear(cooc);
	int maxSeedErrors = (int)(options.errorRate * options.artSeedLength) + 1;
	resize(cooc,maxSeedErrors+1,0);
	for (int i = 0; i < maxSeedErrors+1; ++i)
		cooc[i] = 1;
	
	int count = 0;
	unsigned readNo = -1;
	int preEditDist = -1;
	typename Iterator<TMatches>::Type it = begin(matches,Standard());
	typename Iterator<TMatches>::Type itEnd = end(matches,Standard());
	
	for(; it != itEnd; ++it)
	{
		if ((*it).rseqNo == readNo)
		{
			if(preEditDist > 1) continue;// || dist > options.errorRate * maxReadLength + 1) continue;
			int dist = (*it).seedEditDist - preEditDist;
			if(dist > maxSeedErrors) continue;
			if(dist < 0) ++cooc[0];
			else ++cooc[dist];
		}
		else
		{
			readNo = (*it).rseqNo;
			preEditDist = (*it).seedEditDist;
			if(preEditDist <= 1) ++count;
		}
	}
	for (unsigned i = 0; i < length(cooc); ++i)
	{
		cooc[i] = (int)(-4.343 * log((double)cooc[i]/count) );
		if (cooc[i] < 0) cooc[i] = 0;
	}
	if(options._debugLevel > 1)
	{
		::std::cerr << "[mapping_count] ";
		for(unsigned j = 0; j < length(cooc); ++j)
			::std::cerr << cooc[j] << " ";
		::std::cerr << ::std::endl;
	}

}


template<typename TAlign, typename TString, typename TIter>
void
getCigarLine(TAlign & align, TString & cigar, TString & mutations, TIter ali_it0, TIter ali_it0_stop, TIter ali_it1, TIter ali_it1_stop)
{
	
	typedef typename Source<TAlign>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

// 	typedef typename Row<TAlign>::Type TRow;
// 	typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

	TStringIterator readBase = begin(source(row(align,0))); 

	int readPos = 0;
	bool first = true;
	int inserted = 0;
	int deleted = 0;
	while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
	{
		inserted = 0;
		deleted = 0;
		int matched = 0;
		while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1))
		{
			++readPos;
			if(*ali_it1 != *ali_it0)
			{
				if(first) first = false;
				else mutations << ",";
				mutations << readPos <<*readBase;
			}
			++readBase;
			++ali_it0;
			++ali_it1;
			++matched;
		}
		if(matched>0) cigar << matched<< "M" ;
		while(ali_it0!=ali_it0_stop && isGap(ali_it0))
		{
			++ali_it0;
			++ali_it1;
			++deleted;
		}
		if(deleted>0) cigar << deleted << "D";
		while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
		{
			++ali_it0;
			++ali_it1;
			++readPos;
			if(first) first = false;
			else mutations << ",";
			mutations << readPos << *readBase;
			++readBase;
			++inserted;
		}
		if(inserted>0) cigar << inserted << "I";
	}
	// end gaps can happen in split mapping
	while(ali_it0!=ali_it0_stop)
	{
		++ali_it0;
		++deleted;
	}
	if(deleted>0) cigar << deleted << "D";
	while(ali_it1 != ali_it1_stop)
	{
		++ali_it1;
		++readPos;
		if(first) first = false;
		else mutations << ",";
		mutations << readPos << *readBase;
		++inserted;
	}
	if(inserted>0) cigar << inserted << "I";
	
}

template<typename TAlign, typename TIter>
void
getCigarLine(TAlign & align, 
		String<Pair<char,int> > & cigar, 
		String<Pair<Dna5,int> > & mutations, 
		int offset,
		TIter ali_it0, TIter ali_it0_stop, TIter ali_it1, TIter ali_it1_stop)
{
	
	typedef typename Source<TAlign>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

// 	typedef typename Row<TAlign>::Type TRow;
// 	typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;
 
	TStringIterator readBase = begin(source(row(align,0))); 

	int readPos = 0;
	int inserted = 0;
	int deleted = 0;
 	while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
	{
		inserted = 0;
		deleted = 0;
		int matched = 0;
		while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1))
		{
			++readPos;
			if(*ali_it1 != *ali_it0)
				appendValue(mutations,Pair<Dna5,int>((Dna5)(*readBase),readPos + offset));
			++readBase;
			++ali_it0;
			++ali_it1;
			++matched;
		}
		if(matched>0) appendValue(cigar,Pair<char,int>('M',matched));
		while(ali_it0!=ali_it0_stop && isGap(ali_it0))
		{
			++ali_it0;
			++ali_it1;
			++deleted;
		}
		if(deleted>0) appendValue(cigar,Pair<char,int>('D',deleted));
		while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
		{
			++ali_it0;
			++ali_it1;
			++readPos;
			appendValue(mutations,Pair<Dna5,int>((Dna5)(*readBase),readPos + offset));
			++readBase;
			++inserted;
		}
		if(inserted>0) appendValue(cigar,Pair<char,int>('I',inserted));
	}
	// end gaps can happen in split mapping
	while(ali_it0!=ali_it0_stop)
	{
		++ali_it0;
		++deleted;
	}
	if(deleted>0) appendValue(cigar,Pair<char,int>('D',deleted));
	while(ali_it1 != ali_it1_stop)
	{
		++ali_it1;
		++readPos;
		appendValue(mutations,Pair<Dna5,int>((Dna5)(*readBase),readPos + offset));
		++inserted;
	}
	if(inserted>0) appendValue(cigar,Pair<char,int>('I',inserted + offset));
	
}



#ifdef RAZERS_DIRECT_MAQ_MAPPING
//////////////////////////////////////////////////////////////////////////////
// Assign mapping quality and remove suboptimal matches
template < typename TMatches, typename TReads, typename TCooc, typename TCounts, typename TSpec >
void assignMappingQuality(TMatches &matches, TReads & reads, TCooc & cooc, TCounts &cnts, RazerSOptions<TSpec> & options)
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;

	//matches are already sorted	
	//::std::sort(
	//	begin(matches, Standard()),
	//	end(matches, Standard()), 
	//	LessRNoMQ<TMatch>());
	
	
	int maxSeedErrors = (int)(options.errorRate*options.artSeedLength)+1;
	unsigned readNo = -1;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;

	int bestQualSum, secondBestQualSum;
	int secondBestDist = -1 ,secondBestMatches = -1;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		bool mappingQualityFound = false;
		int mappingQuality = 0;
		int qualTerm1,qualTerm2;

		readNo = (*it).rseqNo;
		bestQualSum = (*it).mScore;
		
		if(++it!=itEnd && (*it).rseqNo==readNo)
		{
			secondBestQualSum = (*it).mScore;
			secondBestDist = (*it).editDist;
			secondBestDist = (*it).editDist;
			secondBestMatches = cnts[1][readNo] >> 5;
//CHECKcnts		secondBestMatches = cnts[secondBestDist][readNo];
//			secondBestMatches = cnts[secondBestDist][readNo];
			(*it).orientation = '-';
		//	if(secondBestDist<=bestDist) unique=0;
		}
		else secondBestQualSum = -1000;
		--it; //it is now at best match of current rseqNo

		int bestDist = (*it).editDist;
		int kPrime = (*it).seedEditDist;
		if((bestQualSum==secondBestQualSum) || (kPrime>maxSeedErrors))
			mappingQualityFound = true;   //mq=0
		else{
			if(secondBestQualSum == -1000) qualTerm1 = 99;
			else
			{
				qualTerm1 = (int)(secondBestQualSum - bestQualSum - 4.343 * log((double)secondBestMatches));
				//if (secondBestKPrime - kPrime <= 1 && qualTerm1 > options.mutationRateQual) qualTerm1 = options.mutationRateQual; //TODO abchecken was mehr sinn macht
				if (secondBestDist - bestDist <= 1 && qualTerm1 > options.mutationRateQual) qualTerm1 = options.mutationRateQual;
			}
			float avgSeedQual = 0.0;
			if(!mappingQualityFound)
			{
				//TODO!!! generalize and adapt to razers lossrates
				// lossrate 0.42 => -10 log10(0.42) = 4
				int kPrimeLoss = 4; // options.kPrimeLoss; // bezieht sich auf 3 fehler in 28bp
				qualTerm2 = kPrimeLoss + cooc[maxSeedErrors-kPrime];
				
				for(unsigned j = 0; j<options.artSeedLength; ++j)
				{
					int q = getQualityValue(reads[readNo][j]);//(int)((unsigned char)(reads[readNo][j])>>3);
					if(q>options.mutationRateQual) q = options.mutationRateQual;
					avgSeedQual+=q;
				}
				avgSeedQual/=options.artSeedLength;
				//-10 log10(28-2) = 14;
				//generalize to -10 log10(artSeedLength - maxSeedErrors +1 ) // 14 fits for seedlength 28 to 32 with 2 errors
				if(avgSeedQual>14) qualTerm2 += (int)((maxSeedErrors-kPrime)*(avgSeedQual-14));
			}
		}
		if (!mappingQualityFound) mappingQuality = (qualTerm1<qualTerm2) ? qualTerm1:qualTerm2;
		if (mappingQuality < 0) mappingQuality = 0;
		(*it).mScore = mappingQuality;
		
		*dit = *it;
	//	if(secondBestQualSum != -1000) ++it;
		++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
}
#endif


//////////////////////////////////////////////////////////////////////////////
// Dump an alignment
template <typename TFile, typename TSource, typename TSpec>
inline void
dumpAlignment(TFile & target, Align<TSource, TSpec> const & source)
{
	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	TRowsPosition row_count = length(rows(source));
	TPosition begin_ = beginPosition(cols(source));
	TPosition end_ = endPosition(cols(source));
	
	// Print sequences
	for(TRowsPosition i=0;i<row_count;++i) {
		if (i == 0)
			streamPut(target, "#Read:   ");
		else
			streamPut(target, "#Genome: ");
		TRow& row_ = row(source, i);
		typedef typename Iterator<typename Row<TAlign>::Type const>::Type TIter;
		TIter begin1_ = iter(row_, begin_);
		TIter end1_ = iter(row_, end_);
		for (; begin1_ != end1_; ++begin1_) {
			if (isGap(begin1_))
			    streamPut(target, gapValue<char>());
			else
			    streamPut(target, *begin1_);
		}
		streamPut(target, '\n');
	}
}

template<typename TMatches, typename TReads, typename TReadNames, typename TSpec>
void
dumpUnmappedReads(TMatches				&matches, 
					TReads const		&reads, 
					TReadNames const	&readIDs,
					RazerSOptions<TSpec> &options)
{
	typedef typename Iterator<TMatches>::Type TIterator;

	::std::ofstream file;
	CharString fileName = options.outputUnmapped;
	if (empty(fileName))
		return;

	file.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) {
		::std::cerr << "Failed to open output file" << ::std::endl;
		return;
	}

	char headerChar = '>';
	if(options.readsWithQualities || options.fastaIdQual) 
		headerChar = '@';

	// sort according to read number
	::std::sort(begin(matches, Standard()), end(matches, Standard()), LessRNoGPos<TMatch>());

	TIterator mIt = begin(matches,Standard());
	TIterator mEnd = end(matches,Standard());

	unsigned readNo = 0;
	// output reads without any matches
	while(mIt != mEnd && readNo < length(reads))
	{
		while((*mIt).rseqNo > readNo)
		{
//			std::cout << "readNo=" << readNo << " + ";
			file << headerChar << readIDs[readNo] << std::endl << reads[readNo] << std::endl;
			if(options.readsWithQualities || options.fastaIdQual) 
			{
				file << '+' << std::endl;
				for(unsigned j=0;j<length(reads[readNo]);++j)
					file << (char)(getQualityValue(reads[readNo][j])+ 33);
				file << std::endl;
			}
			++readNo;
		}
		while((*mIt).rseqNo == readNo && mIt != mEnd)
		{
			++mIt;
		}
		++readNo;
	}
//	std::cout << std::endl << "readNoInBetween=" << readNo << " + " << std::endl;

	while(readNo < length(reads))
	{
//		std::cout << "readNo=" << readNo << " + ";
		file << headerChar << readIDs[readNo] << std::endl << reads[readNo] << std::endl;
		if(options.readsWithQualities || options.fastaIdQual) 
		{
			file << '+' << std::endl;
			for(unsigned j=0;j<length(reads[readNo]);++j)
				file << (char)(getQualityValue(reads[readNo][j])+ 33);
			file << std::endl;
		}
		++readNo;
	}
	return;


}


//////////////////////////////////////////////////////////////////////////////
// Output matches
template <
	typename TMatches,
	typename TGenomeNames,
	typename TReads,
	typename TReadRegions,
	typename TReadNames,
	typename TCounts,
	typename TSpec
>
void dumpMatches(
	TMatches &matches,							// forward/reverse matches
	TGenomeNames const &genomeIDs,				// Genome names (read from Fasta file, currently unused)
	StringSet<CharString> &genomeFileNameList,	// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > &gnoToFileMap, //map to retrieve genome filename and sequence number within that file
	TReads const &reads,
	TReadRegions const &readRegions,
	TCounts & stats,							// Match statistics (possibly empty)
	TReadNames const &readIDs,					// Read names (read from Fasta file, currently unused)
	CharString readFName,					// read name (e.g. "reads.fa"), used for file/read naming
	CharString errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type		TMatch;
	//typedef typename Value<TReads>::Type		TRead;
	//typedef typename Value<TGenomeSet>::Type	TGenome;
	typedef typename TMatch::TGPos				TGPos;

	//dump umapped reads in a separate file if filename was specified
	dumpUnmappedReads(matches, reads, readIDs, options);

	if (options.minMatchLen > 0 && options.outputFormat != 4) options.outputFormat = 33;

	if (options.outputFormat == 2)
	{
		options.sortOrder = 0;		// ... to count multiple matches
	}

	if (options.outputFormat == 2)
	{
		options.maxHits = 1;		// Eland outputs at most one match
		options.sortOrder = 0;		// read numbers are increasing
		options.positionFormat = 1;	// bases in file are numbered starting at 1
		options.dumpAlignment = options.hammingOnly;
	}
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if (options.maqMapping) options.outputFormat = 3;
	int maxSeedErrors = (int)(options.errorRate * options.artSeedLength); //without + 1 here, used to check whether match should be supported if noBelowIdentity option is switched on
#endif

	if (options.outputFormat == 3)
	{
		options.sortOrder = 1;		//  sort according to gPos
		options.positionFormat = 1;	// bases in file are numbered starting at 1
		//options.dumpAlignment = false;
	}
	if (options.outputFormat == 33 || options.outputFormat == 4)
	{
		options.positionFormat = 1;	// bases in file are numbered starting at 1
		options.dumpAlignment = false;
		options.sortOrder = 1;
	}


	// error profile
	unsigned maxReadLength = options.maxReadLength;

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	TGenomeSet genomes;
	if ((options.outputFormat == 0 && !options.hammingOnly) || options.dumpAlignment || !empty(errorPrbFileName) || options.outputFormat == 33 || options.outputFormat == 4)
		if (!loadGenomes(genomes, genomeFileNameList)) {
			::std::cerr << "Failed to load genomes" << ::std::endl;
			options.dumpAlignment = false;
		}

	// how many 0's should be padded?
	int pzeros = 0;
	for (unsigned l = length(reads); l > 9; l = l / 10)
		++pzeros;

	int gzeros = 0;
	for (unsigned l = length(genomes); l > 9; l = l / 10)
		++gzeros;

	// remove the directory prefix of readFName
	::std::string _readName;
	assign(_readName, readFName);
	size_t lastPos = _readName.find_last_of('/') + 1;
	if (lastPos == _readName.npos) lastPos = _readName.find_last_of('\\') + 1;
	if (lastPos == _readName.npos) lastPos = 0;
	CharString readName = _readName.substr(lastPos);
	

	typedef Align<String<Dna5>, ArrayGaps> TAlign;
	TAlign align, alignL, alignR;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
//	Score<int> scoreType(0,-1,-1,-1);

	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);
	resize(rows(alignL), 2);
	resize(rows(alignR), 2);

    VirtualStream<char, Output> file;
    bool success;
    if (!isEqual(options.output, "-"))
        success = open(file, toCString(options.output));
    else
        success = open(file, std::cout, Nothing());

    if (!success) {
        ::std::cerr << "Failed to open output file" << ::std::endl;
        return;
    }

	
	String<short> ambiStates;
	if(options.minMatchLen == 0)
		maskDuplicates(matches);
	if (options.outputFormat > 0
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	 && !options.maqMapping
#endif
	)
	{
		// match statistics
		unsigned maxErrors = (int)(options.errorRate * maxReadLength);
		//if (maxErrors > 10) maxErrors = 10;
		resize(stats, maxErrors + 1);
		for (unsigned i = 0; i <= maxErrors; ++i)
			resize(stats[i], length(reads), 0);
		countMatches(matches, stats);
	}



	Nothing nothing;
	unsigned currSeqNo = 0;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping)
	{
		String<int> cooc;
		compactMatches(matches, stats, options, true, nothing, false); //only best two matches per read are kept
		countCoocurrences(matches,cooc,options);	//coocurrence statistics are filled
		assignMappingQuality(matches,reads,cooc,stats,options);//mapping qualities are assigned and only one match per read is kept
	}
	else	 
#endif
	
#ifdef RAZERS_MICRO_RNA
	if(options.microRNA)purgeAmbiguousRnaMatches(matches,options);
	else
#endif
	{
		if(options.minMatchLen>0)
			compactAndCountSplicedMatches(matches, ambiStates, options, true);
		else
			compactMatches(matches, stats, options, true, nothing);
	}
        int jj = 0;
        if(!(options.minMatchLen > 0 ))
	  switch (options.sortOrder) {
		case 0:
			::std::sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessRNoGPos<TMatch>());
			break;

		case 1:
			::std::sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessGPosRNo<TMatch>());
			break;
			
	}
	
	typename Iterator<TMatches, Standard>::Type	it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type	itEnd = end(matches, Standard());
	
	
	Dna5String gInf, gInfL, gInfR;
	char _sep_ = '\t';
	unsigned viewPosReadFirst = 0;
	unsigned viewPosReadLast  = 0;

	switch (options.outputFormat) 
	{
		case 0:	// Razer Format
//			_sep_ = ',';
			for(; it != itEnd; ++it) 
			{
				unsigned	readLen = length(reads[(*it).rseqNo]);
				double		percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);
#ifdef RAZERS_MICRO_RNA
				percId = 100.0 * (1.0 - (double)(*it).editDist / (double) ((*it).mScore));
#endif
				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << readIDs[(*it).rseqNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << ::std::setw(pzeros) << (*it).rseqNo + 1;
						break;

					// 2..filename is the read sequence itself
					case 2:
						file << reads[(*it).rseqNo];
				}

				file << _sep_ << options.positionFormat << _sep_ << readLen << _sep_ << (*it).orientation << _sep_;

				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[(*it).gseqNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1;
				}
				
				// get alignment to dump or to fix end coordinates
				if (options.dumpAlignment || !options.hammingOnly)
				{
#ifdef RAZERS_MICRO_RNA
					if(options.microRNA)
						assignSource(row(align, 0), prefix(reads[(*it).rseqNo],(*it).mScore));
					else
#endif
					assignSource(row(align, 0), reads[(*it).rseqNo]);
					assignSource(row(align, 1), infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd));
#ifdef RAZERS_MATEPAIRS
					if ((*it).pairId != 0 && ((*it).rseqNo & 1))
						reverseComplement(source(row(align, 0)));
#endif
					if ((*it).orientation == 'R')
						reverseComplement(source(row(align, 1)));
				
					globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());

					// transform first and last read character to genomic positions
					viewPosReadFirst = toViewPosition(row(align, 0), 0);
					viewPosReadLast  = toViewPosition(row(align, 0), readLen - 1);
					unsigned genomePosReadFirst = toSourcePosition(row(align, 1), viewPosReadFirst);
					unsigned genomePosReadLast  = toSourcePosition(row(align, 1), viewPosReadLast);

					if ((*it).orientation == 'R')
					{
						(*it).gBegin = (*it).gEnd - (genomePosReadLast + 1);
						(*it).gEnd -= genomePosReadFirst;
					}
					else
					{
						(*it).gEnd = (*it).gBegin + (genomePosReadLast + 1);
						(*it).gBegin += genomePosReadFirst;
					}
				}

				file << _sep_ << ((*it).gBegin + options.positionFormat) << _sep_ << (*it).gEnd << _sep_ << ::std::setprecision(5) << percId;
#ifdef RAZERS_MICRO_RNA
				if(options.microRNA) file << _sep_ << (*it).mScore;
#endif

#ifdef RAZERS_MATEPAIRS
				if ((*it).pairId != 0)
					file << _sep_ << (*it).pairId << _sep_ << (*it).pairScore << _sep_ << (*it).mateDelta;
#endif
				file << ::std::endl;

				if (options.dumpAlignment)
					dumpAlignment(file, align, viewPosReadFirst, viewPosReadLast + 1);
			}
			break;


		case 1:	// Enhanced Fasta Format
			_sep_ = ',';
			for(unsigned matchReadNo = -1, matchReadCount = 0; it != itEnd; ++it) 
			{
				unsigned	readLen = length(reads[(*it).rseqNo]);
				double		percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

				if (matchReadNo != (*it).rseqNo)
				{
					matchReadNo = (*it).rseqNo;
					matchReadCount = 0;
				} else
					++matchReadCount;

				::std::string fastaID;
				assign(fastaID, readIDs[(*it).rseqNo]);

				::std::string id = fastaID;
				int fragId = (*it).rseqNo;
				bool appendMatchId = options.maxHits > 1;

				size_t left = fastaID.find_first_of('[');
				size_t right = fastaID.find_last_of(']');
				if (left != fastaID.npos && right != fastaID.npos && left < right) 
				{
					fastaID.erase(right);
					fastaID.erase(0, left + 1);
					replace(fastaID.begin(), fastaID.end(), ',', ' ');
					size_t pos = fastaID.find("id=");
					if (pos != fastaID.npos) {
						::std::istringstream iss(fastaID.substr(pos + 3));
						iss >> id;
//						appendMatchId = false;
					}
					pos = fastaID.find("fragId=");
					if (pos != fastaID.npos) {
						::std::istringstream iss(fastaID.substr(pos + 7));
						iss >> fragId;
					}
				}

				if ((*it).orientation == 'F')
					// forward strand
					file << '>' << ((*it).gBegin + options.positionFormat) << _sep_ << (*it).gEnd;
				else
					// reverse strand (switch begin and end)
					file << '>' << (*it).gEnd << _sep_ << ((*it).gBegin + options.positionFormat);
					
				unsigned ambig = 0;
				for (unsigned i = 0; i <= (*it).editDist && i < length(stats); ++i)
					ambig += stats[i][(*it).rseqNo];
				
				file << "[id=" << id;
				if (appendMatchId) file << "_" << matchReadCount;
				file << ",fragId=" << fragId;
				file << ",contigId=" << genomeIDs[(*it).gseqNo];
				file << ",errors=" << (*it).editDist << ",percId=" << ::std::setprecision(5) << percId;
				file << ",ambiguity=" << ambig << ']' << ::std::endl;

				file << reads[(*it).rseqNo] << ::std::endl;
			}
			break;


		case 2:	// Eland Format
			_sep_ = '\t';
			for(unsigned readNo = 0; readNo < length(reads); ++readNo)
			{
				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << '>' << readIDs[readNo] << _sep_;
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << ::std::setw(pzeros) << readNo + 1  << _sep_;
						break;
				}

				if (it == itEnd || readNo < (*it).rseqNo)
				{
					if (!empty(reads[readNo]))
						file << reads[readNo] << _sep_ << "NM" << _sep_ << '0' << _sep_ << '0' << _sep_ << '0' << ::std::endl;
					else
					{
						for (unsigned i = 0; i < maxReadLength; ++i)
							file << '.';
						file << _sep_ << "QC" << _sep_ << '0' << _sep_ << '0' << _sep_ << '0' << ::std::endl;
					}
				} else
				{
					file << reads[readNo] << _sep_;
					unsigned bestMatches = 1;
					if ((unsigned)(*it).editDist < length(stats))
						bestMatches = stats[(*it).editDist][readNo];
					
					if (bestMatches == 0) file << '?';	// impossible
					if (bestMatches == 1) file << 'U';	// unique best match
					if (bestMatches >  1) file << 'R';	// non-unique best matches
					
					file << (*it).editDist << _sep_ << stats[0][readNo] << _sep_ << stats[1][readNo] << _sep_ << stats[2][readNo];
					
					if (bestMatches == 1)
					{
						file << _sep_;
						switch (options.genomeNaming)
						{
							// 0..filename is the read's Fasta id
							case 0:
								file << genomeIDs[(*it).gseqNo];
								break;

							// 1..filename is the read filename + seqNo
							case 1:
								file.fill('0');
								file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1;
						}
						
						if ((*it).orientation == 'F')
							file << _sep_ << ((*it).gBegin + options.positionFormat) << _sep_ << 'F' << _sep_ << "..";
						else
							file << _sep_ << (*it).gEnd << _sep_ << 'R' << _sep_ << "..";

						if ((*it).editDist > 0 && options.dumpAlignment && options.hammingOnly) 
						{
							gInf = infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd);
							if ((*it).orientation == 'R')
								reverseComplement(gInf);
							for (unsigned i = 0; i < length(gInf); ++i)
								if ((options.compMask[ordValue(reads[readNo][i])] & 
									options.compMask[ordValue(gInf[i])]) == 0)
									file << _sep_ << i + 1 << gInf[i];
						}
					}
					file << ::std::endl;
					++it;
				}
			}
			break;
		case 3: // Gff:  printf "$chr $name_$format read $pos %ld . $dir . ID=$col[0]$unique$rest\n",$pos+$len-1;
			// NOTE(marehr): filecount+=2 might be a potential bug [https://github.com/seqan/seqan/issues/2165]
			// In revision 4dbf27b55 and before, filecount was incremented twice at the
			// end of the for loop, which caused a compiler warning (once in the body
			// and once in the iteration_expression of the for loop). We kept this
			// behaviour, because we have no active maintainer for this app.
			for (unsigned filecount = 0; filecount < length(genomeFileNameList); filecount+=2)
			{
				// open genome file	
				SeqFileIn gFile;
				if (!open(gFile, toCString(genomeFileNameList[filecount])))
				{
					std::cerr << "Couldn't open genome file." << std::endl;
					break;
				}

                CharString currId;
				Dna5String currGenome;
				
				// iterate over genome sequences
				for(; !atEnd(gFile); ++currSeqNo)
				{
					readRecord(currId, currGenome, gFile);			// read Fasta sequence
					while(it != itEnd && (*it).gseqNo == currSeqNo)
					{
#ifdef RAZERS_DIRECT_MAQ_MAPPING
						if(options.maqMapping && options.noBelowIdentity && (*it).seedEditDist > maxSeedErrors)
						{
							++it;
							continue;
						}
#endif

						unsigned currReadNo = (*it).rseqNo;
						int unique = 1;
						unsigned bestMatches = 0;
					
#ifdef RAZERS_DIRECT_MAQ_MAPPING
						if(options.maqMapping)
							bestMatches = stats[0][currReadNo] >> 5;
						else
#endif
						{
							if (bestMatches == 0 && (unsigned)(*it).editDist < length(stats))
							bestMatches = stats[(*it).editDist][currReadNo];

						}
						
						bool suboptimal = false;
						if (
#ifdef RAZERS_DIRECT_MAQ_MAPPING
							!options.maqMapping && 
#endif
							(unsigned)(*it).editDist > 0 )
						{
							for(unsigned d = 0; d < (unsigned)(*it).editDist; ++d)
								if (stats[d][currReadNo]>0) suboptimal=true;
						}
	
						if (bestMatches !=  1)
						{
							unique = 0;
							if(options.purgeAmbiguous)
							{
								++it;
								continue;
							}
							
//							if((*it).mScore > 0) std::cout << (*it).mScore << "<-non uniq but score > 0\n";
//							++it;
//							continue; // TODO: output non-unique matches (concerns maq mapping only)
						}
						unsigned readLen = length(reads[currReadNo]);

						String<Dna5Q> readInf = infix(reads[currReadNo],0,readLen);
						switch (options.genomeNaming)
						{
							// 0..filename is the read's Fasta id
							case 0:
								file << genomeIDs[(*it).gseqNo] <<'\t';
								break;
							// 1..filename is the read filename + seqNo
							case 1:
								file.fill('0');
								file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1 << '\t';
								break;
						}
						
						// get alignment to dump or to fix end coordinates
						if (!options.hammingOnly)
						{
							assignSource(row(align, 0), readInf);
							//std::cout << "begin=" << (*it).gBegin <<" end="<< (*it).gEnd << std::endl;
							assignSource(row(align, 1), infix(currGenome, (*it).gBegin, (*it).gEnd));
							if ((*it).orientation == 'R')
								reverseComplement(source(row(align, 1)));

							globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());
								
							//std::cout << align << std::endl; 
							//dumpAlignment(::std::cout, align);

							// transform first and last read character to genomic positions
							viewPosReadFirst  = toViewPosition(row(align, 0), 0);
							viewPosReadLast   = toViewPosition(row(align, 0), readLen - 1);
							TGPos genomePosReadFirst = toSourcePosition(row(align, 1), viewPosReadFirst);
							TGPos genomePosReadLast  = toSourcePosition(row(align, 1), viewPosReadLast);
							if ((*it).orientation == 'R')
							{
								// watch out at chromosome borders
								(*it).gBegin = (*it).gEnd - _min((TGPos)(genomePosReadLast + 1),(*it).gEnd); 
								(*it).gEnd -= genomePosReadFirst;
							}
							else
							{
								(*it).gEnd = (*it).gBegin + _min(genomePosReadLast + 1, static_cast<TGPos>(length(currGenome) - (*it).gBegin));
								(*it).gBegin += genomePosReadFirst;
							}
							
							clear(row(align,1));
							clear(row(align,0));
							assignSource(row(align, 0), readInf);
							//std::cout << "begin=" << (*it).gBegin <<" end="<< (*it).gEnd << std::endl;
							assignSource(row(align, 1), infix(currGenome, (*it).gBegin, (*it).gEnd));
							if ((*it).orientation == 'R')
								reverseComplement(source(row(align, 1)));

							globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());
								
							//std::cout << align << std::endl; 
							//dumpAlignment(::std::cout, align);

							// transform first and last read character to genomic positions
							viewPosReadFirst  = toViewPosition(row(align, 0), 0);
							viewPosReadLast   = toViewPosition(row(align, 0), readLen - 1);
							
							
						}
						
						//file <<  options.runID << "_razers\tread";
						file << "razers\tread";
						file << '\t' << ((*it).gBegin + options.positionFormat) << '\t' << (*it).gEnd << '\t';
						double percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);
#ifdef RAZERS_DIRECT_MAQ_MAPPING
						if(options.maqMapping)
							file << (*it).mScore << "\t";
						else
#endif
							file << percId << "\t";

						if ((*it).orientation == 'F')
							file << '+' << '\t' << '.' <<'\t';
						else
							file << '-' << '\t' << '.' <<'\t';
		
						switch (options.readNaming)
						{
							// 0..filename is the read's Fasta id
							case 0:
								file << "ID=" <<readIDs[currReadNo];
								break;
							
							// 1..filename is the read filename + seqNo
							case 1:
								file.fill('0');
								file << "ID=" << readName << '#' << ::std::setw(pzeros) << currReadNo + 1;
								break;
						}
						if(suboptimal) file << ";suboptimal";
						else 
						{
							if(unique==1) file << ";unique";
							if(unique==0) file << ";multi";
						}
						
					
						if ((*it).editDist > 0)
						{
							if (options.hammingOnly)
							{
								gInf = infix(currGenome, (*it).gBegin, (*it).gEnd);
								if ((*it).orientation == 'R')
									reverseComplement(gInf);
								bool first = true;
								file << ";cigar=" << readLen << "M";
								file << ";mutations=";
								for (unsigned i = 0; i < length(gInf); ++i)
									if ((options.compMask[ordValue(readInf[i])] & 
										options.compMask[ordValue(gInf[i])]) == 0)
									{
										if(first){ file << i + 1 << (Dna5)readInf[i]; first = false;}
										else file <<','<< i + 1 << (Dna5)readInf[i];
									}
								
							}
							else
							{
								std::stringstream cigar, mutations;
								typedef typename Row<TAlign>::Type TRow;
								typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;
								TAlignIterator ali_it0 = iter(row(align,0),viewPosReadFirst);
								TAlignIterator ali_it1 = iter(row(align,1),viewPosReadFirst);
								TAlignIterator ali_it0stop = iter(row(align,0),viewPosReadLast + 1);
								TAlignIterator ali_it1stop = iter(row(align,1),viewPosReadLast + 1);

								getCigarLine(align,cigar,mutations,ali_it0,ali_it0stop,ali_it1,ali_it1stop);
								file << ";cigar="<<cigar.str();
								if(length(mutations.str())>0)
									file << ";mutations=" << mutations.str();								
							}
						}
						
						if(
#ifdef RAZERS_DIRECT_MAQ_MAPPING
						options.maqMapping || 
#endif
						options.fastaIdQual || options.readsWithQualities)
						{
		//					file << ";read=";
		//					for(unsigned j=0;j<length(reads[currReadNo]);++j)
		//					{
		//						file << (Dna5)reads[currReadNo][j];
		//					}
							file << ";quality=";
							for(unsigned j=0;j<readLen;++j)
							{
								file << (char)(getQualityValue(readInf[j])+ 33);
							}
						}

						file << ::std::endl;
						if(options.dumpAlignment)
						{
							if((*it).editDist > 0 && !options.hammingOnly)
								file << align;
							else
							{
								clear(row(align,1));
								clear(row(align,0));
								assignSource(row(align, 0), readInf);
								assignSource(row(align, 1), infix(currGenome, (*it).gBegin, (*it).gEnd));
								if ((*it).orientation == 'R')
									reverseComplement(source(row(align, 1)));
								globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());
								file << align;

							}
						}
						++it;
					}
				}
				close(gFile);
			}
			break;
		case 4: // SAM output for splitMatches
			while(it != itEnd)// && (*it).gseqNo == currSeqNo)
			{
				unsigned currReadNo = (*it).rseqNo;
				Dna5String &currGenome = genomes[(*it).gseqNo];
				unsigned readLen = length(reads[currReadNo]);
				TMatch& mL = *it;
				++it;
				TMatch& mR = *it;
				unsigned readLenL = mL.mScore;
				unsigned readLenR = mR.mScore;
				String<Dna5Q> readInfL = infix(reads[currReadNo],0,readLenL);
				String<Dna5Q> readInfR = infix(reads[currReadNo],length(reads[currReadNo])-readLenR,length(reads[currReadNo]));
				int offsetL = 0;
				int offsetR = readLen - readLenR;
				if(mL.orientation == 'R')
				{
					offsetR = 0;
					offsetL = readLen - readLenL;
					readInfL = infix(reads[currReadNo],length(reads[currReadNo])-readLenL,length(reads[currReadNo]));
					readInfR = infix(reads[currReadNo],0,readLenR);
				}

				// get alignment to dump and get cigar string
				if (!options.hammingOnly)
				{
					assignSource(row(alignL, 0), readInfL);
					assignSource(row(alignL, 1), infix(currGenome, mL.gBegin, mL.gEnd));
					if (mL.orientation == 'R')
						reverseComplement(source(row(alignL, 1)));

					globalAlignment(alignL, scoreType, AlignConfig<false,false,false,false>(), Gotoh());

					assignSource(row(alignR, 0), readInfR);
					assignSource(row(alignR, 1), infix(currGenome, mR.gBegin, mR.gEnd));
					if (mR.orientation == 'R')
						reverseComplement(source(row(alignR, 1)));

					globalAlignment(alignR, scoreType, AlignConfig<false,false,false,false>(), Gotoh());
				}

				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << "" <<readIDs[currReadNo];
						break;
					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << "" << readName << '#' << ::std::setw(pzeros) << currReadNo + 1;
						break;
				}
				int samFlag = 0;
				if(ambiStates[jj] == 2) samFlag |= 0x100; //suboptimal
				if (mR.orientation == 'R') samFlag |= 0x10;
				if(!empty(readRegions) && options.anchored)
				{
					samFlag |= 1;
					if(readRegions[currReadNo].i2.i1 > 1) samFlag |= 0x10;
					else samFlag |= 0x20;
					if((readRegions[currReadNo].i2.i1 & 1) == 1) samFlag |= 0x80;
					else samFlag |= 0x40;
				}
				file << '\t' << samFlag << '\t';

				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[(*it).gseqNo] <<'\t';
						break;
					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1 << '\t';
						break;
				}
				
				file << (mL.gBegin + 1) << '\t' ;

				//double percId = 100.0 * (1.0 - (double)(mL.editDist + mR.editDist) /(double)(mL.mScore+mR.mScore));

				if(ambiStates[jj] == 0) file << 255 << '\t'; //unique
				if(ambiStates[jj] > 0) file << 0 << '\t'; //multi or suboptimal
				++jj;
					
				SEQAN_ASSERT(mL.pairId == mR.pairId); 
		
				int indelLen = mR.mScore + mL.mScore - readLen;
				if(mL.gEnd != mR.gBegin) indelLen = mR.gBegin - mL.gEnd;
				
				String<Pair<Dna5,int> > mutationsStr, mutStrMid, mutStrL, mutStrR;
				String<Pair<char,int> > cigarStr, cigarMid, cigarL, cigarR;
				if(indelLen < 0) 
				{
					int offset = readLenL;
					if(mL.orientation == 'R') 
						offset = readLenR;
					String<Dna5Q> readInfMid = infix(reads[currReadNo],offset,offset - indelLen);
					for (int i = 0; i < -indelLen; ++i)
						appendValue(mutStrMid, Pair<Dna5,int>(readInfMid[i],i+offset+1));
					appendValue(cigarMid, Pair<char,int>('I',-indelLen));
					
				}
				if(indelLen > 0) appendValue(cigarMid, Pair<char,int>('N',indelLen));
						
				if (options.hammingOnly)
				{
					if(indelLen == 0)
						appendValue(cigarL, Pair<char,int>('M',readLen));
					else
					{
						appendValue(cigarL, Pair<char,int>('M',readLenL));
						appendValue(cigarR, Pair<char,int>('M',readLenR));
					}

					//make left part of mutation string
					if (mL.editDist > 0)
					{
						gInfL = infix(currGenome, mL.gBegin, mL.gEnd);
						if (mL.orientation == 'R')
							reverseComplement(gInfL);
						for (unsigned i = 0; i < length(gInfL); ++i)
							if ((options.compMask[ordValue(readInfL[i])] & options.compMask[ordValue(gInfL[i])]) == 0)
								appendValue(mutStrL, Pair<Dna5,int>(readInfL[i],i+offsetL+1));
					}
							
					//make right part of mutation string
					if (mR.editDist > 0)
					{
						gInfR = infix(currGenome, mR.gBegin, mR.gEnd);
						if (mR.orientation == 'R')
							reverseComplement(gInfR);
						for (unsigned i = 0; i < length(gInfR); ++i)
							if ((options.compMask[ordValue(readInfR[i])] & options.compMask[ordValue(gInfR[i])]) == 0)
								appendValue(mutStrR, Pair<Dna5,int>(readInfR[i],i+offsetR+1));
					}
				}
				else
				{
					typedef typename Row<TAlign>::Type TRow;
					typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

					//left business
					TAlignIterator aliL_it0 = begin(row(alignL,0));
					TAlignIterator aliL_it1 = begin(row(alignL,1));
					TAlignIterator aliL_it0stop = end(row(alignL,0));
					TAlignIterator aliL_it1stop = end(row(alignL,1));
					getCigarLine(alignL,cigarL,mutStrL,offsetL,aliL_it0,aliL_it0stop,aliL_it1,aliL_it1stop);
						
					//right business
					TAlignIterator aliR_it0 = begin(row(alignR,0));
					TAlignIterator aliR_it1 = begin(row(alignR,1));
					TAlignIterator aliR_it0stop = end(row(alignR,0));
					TAlignIterator aliR_it1stop = end(row(alignR,1));
					getCigarLine(alignR,cigarR,mutStrR,offsetR,aliR_it0,aliR_it0stop,aliR_it1,aliR_it1stop);
				}			
				//now plug together the parts: cigar
				if (mL.orientation == 'F')
				{
					append(cigarStr,cigarL);
					append(cigarStr,cigarMid);
					append(cigarStr,cigarR);
					append(mutationsStr,mutStrL);
					append(mutationsStr,mutStrMid);
					append(mutationsStr,mutStrR);
					
					// print cigar string
					for (unsigned i = 0; i < length(cigarStr); ++i)
						file << cigarStr[i].i2 << cigarStr[i].i1;
				}
				else
				{
					append(cigarStr,cigarR);
					append(cigarStr,cigarMid);
					append(cigarStr,cigarL);
					append(mutationsStr,mutStrR);
					append(mutationsStr,mutStrMid);
					append(mutationsStr,mutStrL);
					
					// print cigar string
					for (unsigned i = length(cigarStr); i != 0; )
					{
						--i;
						file << cigarStr[i].i2 << cigarStr[i].i1;
					}
				}

				file << '\t';

				if(!empty(readRegions) && options.anchored)
				{
				
					switch (options.genomeNaming)
					{
						// 0..filename is the read's Fasta id
						case 0:
							file << genomeIDs[(*it).gseqNo] <<'\t';
							break;
						// 1..filename is the read filename + seqNo
						case 1:
							file.fill('0');
							file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1 << '\t';
							break;
					}
				
					file << readRegions[currReadNo].i2.i2 << '\t' ;
					//int libraryError;
					int templateLen = 0;
					if(readRegions[currReadNo].i2.i1 < 2) 
						templateLen = readRegions[currReadNo].i2.i2 + readLen -(int) mL.gBegin;
						//libraryError = (int) mR.gEnd + readRegions[currReadNo].i2.i2 - options.libraryLength + 2*readLen;
					else 
						templateLen = (int) mR.gEnd - readRegions[currReadNo].i2.i2;
						//libraryError = (int) mL.gBegin - readRegions[currReadNo].i2.i2 + options.libraryLength - readLen ;
					file << templateLen << '\t';
				}
				else file << "*\t*\t*\t";

                if (mL.orientation == 'F')
				    file << reads[currReadNo] << "\t*\t";
                else
				    file << reverseComplementString(reads[currReadNo]) << "\t*\t";
				// quality
				if(options.fastaIdQual)
				{
					for(unsigned j=0;j<length(reads[currReadNo]);++j)
					{
						file << (char)(getQualityValue(reads[currReadNo][j])+ 33);
					}
				}
				file << "AS:i:" << mR.pairScore;

				file << ::std::endl;
				++it;
			}
			break;
		case 33: // special Gff for split reads
			while(it != itEnd)// && (*it).gseqNo == currSeqNo)
			{
				unsigned currReadNo = (*it).rseqNo;
				Dna5String &currGenome = genomes[(*it).gseqNo];
				unsigned readLen = length(reads[currReadNo]);
				TMatch& mL = *it;
				++it;
				TMatch& mR = *it;
				unsigned readLenL = mL.mScore;
				unsigned readLenR = mR.mScore;
				SEQAN_ASSERT_LEQ(mL.mScore+mR.mScore,(int)readLen);
				String<Dna5Q> readInfL = infix(reads[currReadNo],0,readLenL);
				String<Dna5Q> readInfR = infix(reads[currReadNo],length(reads[currReadNo])-readLenR,length(reads[currReadNo]));
				int offsetL = 0;
				int offsetR = readLen - readLenR;
				if(mL.orientation == 'R')
				{
					offsetR = 0;
					offsetL = readLen - readLenL;
					readInfL = infix(reads[currReadNo],length(reads[currReadNo])-readLenL,length(reads[currReadNo]));
					readInfR = infix(reads[currReadNo],0,readLenR);
				}
				
				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[(*it).gseqNo] <<'\t';
						break;
					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1 << '\t';
						break;
				}
				
				// get alignment to dump and get cigar string
				if (!options.hammingOnly)
				{
					assignSource(row(alignL, 0), readInfL);
					assignSource(row(alignL, 1), infix(currGenome, mL.gBegin, mL.gEnd));
					if (mL.orientation == 'R')
						reverseComplement(source(row(alignL, 1)));

					globalAlignment(alignL, scoreType, AlignConfig<false,false,false,false>(), Gotoh());
						
					assignSource(row(alignR, 0), readInfR);
					assignSource(row(alignR, 1), infix(currGenome, mR.gBegin, mR.gEnd));
					if (mR.orientation == 'R')
						reverseComplement(source(row(alignR, 1)));

					globalAlignment(alignR, scoreType, AlignConfig<false,false,false,false>(), Gotoh());
				}
						
				//file <<  options.runID << "_razers\tread";
				file << "razers\tread";
				file << '\t' << (mL.gBegin + options.positionFormat) << '\t' << mR.gEnd << '\t';
//?//			double percId = 100.0 * (1.0 + (double)(mL.pairScore-mL.mScore-mR.mScore) / (double)readLen);
				double percId = 100.0 * (1.0 - (double)(mL.editDist + mR.editDist) /(double)(mL.mScore+mR.mScore));
				file << percId << "\t";
//				if(!empty(readRegions) && options.anchored)
//				{
//					if(readRegions[currReadNo].i2 < 0) 
//						file << '-' << '\t' << '.' <<'\t';
//					else
//						file << '+' << '\t' << '.' <<'\t';
//				}			
//				else
				{
					if (mL.orientation == 'F')
						file << '+' << '\t' << '.' <<'\t';
					else
						file << '-' << '\t' << '.' <<'\t';
				}
		
				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << "ID=" <<readIDs[currReadNo];
						break;
					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << "ID=" << readName << '#' << ::std::setw(pzeros) << currReadNo + 1;
						break;
				}
				if(ambiStates[jj] == 0) file << ";unique";
				if(ambiStates[jj] == 1) file << ";multi";
				if(ambiStates[jj] == 2) file << ";suboptimal";
				++jj;
					
				SEQAN_ASSERT(mL.pairId == mR.pairId); 
		
				int indelLen = mR.mScore + mL.mScore - readLen;
				if(mL.gEnd != mR.gBegin) indelLen = mR.gBegin - mL.gEnd;
				
				file << ";pairScore=" << (unsigned int) mR.pairScore;
				if(!empty(readRegions) && options.anchored)
				{
					if(readRegions[currReadNo].i2.i1 < 2) 
						file << ";libraryError=" << (int) readRegions[currReadNo].i2.i2 - mR.gEnd - (int) options.libraryLength + 2*readLen;
                    else {
                        int tmp = (int)mL.gBegin + readLen - (int)options.libraryLength - readRegions[currReadNo].i2.i2;
						file << ";libraryError=" <<  tmp ;
                    }
					
				}

                if(mL.gBegin > mR.gBegin) std::cout << "falsch\n";
//				if(mR.traceExtension >= abs(indelLen))
//					file << ";traceExt=" << (unsigned int) mR.traceExtension;

				String<Pair<Dna5,int> > mutationsStr, mutStrMid, mutStrL, mutStrR;
				String<Pair<char,int> > cigarStr, cigarMid, cigarL, cigarR;
				if(indelLen < 0) 
				{
					int offset = readLenL;
					if(mL.orientation == 'R') 
						offset = readLenR;
					String<Dna5Q> readInfMid = infix(reads[currReadNo],offset,offset - indelLen);
					for (int i = 0; i < -indelLen; ++i)
						appendValue(mutStrMid, Pair<Dna5,int>(readInfMid[i],i+offset+1));
					appendValue(cigarMid, Pair<char,int>('I',-indelLen));
					
				}
				if(indelLen > 0) appendValue(cigarMid, Pair<char,int>('D',indelLen));
				if(indelLen != 0) 
					file << ";split";
						
				if (options.hammingOnly)
				{
					if(indelLen == 0)
						appendValue(cigarL, Pair<char,int>('M',readLen));
					else
					{
						appendValue(cigarL, Pair<char,int>('M',readLenL));
						appendValue(cigarR, Pair<char,int>('M',readLenR));
					}

					//make left part of mutation string
					if (mL.editDist > 0)
					{
						gInfL = infix(currGenome, mL.gBegin, mL.gEnd);
						if (mL.orientation == 'R')
							reverseComplement(gInfL);
						for (unsigned i = 0; i < length(gInfL); ++i)
							if ((options.compMask[ordValue(readInfL[i])] & options.compMask[ordValue(gInfL[i])]) == 0)
								appendValue(mutStrL, Pair<Dna5,int>(readInfL[i],i+offsetL+1));
					}
							
					//make right part of mutation string
					if (mR.editDist > 0)
					{
						gInfR = infix(currGenome, mR.gBegin, mR.gEnd);
						if (mR.orientation == 'R')
							reverseComplement(gInfR);
						for (unsigned i = 0; i < length(gInfR); ++i)
							if ((options.compMask[ordValue(readInfR[i])] & options.compMask[ordValue(gInfR[i])]) == 0)
								appendValue(mutStrR, Pair<Dna5,int>(readInfR[i],i+offsetR+1));
					}
				}
				else
				{
					typedef typename Row<TAlign>::Type TRow;
					typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

					//left business
					TAlignIterator aliL_it0 = begin(row(alignL,0));
					TAlignIterator aliL_it1 = begin(row(alignL,1));
					TAlignIterator aliL_it0stop = end(row(alignL,0));
					TAlignIterator aliL_it1stop = end(row(alignL,1));
					getCigarLine(alignL,cigarL,mutStrL,offsetL,aliL_it0,aliL_it0stop,aliL_it1,aliL_it1stop);
						
					//right business
					TAlignIterator aliR_it0 = begin(row(alignR,0));
					TAlignIterator aliR_it1 = begin(row(alignR,1));
					TAlignIterator aliR_it0stop = end(row(alignR,0));
					TAlignIterator aliR_it1stop = end(row(alignR,1));
					getCigarLine(alignR,cigarR,mutStrR,offsetR,aliR_it0,aliR_it0stop,aliR_it1,aliR_it1stop);
				}			
				//now plug together the parts: cigar
				if(mL.orientation=='F')
				{
					append(cigarStr,cigarL);
					append(cigarStr,cigarMid);
					append(cigarStr,cigarR);
					append(mutationsStr,mutStrL);
					append(mutationsStr,mutStrMid);
					append(mutationsStr,mutStrR);
				}
				else
				{
					append(cigarStr,cigarR);
					append(cigarStr,cigarMid);
					append(cigarStr,cigarL);
					append(mutationsStr,mutStrR);
					append(mutationsStr,mutStrMid);
					append(mutationsStr,mutStrL);
				}
				bool first = true;
				file << ";cigar=";// << cigarStr.str();
				for (unsigned i = 0; i < length(cigarStr); ++i)
					file << cigarStr[i].i2 << cigarStr[i].i1;
 					
				if(length(mutationsStr)>0)file << ";mutations=";
				for (unsigned i = 0; i < length(mutationsStr); ++i)
				{
					if(first)
					{
						first = false;
						file << mutationsStr[i].i2 << mutationsStr[i].i1;
					}
					else
						file << "," << mutationsStr[i].i2 << mutationsStr[i].i1;
				}
						
				if(options.fastaIdQual)
				{
					file << ";quality=";
					for(unsigned j=0;j<length(reads[currReadNo]);++j)
					{
						file << (char)(getQualityValue(reads[currReadNo][j])+ 33);
					}
				}
				file << ::std::endl;
				++it;
			}
			break;


	}

	close(file);

	// get empirical error distribution
	if (!empty(errorPrbFileName) && maxReadLength > 0)
	{
        std::ofstream file;
		file.open(toCString(errorPrbFileName), ::std::ios_base::out | ::std::ios_base::trunc);
		if (file.is_open())
		{
			String<long double> posError;
			unsigned unique = 0;
			unsigned insertions = 0;
			unsigned deletions = 0;
			resize(posError, maxReadLength, 0);
			
			if (options.hammingOnly)
				unique = getErrorDistribution(posError, matches, reads, genomes, options);
			else
			{
				unique = getErrorDistribution(posError, insertions, deletions, matches, reads, genomes, options);
				::std::cerr << "insertProb: " << (double)insertions / ((double)length(posError) * (double)unique) << ::std::endl;
				::std::cerr << "deleteProb: " << (double)deletions / ((double)length(posError) * (double)unique) << ::std::endl;
			}

			file << (double)posError[0] / (double)unique;
			for (unsigned i = 1; i < length(posError); ++i)
				file << '\t' << (double)posError[i] / (double)unique;
			file << ::std::endl;
			file.close();
		} else
			::std::cerr << "Failed to open error distribution file" << ::std::endl;
	}

	options.timeDumpResults = SEQAN_PROTIMEDIFF(dump_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Dumping results took             \t" << options.timeDumpResults << " seconds" << ::std::endl;
}


}

#endif

