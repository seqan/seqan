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

namespace SEQAN_NAMESPACE_MAIN
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

	// ... to sort matches and remove duplicates with equal beginPos
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore>
	struct LessGPosRNo : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TAlignedReadQualityStore &qualStore;
		
		LessGPosRNo(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (TAlignedRead const &a, TAlignedRead const &b) const 
		{
			// contig
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// beginning position
			typename TAlignedRead::TPos ba = _min(a.beginPos, a.endPos);
			typename TAlignedRead::TPos bb = _min(b.beginPos, b.endPos);
			if (ba < bb) return true;
			if (ba > bb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;

			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// qualities
			SEQAN_ASSERT_NEQ(a.id, TAlignedRead::INVALID_ID);
			SEQAN_ASSERT_NEQ(b.id, TAlignedRead::INVALID_ID);
			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			return qa.score > qb.score;
		}
	};

    template <typename TAlignedReadStore, typename TReadNameStore, typename TAlignedReadQualityStore>
	struct LessGPosRName : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TReadNameStore & nameStore;
		TAlignedReadQualityStore &qualStore;
		
		LessGPosRName(TReadNameStore & nameStore, TAlignedReadQualityStore &_qualStore):
            nameStore(nameStore), qualStore(_qualStore) {}
		
		inline bool operator() (TAlignedRead const &a, TAlignedRead const &b) const 
		{
			// contig
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// beginning position
			typename TAlignedRead::TPos ba = _min(a.beginPos, a.endPos);
			typename TAlignedRead::TPos bb = _min(b.beginPos, b.endPos);
			if (ba < bb) return true;
			if (ba > bb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;

			// read name
			if (nameStore[a.readId] < nameStore[b.readId]) return true;
			if (nameStore[a.readId] > nameStore[b.readId]) return false;

			// qualities
			SEQAN_ASSERT_NEQ(a.id, TAlignedRead::INVALID_ID);
			SEQAN_ASSERT_NEQ(b.id, TAlignedRead::INVALID_ID);
			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			return qa.score > qb.score;
		}
	};

//////////////////////////////////////////////////////////////////////////////
// Determine error distribution
template <typename TErrDistr, typename TFragmentStore, typename TOptions>
inline unsigned
getErrorDistribution(
	TErrDistr &posError, 
	TFragmentStore &store, 
	TOptions &options)
{
	typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;
	typedef typename Value<TAlignedReadStore>::Type		TAlignedRead;
	typedef typename TFragmentStore::TContigPos			TContigPos;
	
	typename Iterator<TAlignedReadStore, Standard>::Type	it = begin(store.alignedReadStore, Standard());
	typename Iterator<TAlignedReadStore, Standard>::Type	itEnd = end(store.alignedReadStore, Standard());

	Dna5String genome;
	TContigPos left, right;
	unsigned unique = 0;
	
	for (; it != itEnd; ++it) 
	{
		if ((*it).contigId == TAlignedRead::INVALID_ID) continue;

		Dna5String const &read = store.readSeqStore[(*it).readId];
		left = (*it).beginPos;
		right = (*it).endPos;
		
		if (left < right)
			genome = infix(store.contigStore[(*it).contigId].seq, left, right);
		else
		{
			genome = infix(store.contigStore[(*it).contigId].seq, right, left);
			reverseComplement(genome);
		}
		for (unsigned i = 0; i < length(posError) && i < length(read); ++i)
			if ((options.compMask[ordValue(genome[i])] & options.compMask[ordValue(read[i])]) == 0)
				++posError[i]; 
		++unique;
	}
	return unique;
}

template <typename TErrDistr, typename TCount1, typename TCount2, typename TFragmentStore, typename TSpec>
inline unsigned
getErrorDistribution(
	TErrDistr &posError,
	TCount1 &insertions,
	TCount2 &deletions,
	TFragmentStore &store, 
	RazerSOptions<TSpec> &options)
{
	typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;
	typedef typename Value<TAlignedReadStore>::Type		TAlignedRead;
	typedef typename TFragmentStore::TContigPos			TContigPos;

	typedef Align<String<Dna5>, ArrayGaps> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TIter;

	//typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	typename Iterator<TAlignedReadStore, Standard>::Type	it = begin(store.alignedReadStore, Standard());
	typename Iterator<TAlignedReadStore, Standard>::Type	itEnd = end(store.alignedReadStore, Standard());

	Align<Dna5String, ArrayGaps> align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	unsigned unique = 0;
	for (; it != itEnd; ++it) 
	{
		if ((*it).contigId == TAlignedRead::INVALID_ID) continue;

		assignSource(row(align, 0), store.readSeqStore[(*it).readId]);
		TContigPos left = (*it).beginPos;
		TContigPos right = (*it).endPos;
		
		if (left < right)
			assignSource(row(align, 1), infix(store.contigStore[(*it).contigId].seq, left, right));
		else
		{
			assignSource(row(align, 1), infix(store.contigStore[(*it).contigId].seq, right, left));
			reverseComplement(source(row(align, 1)));
		}
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
		if ((*it).readId == readNo)
		{
			if(preEditDist > 1) continue;// || dist > options.errorRate * maxReadLength + 1) continue;
			int dist = (*it).seedEditDist - preEditDist;
			if(dist > maxSeedErrors) continue;
			if(dist < 0) ++cooc[0];
			else ++cooc[dist];
		}
		else
		{
			readNo = (*it).readId;
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
		std::cerr << "[mapping_count] ";
		for(unsigned j = 0; j < length(cooc); ++j)
			std::cerr << cooc[j] << " ";
		std::cerr << std::endl;
	}

}


template<typename TAlign, typename TString>
void
getCigarLine(TAlign & align, TString & cigar, TString & mutations)
{
	
	typedef typename Source<TAlign>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

	TAlignIterator ali_it0_stop = end(row(align, 0));
	TAlignIterator ali_it1_stop = end(row(align, 1));
	TAlignIterator ali_it0 = begin(row(align, 0));
	TAlignIterator ali_it1 = begin(row(align, 1));
	TStringIterator readBase = begin(source(row(align,0))); 
	// std::cout << "getting cigar line\nali0 len = " <<ali_it0_stop-ali_it0 << " \t ali1 len = "<<ali_it1_stop-ali_it1<<"\n";
	int readPos = 0;
	bool first = true;
	while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
	{
		int matched = 0;
		int inserted = 0;
		int deleted = 0;
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
}

//////////////////////////////////////////////////////////////////////////////
// Output matches
template <
	typename TFSSpec,
	typename TFSConfig,
	typename TGNoToFile,
	typename TSpec
>
void dumpMatches(
	FragmentStore<TFSSpec, TFSConfig> &store,		// forward/reverse matches
	StringSet<CharString> &genomeFileNameList,		// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	String<TGNoToFile> &gnoToFileMap,				// map to retrieve genome filename and sequence number within that file
	CharString readFName,							// read name (e.g. "reads.fa"), used for file/read naming
	CharString errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	typedef typename TFragmentStore::TReadNameStore                 TReadNameStore;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	//typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename GetValue<TAlignQualityStore>::Type				TQuality;
	//typedef typename Value<TGenomeSet>::Type						TGenome;
	//typedef typename TFragmentStore::TContigPos						TGPos;

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
	if (options.outputFormat == 3)
	{
		options.sortOrder = 1;		//  sort according to gPos
		options.positionFormat = 1;	// bases in file are numbered starting at 1
		options.dumpAlignment = false;
	}


	// error profile
	unsigned maxReadLength = 0;
	for (unsigned i = 0; i < length(store.readSeqStore); ++i)
		if (maxReadLength < length(store.readSeqStore[i]))
			maxReadLength = length(store.readSeqStore[i]);

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	TGenomeSet genomes;
	if (options.dumpAlignment || options.outputFormat == 4 || options.outputFormat == 5 || !empty(errorPrbFileName))
		if (!loadContigs(store, genomeFileNameList)) {
			std::cerr << "Failed to load genomes" << std::endl;
			options.dumpAlignment = false;
		}

	// how many 0's should be padded?
	int pzeros = 0;
	for (unsigned l = length(store.readSeqStore); l > 9; l = l / 10)
		++pzeros;

	int gzeros = 0;
	for (unsigned l = length(store.contigStore); l > 9; l = l / 10)
		++gzeros;

	// remove the directory prefix of readFName
	std::string _readName;
	assign(_readName, readFName);
	size_t lastPos = _readName.find_last_of('/') + 1;
	if (lastPos == _readName.npos) lastPos = _readName.find_last_of('\\') + 1;
	if (lastPos == _readName.npos) lastPos = 0;
	CharString readName = _readName.substr(lastPos);
	

	Align<String<Dna5>, ArrayGaps> align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)

	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	std::ofstream file;
	CharString fileName = options.output;
	if (empty(fileName))
	{
		fileName = readFName;
		append(fileName, ".result");
	}

	file.open(toCString(fileName), std::ios_base::out | std::ios_base::trunc);
	if (!file.is_open()) {
		std::cerr << "Failed to open output file" << std::endl;
		return;
	}
	
#ifndef RAZERS_DONTMASKDUPLICATES
    maskDuplicates(store);
#endif

	String< String<unsigned short> > stats;
	if (options.outputFormat >= 1 && options.outputFormat <= 3)
	{
		// match statistics
		unsigned maxErrors = (int)(options.errorRate * maxReadLength);
		if (maxErrors > 10) maxErrors = 10;
		resize(stats, maxErrors + 1);
		for (unsigned i = 0; i <= maxErrors; ++i)
			resize(stats[i], length(store.readStore), 0);
		countMatches(store, stats);
	}
	
	
	String<int> libSize;	// store outer library size for each pair match (indexed by pairMatchId)
	Nothing nothing;
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
	{
		compactPairMatches(store, options, nothing, nothing);
		if (options.outputFormat == 0)
			calculateInsertSizes(libSize, store);
	} else
#endif
		compactMatches(store, options, nothing);


	switch (options.sortOrder) {
		case 0: 
			sortAlignedReads(
				store.alignedReadStore, 
				LessRNoGPos<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));
			break;

		case 1:
			sortAlignedReads(
				store.alignedReadStore, 
				LessGPosRNo<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));
			break;
		case 2:
			sortAlignedReads(
				store.alignedReadStore, 
				LessRNameGPos<TAlignedReadStore, TReadNameStore, TAlignQualityStore>(store.readNameStore, store.alignQualityStore));
			break;
		case 3:
			sortAlignedReads(
				store.alignedReadStore, 
				LessGPosRName<TAlignedReadStore, TReadNameStore, TAlignQualityStore>(store.readNameStore, store.alignQualityStore));
			break;
	}
	
	TAlignedReadIter it = begin(store.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
	
	Dna5String gInf;
	char _sep_ = '\t';
	unsigned curReadId = 0;

	switch (options.outputFormat) 
	{
		case 0:	// Razer Format
//			_sep_ = ',';
			for(; it != itEnd; ++it) 
			{
				TQuality	qual = getValue(store.alignQualityStore, (*it).id);
				unsigned	readLen = length(store.readSeqStore[(*it).readId]);
				double		percId = 100.0 * (1.0 - (double)qual.errors / (double)readLen);

				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
					case 3:  // same as 0 if non-paired
						file << store.readNameStore[(*it).readId];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << std::setw(pzeros) << (*it).readId + 1;
						break;

					// 2..filename is the read sequence itself
					case 2:
						file << store.readSeqStore[(*it).readId];
				}

				file << _sep_ << options.positionFormat << _sep_ << readLen << _sep_ << (((*it).beginPos < (*it).endPos)? 'F': 'R') << _sep_;

				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << store.contigNameStore[(*it).contigId];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << gnoToFileMap[(*it).contigId].i1 << '#' << std::setw(gzeros) << gnoToFileMap[(*it).contigId].i2 + 1;
				}

				if ((*it).beginPos < (*it).endPos)
					file << _sep_ << ((*it).beginPos + options.positionFormat) << _sep_ << (*it).endPos << _sep_ << std::setprecision(5) << percId;
				else
					file << _sep_ << ((*it).endPos + options.positionFormat) << _sep_ << (*it).beginPos << _sep_ << std::setprecision(5) << percId;

				if ((*it).pairMatchId != TAlignedRead::INVALID_ID)
				{
					file << _sep_ << (*it).pairMatchId << _sep_ << (int)store.alignQualityStore[(*it).id].pairScore << _sep_;
					unsigned char no = getMateNo(store, (*it).readId);
					if (no == 0) file << libSize[(*it).pairMatchId];
					if (no == 1) file << -libSize[(*it).pairMatchId];
				}
				file << std::endl;

				if (options.dumpAlignment) 
				{
					assignSource(row(align, 0), store.readSeqStore[(*it).readId]);
					
					TContigPos left = (*it).beginPos;
					TContigPos right = (*it).endPos;
					
					if (left < right)
						assignSource(row(align, 1), infix(store.contigStore[(*it).contigId].seq, left, right));
					else
					{
						assignSource(row(align, 1), infix(store.contigStore[(*it).contigId].seq, right, left));
						reverseComplement(source(row(align, 1)));
					}
					globalAlignment(align, scoreType);
					dumpAlignment(file, align);
				}

			}
			break;


		case 1:	// Enhanced Fasta Format
			_sep_ = ',';
			for(unsigned matchReadNo = -1, matchReadCount = 0; it != itEnd; ++it) 
			{
				TQuality	qual = getValue(store.alignQualityStore, (*it).id);
				unsigned	readLen = length(store.readSeqStore[(*it).readId]);
				double		percId = 100.0 * (1.0 - (double)qual.errors / (double)readLen);

				if (matchReadNo != (*it).readId)
				{
					matchReadNo = (*it).readId;
					matchReadCount = 0;
				} else
					++matchReadCount;

				std::string fastaID;
				assign(fastaID, store.readNameStore[(*it).readId]);

				std::string id = fastaID;
				int fragId = (*it).readId;
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
						std::istringstream iss(fastaID.substr(pos + 3));
						iss >> id;
//						appendMatchId = false;
					}
					pos = fastaID.find("fragId=");
					if (pos != fastaID.npos) {
						std::istringstream iss(fastaID.substr(pos + 7));
						iss >> fragId;
					}
				}

				if ((*it).beginPos < (*it).endPos)
					// forward strand
					file << '>' << ((*it).beginPos + options.positionFormat) << _sep_ << (*it).endPos;
				else
					// reverse strand (switch begin and end)
					file << '>' << (*it).beginPos << _sep_ << ((*it).endPos + options.positionFormat);
					
				unsigned ambig = 0;
				for (unsigned i = 0; i <= qual.errors && i < length(stats); ++i)
					ambig += stats[i][(*it).readId];
				
				file << "[id=" << id;
				if (appendMatchId) file << "_" << matchReadCount;
				file << ",fragId=" << fragId;
				file << ",contigId=" << store.contigNameStore[(*it).contigId];
				file << ",errors=" << (unsigned)qual.errors << ",percId=" << std::setprecision(5) << percId;
				file << ",ambiguity=" << ambig << ']' << std::endl;

				file << store.readSeqStore[(*it).readId] << std::endl;
			}
			break;


		case 2:	// Eland Format
			_sep_ = '\t';
			for(unsigned readNo = 0; readNo < length(store.readSeqStore); ++readNo)
			{
				TQuality	qual = getValue(store.alignQualityStore, (*it).id);
				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
					case 3:  // same as 0 if non-paired
						file << '>' << store.readNameStore[readNo] << _sep_;
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << std::setw(pzeros) << readNo + 1  << _sep_;
						break;
				}

				if (it == itEnd || readNo < (*it).readId)
				{
					if (!empty(store.readSeqStore[readNo]))
						file << store.readSeqStore[readNo] << _sep_ << "NM" << _sep_ << '0' << _sep_ << '0' << _sep_ << '0' << std::endl;
					else
					{
						for (unsigned i = 0; i < maxReadLength; ++i)
							file << '.';
						file << _sep_ << "QC" << _sep_ << '0' << _sep_ << '0' << _sep_ << '0' << std::endl;
					}
				} else
				{
					SEQAN_ASSERT_EQ(readNo, (*it).readId);
					file << store.readSeqStore[readNo] << _sep_;
					unsigned bestMatches = 1;
					if ((unsigned)qual.errors < length(stats))
						bestMatches = stats[qual.errors][readNo];
					
					if (bestMatches == 0) file << '?';	// impossible
					if (bestMatches == 1) file << 'U';	// unique best match
					if (bestMatches >  1) file << 'R';	// non-unique best matches
					
					file << (unsigned)qual.errors << _sep_ << stats[0][readNo] << _sep_ << stats[1][readNo] << _sep_ << stats[2][readNo];
					
					if (bestMatches == 1)
					{
						file << _sep_;
						switch (options.genomeNaming)
						{
							// 0..filename is the read's Fasta id
							case 0:
								file << store.contigNameStore[(*it).contigId];
								break;

							// 1..filename is the read filename + seqNo
							case 1:
								file.fill('0');
								file << gnoToFileMap[(*it).contigId].i1 << '#' << std::setw(gzeros) << gnoToFileMap[(*it).contigId].i2 + 1;
						}
						
						if ((*it).beginPos < (*it).endPos)
							file << _sep_ << ((*it).beginPos + options.positionFormat) << _sep_ << 'F' << _sep_ << "..";
						else
							file << _sep_ << (*it).beginPos << _sep_ << 'R' << _sep_ << "..";

						if (qual.errors > 0 && options.dumpAlignment && options.hammingOnly) 
						{
							TContigPos left = (*it).beginPos;
							TContigPos right = (*it).endPos;
							
							if (left < right)
								gInf = infix(store.contigStore[(*it).contigId].seq, left, right);
							else
							{
								gInf = infix(store.contigStore[(*it).contigId].seq, right, left);
								reverseComplement(gInf);
							}
							for (unsigned i = 0; i < length(gInf); ++i)
								if ((options.compMask[ordValue(store.readSeqStore[readNo][i])] & 
									options.compMask[ordValue(gInf[i])]) == 0)
									file << _sep_ << i + 1 << gInf[i];
						}
					}
					file << std::endl;
					++it;
				}
			}
			break;
		case 3: // Gff:  printf "$chr $name_$format read $pos %ld . $dir . ID=$col[0]$unique$rest\n",$pos+$len-1;
			for (unsigned filecount = 0; filecount < length(genomeFileNameList); ++filecount)
			{
				TQuality	qual = getValue(store.alignQualityStore, (*it).id);

				// open genome file	
				std::ifstream gFile;
				gFile.open(toCString(genomeFileNameList[filecount]), std::ios_base::in | std::ios_base::binary);
				if (!gFile.is_open())
				{
					std::cerr << "Couldn't open genome file." << std::endl;
					break;
				}

				Dna5String	currGenome;
				
				// iterate over genome sequences
				for(; !streamEof(gFile); ++curReadId)
				{
					read(gFile, currGenome, Fasta());			// read Fasta sequence
					while(it != itEnd && (*it).contigId == curReadId)
					{
						unsigned currReadNo = (*it).readId;
						int unique = 1;
						unsigned bestMatches = 0;
						//would seedEditDist make more sense here?
//CHECKcnts					if ((unsigned)qual.errors < length(r))
//							bestMatches = stats[qual.errors][currReadNo];
							if (bestMatches == 0 && (unsigned)qual.errors < length(stats))
								bestMatches = stats[qual.errors][currReadNo];

						bool suboptimal = false;
						if ((unsigned)qual.errors > 0)
						{
							for(unsigned d = 0; d < (unsigned)qual.errors; ++d)
								if (stats[d][currReadNo]>0) suboptimal=true;
						}
						//std::cout << (stats[0][currReadNo] & 31) <<"<-dist "<< (stats[0][currReadNo] >> 5) <<"<-count\n";
					//	std::cout << "hier1\n";
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
//							continue; // TODO: output non-unique matches
						}
					//	std::cout << "hier2\n";
						unsigned readLen = length(store.readSeqStore[currReadNo]);
		
						switch (options.genomeNaming)
						{
							// 0..filename is the read's Fasta id
							case 0:
								file << store.contigNameStore[(*it).contigId] <<'\t';
								break;
							// 1..filename is the read filename + seqNo
							case 1:
								file.fill('0');
								file << gnoToFileMap[(*it).contigId].i1 << '#' << std::setw(gzeros) << gnoToFileMap[(*it).contigId].i2 + 1 << '\t';
								break;
						}
					//	std::cout << "hier3\n";
						//file <<  options.runID << "_razers\tread";
						file << "razers\tread\t";
						if ((*it).beginPos < (*it).endPos)
							file << ((*it).beginPos + options.positionFormat) << '\t' << (*it).endPos << '\t';
						else
							file << ((*it).endPos + options.positionFormat) << '\t' << (*it).beginPos << '\t';
		//				if ((*it).orientation == 'F')
		//					file << '\t' << ((*it).beginPos + options.positionFormat) << '\t' << (*it).endPos <<'\t';
		//				else
		//					file << '\t' << (*it).endPos << '\t'<<((*it).beginPos + options.positionFormat)<< '\t';
						double percId = 100.0 * (1.0 - (double)qual.errors / (double)readLen);
						file << percId << "\t";
					//	std::cout << "hier4\n";

						if ((*it).beginPos < (*it).endPos)
							file << '+' << '\t' << '.' <<'\t';
						else
							file << '-' << '\t' << '.' <<'\t';
		
						switch (options.readNaming)
						{
							// 0..filename is the read's Fasta id
							case 0:
							case 3:  // same as 0 if non-paired
								file << "ID=" <<store.readNameStore[currReadNo];
								break;
							
							// 1..filename is the read filename + seqNo
							case 1:
								file.fill('0');
								file << "ID=" << readName << '#' << std::setw(pzeros) << currReadNo + 1;
								break;
						}
					//	std::cout << "hier5\n";
						if(suboptimal) file << ";suboptimal";
						else 
						{
							if(unique==1) file << ";unique";
							if(unique==0) file << ";multi";
						}
						if (qual.errors > 0)
						{
							if (options.hammingOnly)
							{
								TContigPos left = (*it).beginPos;
								TContigPos right = (*it).endPos;
								
								if (left < right)
									gInf = infix(store.contigStore[(*it).contigId].seq, left, right);
								else
								{
									gInf = infix(store.contigStore[(*it).contigId].seq, right, left);
									reverseComplement(gInf);
								}
								bool first = true;
								file << ";cigar=" << length(store.readSeqStore[currReadNo]) << "M";
								file << ";mutations=";
								unsigned i = 0;
//								while ((*it).beginPos == 0 && i < length(store.readSeqStore[currReadNo])-length(gInf) )
//								{
//									if(first){ file << i + 1 << (Dna5)store.readSeqStore[currReadNo][i]; first = false;}
//									else file <<','<< i + 1 << (Dna5)store.readSeqStore[currReadNo][i];
//									++i;
//								}
								for (; i < length(gInf); ++i)
									if ((options.compMask[ordValue(store.readSeqStore[currReadNo][i])] & 
										options.compMask[ordValue(gInf[i])]) == 0)
									{
				//						if(first){ file << i + 1 << gInf[i]; first = false;}
				//						else file <<','<< i + 1 << gInf[i];
										if(first){ file << i + 1 << (Dna5)store.readSeqStore[currReadNo][i]; first = false;}
										else file <<','<< i + 1 << (Dna5)store.readSeqStore[currReadNo][i];
									}
//								while ((*it).endPos == length(currGenome) && i < length(store.readSeqStore[currReadNo]) )
//								{
//									if(first){ file << i + 1 << (Dna5)store.readSeqStore[currReadNo][i]; first = false;}
//									else file <<','<< i + 1 << (Dna5)store.readSeqStore[currReadNo][i];
//									++i;
//								}
								
							}
							else
							{
								assignSource(row(align, 0), store.readSeqStore[currReadNo]);
								TContigPos left = (*it).beginPos;
								TContigPos right = (*it).endPos;

								if (left < right)
									assignSource(row(align, 1), infix(currGenome, left, right));
								else
								{
									assignSource(row(align, 1), infix(currGenome, right, left));
									reverseComplement(source(row(align, 1)));
								}
								globalAlignment(align, scoreType);

								std::stringstream cigar, mutations;
								getCigarLine(align,cigar,mutations);
								file << ";cigar="<<cigar.str();

								if(length(mutations.str())>0)
									file << ";mutations=" << mutations.str();
								
							}
						}
						file << std::endl;
						++it;
					}
				}
				gFile.close();
				++filecount;
			}
			break;

		case 4: // Sam
            if (options.dontShrinkAlignments)
                convertMatchesToGlobalAlignment(store, scoreType, False());
            else
                convertMatchesToGlobalAlignment(store, scoreType, True());

            switch (options.sortOrder) {
                case 0: 
                    sortAlignedReads(
                        store.alignedReadStore, 
                        LessRNoGPos<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));
                    break;

                case 1:
                    sortAlignedReads(
                        store.alignedReadStore, 
                        LessGPosRNo<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));
                    break;
                case 2:
                    sortAlignedReads(
                            store.alignedReadStore, 
                            LessRNameGPos<TAlignedReadStore, TReadNameStore, TAlignQualityStore>(store.readNameStore, store.alignQualityStore));
                    break;
                case 3:
                    sortAlignedReads(
                            store.alignedReadStore, 
                            LessGPosRName<TAlignedReadStore, TReadNameStore, TAlignQualityStore>(store.readNameStore, store.alignQualityStore));
                    break;
            }
//            {
//                std::cerr << "AFTER CONVERSION" << std::endl << std::endl;
//                AlignedReadLayout layout;
//                layoutAlignment(layout, store);
//                typedef typename TFragmentStore::TContigStore                                        TContigStore;
//                typedef typename Value<TContigStore>::Type                                       TContig;
//                typedef Gaps<Nothing, AnchorGaps<typename TContig::TGapAnchors> >               TContigGaps;
//                TContigGaps contigGaps(store.contigStore[0].gaps);
//                unsigned l = positionSeqToGap(contigGaps, 5349500);
//                unsigned r = positionSeqToGap(contigGaps, 5349800);
//                printAlignment(std::cerr, Raw(), layout, store, 0, l, r, 0, 1000);
//            }
			write(file, store, Sam());
			break;

		case 5: // AFG
			convertMatchesToGlobalAlignment(store, scoreType, True());
			write(file, store, Amos());
			break;
	}

	file.close();

	// get empirical error distribution
	if (!empty(errorPrbFileName) && maxReadLength > 0)
	{
		file.open(toCString(errorPrbFileName), std::ios_base::out | std::ios_base::trunc);
		if (file.is_open())
		{
			String<long double> posError;
			unsigned unique = 0;
			unsigned insertions = 0;
			unsigned deletions = 0;
			resize(posError, maxReadLength, 0);
			
			if (options.hammingOnly)
				unique = getErrorDistribution(posError, store, options);
			else
			{
				unique = getErrorDistribution(posError, insertions, deletions, store, options);
				std::cerr << "insertProb: " << (double)insertions / ((double)length(posError) * (double)unique) << std::endl;
				std::cerr << "deleteProb: " << (double)deletions / ((double)length(posError) * (double)unique) << std::endl;
			}

			file << (double)posError[0] / (double)unique;
			for (unsigned i = 1; i < length(posError); ++i)
				file << '\t' << (double)posError[i] / (double)unique;
			file << std::endl;
			file.close();
		} else
			std::cerr << "Failed to open error distribution file" << std::endl;
	}

	options.timeDumpResults = SEQAN_PROTIMEDIFF(dump_time);

	if (options._debugLevel >= 1)
		std::cerr << "Dumping results took             \t" << options.timeDumpResults << " seconds" << std::endl;
}


}

#endif

