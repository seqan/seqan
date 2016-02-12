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




//////////////////////////////////////////////////////////////////////////////
// Output matches
template <
	typename TMatches,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TCounts,
	typename TSpec
>
void dumpMatches(
	TMatches &matches,							// forward/reverse matches
	TGenomeNames const &genomeIDs,				// Genome names (read from Fasta file, currently unused)
	StringSet<CharString> &genomeFileNameList,	// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > &gnoToFileMap, //map to retrieve genome filename and sequence number within that file
	TReads  &reads,
	TCounts & stats,							// Match statistics (possibly empty)
	TReadNames const &readIDs,					// Read names (read from Fasta file, currently unused)
	CharString readFName,                       // read name (e.g. "reads.fa"), used for file/read naming
	CharString errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type		TMatch;
	//typedef typename Value<TReads>::Type		TRead;
	//typedef typename Value<TGenomeSet>::Type	TGenome;
	//typedef typename TMatch::TGPos				TGPos;

	//if (options.outputFormat == 4)
	//{
	//	options.sortOrder = 2;
	//}


	// error profile
	unsigned maxReadLength = 0;
	for (unsigned i = 0; i < length(reads); ++i)
		if (maxReadLength < length(reads[i]))
			maxReadLength = length(reads[i]);

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	TGenomeSet genomes;
	if ((options.outputFormat == 0 && !options.hammingOnly) || options.dumpAlignment || !empty(errorPrbFileName))
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
	TAlign align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)

	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

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

	maskDuplicates(matches, options);

	resize(stats, length(reads),0);
	purgeAmbiguousRnaMatches(matches,stats,options);


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
        case 2:
            ::std::sort(
                begin(matches, Standard()),
                end(matches, Standard()),
                LessRNoEdistHLen<TMatch>());
            break;
	}
	
	typename Iterator<TMatches, Standard>::Type	it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type	itEnd = end(matches, Standard());
	
	
	Dna5String gInf;
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
				percId = 100.0 * (1.0 - (double)(*it).editDist / (double) ((*it).mScore));
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
					if(options.microRNA)
						assignSource(row(align, 0), prefix(reads[(*it).rseqNo],(*it).mScore));
					else
    					assignSource(row(align, 0), reads[(*it).rseqNo]);
    					assignSource(row(align, 1), infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd));
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
				if(options.microRNA) file << _sep_ << (*it).mScore;

				file << ::std::endl;

				if (options.dumpAlignment)
					dumpAlignment(file, align, viewPosReadFirst, viewPosReadLast + 1);
			}
			break;
		case 4:	// SAM Format
//			_sep_ = ',';
			for(; it != itEnd; ++it) 
			{
				int     	readLen = (int) length(reads[(*it).rseqNo]);
				double		percId = 100.0 * (1.0 - (double)(*it).editDist / (double) ((*it).mScore));
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
                int samFlag = 0;
                if ((*it).orientation == 'R') samFlag |= 16;
                // add sth with stats? 
                // if this is not the single best alignment (or through mapping quality)
				file << '\t' << samFlag << '\t';

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
                int mapQ = 255;
                if(stats[(*it).rseqNo] > 1) mapQ = 0; 
                //if() sth with stats
				file << '\t' << (*it).gBegin+1 << '\t' << mapQ << '\t';
                if((*it).orientation == 'R' && (*it).mScore < readLen) 
                    file << readLen - (*it).mScore << "S";
                file << (*it).mScore << "M";
                if((*it).orientation == 'F' && (*it).mScore < readLen) 
                    file << readLen - (*it).mScore << "S";
                file << "\t*\t0\t0\t";
                if((*it).orientation == 'R')
                    reverseComplement(reads[(*it).rseqNo]);
				
                file << reads[(*it).rseqNo] << '\t';
				for(unsigned j = 0; j < length(reads[(*it).rseqNo]); ++j)
				{
					//file << (char) ((int)((unsigned char)store.readSeqStore[currReadNo][j] >> 3) + 33);
					file << (char)(getQualityValue(reads[(*it).rseqNo][j])+ 33);
				}
                if((*it).orientation == 'R')
                    reverseComplement(reads[(*it).rseqNo]);
				
                file << '\t' << "AS:i:" << (int)percId << std::endl;
	
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

