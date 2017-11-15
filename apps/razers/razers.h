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

#ifndef SEQAN_HEADER_RAZERS_H
#define SEQAN_HEADER_RAZERS_H

#include <iostream>
#include <fstream>

#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Default options

#ifdef RAZERS_OPENADDRESSING
	typedef OpenAddressing	TQGramIndexSpec;
#else
	typedef Default	TQGramIndexSpec;
#endif

	template < bool DONT_VERIFY_ = false, bool DONT_DUMP_RESULTS_ = false >
	struct RazerSSpec 
	{
		enum { DONT_VERIFY = DONT_VERIFY_ };				// omit verifying potential matches
		enum { DONT_DUMP_RESULTS = DONT_DUMP_RESULTS_ };	// omit dumping results
		enum { DUMP_VERIFICATION_TASKS = false };
	};

    struct VerificationTask_
    {
        unsigned genomeBegin;
        unsigned genomeEnd;
        unsigned genomeId;
        unsigned readId;
        unsigned leftClip;
    };

	template < typename TSpec = RazerSSpec<> >
	struct RazerSOptions
	{
	// main options
		TSpec		spec;
		bool		forward;			// compute forward oriented read matches
		bool		reverse;			// compute reverse oriented read matches
		double		errorRate;			// Criteria 1 threshold
		unsigned	maxHits;			// output at most maxHits many matches
		unsigned	distanceRange;		// output only the best, second best, ..., distanceRange best matches
		bool		purgeAmbiguous;		// true..remove reads with more than maxHits best matches, false..keep them
		CharString	output;				// name of result file
		int			_debugLevel;		// level of verbosity
		bool		printVersion;		// print version number
		bool		hammingOnly;		// no indels
		int			trimLength;			// if >0, cut reads to #trimLength characters
		
	// output format options
		unsigned	outputFormat;		// 0..Razer format
										// 1..enhanced Fasta
										// 2..ELAND format
		bool		dumpAlignment;		// compute and dump the match alignments in the result files
		unsigned	genomeNaming;		// 0..use Fasta id
										// 1..enumerate reads beginning with 1
		unsigned	readNaming;			// 0..use Fasta id
										// 1..enumerate reads beginning with 1
										// 2..use the read sequence (only for short reads!)
		unsigned	sortOrder;			// 0..sort keys: 1. read number, 2. genome position
										// 1..           1. genome pos50ition, 2. read number
		int			positionFormat;		// 0..gap space
										// 1..position space
		const char	*runID;				// runID needed for gff output

	// filtration parameters
		::std::string shape;			// shape (e.g. 11111111111)
		int			threshold;			// threshold
		int			tabooLength;		// taboo length
		int			repeatLength;		// repeat length threshold
		double		abundanceCut;		// abundance threshold

	// mate-pair parameters
		int			libraryLength;		// offset between two mates
		int			libraryError;		// offset tolerance
		unsigned	nextMatePairId;		// use this id for the next mate-pair

	// verification parameters
		bool		matchN;				// false..N is always a mismatch, true..N matches with all
		unsigned char compMask[5];

	// statistics
		int64_t		FP;					// false positives (threshold reached, no match)
		int64_t		TP;					// true positives (threshold reached, match)
        double      timeCompactMatches;     // time for compacting reads
        double      timeMaskDuplicates; // time spent masking duplicates
		double		timeLoadFiles;		// time for loading input files
		double		timeMapReads;		// time for mapping reads
		double		timeDumpResults;	// time for dumping the results
        double      timeBuildQGramIndex;  // time for q-gram index building.
		
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		bool		maqMapping;
		int		maxMismatchQualSum;
		int		mutationRateQual;
		int		absMaxQualSumErrors;
		unsigned	artSeedLength;
		bool		noBelowIdentity;
#endif

#ifdef RAZERS_MICRO_RNA
		bool		microRNA;
		unsigned	rnaSeedLength;
		bool 		exactSeed;
#endif			

		bool		lowMemory;		// set maximum shape weight to 13 to limit size of q-gram index
		bool		fastaIdQual;		// hidden option for special fasta+quality format we use
		int			minClippedLen;

	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh

#ifdef RAZERS_SPLICED
		unsigned        minMatchLen;
		unsigned        maxDistance;
		unsigned        minDistance;
#endif
		
	// multi-threading

#ifdef RAZERS_PARALLEL
		typedef ::tbb::spin_mutex	TMutex;

		TMutex		*patternMutex;
		TMutex		optionsMutex;
		TMutex		matchMutex;
#endif

		String<VerificationTask_, MMap<> > verifications;

        SeqFileIn   readFile;           // left read's SeqFile (we have to keep it open and store it here to stream it only once)

		RazerSOptions() 
		{
			forward = true;
			reverse = true;
			errorRate = 0.08;
			maxHits = 100;
			distanceRange = 0;	// disabled
			purgeAmbiguous = false;
			output = "";
			_debugLevel = 0;
			printVersion = false;
			hammingOnly = false;
			trimLength = 0;
			
			outputFormat = 0;
			dumpAlignment = false;
			genomeNaming = 0;
			readNaming = 0;
			sortOrder = 0;
			positionFormat = 0;
			runID = "s"; 	//

			matchN = false;

			shape = "11111111111";
			threshold = 1;
			tabooLength = 1;
			repeatLength = 1000;
			abundanceCut = 1;

			libraryLength = 220;
			libraryError = 50;
			nextMatePairId = 1;
			
			for (unsigned i = 0; i < 4; ++i)
				compMask[i] = 1 << i;
			compMask[4] = 0;

			compactThresh = 1024;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			maqMapping = false;
			maxMismatchQualSum = 70;
			mutationRateQual = 30;
			artSeedLength = 28;	// the "artificial" seed length that is used for mapping quality assignment 
						// (28bp is maq default)
			absMaxQualSumErrors = 100; // maximum for sum of mism qualities in total readlength
			noBelowIdentity = false;
#endif

#ifdef RAZERS_MICRO_RNA
			microRNA = false;
			rnaSeedLength = 16;
			exactSeed = true;
#endif			
#ifdef RAZERS_SPLICED
			minMatchLen = 0;
			maxDistance = 4000;
			minDistance = 0;
#endif

			lowMemory = false;		// set maximum shape weight to 13 to limit size of q-gram index
			fastaIdQual = false;
			minClippedLen = 0;

            if (TSpec::DUMP_VERIFICATION_TASKS)
            {
                open(verifications, "verification_tasks.bin", OPEN_WRONLY|OPEN_CREATE);
            }
		}
	};

struct MicroRNA{};	

#ifdef RAZERS_MICRO_RNA
#define RAZERS_EXTENDED_MATCH
#endif

#ifdef RAZERS_DIRECT_MAQ_MAPPING 
#define RAZERS_EXTENDED_MATCH
#endif

#ifdef RAZERS_SPLICED
#define RAZERS_EXTENDED_MATCH
#endif

//////////////////////////////////////////////////////////////////////////////
// Typedefs

	// definition of a Read match
	template <typename TGPos_>
	struct ReadMatch 
	{
		typedef typename MakeSigned_<TGPos_>::Type TGPos;

		unsigned		gseqNo;			// genome seqNo
		unsigned		rseqNo;			// read seqNo
		TGPos			gBegin;			// begin position of the match in the genome
		TGPos			gEnd;			// end position of the match in the genome
#ifdef RAZERS_MATEPAIRS
		unsigned		pairId;			// unique id for the two mate-pair matches (0 if unpaired)
		int				mateDelta:24;	// outer coordinate delta to the other mate 
		int				pairScore:8;	// combined score of both mates
#endif
		unsigned short	editDist;		// Levenshtein distance
#ifdef RAZERS_EXTENDED_MATCH
		short	 		mScore;
		short			seedEditDist:8;	
#endif
#ifdef RAZERS_SPLICED
		short			gSeedLen:8;// used as gMinMatchLen to store genomic end position of seed match
#endif
		char			orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};
	
	enum RAZERS_ERROR 
	{
		RAZERS_INVALID_OPTIONS = -1,
		RAZERS_READS_FAILED    = -2,
		RAZERS_GENOME_FAILED   = -3,
		RAZERS_INVALID_SHAPE   = -4
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions

	typedef Dna5String									TGenome;
	typedef StringSet<TGenome>							TGenomeSet;
//	typedef Dna5String									TRead;
	typedef String<Dna5Q>								TRead;
#ifdef RAZERS_CONCATREADS
	typedef StringSet<TRead, Owner<ConcatDirect<> > >	TReadSet;
#else
	typedef StringSet<TRead>							TReadSet;
#endif

	typedef ReadMatch<Difference<TGenome>::Type>		TMatch;		// a single match
	typedef String<TMatch/*, MMap<>*/ >					TMatches;	// array of matches


	template <typename TSpec>
	struct Cargo< Index<TReadSet, TSpec> > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};
	
	template <typename TAnyReadSet, typename TSpec>
	struct Cargo< Index<TAnyReadSet, TSpec> > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};
//////////////////////////////////////////////////////////////////////////////
// Memory tuning

#ifdef RAZERS_MEMOPT

	template <typename TSpec>
	struct SAValue< Index<TReadSet, TSpec> > {
		typedef Pair<
			unsigned,				
			unsigned,
			BitPacked<24, 8>	// max. 16M reads of length < 256
		> Type;
	};
	//454 
//	template <typename TSpec>
//	struct SAValue< Index<TReadSet, TSpec> > {
//		typedef Pair<
//			unsigned,				
//			unsigned,
//			BitPacked<22, 10>	// max. 4M reads of length < 1024
//		> Type;
//	};
	
#else

	template <typename TSpec>
	struct SAValue< Index<TReadSet, TSpec> > {
		typedef Pair<
			unsigned,			// many reads
			unsigned,			// of arbitrary length
			Pack
		> Type;
	};

#endif

	template <>
	struct Size<Dna5String>
	{
		typedef unsigned Type;
	};

	template <typename TShape, typename TSpec>
	struct Size< Index<TReadSet, IndexQGram<TShape, TSpec> > >
	{
		typedef unsigned Type;
	};
	

#ifdef RAZERS_PRUNE_QGRAM_INDEX

	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <typename TShape, typename TSpec>
	inline bool _qgramDisableBuckets(Index<TReadSet, IndexQGram<TShape, TSpec> > &index) 
	{
		typedef Index<TReadSet, IndexQGram<TShape, TSpec>	>	TReadIndex;
		typedef typename Fibre<TReadIndex, QGramDir>::Type		TDir;
		typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
		typedef typename Value<TDir>::Type						TSize;

		TDir &dir    = indexDir(index);
		bool result  = false;
		unsigned counter = 0;
		TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
		if (thresh < 100) thresh = 100;

		TDirIterator it = begin(dir, Standard());
		TDirIterator itEnd = end(dir, Standard());
		for (; it != itEnd; ++it)
			if (*it > thresh) 
			{
				*it = (TSize)-1;
				result = true;
				++counter;
			}

		if (counter > 0 && cargo(index)._debugLevel >= 1)
			::std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

		return result;
	}

#endif


template<typename TIdString, typename TSeqString, typename TQString, typename TOptions>
bool
_clipReads(TIdString & fastaID, TSeqString & seq, TQString & qual, TOptions & options)
{
	typedef typename Value<TIdString>::Type TChar;
        
	int tagStart = length(fastaID);
	int clipFront = -1;
	int clipBack = -1;
	for(unsigned i = 0; i < length(fastaID); ++i)
	{
		TChar c = fastaID[i];
		if(c == 'c')
		{
			tagStart = i-1;
			if(infix(fastaID,i,i+5)=="clip=")
			{
				i += 5;
				String<TChar> clipStr = "";
				while(i < length(fastaID) && _parseIsDigit(fastaID[i]))
				{
					append(clipStr,fastaID[i]);
					++i;
				}
				std::istringstream istr(toCString(clipStr));
				istr >> clipFront;
				if(i < length(fastaID) && fastaID[i] == ',') ++i;
				clipStr = "";
				while(i < length(fastaID) && _parseIsDigit(fastaID[i]))
				{
					append(clipStr,fastaID[i]);
					++i;
				}
				std::istringstream istr2(toCString(clipStr));
				istr2 >> clipBack;
				break;
			}
		}
	}
	if(clipFront<0 && clipBack<0) return true;	//no clip tag found
	if(clipFront<0) clipFront = 0;
	if(clipBack<0) clipBack = 0;
	resize(fastaID,tagStart);	// only resize fastaID, as it might contain the quality string

	if(options.minClippedLen == 0) // clipping tag was found, but clipping option is turned off
		return true;

	if(((int)length(seq)-clipBack-clipFront) < options.minClippedLen)  // sequence is too short after clipping
		return false;

//	int newLen = (int)length(seq)-clipBack-clipFront;
	
	//clip
	if(length(qual) == 0) // meaning that quality is in fasta header
	{
		// first adapt the fasta header
		TIdString tmp = fastaID;
		fastaID = infix(tmp,0,length(tmp)-length(seq)); 
		append(fastaID,infix(tmp,length(tmp)-length(seq)+clipFront,length(tmp)-clipBack)); 
	}
	else
		qual = infix(qual,clipFront,length(qual)-clipBack);

	// now adapt the sequence
	seq = infix(seq,clipFront,length(seq)-clipBack);

	return true;	
	
}


//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences with or w/o quality values
template <typename TReadSet, typename TNameSet, typename TRazerSOptions>
bool loadReads(
	TReadSet &reads, 
	TNameSet &fastaIDs, 
	SeqFileIn &seqFile,
	TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);
#ifdef RAZERS_MICRO_RNA
	if(options.microRNA) countN = false;
#endif

	CharString fastaId;
	String<Dna5Q> seq;
	CharString qual;
	
	unsigned kickoutcount = 0;
	unsigned maxReadLength = 0;
	while (!atEnd(seqFile))
	{
        readRecord(fastaId, seq, qual, seqFile);
		if (options.readNaming == 0
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			|| options.fastaIdQual
#endif
			)
			appendValue(fastaIDs, fastaId); // append Fasta id
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		//check if sequence has a clip tag
		if(options.minClippedLen > 0)
		{
			if(!_clipReads(back(fastaIDs),seq,qual,options))
			{
				clear(seq);
				clear(qual);
				++kickoutcount;
			}
		}
		if(options.fastaIdQual && !empty(seq))
		{
			if(options.minClippedLen == 0)_clipReads(back(fastaIDs),seq,qual,options); // if the header wasnt clipped before, then clip now!! necessary for quality in header
			qual = suffix(back(fastaIDs),length(back(fastaIDs))-length(seq));
			if(options.readNaming == 0)
				back(fastaIDs) = prefix(back(fastaIDs),length(back(fastaIDs))-length(seq));
			else clear(back(fastaIDs));
		}
#endif
		if (countN)
		{
			int count = 0;
			int cutoffCount = (int)(options.errorRate * length(seq));
			for (unsigned j = 0; j < length(seq); ++j)
				if (getValue(seq, j) == 'N')
					if (++count > cutoffCount)
					{
						clear(seq);
						clear(qual);
						++kickoutcount;
						break;
					}
// low qual. reads are empty to output them and their id later as LQ reads
//			if (count > cutoffCount) continue;
		}
#ifdef RAZERS_MICRO_RNA
		if(options.microRNA && length(seq)<options.rnaSeedLength) 
		{
			clear(seq);
			clear(qual);
		}			
#endif

		// store dna and quality together
		assignQualities(seq, qual); 
		if (options.trimLength > 0 && length(seq) > (unsigned)options.trimLength)
			resize(seq, options.trimLength);
		appendValue(reads, seq, Generous());
		if (maxReadLength < length(seq))
			maxReadLength = length(seq);
	}

	typedef Shape<Dna, SimpleShape> TShape;
	typedef typename SAValue< Index<TReadSet, IndexQGram<TShape, OpenAddressing> > >::Type TSAValue;
	TSAValue sa(0, 0);
	sa.i1 = ~sa.i1;
	sa.i2 = ~sa.i2;

    bool success = true;
	if ((unsigned)sa.i1 < length(reads) - 1)
	{
		::std::cerr << "Maximal read number of " << (unsigned)sa.i1 + 1 << " exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		success = false;
	}
	if ((unsigned)sa.i2 < maxReadLength - 1)
	{
		::std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		success = false;
	}
	
	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality reads.\n";
	return success;
}

//////////////////////////////////////////////////////////////////////////////
// Read the first sequence of a multi-sequence file
// and return its length
inline int estimateReadLength(SeqFileIn &seqFile)
{
	if (atEnd(seqFile))
		return RAZERS_READS_FAILED;

    typedef String<char, Array<1000> > TBuffer;

    // read chunk into buffer
    TBuffer buffer;
    resize(buffer, capacity(buffer));
    size_t len = seqFile.stream.readsome(&buffer[0], length(buffer));
    for (size_t i = 0; i < len; ++i)
        seqFile.stream.unget();
    resize(buffer, len);

    // parse record from buffer
    DirectionIterator<TBuffer, Input>::Type iter = directionIterator(buffer, Input());
    CharString fastaId, seq;
    readRecord(fastaId, seq, iter, seqFile.format);
    return length(seq);
}


/*
//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TGenomeSet>
bool loadGenomes(TGenomeSet &genomes, const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
	resize(genomes, seqCount, Exact());
	for(unsigned i = 0; i < seqCount; ++i)
		assignSeq(genomes[i], multiFasta[i], Fasta());		// read Genome sequence

	return (seqCount > 0);
}*/

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences from multiple files
template <typename TGenomeSet>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList)
{
    SeqFileIn seqFile;
    StringSet<CharString> ids;
	for (size_t i = 0; i != length(fileNameList); ++i)
	{
		if (!open(seqFile, toCString(fileNameList[i])))
            return false;

        readRecords(ids, genomes, seqFile);
        clear(ids);
        close(seqFile);
	}
	return length(genomes);
}

#ifdef RAZERS_MICRO_RNA
	template <typename TReadMatch>
	struct LessRNoGPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

#ifdef RAZERS_MATEPAIRS
			// pair match id
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;
#endif

			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			if(a.editDist < b.editDist) return true;
			if(a.editDist > b.editDist) return false;
			return (a.mScore > b.mScore);
		}
	};

	template <typename TReadMatch>
	struct LessRNoEdistHLen : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

#ifdef RAZERS_MATEPAIRS
			// pair match id
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;
#endif

			if(a.editDist < b.editDist) return true;
			if(a.editDist > b.editDist) return false;
			return (a.mScore > b.mScore);

		}
	};
#else
	
	
	template <typename TReadMatch>
	struct LessRNoGPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// pair match id
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;

			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
#ifdef RAZERS_MATEPAIRS
			return a.pairScore > b.pairScore;
#else
			return a.editDist < b.editDist;
#endif
		}
	};
#endif

	// ... to sort matches and remove duplicates with equal gEnd
	template <typename TReadMatch>
	struct LessRNoGEndPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

#ifdef RAZERS_MATEPAIRS
			// pair match id
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;
#endif

			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gEnd   < b.gEnd) return true;
			if (a.gEnd   > b.gEnd) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
#ifdef RAZERS_MATEPAIRS
			return a.pairScore > b.pairScore;
#else
			return a.editDist < b.editDist;
#endif
		}
	};

	template <typename TReadMatch>
	struct LessErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

#ifdef RAZERS_MATEPAIRS
			// pair match id
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;
#endif

			// quality
#ifdef RAZERS_MATEPAIRS
			return a.pairScore > b.pairScore;
#else
			return a.editDist < b.editDist;
#endif
		}
	};
	
#ifdef RAZERS_SPLICED
template <typename TReadMatch>
struct LessSplicedErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
{
	inline bool operator() (TReadMatch const &a, TReadMatch const &b) const
	{
		// read number
		if (a.rseqNo < b.rseqNo ) return true;
		if (a.rseqNo > b.rseqNo ) return false;
		//if ((a.rseqNo >> 1) < (b.rseqNo >> 1)) return true;
		//if ((a.rseqNo >> 1) > (b.rseqNo >> 1)) return false;

			// pair match id
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;

		// quality
		if (a.pairScore > b.pairScore) return true;
		if (a.pairScore < b.pairScore) return false;
		return a.pairId < b.pairId;
	}
};
#endif	
	
	
#ifdef RAZERS_DIRECT_MAQ_MAPPING

	struct QualityBasedScoring{};

	template <typename TReadMatch>
	struct LessRNoMQ : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;
			
			// pair match id
			if (a.pairId < b.pairId) return true;
			if (a.pairId > b.pairId) return false;

			// quality
			if (a.mScore < b.mScore) return true; // sum of quality values of mismatches (the smaller the better)
			if (a.mScore > b.mScore) return false;
			
			return (a.editDist < b.editDist); // seedEditDist?
			// genome position and orientation
	/*		if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;
	*/		
		}
	};
#endif

//////////////////////////////////////////////////////////////////////////////
// Remove duplicate matches and leave at most maxHits many distanceRange
// best matches per read
template < typename TMatches, typename TOptions >
void maskDuplicates(TMatches &matches, TOptions & options)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;

    double beginTime = sysTime();
    
	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal ends

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoGEndPos<TMatch>());

	typename	TMatch::TGPos gBegin = -1;
	typename	TMatch::TGPos gEnd = -1;
	unsigned	gseqNo = -1;
	unsigned	readNo = -1;
	char		orientation = '-';
    unsigned    masked = 0;

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	for (; it != itEnd; ++it) 
	{
#ifdef RAZERS_MATEPAIRS
		if ((*it).pairId != 0) continue;
#endif
		if (gEnd == (*it).gEnd && orientation == (*it).orientation &&
			gseqNo == (*it).gseqNo && readNo == (*it).rseqNo) 
		{
			(*it).orientation = '-';
            masked += 1;
			continue;
		}
		readNo = (*it).rseqNo;
		gseqNo = (*it).gseqNo;
		gEnd = (*it).gEnd;
		orientation = (*it).orientation;
	}

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal begins

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoGPos<TMatch>());

	orientation = '-';

	it = begin(matches, Standard());
	itEnd = end(matches, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-'
#ifdef RAZERS_MATEPAIRS
			|| ((*it).pairId != 0)
#endif
			) continue;
		if (gBegin == (*it).gBegin && readNo == (*it).rseqNo &&
			gseqNo == (*it).gseqNo && orientation == (*it).orientation) 
		{
			(*it).orientation = '-';
            masked += 1;
			continue;
		}
		readNo = (*it).rseqNo;
		gseqNo = (*it).gseqNo;
		gBegin = (*it).gBegin;
		orientation = (*it).orientation;
	}

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessErrors<TMatch>());

    options.timeMaskDuplicates = sysTime() - beginTime;
    if (options._debugLevel >= 2)
        fprintf(stderr, " [%u matches masked]", masked);
}

//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template < typename TMatches, typename TCounts >
void countMatches(TMatches &matches, TCounts &cnt)
{
	//typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	typedef typename Value<TCounts>::Type					TRow;
	typedef typename Value<TRow>::Type						TValue;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	
	unsigned readNo = -1;
	short editDist = -1;
	int64_t count = 0;
	int64_t maxVal = std::numeric_limits<TValue>::max();

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo && editDist == (*it).editDist)
			++count;
		else
		{
			if (readNo != (unsigned)-1 && (unsigned)editDist < length(cnt))
				cnt[editDist][readNo] = (maxVal < count)? maxVal : count;
			readNo = (*it).rseqNo;
			editDist = (*it).editDist;
			count = 1;
		}
	}
	if (readNo != (unsigned)-1 && (unsigned)editDist < length(cnt))
		cnt[editDist][readNo] = count;


}

//////////////////////////////////////////////////////////////////////////////

template < typename TReadNo, typename TMaxErrors >
inline void 
setMaxErrors(Nothing &, TReadNo, TMaxErrors)
{
}

template < typename TSwift, typename TReadNo, typename TMaxErrors >
inline void 
setMaxErrors(TSwift &swift, TReadNo readNo, TMaxErrors maxErrors)
{
	int minT = _qgramLemma(swift, readNo, maxErrors);
	if (minT > 1)
	{
		if (maxErrors < 0) minT = std::numeric_limits<int>::max();
//		::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
		setMinThreshold(swift, readNo, (unsigned)minT);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TCounts, typename TSpec, typename TSwift >
void compactMatches(TMatches &matches, TCounts & 
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		cnts
#endif
	, RazerSOptions<TSpec> &options
	, bool compactFinal ,
	TSwift &
#if defined RAZERS_DIRECT_MAQ_MAPPING || defined RAZERS_MASK_READS
		swift
#endif
	)
{
    double beginTime = sysTime();
    
    // Get rid of "unused variable" warnings.  This is hard to read
    // and should not be done anywhere else. Better not use ifdefs.
    (void)compactFinal;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping) compactMatches(matches, cnts,options,compactFinal,swift,true);
#endif
	//typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
#ifdef RAZERS_MICRO_RNA
	if(options.microRNA)
		::std::sort(
			begin(matches, Standard()),
			end(matches, Standard()), 
			LessRNoEdistHLen<TMatch>());
	int bestMScore = 0;
#endif
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int64_t disabled = 0;
#ifdef RAZERS_MICRO_RNA
	if(options.microRNA && options.purgeAmbiguous)
		++hitCountCutOff;	// we keep one more match than we actually want, so we can later decide
							// whether the read mapped more than maxhits times 
#endif
	int editDistCutOff = std::numeric_limits<int>::max();

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
#ifdef RAZERS_MASK_READS
	TIterator ditBeg = it;
#endif

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo
#ifdef RAZERS_MATEPAIRS
			&& (*it).pairId == 0
#endif
			)
		{ 
			if ((*it).editDist >= editDistCutOff) continue;
#ifdef RAZERS_MICRO_RNA
			if ( (*it).mScore < bestMScore) continue;
#endif
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = (*it).editDist - 1;
					if (options.purgeAmbiguous && (options.distanceRange == 0 || (*it).editDist < options.distanceRange))
						maxErrors = -1;

					setMaxErrors(swift, readNo, maxErrors);

					if (maxErrors == -1)
						++disabled;
//						::std::cerr << "(read #" << readNo << " disabled)";

					if (options.purgeAmbiguous)
					{
						if (options.distanceRange == 0 || (*it).editDist < options.distanceRange || compactFinal)
							dit = ditBeg;
						else {
							*dit = *it;
							++dit;
						}
					}
				}
#endif
				continue;
			}
		}
		else
		{
			readNo = (*it).rseqNo;
			hitCount = 0;
			if (options.distanceRange > 0)
				editDistCutOff = (*it).editDist + options.distanceRange;
#ifdef RAZERS_MICRO_RNA
			bestMScore = (*it).mScore;
#endif
#ifdef RAZERS_MASK_READS
			ditBeg = dit;
#endif
		}
		*dit = *it;
		++dit;
	}
	if (options._debugLevel >= 2)
	{
		std::cerr << " [" << length(matches) - (dit - begin(matches, Standard())) << " matches removed]";
		std::cerr << " [" << disabled << " reads disabled]";
	}
	
	resize(matches, dit - begin(matches, Standard()));

    options.timeCompactMatches += sysTime() - beginTime;
}


#ifdef RAZERS_DIRECT_MAQ_MAPPING
//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TCounts, typename TSpec, typename TSwift >
void compactMatches(TMatches &matches, TCounts &cnts, RazerSOptions<TSpec> &, bool, TSwift &, bool dontCountFirstTwo)
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;

    double beginTime = sysTime();
	
	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoMQ<TMatch>());
	
	unsigned readNo = -1;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;

	//number of errors may not exceed 31!
	bool second = true;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo)
		{
			//second best match
			if (second)
			{
				second = false;
				if((cnts[1][(*it).rseqNo] & 31)  > (*it).editDist)
				{
					//this second best match is better than any second best match before
					cnts[1][(*it).rseqNo] = (*it).editDist; // second best dist is this editDist
										// count is 0 (will be updated if countFirstTwo)
				}
				if(!dontCountFirstTwo) 
					if((cnts[1][(*it).rseqNo]>>5) != 2047) cnts[1][(*it).rseqNo] += 32;
			}
			else
			{
				if ((*it).editDist <= (cnts[0][(*it).rseqNo] & 31) )
					if(cnts[0][(*it).rseqNo]>>5 != 2047)
						cnts[0][(*it).rseqNo] +=32;
				if ((*it).editDist <= (cnts[1][(*it).rseqNo] & 31) )
					if((cnts[1][(*it).rseqNo]>>5) != 2047)
						cnts[1][(*it).rseqNo] +=32;
				continue;
			}
		} else
		{	//best match
			second = true;
			readNo = (*it).rseqNo;
			//cnts has 16bits, 11:5 for count:dist
			if((cnts[0][(*it).rseqNo] & 31)  > (*it).editDist)
			{
				//this match is better than any match before
				cnts[1][(*it).rseqNo] = cnts[0][(*it).rseqNo]; // best before is now second best 
									       // (count will be updated when match is removed)
				cnts[0][(*it).rseqNo] = (*it).editDist; // best dist is this editDist
									// count is 0 (will be updated if countFirstTwo)
			}
			if(!dontCountFirstTwo) 
				if((cnts[0][(*it).rseqNo]>>5) != 2047) cnts[0][(*it).rseqNo] += 32;	// shift 5 to the right, add 1, shift 5 to the left, and keep count
		}
		*dit = *it;
		++dit;
	}

	resize(matches, dit - begin(matches, Standard()));
    options.timeCompactMatches += sysTime() - beginTime;
}
#endif


#ifdef RAZERS_MICRO_RNA

template < typename TMatches, typename TStats, typename TSpec >
void purgeAmbiguousRnaMatches(TMatches &matches, TStats &stats, RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type                          TMatch;
	typedef typename Iterator<TMatches, Standard>::Type             TIterator;

	::std::sort(begin(matches, Standard()),end(matches, Standard()),LessRNoEdistHLen<TMatch>());
	int bestMScore = 0;

	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int editDistCutOff = std::numeric_limits<int>::max();

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	for (; it != itEnd; ++it)
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo)
		{
            if ((*it).editDist >= editDistCutOff) continue;
			if ( (*it).mScore < bestMScore) continue;
            stats[readNo] += 1; 
			if (++hitCount >= hitCountCutOff)
			{
				if (hitCount == hitCountCutOff)
					dit = ditBeg;
				continue;
			}
		}
		else
		{
			readNo = (*it).rseqNo;
            stats[readNo] = 1;
			hitCount = 0;
			if (options.distanceRange > 0)
				editDistCutOff = (*it).editDist + options.distanceRange;
			bestMScore = (*it).mScore;
			ditBeg = dit;
		}
		*dit = *it;
		++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
}
 


//////////////////////////////////////////////////////////////////////////////
// Hamming verification recording sum of mismatching qualities in m.mScore
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,					// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,				// read number
	TReadSet &readSet,				// reads
	RazerSOptions<TSpec> const &options,		// RazerS options
	MicroRNA)					// MaqMapping
{

	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

	unsigned ndlLength = sequenceLength(rseqNo, readSet);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read		= readSet[rseqNo];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git	= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);

	// this is max number of errors the seed should have
	unsigned maxErrorsSeed = (unsigned)(options.rnaSeedLength * options.errorRate);	
	unsigned minSeedErrors = maxErrorsSeed + 1;
	unsigned bestHitLength = 0;

	for (; git < gitEnd; ++git)
	{
		bool hit = true;
		unsigned hitLength = 0;
		unsigned count = 0;
		unsigned seedErrors = 0;
		TGenomeIterator g = git;	//maq would count errors in the first 28bp only (at least initially. not for output)
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
		{
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
			{
				if (count < options.rnaSeedLength)		// the maq (28bp-)seed
				{
					if(++seedErrors > maxErrorsSeed)
					{
						hit = false;
						break;
					}
				}
				else 
					break;
			}
			++count;
		}
		if (hit) hitLength = count;
		if (hitLength > bestHitLength ) //simply take the longest hit
		{
			minSeedErrors = seedErrors;
			bestHitLength = hitLength;
			m.gBegin = git - begin(host(inf), Standard());
		}
	}

//	std::cout  << "options.absMaxQualSumErrors" << options.absMaxQualSumErrors << std::endl;
//	std::cout  << "maxSeedErrors" << maxErrorsSeed << std::endl;
//	std::cout  << "minErrors" << minSeedErrors << std::endl;
//	if(derBesgte) ::std::cout << minErrors <<"minErrors\n";
	if (minSeedErrors > maxErrorsSeed) return false;
	
	m.gEnd = m.gBegin + bestHitLength;
	m.editDist = minSeedErrors;			// errors in seed or total number of errors?
	m.mScore = bestHitLength;
	m.seedEditDist = minSeedErrors;
	return true;
}	

#endif



//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadSet &readSet,						// reads
	TMyersPatterns const &,					// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,	// RazerS options
	SwiftSemiGlobalHamming)					// Hamming only
{
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping) 
		return matchVerify(m,inf,rseqNo,readSet,options,QualityBasedScoring());
#endif

#ifdef RAZERS_MICRO_RNA
	if(options.microRNA) 
		return matchVerify(m,inf,rseqNo,readSet,options,MicroRNA());
#endif

	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[rseqNo] << ::std::endl;
#endif

	unsigned ndlLength = sequenceLength(rseqNo, readSet);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read				= readSet[rseqNo];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git		= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);

	unsigned maxErrors = (unsigned)(ndlLength * options.errorRate);
	unsigned minErrors = maxErrors + 1;

	for (; git < gitEnd; ++git)
	{
		unsigned errors = 0;
		TGenomeIterator g = git;
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
				if (++errors > maxErrors)
					break;
		if (minErrors > errors)
		{
			minErrors = errors;
			m.gBegin = git - begin(host(inf), Standard());
		}
	}

	if (minErrors > maxErrors) return false;

	m.gEnd = m.gBegin + ndlLength;
	m.editDist = minErrors;
	return true;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadSet &readSet,	    				// reads
	TMyersPatterns &forwardPatterns,		// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,	// RazerS options
	SwiftSemiGlobal)						// Swift specialization
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;

	// find read match end
	typedef Finder<TGenomeInfix>							TMyersFinder;
	typedef typename Value<TMyersPatterns>::Type			TMyersPattern;

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>		TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>				TReadRev;
	typedef Finder<TGenomeInfixRev>							TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>			TMyersPatternRev;

	TMyersFinder myersFinder(inf);
	TMyersPattern &myersPattern = forwardPatterns[rseqNo];

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[rseqNo]<<::std::endl;
#endif

    unsigned ndlLength = sequenceLength(rseqNo, readSet);
	int maxScore = std::numeric_limits<int>::min();
	int minScore = -(int)(ndlLength * options.errorRate);
	TMyersFinder maxPos;

	// find end of best semi-global alignment
	while (find(myersFinder, myersPattern, minScore))
		if (maxScore <= getScore(myersPattern)) 
		{
			maxScore = getScore(myersPattern);
			maxPos = myersFinder;
		}
	
	if (maxScore < minScore) return false;
	m.editDist	= (unsigned)-maxScore;
	setEndPosition(inf, m.gEnd = (beginPosition(inf) + position(maxPos) + 1));

	// limit the beginning to needle length plus errors (== -maxScore)
	if (length(inf) > ndlLength - maxScore)
		setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
	
	// find beginning of best semi-global alignment
	TGenomeInfixRev		infRev(inf);
	TReadRev			readRev(readSet[rseqNo]);
	TMyersFinderRev		myersFinderRev(infRev);
	TMyersPatternRev	myersPatternRev(readRev);

	_patternMatchNOfPattern(myersPatternRev, options.matchN);
	_patternMatchNOfFinder(myersPatternRev, options.matchN);
	while (find(myersFinderRev, myersPatternRev, maxScore))
		m.gBegin = m.gEnd - (position(myersFinderRev) + 1);

	return true;
}


#ifdef RAZERS_DIRECT_MAQ_MAPPING
//////////////////////////////////////////////////////////////////////////////
// Hamming verification recording sum of mismatching qualities in m.mScore
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,					// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,				// read number
	TReadSet &readSet,				// reads
	RazerSOptions<TSpec> const &options,		// RazerS options
	QualityBasedScoring)					// MaqMapping
{
	
	typedef Segment<TGenome, InfixSegment>				TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

// #ifdef RAZERS_DEBUG
// 	cout<<"Verify: "<<::std::endl;
// 	cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
// 	cout<<"Read:   "<<host(myersPattern)<<::std::endl;
// #endif

//	bool derBesgte = false;
	//if(rseqNo == 2) derBesgte = true;
//	if(derBesgte) ::std::cout << "der besagte\n";
	unsigned ndlLength = sequenceLength(rseqNo, readSet);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read		= readSet[rseqNo];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git	= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	TGenomeIterator bestIt	= begin(inf, Standard());

	// this is max number of errors the 28bp 'seed' should have
	//assuming that maxErrors-1 matches can be found with 100% SN 
	unsigned maxErrorsSeed = (unsigned)(options.artSeedLength * options.errorRate) + 1;	
	unsigned maxErrorsTotal = (unsigned)(ndlLength * 0.25); //options.maxErrorRate);
	unsigned minErrors = maxErrorsTotal + 1;
	int minQualSumErrors = options.absMaxQualSumErrors + 10;
	unsigned minSeedErrors = maxErrorsSeed + 1;

	for (; git < gitEnd; ++git)
	{
		bool hit = true;
		unsigned errors = 0;
		unsigned count = 0;
		unsigned seedErrors = 0;
		int qualSumErrors = 0;
		TGenomeIterator g = git;	//maq would count errors in the first 28bp only (at least initially. not for output)
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
		{
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
			{
			//	::std::cout << count << "<-";
				if (++errors > maxErrorsTotal)
				{
					hit = false;
					break;
				}
				//int qualityValue = (int)((unsigned char)*r >> 3);
				//qualSumErrors += (qualityValue < options.mutationRateQual) ? qualityValue : options.mutationRateQual;
				qualSumErrors += (getQualityValue(*r) < options.mutationRateQual) ? getQualityValue(*r) : options.mutationRateQual;
				if(qualSumErrors > options.absMaxQualSumErrors || qualSumErrors > minQualSumErrors)
				{
					hit = false;
					break;
				}
				if (count < options.artSeedLength)		// the maq (28bp-)seed
				{
					if(++seedErrors > maxErrorsSeed)
					{
						hit = false;
						break;
					}
					if(qualSumErrors > options.maxMismatchQualSum)
					{
						hit = false;
						break;							
					}// discard match, if 'seed' is bad (later calculations are done with the quality sum over the whole read)
				}
			}
			++count;
		}
		if (hit && (qualSumErrors < minQualSumErrors /*&& seedErrors <=maxErrorsSeed*/) ) //oder (seedErrors < minSeedErrors)
		{
			minErrors = errors;
			minSeedErrors = seedErrors;
			minQualSumErrors = qualSumErrors;
			m.gBegin = git - begin(host(inf), Standard());
			bestIt = git;
		}
	}
	if (minQualSumErrors > options.absMaxQualSumErrors || minSeedErrors > maxErrorsSeed || minErrors > maxErrorsTotal) return false;
	
	m.gEnd = m.gBegin + ndlLength;
	m.editDist = minErrors;			// errors in seed or total number of errors?
	m.mScore = minQualSumErrors;
	m.seedEditDist = minSeedErrors;
	return true;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet,
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadSet &readSet,	    				// reads
	TMyersPatterns const & pat,				// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,		// RazerS options
	SwiftSemiGlobal const &swiftsemi,				// Hamming only
	QualityBasedScoring)						// Swift specialization
{
	//if(!options.maqMapping) 
		return matchVerify(m,inf,rseqNo,readSet,pat,options,swiftsemi);
	//else
	//	return matchVerify(m,inf,rseqNo,readSet,pat,options,swiftsemi); // todo!
}
#endif




#ifndef RAZERS_PARALLEL
//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TCounts,
	typename TSpec >
void mapSingleReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPattern,
	TVerifier &forwardPatterns,
	TCounts & cnts,
	char orientation,				// q-gram index of reads
	RazerSOptions<TSpec> &options)
{
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Value<TMatches>::Type					TMatch;

	
	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinder;
	//typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	TReadSet &readSet = host(host(swiftPattern));
	TSwiftFinder swiftFinder(genome, options.repeatLength, 1);
	
	TMatch m = {	// to supress uninitialized warnings
		0, 0, 0, 0,
#ifdef RAZERS_MATEPAIRS
		0, 0, 0,
#endif
		0,
#ifdef RAZERS_EXTENDED_MATCH
		0, 0,
#endif
#ifdef RAZERS_SPLICED
		0,
#endif
		0
	};

	SEQAN_PROTIMESTART(map_contig_time);

	// iterate all verification regions returned by SWIFT
	TSize gLength = length(genome);
	int64_t localTP = 0;
	int64_t localFP = 0;

    double beginTime = sysTime();
    // Build q-gram index separately, so we can better compute the time for it.
    indexRequire(host(swiftPattern), QGramSADir());
    options.timeBuildQGramIndex += sysTime() - beginTime;

#ifdef RAZERS_MICRO_RNA
	while (find(swiftFinder, swiftPattern, 0.2)) 
#else
	while (find(swiftFinder, swiftPattern, options.errorRate)) 
#endif
	{
        // std::cout << "read id = " << (*swiftFinder.curHit).ndlSeqNo << ", " << beginPosition(swiftFinder) << std::endl;

		unsigned rseqNo = (*swiftFinder.curHit).ndlSeqNo;
        if (options.spec.DUMP_VERIFICATION_TASKS)
        {
            VerificationTask_ vt;
            vt.genomeBegin = beginPosition(infix(swiftFinder));
            vt.genomeEnd = beginPosition(infix(swiftFinder));
            vt.readId = rseqNo;
            vt.genomeId = gseqNo*2;
            if (orientation == 'R') ++vt.genomeId;
			vt.leftClip = (beginPosition(swiftFinder) >= 0)? 0: -beginPosition(swiftFinder);	// left clip if match begins left of the genome
            appendValue(options.verifications, vt);
        }
		if (!options.spec.DONT_VERIFY && 
			matchVerify(m, infix(swiftFinder), rseqNo, readSet, forwardPatterns, options, TSwiftSpec()))
		{
			// transform coordinates to the forward strand
			if (orientation == 'R') 
			{
				TSize temp = m.gBegin;
				m.gBegin = gLength - m.gEnd;
				m.gEnd = gLength - temp;
			}
			m.gseqNo = gseqNo;
			m.rseqNo = rseqNo;
			m.orientation = orientation;
#ifdef RAZERS_MATEPAIRS
			m.pairId = 0;
			m.pairScore = 0 - m.editDist;
#endif

			if (!options.spec.DONT_DUMP_RESULTS)
			{
				appendValue(matches, m, Generous());
				if (length(matches) > options.compactThresh)
				{
#ifndef RAZERS_MICRO_RNA
					typename Size<TMatches>::Type oldSize = length(matches);
#endif
                    maskDuplicates(matches, options);	// overlapping parallelograms cause duplicates
#ifdef RAZERS_DIRECT_MAQ_MAPPING
					if(options.maqMapping)
						compactMatches(matches, cnts, options, false, swiftPattern, true);
					else	
#endif
						compactMatches(matches, cnts, options, false, swiftPattern);
#ifdef RAZERS_MICRO_RNA
					options.compactThresh = 2 * length(matches);
#else
					if (length(matches) * 4 > oldSize)			// the threshold should not be raised
						options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed
#endif
				}
			}

			++localTP;
//			::std::cerr << "\"" << infix(swiftFinder) << "\"  ";
//			::std::cerr << hstkPos << " + ";
//			::std::cerr << ::std::endl;
		} else {
			++localFP;
//			::std::cerr << "\"" << infx(swiftFinder) << "\"   \"" << infix(swiftPattern) << "\"  ";
//			::std::cerr << rseqNo << " : ";
//			::std::cerr << hstkPos << " + ";
//			::std::cerr << bucketWidth << "  " << TP << ::std::endl;
		}
	}
	options.TP += localTP;
	options.FP += localFP;
	if (options._debugLevel >= 2)
	{
		double spec = 100;
		double time = SEQAN_PROTIMEDIFF(map_contig_time);
		if (localFP+localTP != 0)
			spec = (1000*localTP/(localFP+localTP))/10.0;
		::std::cerr << " [" << (int)((gLength / time)/100.0)/10.0 << "kbp/s " << time << "s] [" << spec << "%SP " << localFP+localTP << "V]";
	}
}
#endif



#ifdef RAZERS_MICRO_RNA

	// multiple sequences
	template <
		typename TSA, 
		typename TStringSet, 
		typename TShape, 
		typename TDir, 
		typename TValue, 
		typename TWithConstraints
	>
	inline void
	_qgramFillSuffixArray(
		TSA &sa, 
		TStringSet &stringSet,
		TShape &shape, 
		TDir &dir,
		Nothing,
		TWithConstraints const,
		TValue prefixLen,
		MicroRNA)
	{
		typedef typename Value<TStringSet>::Type					TString;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type				TSize;
		typedef typename Value<TShape>::Type				THash;

		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			int num_qgrams = prefixLen - (int)length(shape) + 1;

			typename Value<TSA>::Type localPos;
			assignValueI1(localPos, seqNo);
			assignValueI2(localPos, 0);

			TIterator itText = begin(sequence, Standard());
			if (TWithConstraints::VALUE) {
				THash h = hash(shape, itText) + 1;						// first hash
				if (dir[h] != (TSize)-1) sa[dir[h]++] = localPos;		// if bucket is enabled
			} else
				sa[dir[hash(shape, itText) + 1]++] = localPos;			// first hash

			for(int i = 1; i < num_qgrams; ++i)
			{
				++itText;
				assignValueI2(localPos, i);
				if (TWithConstraints::VALUE) {
					THash h = hashNext(shape, itText) + 1;				// next hash
					if (dir[h] != (TSize)-1) sa[dir[h]++] = localPos;	// if bucket is enabled
				} else
					sa[dir[hashNext(shape, itText) + 1]++] = localPos;	// next hash
			}
		}
	}

	template < typename TDir, typename TStringSet, typename TShape, typename TValue >
	inline void
	_qgramCountQGrams(TDir &dir, TStringSet &stringSet, TShape &shape, TValue prefixLen, MicroRNA)
	{
		typedef typename Value<TStringSet>::Type					TString;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type							TSize;
	
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = prefixLen - length(shape) + 1;

			TIterator itText = begin(sequence, Standard());
			++dir[hash(shape, itText)];
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				++dir[hashNext(shape, itText)];
			}
		}
	}
	

	template < typename TIndex, typename TValue>
	void createQGramIndex(TIndex &index, TValue prefixLen, MicroRNA)
	{
		typename Fibre<TIndex, QGramText>::Type	   &text  = indexText(index);
		typename Fibre<TIndex, QGramSA>::Type         &sa    = indexSA(index);
		typename Fibre<TIndex, QGramDir>::Type        &dir   = indexDir(index);
		typename Fibre<TIndex, QGramShape>::Type      &shape = indexShape(index);

		Nothing nothing;
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		_qgramCountQGrams(dir, text, shape, prefixLen,MicroRNA());

		if (_qgramDisableBuckets(index))
		{
			// 3. cumulative sum
			_qgramCummulativeSum(dir, True());

			// 4. fill suffix array
			_qgramFillSuffixArray(sa, text, shape, dir, nothing, True(), prefixLen, MicroRNA());

			// 5. correct disabled buckets
			_qgramPostprocessBuckets(dir);
		}
		else
		{
			// 3. cumulative sum
			_qgramCummulativeSum(dir, False());
			
			// 4. fill suffix array
			_qgramFillSuffixArray(sa, text, shape, dir, nothing, False(),prefixLen, MicroRNA());
		} 
	}


#endif



//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TMatches, 
	typename TReadSet, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet const &		readSet,
	TCounts & cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type							TRead;
	typedef Index<TReadSet, IndexQGram<TShape, TQGramIndexSpec> >	TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >						TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>							TMyersPattern;	// verifier

/*	// try opening each genome file once before running the whole mapping procedure
	int filecount = 0;
	int numFiles = length(genomeFileNameList);
	while(filecount < numFiles)
	{
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;
		file.close();
		++filecount;
	}
	*/

	// configure q-gram index
	TIndex swiftIndex(readSet, shape);
#ifdef RAZERS_OPENADDRESSING
	swiftIndex.alpha = 2;
#endif
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPattern(swiftIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;
	swiftPattern.params.printDots = options._debugLevel > 0;

	// init edit distance verifiers
	unsigned readCount = countSequences(swiftIndex);
	String<TMyersPattern> forwardPatterns;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}
#ifdef RAZERS_MICRO_RNA
	typename Size<TIndex>::Type qgram_count = 0;
	if(options.microRNA)
	{
		for(unsigned i = 0; i < countSequences(swiftIndex); ++i)
			if (sequenceLength(i, swiftIndex) >= options.rnaSeedLength)
				qgram_count += options.rnaSeedLength - (length(shape) - 1);
		resize(indexSA(swiftIndex), qgram_count, Exact());
		resize(indexDir(swiftIndex), _fullDirLength(swiftIndex), Exact());
		createQGramIndex(swiftIndex,options.rnaSeedLength,MicroRNA());
	}
#endif

	
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping)
	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			resize(cnts[i], readCount, 31); //initialize with maxeditDist, 11:5 for count:dist
	}
#endif

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
    options.timeCompactMatches = 0;
	options.timeMaskDuplicates = 0;

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
		Dna5String	genome;
		unsigned gseqNoWithinFile = 0;
		// iterate over genome sequences
		SEQAN_PROTIMESTART(find_time);
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
				mapSingleReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, cnts, 'F', options);

			if (options.reverse)
			{
				reverseComplement(genome);
				mapSingleReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, cnts, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
	}

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << "Masking duplicates took          \t" << options.timeMaskDuplicates << " seconds" << ::std::endl;
		::std::cerr << "Compacting matches took          \t" << options.timeCompactMatches << " seconds" << ::std::endl;
		::std::cerr << "Building q-gram index took       \t" << options.timeBuildQGramIndex << " seconds" << ::std::endl;
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:              \t" << options.FP + options.TP << ::std::endl;
		::std::cerr << "Verification counter:            \t" << options.TP << ::std::endl;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (given as StringSet)
template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet,
	typename TCounts, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type							TRead;
	typedef Index<TReadSet, IndexQGram<TShape, TQGramIndexSpec> >	TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >						TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>							TMyersPattern;	// verifier

	// configure q-gram index
	TIndex swiftIndex(readSet, shape);
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPattern(swiftIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;

	// init edit distance verifiers
	String<TMyersPattern> forwardPatterns;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		unsigned readCount = countSequences(swiftIndex);
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
    options.timeMaskDuplicates = 0;
	options.timeCompactMatches = 0;

	CharString	id;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
 	if(options.maqMapping)
 	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			resize(cnts[i], length(readSet), 31);
	}
#endif
	
	
	
	// iterate over genome sequences
	SEQAN_PROTIMESTART(find_time);
	for(unsigned gseqNo = 0; gseqNo < length(genomeSet); ++gseqNo)
	{
		if (options.forward)
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, cnts, 'F', options);

		if (options.reverse)
		{
			reverseComplement(genomeSet[gseqNo]);
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, cnts, 'R', options);
			reverseComplement(genomeSet[gseqNo]);
		}

	}
	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << "Masking duplicates took          \t" << options.timeMaskDuplicates << " seconds" << ::std::endl;
		::std::cerr << "Compacting matches took          \t" << options.timeCompactMatches << " seconds" << ::std::endl;
		::std::cerr << "Building q-gram index took       \t" << options.timeBuildQGramIndex << " seconds" << ::std::endl;
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:              \t" << options.FP + options.TP << ::std::endl;
		::std::cerr << "Verification counter:            \t" << options.TP << ::std::endl;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for single/mate-pair mapping
template <
	typename TMatches, 
	typename TReadSet,
	typename TCounts,
	typename TSpec,
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet &		readSet, 
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{

#ifdef RAZERS_SPLICED
	if (options.minMatchLen > 0)
	{
		//std::cout << "Spliced mapping\n";
		return mapSplicedReads(matches,genomeFileNameList,  genomeNames,gnoToFileMap, readSet, cnts, options, shape, Swift<TSwiftSpec>());
	}
#endif


#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, shape, Swift<TSwiftSpec>());
}

template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_SPLICED
	if (options.minMatchLen > 0)
	{
	//            std::cout << "Spliced mapping\n";
		return mapSplicedReads(matches, genomeSet, readSet, cnts, options, shape, Swift<TSwiftSpec>());
	
	}
#endif

#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(matches, genomeSet, readSet, cnts, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(matches, genomeSet, readSet, cnts, options, shape, Swift<TSwiftSpec>());
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different template specializations
template <typename TMatches, typename TReadSet, typename TCounts, typename TSpec>
int mapReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet &		readSet, 
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GenericShape>	gapped;

	// 2x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting shape
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

template <typename TMatches, typename TGenomeSet, typename TReadSet, typename TCounts, typename TSpec>
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet &		readSet, 
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GenericShape>	gapped;

	// 2x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting shape
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

}

#endif
