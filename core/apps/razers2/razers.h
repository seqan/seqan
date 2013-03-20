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
#include <seqan/index.h>
#include <seqan/store.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default options

	template < bool DONT_VERIFY_ = false, bool DONT_DUMP_RESULTS_ = false >
	struct RazerSSpec 
	{
		enum { DONT_VERIFY = DONT_VERIFY_ };				// omit verifying potential matches
		enum { DONT_DUMP_RESULTS = DONT_DUMP_RESULTS_ };	// omit dumping results
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
										// 3..use Fasta id, do not append /L and /R for mate pairs.
		unsigned	sortOrder;			// 0..sort keys: 1. read number, 2. genome position
										// 1..           1. genome pos50ition, 2. read number
		                       			// 2..           1. read name, 2. genome position
										// 3..           1. genome position, 2. read name
		int			positionFormat;		// 0..gap space
										// 1..position space
		const char	*runID;				// runID needed for gff output
        bool dontShrinkAlignments;

	// filtration parameters
		::std::string shape;			// shape (e.g. 11111111111)
		int			threshold;			// threshold
		int			tabooLength;		// taboo length
		int			repeatLength;		// repeat length threshold
		double		abundanceCut;		// abundance threshold

	// mate-pair parameters
		int			libraryLength;		// offset between two mates
		int			libraryError;		// offset tolerance
		unsigned	nextPairMatchId;	// use this id for the next mate-pair

	// verification parameters
		bool		matchN;				// false..N is always a mismatch, true..N matches with all
		unsigned char compMask[5];

	// statistics
		__int64		countFiltration;	// matches returned by the filter
		__int64		countVerification;	// matches returned by the verifier
		double		timeLoadFiles;		// time for loading input files
		double		timeMapReads;		// time for mapping reads
		double		timeDumpResults;	// time for dumping the results
		
		bool		lowMemory;		// set maximum shape weight to 13 to limit size of q-gram index
		bool		fastaIdQual;		// hidden option for special fasta+quality format we use


	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh

	// multi-threading

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
            dontShrinkAlignments = false;

			matchN = false;

			shape = "11111111111";
			threshold = 1;
			tabooLength = 1;
			repeatLength = 1000;
			abundanceCut = 1;

			libraryLength = 220;
			libraryError = 50;
			nextPairMatchId = 0;
			
			for (unsigned i = 0; i < 4; ++i)
				compMask[i] = 1 << i;
			compMask[4] = 0;

			compactThresh = 1024;

			lowMemory = false;		// set maximum shape weight to 13 to limit size of q-gram index
			fastaIdQual = false;

		}
	};
	
//////////////////////////////////////////////////////////////////////////////
// Typedefs
/*
	// definition of a Read match
	template <typename TGPos_>
	struct ReadMatch 
	{
		typedef typename MakeSigned_<TGPos_>::Type TGPos;

		unsigned		gseqNo;			// genome seqNo
		unsigned		rseqNo;			// read seqNo
		TGPos			beginPos;			// begin position of the match in the genome
		TGPos			gEnd;			// end position of the match in the genome
#ifdef RAZERS_MATEPAIRS
		unsigned		pairId;			// unique id for the two mate-pair matches (0 if unpaired)
		int				mateDelta:24;	// outer coordinate delta to the other mate 
		int				pairScore:8;	// combined score of both mates
#endif
		unsigned short	editDist;		// Levenshtein distance
#ifdef RAZERS_EXTENDED_MATCH
		short	 		mScore;
		short			seedEditDist;
#endif
		char			orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};
*/	
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
/*#ifdef RAZERS_CONCATREADS
	typedef StringSet<TRead, Owner<ConcatDirect<> > >	TReadSet;
#else
	typedef StringSet<TRead>							TReadSet;
#endif
*/
/*	typedef ReadMatch<Difference<TGenome>::Type>		TMatch;		// a single match
	typedef String<TMatch>								TMatches;	// array of matches
*/

	template <typename TReadSet, typename TShape, typename TSpec>
	struct Cargo< Index<TReadSet, IndexQGram<TShape, TSpec> > > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};

//////////////////////////////////////////////////////////////////////////////
// Memory tuning

#ifdef RAZERS_MEMOPT

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, IndexQGram<TShape, TSpec> > > 
	{
		typedef Pair<
			unsigned,				
			unsigned,
			BitPacked<22, 10>	// max. 4M reads of length < 1024
		> Type;
	};
	
#else

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, IndexQGram<TShape, TSpec> > > 
	{
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

	template <typename TReadSet, typename TShape>
	struct Size< Index<TReadSet, IndexQGram<TShape> > >
	{
		typedef unsigned Type;
	};

	template <typename TReadSet, typename TShape>
	struct Size< Index<TReadSet, IndexQGram<TShape, OpenAddressing> > >
	{
		typedef unsigned Type;
	};
	

#ifdef RAZERS_PRUNE_QGRAM_INDEX

	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <typename TReadSet, typename TShape, typename TSpec>
	inline bool _qgramDisableBuckets(Index<TReadSet, IndexQGram<TShape, TSpec> > &index) 
	{
		typedef Index<TReadSet, IndexQGram<TShape>	>		TReadIndex;
		typedef typename Fibre<TReadIndex, QGramDir>::Type	TDir;
		typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
		typedef typename Value<TDir>::Type					TSize;

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


	template <
		typename TFragmentStore_, 
		typename TRazerSOptions_,
		typename TPreprocessing_,
		typename TSwiftPattern_
	>
	struct MatchVerifier
	{
		typedef TFragmentStore_									TFragmentStore;
		typedef TRazerSOptions_									TOptions;
		typedef TPreprocessing_									TPreprocessing;
		typedef TSwiftPattern_									TSwiftPattern;
		
		typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
		typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
		typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
		typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
		typedef typename Size<TGenome>::Type					TSize;

		TFragmentStore	&store;
		TOptions		&options;			// RazerS options
		TPreprocessing	&preprocessing;
		TSwiftPattern	&swiftPattern;
		
		TAlignedRead	m;
		TAlignQuality	q;
		bool			onReverseComplement;
		TSize			genomeLength;
		bool			oneMatchPerBucket;
		
		MatchVerifier(TFragmentStore_ &_store, TOptions &_options, TPreprocessing &_preprocessing, TSwiftPattern &_swiftPattern):
			store(_store),
			options(_options),
			preprocessing(_preprocessing),
			swiftPattern(_swiftPattern)
		{
			onReverseComplement = false;
			genomeLength = 0;
			oneMatchPerBucket = false;
		}
		
		inline void push()
		{
			if (onReverseComplement) 
			{
				// transform coordinates to the forward strand
				m.beginPos = genomeLength - m.beginPos;
				m.endPos = genomeLength - m.endPos;
			}
			if (!options.spec.DONT_DUMP_RESULTS)
			{
				m.id = length(store.alignedReadStore);
				appendValue(store.alignedReadStore, m, Generous());
//				if (infix(store.readNameStore[m.readId], 0, length("SRR049254.14375884")) == "SRR049254.14375884")
//				    std::cerr << "append value " << m.beginPos << ", " << m.endPos << std::endl;
				appendValue(store.alignQualityStore, q, Generous());
				if (length(store.alignedReadStore) > options.compactThresh)
				{
					typename Size<TAlignedReadStore>::Type oldSize = length(store.alignedReadStore);
#ifndef RAZERS_DONTMASKDUPLICATES
					maskDuplicates(store);	// overlapping parallelograms cause duplicates
#endif
					compactMatches(store, options, swiftPattern);
					if (length(store.alignedReadStore) * 4 > oldSize)			// the threshold should not be raised
						options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed
					if (options._debugLevel >= 2)
						::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
				}
			}
			++options.countVerification;
		}
	};
 


//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences with or w/o quality values
template <typename TFSSpec, typename TFSConfig, typename TRazerSOptions>
bool loadReads(
	FragmentStore<TFSSpec, TFSConfig> &store,
	const char *fileName, 
	TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);

	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);

	unsigned seqCount = length(multiFasta);

	String<Dna5Q>	seq;
	CharString		qual;
	CharString		id;
	
	unsigned kickoutcount = 0;
	unsigned maxReadLength = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0 || options.readNaming == 3)
			assignSeqId(id, multiFasta[i], format);	// read Fasta id
		assignSeq(seq, multiFasta[i], format);					// read Read sequence
		assignQual(qual, multiFasta[i], format);				// read ascii quality values  
		if (countN)
		{
			int count = 0;
			int cutoffCount = (int)(options.errorRate * length(seq));
			for (unsigned j = 0; j < length(seq); ++j)
				if (getValue(seq, j) == 'N')
					if (++count > cutoffCount)
					{
						clear(seq);
                        clear(qual);  // So no qualities are assigned below.
						clear(id);
						++kickoutcount;
						break;
					}
// low qual. reads are empty to output them and their id later as LQ reads
//			if (count > cutoffCount) continue;
		}

		// store dna and quality together
		assignQualities(seq, qual); 
		if (options.trimLength > 0 && length(seq) > (unsigned)options.trimLength)
			resize(seq, options.trimLength);

		appendRead(store, seq, id);
		if (maxReadLength < length(seq))
			maxReadLength = length(seq);
	}
	// memory optimization
	reserve(store.readSeqStore.concat, length(store.readSeqStore.concat), Exact());
//	reserve(store.readNameStore.concat, length(store.readNameStore.concat), Exact());

	typedef Shape<Dna, SimpleShape> TShape;
	typedef typename SAValue< Index<StringSet<TRead>, IndexQGram<TShape, OpenAddressing> > >::Type TSAValue;
	TSAValue sa(0, 0);
	sa.i1 = ~sa.i1;
	sa.i2 = ~sa.i2;
	
	if ((unsigned)sa.i1 < length(store.readSeqStore) - 1)
	{
		::std::cerr << "Maximal read number of " << (unsigned)sa.i1 + 1 << " exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}
	if ((unsigned)sa.i2 < maxReadLength - 1)
	{
		::std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality reads.\n";
	return (seqCount > 0);
}

//////////////////////////////////////////////////////////////////////////////
// Read the first sequence of a multi-sequence file
// and return its length
inline int estimateReadLength(char const *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY))	// open the whole file
		return RAZERS_READS_FAILED;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);					// guess file format
	split(multiFasta, format);								// divide into single sequences

	if (length(multiFasta) == 0)
		return 0;

	Dna5String firstRead;
	assignSeq(firstRead, multiFasta[0], format);			// read the first sequence
	return length(firstRead);
}
	
	
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore>
	struct LessRNoGPos : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TAlignedReadQualityStore &qualStore;
		
		LessRNoGPos(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (TAlignedRead const &a, TAlignedRead const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// contig number
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

			// qualities
			SEQAN_ASSERT_NEQ(a.id, TAlignedRead::INVALID_ID);
			SEQAN_ASSERT_NEQ(b.id, TAlignedRead::INVALID_ID);
			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			if (qa.score > qb.score) return true;
			if (qb.score > qa.score) return false;
			
			// prefer reads that support more of the reference
			return _max(a.beginPos, a.endPos) > _max(b.beginPos, b.endPos);
		}
	};

	// ... to sort matches and remove duplicates with equal gEnd
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore>
	struct LessRNoGEndPos : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TAlignedReadQualityStore &qualStore;
		
		LessRNoGEndPos(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// contig number
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// end position
			typename TAlignedRead::TPos ea = _max(a.beginPos, a.endPos);
			typename TAlignedRead::TPos eb = _max(b.beginPos, b.endPos);
			if (ea < eb) return true;
			if (ea > eb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;

			// qualities
			SEQAN_ASSERT_NEQ(a.id, TAlignedRead::INVALID_ID);
			SEQAN_ASSERT_NEQ(b.id, TAlignedRead::INVALID_ID);
			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			if (qa.score > qb.score) return true;
			if (qb.score > qa.score) return false;
			
			// prefer reads that support more of the reference
			return _min(a.beginPos, a.endPos) < _min(b.beginPos, b.endPos);
		}
	};
	
    template <typename TAlignedReadStore, typename TReadNameStore, typename TAlignedReadQualityStore>
	struct LessRNameGPos : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TReadNameStore &nameStore;
		TAlignedReadQualityStore &qualStore;
		
		LessRNameGPos(TReadNameStore &nameStore, TAlignedReadQualityStore &_qualStore):
            nameStore(nameStore), qualStore(_qualStore) {}
		
		inline bool operator() (TAlignedRead const &a, TAlignedRead const &b) const 
		{
			// read name
			if (nameStore[a.readId] < nameStore[b.readId]) return true;
			if (nameStore[a.readId] > nameStore[b.readId]) return false;

			// contig number
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

			// qualities
			SEQAN_ASSERT_NEQ(a.id, TAlignedRead::INVALID_ID);
			SEQAN_ASSERT_NEQ(b.id, TAlignedRead::INVALID_ID);
			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			if (qa.score > qb.score) return true;
			if (qb.score > qa.score) return false;
			
			// prefer reads that support more of the reference
			return _max(a.beginPos, a.endPos) > _max(b.beginPos, b.endPos);
		}
	};
	
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore>
	struct LessErrors : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		TAlignedReadQualityStore &qualStore;
		
		LessErrors(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// quality
			SEQAN_ASSERT_NEQ(a.id, TAlignedRead::INVALID_ID);
			SEQAN_ASSERT_NEQ(b.id, TAlignedRead::INVALID_ID);
			if (a.contigId == TAlignedRead::INVALID_ID) return false;
			if (b.contigId == TAlignedRead::INVALID_ID) return true;
			return qualStore[a.id].errors < qualStore[b.id].errors;
		}
	};
	

//////////////////////////////////////////////////////////////////////////////
// Mark duplicate matches for deletion
template <typename TFragmentStore>
void maskDuplicates(TFragmentStore &store)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal ends

	sortAlignedReads(store.alignedReadStore, LessRNoGEndPos<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));

	TContigPos	beginPos = -1;
	TContigPos	endPos = -1;
	unsigned	contigId = TAlignedRead::INVALID_ID;
	unsigned	readId = TAlignedRead::INVALID_ID;
	bool		orientation = false;

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).pairMatchId != TAlignedRead::INVALID_ID && getMateNo(store, it->readId) != 0) continue;	// remove only single reads or left mates
		TContigPos itEndPos = _max((*it).beginPos, (*it).endPos);
		if (endPos == itEndPos && orientation == ((*it).beginPos < (*it).endPos) &&
			contigId == (*it).contigId && readId == (*it).readId) 
		{
			(*it).contigId = TAlignedRead::INVALID_ID;	// mark this alignment for deletion
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		endPos = itEndPos;
		orientation = (*it).beginPos < (*it).endPos;
	}

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal begins

	sortAlignedReads(store.alignedReadStore, LessRNoGPos<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));

	contigId = TAlignedRead::INVALID_ID;
	it = begin(store.alignedReadStore, Standard());
	itEnd = end(store.alignedReadStore, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).contigId == TAlignedRead::INVALID_ID) continue;
		if ((*it).pairMatchId != TAlignedRead::INVALID_ID && getMateNo(store, it->readId) != 0) continue;	// remove only single reads or left mates

		TContigPos itBeginPos = _min((*it).beginPos, (*it).endPos);
		if (beginPos == itBeginPos && readId == (*it).readId &&
			contigId == (*it).contigId && orientation == ((*it).beginPos < (*it).endPos))
		{
			(*it).contigId = TAlignedRead::INVALID_ID;	// mark this alignment for deletion
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		beginPos = itBeginPos;
		orientation = (*it).beginPos < (*it).endPos;
	}

	sortAlignedReads(store.alignedReadStore, LessErrors<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));
}

//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template <typename TFragmentStore, typename TCounts>
void countMatches(TFragmentStore &store, TCounts &cnt)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	//typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef typename Value<TCounts>::Type							TRow;
	typedef typename Value<TRow>::Type								TValue;
	
	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	
	TValue maxVal = MaxValue<TValue>::VALUE;
	unsigned maxError = length(cnt);
/*	
	for (; it != itEnd; ++it)
	{
		if (it->contigId == TAlignedRead::INVALID_ID) continue;
		unsigned errors = store.alignQualityStore[it->id].errors;
		if (errors >= maxError) continue;		
		if (cnt[errors][it->readId] < maxVal)
			++cnt[errors][it->readId];
	}
*/
	short errors = -1;
	__int64 count = 0;
	unsigned readId = TAlignedRead::INVALID_ID;
	for (; it != itEnd; ++it)
	{
		if (it->contigId == TAlignedRead::INVALID_ID) continue;
		if (readId == it->readId && errors == store.alignQualityStore[it->id].errors)
			++count;
		else
		{
			if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < maxError)
				cnt[errors][readId] = (maxVal < count)? (TValue)maxVal : (TValue)count;
			readId = it->readId;
			errors = store.alignQualityStore[it->id].errors;
			count = 1;
		}
	}
	if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < maxError)
		cnt[errors][readId] = (TValue)count;
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
		if (maxErrors < 0) minT = MaxValue<int>::VALUE;
//		::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
		setMinThreshold(swift, readNo, (unsigned)minT);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TFragmentStore, typename TSpec, typename TSwift >
void compactMatches(TFragmentStore &store, RazerSOptions<TSpec> &options, TSwift & swift)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	//typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;

	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	unsigned errorsCutOff = MaxValue<unsigned>::VALUE;

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	for (; it != itEnd; ++it) 
	{
		if ((*it).contigId == TAlignedRead::INVALID_ID) continue;
		unsigned errors = store.alignQualityStore[(*it).id].errors;
		if (readNo == (*it).readId && (*it).pairMatchId == TAlignedRead::INVALID_ID)
		{ 
			if (errors >= errorsCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = errors - 1;
					if (options.purgeAmbiguous && (options.distanceRange == 0 || errors < options.distanceRange))
						maxErrors = -1;

					setMaxErrors(swift, readNo, maxErrors);

					if (maxErrors == -1 && options._debugLevel >= 2)
						::std::cerr << "(read #" << readNo << " disabled)";

					if (options.purgeAmbiguous)
					{
						if (options.distanceRange == 0 || errors < options.distanceRange || IsSameType<TSwift, Nothing>::VALUE)
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
			readNo = (*it).readId;
			hitCount = 0;
			if (options.distanceRange > 0)
				errorsCutOff = errors + options.distanceRange;
			ditBeg = dit;
		}
		*dit = *it;
		++dit;
	}
	resize(store.alignedReadStore, dit - begin(store.alignedReadStore, Standard()));
	compactAlignedReads(store);
}



struct SemiGlobalHamming_;
struct SemiGlobalEdit_;


//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned readId,						// read number
	TReadSet &readSet,						// reads
	SwiftSemiGlobalHamming)					// Hamming only
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[readId] << ::std::endl;
#endif

	// verify
	TRead &read				= readSet[readId];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	unsigned ndlLength		= ritEnd - ritBeg;

	if (length(inf) < ndlLength) return false;
	TGenomeIterator git		= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	
	unsigned maxErrors = (unsigned)(ndlLength * verifier.options.errorRate);
	unsigned minErrors = maxErrors + 1;
	unsigned errorThresh = (verifier.oneMatchPerBucket)? MaxValue<unsigned>::VALUE: maxErrors;

	for (; git < gitEnd; ++git)
	{
		unsigned errors = 0;
		TGenomeIterator g = git;
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
			if ((verifier.options.compMask[ordValue(*g)] & verifier.options.compMask[ordValue(*r)]) == 0)
				if (++errors > maxErrors)
					break;
		
		if (errors < minErrors)
		{
			minErrors = errors;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		} else if (errorThresh < errors)
		{
			if (minErrors <= maxErrors)
			{
				verifier.m.endPos = verifier.m.beginPos + ndlLength;
				verifier.q.pairScore = verifier.q.score = -(int)minErrors;
				verifier.q.errors = minErrors;
				verifier.push();
				minErrors = maxErrors + 1;
			}
		}
	}

	if (minErrors <= maxErrors)
	{
		verifier.m.endPos = verifier.m.beginPos + ndlLength;
		verifier.q.pairScore = verifier.q.score = -(int)minErrors;
		verifier.q.errors = minErrors;
		if (!verifier.oneMatchPerBucket)
			verifier.push();
		return true;
	}
	return false;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned readId,						// read number
	TReadSet &readSet,						// reads
	SwiftSemiGlobal)						// Mismatches and Indels
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Prefix<TRead>::Type					TReadPrefix;
	typedef typename Position<TGenomeInfix>::Type			TPosition;

	// find read match end
	typedef Finder<TGenomeInfix>							TMyersFinder;
	typedef typename TMatchVerifier::TPreprocessing			TPreprocessing;
	typedef typename Value<TPreprocessing>::Type			TMyersPattern;

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>		TGenomeInfixRev;
	typedef Finder<TGenomeInfixRev>							TMyersFinderRev;

#ifdef RAZERS_NOOUTERREADGAPS
	typedef ModifiedString<TReadPrefix, ModReverse>			TReadPrefixRev;
	typedef Pattern<TReadPrefixRev, MyersUkkonenGlobal>		TMyersPatternRev;
#else
	typedef ModifiedString<TRead, ModReverse>				TReadRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>			TMyersPatternRev;
#endif

    unsigned ndlLength = sequenceLength(readId, readSet);
	int maxScore = MinValue<int>::VALUE;
	int minScore = -(int)(ndlLength * verifier.options.errorRate);
	TPosition maxPos = 0;
	TPosition lastPos = length(inf);
	unsigned minDistance = (verifier.oneMatchPerBucket)? lastPos: 1;


#ifdef RAZERS_NOOUTERREADGAPS
	TGenomeInfix origInf(inf);
	setEndPosition(inf, endPosition(inf) - 1);
	--ndlLength;
#endif

	TMyersFinder myersFinder(inf);
	TMyersPattern &myersPattern = verifier.preprocessing[readId];

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << "; " << length(host(inf)) - beginPosition(inf) << "," << length(host(inf)) - endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[readId]<<::std::endl;
#endif

	// find end of best semi-global alignment
	while (find(myersFinder, myersPattern, minScore))
	{
		TPosition pos = position(hostIterator(myersFinder));
		int score = getScore(myersPattern);
        // std::cerr << "found " << pos << ", " << length(host(inf)) - pos << ", " << score << std::endl;
		
#ifdef RAZERS_NOOUTERREADGAPS
		// Manually align the last base of the read
		//
		// In this case myersPattern contains the whole read without the
		// last base. We compare the bases and adjust the score.
		// We also have to adjust inf and remove the last base of the
		// genomic region that has to be verified.
		SEQAN_ASSERT_LT(pos + 1, length(origInf));
		if ((verifier.options.compMask[ordValue(origInf[pos + 1])] & verifier.options.compMask[ordValue(back(readSet[readId]))]) == 0)
			if (--score < minScore) continue;
#endif		
		if (lastPos + minDistance < pos)
		{
			if (minScore <= maxScore)
			{
				TPosition infBeginPos = beginPosition(inf);
				TPosition infEndPos = endPosition(inf);
				verifier.q.pairScore = verifier.q.score = maxScore;
				verifier.q.errors = -maxScore;

				// find beginning of best semi-global alignment
				setEndPosition(inf, verifier.m.endPos = (beginPosition(inf) + maxPos + 1));

#ifdef RAZERS_NOOUTERREADGAPS
				// The best score must be corrected to hold the score of the prefix w/o the last read base
				if ((verifier.options.compMask[ordValue(origInf[maxPos + 1])] & verifier.options.compMask[ordValue(back(readSet[readId]))]) == 0)
					++maxScore;

				TReadPrefixRev		readRev(prefix(readSet[readId], ndlLength));
				TMyersPatternRev	myersPatternRev(readRev);
#else
				TReadRev			readRev(readSet[readId]);
				TMyersPatternRev	myersPatternRev(readRev);
#endif

//				// limit the beginning to needle length plus errors (== -maxScore)
//				if (length(inf) > ndlLength - maxScore)
//					setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);

				// we eventually have to search before the beginning of our parallelogram
				// otherwise alignments of an island in the previous parallelogram
				// could be cut and prevent that an island in this parallelgram is found
				if (endPosition(inf) > (unsigned)(ndlLength - maxScore))
					setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
				else
					setBeginPosition(inf, 0);
				
				TGenomeInfixRev		infRev(inf);
				TMyersFinderRev		myersFinderRev(infRev);

				_patternMatchNOfPattern(myersPatternRev, verifier.options.matchN);
				_patternMatchNOfFinder(myersPatternRev, verifier.options.matchN);
				while (find(myersFinderRev, myersPatternRev, /*score*/maxScore))
					verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);

				setBeginPosition(inf, infBeginPos);
				setEndPosition(inf, infEndPos);

#ifdef RAZERS_NOOUTERREADGAPS
				// The match end position must be increased by the omitted base.
				++verifier.m.endPos;
#endif
				verifier.push();
//            std::cerr << "push(" << verifier.m.beginPos << ", " << verifier.m.endPos << ")" << std::endl;
//            std::cerr << "  push(" << length(host(inf)) - verifier.m.beginPos << ", " << length(host(inf)) - verifier.m.endPos << ")" << std::endl;
				maxScore = minScore - 1;
			}
		}
		if (score >= maxScore) 
		{
			maxScore = score;
			maxPos = pos;
		}
		lastPos = pos;
	}

	if (minScore <= maxScore)
	{
		verifier.q.pairScore = verifier.q.score = maxScore;
		verifier.q.errors = -maxScore;
		setEndPosition(inf, verifier.m.endPos = (beginPosition(inf) + maxPos + 1));

#ifdef RAZERS_NOOUTERREADGAPS
		// The best score must be corrected to hold the score of the prefix w/o the last read base
		if ((verifier.options.compMask[ordValue(origInf[maxPos + 1])] & verifier.options.compMask[ordValue(back(readSet[readId]))]) == 0)
			++maxScore;

		TReadPrefixRev		readRev(prefix(readSet[readId], ndlLength));
		TMyersPatternRev	myersPatternRev(readRev);
#else
		TReadRev			readRev(readSet[readId]);
		TMyersPatternRev	myersPatternRev(readRev);
#endif

//		// limit the beginning to needle length plus errors (== -maxScore)
//		if (length(inf) > ndlLength - maxScore)
//			setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);

		// we eventually have to search before the beginning of our parallelogram
		// otherwise alignments of an island in the previous parallelogram
		// could be cut and prevent that an island in this parallelgram is found
		if (endPosition(inf) > (unsigned)(ndlLength - maxScore))
			setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
		else
			setBeginPosition(inf, 0);
		
		// find beginning of best semi-global alignment
		TGenomeInfixRev		infRev(inf);
		TMyersFinderRev		myersFinderRev(infRev);

		_patternMatchNOfPattern(myersPatternRev, verifier.options.matchN);
		_patternMatchNOfFinder(myersPatternRev, verifier.options.matchN);
		while (find(myersFinderRev, myersPatternRev, maxScore))
			verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);

#ifdef RAZERS_NOOUTERREADGAPS
		// The match end position must be increased by the omitted base.
		++verifier.m.endPos;
#endif
		if (!verifier.oneMatchPerBucket)
//        {
			verifier.push();
//            std::cerr << "push(" << verifier.m.beginPos << ", " << verifier.m.endPos << ")" << std::endl;
//            std::cerr << "  push(" << length(host(inf)) - verifier.m.beginPos << ", " << length(host(inf)) - verifier.m.endPos << ")" << std::endl;
//        }
		return true;
	}
	return false;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TFragmentStore, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TPreprocessing,
	typename TRazerSOptions >
void mapSingleReads(
	TFragmentStore							& store,
	TGenome									& genome,				// genome ...
	unsigned								  contigId,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPattern,
	TPreprocessing							& preprocessing,
	char									  orientation,				// q-gram index of reads
	TRazerSOptions							& options)
{
	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >								TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >							TSwiftPattern;

	// VERIFICATION
	typedef MatchVerifier <
		TFragmentStore, 
		TRazerSOptions, 
		TPreprocessing, 
		TSwiftPattern >													TVerifier;
	typedef typename Fibre<TReadIndex, FibreText>::Type					TReadSet;
	
	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F') ::std::cerr << "[fwd]";
		else                    ::std::cerr << "[rev]";
	}

	TReadSet		&readSet = host(host(swiftPattern));
	TSwiftFinder	swiftFinder(genome, options.repeatLength, 1);
	TVerifier		verifier(store, options, preprocessing, swiftPattern);

	verifier.onReverseComplement = (orientation == 'R');
	verifier.genomeLength = length(genome);
	verifier.m.contigId = contigId;

	// iterate all verification regions returned by SWIFT
	while (find(swiftFinder, swiftPattern, options.errorRate))
	{
		verifier.m.readId = (*swiftFinder.curHit).ndlSeqNo;
		if (!options.spec.DONT_VERIFY)
			matchVerify(verifier, infix(swiftFinder), verifier.m.readId, readSet, TSwiftSpec());
		++options.countFiltration;
	}
}

//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TGNoToFile,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	StringSet<CharString>				& genomeFileNameList,
	String<TGNoToFile>					& gnoToFileMap,
	RazerSOptions<TSpec>				& options,
	TShape const						& shape,
	Swift<TSwiftSpec> const)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	typedef Index<TReadSeqStore, IndexQGram<TShape,OpenAddressing> >	TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

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
	TIndex swiftIndex(store.readSeqStore, shape);
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
#ifdef RAZERS_NOOUTERREADGAPS
			if (!empty(indexText(swiftIndex)[i]))
				setHost(forwardPatterns[i], prefix(indexText(swiftIndex)[i], length(indexText(swiftIndex)[i]) - 1));
#else
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
#endif
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}
	
	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
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
		Dna5String	genome;
		unsigned gseqNoWithinFile = 0;
		// iterate over genome sequences
		SEQAN_PROTIMESTART(find_time);
		for(; !_streamEOF(file); ++gseqNo)
		{
			if (options.genomeNaming == 0)
			{
				//readID(file, id, Fasta());			// read Fasta id
				readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
				appendValue(store.contigNameStore, id, Generous());
			}
			read(file, genome, Fasta());			// read Fasta sequence
			
			appendValue(gnoToFileMap, TGNoToFile(genomeName, gseqNoWithinFile));
			
			if (options.forward)
				mapSingleReads(store, genome, gseqNo, swiftPattern, forwardPatterns, 'F', options);

			if (options.reverse)
			{
				reverseComplement(genome);
				mapSingleReads(store, genome, gseqNo, swiftPattern, forwardPatterns, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
		file.close();
		++filecount;
	}

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (given as StringSet)
template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type				TRead;
	typedef Index<TReadSet, IndexQGram<TShape,OpenAddressing> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

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
#ifdef RAZERS_NOOUTERREADGAPS
			if (!empty(indexText(swiftIndex)[i]))
				setHost(forwardPatterns[i], prefix(indexText(swiftIndex)[i], length(indexText(swiftIndex)[i]) - 1));
#else
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
#endif
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}

	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	CharString	id;
	
	// iterate over genome sequences
	SEQAN_PROTIMESTART(find_time);
	for(unsigned gseqNo = 0; gseqNo < length(genomeSet); ++gseqNo)
	{
		if (options.forward)
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, 'F', options);

		if (options.reverse)
		{
			reverseComplement(genomeSet[gseqNo]);
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, 'R', options);
			reverseComplement(genomeSet[gseqNo]);
		}

	}
	
	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for single/mate-pair mapping
template <
	typename TFSSpec, 
	typename TFSConfig,
	typename TGNoToFile,
	typename TSpec,
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	StringSet<CharString> &	genomeFileNameList,
	String<TGNoToFile> & gnoToFileMap,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(store, genomeFileNameList, gnoToFileMap, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(store, genomeFileNameList, gnoToFileMap, options, shape, Swift<TSwiftSpec>());
}

template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(matches, genomeSet, readSet, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(matches, genomeSet, readSet, options, shape, Swift<TSwiftSpec>());
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different template specializations
template <typename TFSSpec, typename TFSConfig, typename TGNoToFile, typename TSpec>
int mapReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	StringSet<CharString> &	genomeFileNameList,
	String<TGNoToFile> & gnoToFileMap,
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
			return mapReads(store, genomeFileNameList, gnoToFileMap, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

template <typename TMatches, typename TGenomeSet, typename TReadSet, typename TSpec>
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet, 
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
			return mapReads(matches, genomeSet, readSet, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

}

#endif
