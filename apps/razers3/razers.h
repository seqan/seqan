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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/index/find_pigeonhole.h>
#include <seqan/store.h>
#include <seqan/pipe.h>
#include <seqan/parallel.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#ifdef RAZERS_PROFILE
#include "profile_timeline.h"
#endif  // #ifdef RAZERS_PROFILE

// No parallelism for less than MIN_PARALLEL_WORK reads.
const unsigned MIN_PARALLEL_WORK = 1; //100/*0*/; // TODO(holtgrew): Set to some useful value after development.

namespace seqan {

// Compact representation of a match.
template <typename TContigPos_>
struct MatchRecord
{
    typedef typename MakeSigned_<TContigPos_>::Type TContigPos;

    unsigned        contigId;       // genome seqNo
    unsigned        readId;         // read seqNo
    TContigPos      beginPos;       // begin position of the match in the genome
    TContigPos      endPos;         // end position of the match in the genome
#ifdef RAZERS_DEFER_COMPACTION
    bool            isRegistered; // registered in masking process.
#endif  // #ifdef RAZERS_DEFER_COMPACTION
    char            orientation;    // 'F', 'R', '-'
    short int       score;          // Levenshtein distance / score.
    unsigned        pairMatchId;            // unique id for the two mate-pair matches (0 if unpaired)
    int             libDiff : 24;     // outer distance difference from librarySize
    int             pairScore : 8;    // combined score of both mates

    static const unsigned INVALID_ID;

    MatchRecord() :
        contigId(std::numeric_limits<unsigned>::max()), readId(std::numeric_limits<unsigned>::max()),
        beginPos(0), endPos(0),
#ifdef RAZERS_DEFER_COMPACTION
        isRegistered(false),
#endif  // #ifndef RAZERS_DEFER_COMPACTION
        orientation('-'), score(0), pairMatchId(std::numeric_limits<unsigned>::max()),
        libDiff(0), pairScore(0)
    {}
};

template <typename TStream, typename TPos>
TStream &
operator<<(TStream & stream, MatchRecord<TPos> & record)
{
    stream << "(contigId=" << record.contigId << ", readId=" << record.readId << ", beginPos=" << record.beginPos << ", endPos = " << record.endPos << ", orientation=" << record.orientation << ", score=" << record.score << ", pairMatchId=" << record.pairMatchId << ", libDiff=" << record.libDiff << ", pairScore=" << record.pairScore << ")";
    return stream;
}

template <typename TGPos_>
const unsigned MatchRecord<TGPos_>::INVALID_ID = std::numeric_limits<unsigned>::max();

#ifdef RAZERS_PROFILE
enum
{
    TASK_WAIT,
    TASK_ON_CONTIG,
    TASK_INIT,
    TASK_REVCOMP,
    TASK_FILTER,
    TASK_VERIFY,
    TASK_WRITEBACK,
    TASK_COMPACT,
    TASK_DUMP_MATCHES,
    TASK_LOAD,
    TASK_SORT,
    TASK_COPY_FINDER
};
#endif  // #ifdef RAZERS_PROFILE

//////////////////////////////////////////////////////////////////////////////
// RazerS modes

// Alignment mode
struct RazerSLocal;
struct RazerSGlobal;
struct RazerSPrefix;

// Gap mode
struct RazerSGapped;
struct RazerSUngapped;

// Score mode
struct RazerSErrors;
struct RazerSScore;
struct RazerSMAQ;

template <typename TSpec = Default>
struct RazerSQuality;

template <typename TAlignMode_, typename TGapMode_, typename TScoreMode_, typename TMatchNPolicy_>
struct RazerSMode
{
    typedef TAlignMode_ TAlignMode;
    typedef TGapMode_   TGapMode;
    typedef TScoreMode_ TScoreMode;
    typedef TMatchNPolicy_  TMatchNPolicy;
};

enum AlignMode          {RAZERS_LOCAL, RAZERS_PREFIX, RAZERS_GLOBAL};
enum GapMode            {RAZERS_GAPPED, RAZERS_UNGAPPED};
enum ScoreMode          {RAZERS_ERRORS, RAZERS_SCORE, RAZERS_QUALITY};
enum CompactMatchesMode {COMPACT, COMPACT_FINAL, COMPACT_FINAL_EXTERNAL};



//////////////////////////////////////////////////////////////////////////////
// Default options

template <bool DONT_VERIFY_ = false, bool DONT_DUMP_RESULTS_ = false>
struct RazerSSpec
{
    enum {DONT_VERIFY = DONT_VERIFY_};                      // omit verifying potential matches
    enum {DONT_DUMP_RESULTS = DONT_DUMP_RESULTS_};          // omit dumping results
};

template <typename TSpec = RazerSSpec<> >
struct RazerSCoreOptions
{
    // major options
    AlignMode  alignMode;
    GapMode    gapMode;
    ScoreMode  scoreMode;

    // main options
    TSpec       spec;
    bool        forward;                // compute forward oriented read matches
    bool        reverse;                // compute reverse oriented read matches
    double      errorRate;              // Criteria 1 threshold
    unsigned    maxHits;                // output at most maxHits many matches
    unsigned    scoreDistanceRange;     // output only the best, second best, ..., scoreDistanceRange best matches
    int         dRange;                 // used in matchVerify
                                        // to a best match with e errors
    bool        purgeAmbiguous;         // true..remove reads with more than maxHits best matches, false..keep them
    CharString  output;                 // name of result file
    int         _debugLevel;            // level of verbosity
    bool        printVersion;           // print version number
    int         trimLength;             // if >0, cut reads to #trimLength characters
    // controlled pigeonhole extensions
    double      mutationRate;           // difference between reference genome and sequenced genome
    double      lossRate;               // 1.0 - sensitivity

    // output format options
    unsigned    outputFormat;           // 0..Razer format
                                        // 1..enhanced Fasta
                                        // 2..ELAND format
    bool        dumpAlignment;          // compute and dump the match alignments in the result files
    unsigned    genomeNaming;           // 0..use Fasta id
                                        // 1..enumerate reads beginning with 1
    // TODO(holtgrew): SAM export should imply --read-naming 3
    unsigned    readNaming;             // 0..use Fasta id
                                        // 1..enumerate reads beginning with 1
                                        // 2..use the read sequence (only for short reads!)
                                        // 3..use Fasta id, do not append /L and /R for mate pairs.
    bool        fullFastaId;            // read full FastaId or clip after first whitespace
    unsigned    sortOrder;              // 0..sort keys: 1. read number, 2. genome position
                                        // 1..           1. genome pos50ition, 2. read number
    int         positionFormat;         // 0..gap space
                                        // 1..position space
    const char  * runID;                // runID needed for gff output
    bool        dontShrinkAlignments;   // Required when used for building gold Rabema mapping.
    bool        computeGlobal;          // compute global alignment in SAM output

    // filtration parameters
    std::string shape;                  // shape (e.g. 11111111111)
    int         threshold;              // threshold (minimal threshold, 0 activates pigeonhole mode)
    int         tabooLength;            // taboo length
    int         repeatLength;           // repeat length threshold
    double      abundanceCut;           // abundance threshold
    int         delta;                  // q-gram delta (in pigeonhole mode), 0=automatic
    int         overlap;                // q-gram overlap (in pigeonhole mode), 0=lossless
    unsigned    maxOverlap;             // limits the overlap in automatic mode

    // mate-pair parameters
    int         libraryLength;          // offset between two mates
    int         libraryError;           // offset tolerance
    unsigned    nextPairMatchId;        // use this id for the next mate-pair

    // verification parameters
    unsigned    prefixSeedLength;       // length of the prefix seed
    bool        matchN;                 // false..N is always a mismatch, true..N matches with all
    unsigned char compMask[5];
    Score<int, Simple> scoringScheme;
    int         minScore;               // minimal alignment score

    // statistics
    typedef LogProb<> TProb;
//		typedef double TProb;
    String<unsigned> readLengths;       // read length histogram (i -> #reads of length i)
    String<double>   avrgQuality;       // average error quality per base
    String<TProb>    errorProb;         // error probability per base
    CharString  errorPrbFileName;
    CharString  mismatchFilename;

    String<double> errorDist;           // error distribution
    int64_t     countFiltration;        // matches returned by the filter
    int64_t     countVerification;      // matches returned by the verifier
    double      timeLoadFiles;          // time for loading input files
    double      timeMapReads;           // time for mapping reads
    double      timeDumpResults;        // time for dumping the results
    double      timeBuildQGramIndex;    // time for q-gram index building.
    double      timeCompactMatches;     // time for compacting reads
    double      timeMaskDuplicates;     // time spent masking duplicates
    double      timeFsCopy;             // time spent copying alignments back into the fragment store
    double      timeFiltration;
    double      timeVerification;

    bool        maqMapping;
    int         absMaxQualSumErrors;

    bool        lowMemory;              // set maximum shape weight to 13 to limit size of q-gram index
    bool        fastaIdQual;            // hidden option for special fasta+quality format we use

    // misc
    double      noCompactFrac;          // If in last noCompactFrac of genome, don't compact.
    double      compactMult;            // Multiplicator for compaction threshold.
    int64_t     compactThresh;          // compact match array if larger than compactThresh

    // multi-threading

    unsigned    threadCount;      // Number of threads to use in the parallel version.
    unsigned    windowSize;      // Collect SWIFT hits in windows of this length.
    unsigned    verificationPackageSize;      // This number of SWIFT hits per verification.
    unsigned    maxVerificationPackageCount;      // Maximum number of verification packages to create.
    int64_t     availableMatchesMemorySize;      // Memory available for matches.  Used for switching to external memory algorithms. -1 for always external, 0 for never.
    int         matchHistoStartThreshold;      // Threshold to use for starting histogram. >= 1

#ifdef RAZERS_OPENADDRESSING
    double      loadFactor;
#endif

    // global preprocessing information and maximal allowed errors

    typedef Infix<String<Dna5Q> >::Type         TRead;
    typedef Pattern<TRead const, MyersUkkonen>  TMyersPattern;      // verifier
    typedef String<TMyersPattern>               TPreprocessing;

    static String<unsigned char>                errorCutOff;        // ignore matches with >=errorCutOff errors
    static TPreprocessing                       forwardPatterns;

    CharString commandLine;
    std::string version;

    RazerSCoreOptions()
    {
        alignMode = RAZERS_GLOBAL;
        gapMode = RAZERS_GAPPED;
        scoreMode = RAZERS_ERRORS;

        forward = true;
        reverse = true;
        errorRate = 0.05;
        maxHits = 100;
        scoreDistanceRange = 0;     // disabled
        dRange = 1 << 30;
        purgeAmbiguous = false;
        output = "";
        _debugLevel = 0;
        printVersion = false;
        trimLength = 0;
        mutationRate = 0.05;

        outputFormat = 0;
        dumpAlignment = false;
        genomeNaming = 0;
        readNaming = 0;
        fullFastaId = false;
        sortOrder = 0;
        positionFormat = 0;
        runID = "s";        //
        dontShrinkAlignments = false;
        computeGlobal = false;

        matchN = false;
        shape = "11111111111";
        threshold = 1;
        tabooLength = 1;
        repeatLength = 1000;
        abundanceCut = 1;
        delta = 0;
        overlap = 0;
        maxOverlap = 10;

        libraryLength = 220;
        libraryError = 50;
        nextPairMatchId = 0;

        prefixSeedLength = 28;      // the "artificial" seed length that is used for mapping quality assignment
        for (unsigned i = 0; i < 4; ++i)
            compMask[i] = 1 << i;
        compMask[4] = 0;

        noCompactFrac = 0.05;
        compactMult = 2.2;
        compactThresh = 1024;
        // compactThresh = 40;

        absMaxQualSumErrors = 100;      // maximum for sum of mism qualities in total readlength
        lowMemory = false;          // set maximum shape weight to 13 to limit size of q-gram index
        fastaIdQual = false;

        threadCount = 1;
        // TODO(holtgrew): Tune this!
        windowSize = 500000;
        verificationPackageSize = 100;
        maxVerificationPackageCount = 100;
        availableMatchesMemorySize = 0;
        matchHistoStartThreshold = 5;

#ifdef RAZERS_OPENADDRESSING
        loadFactor = 1.6;
#endif

        lossRate = 0.0;
        minScore = 0;
        countFiltration = 0;
        countVerification = 0;
        timeLoadFiles = 0.0;
        timeMapReads = 0.0;
        timeDumpResults = 0.0;
        timeBuildQGramIndex = 0.0;
        timeCompactMatches = 0.0;
        timeMaskDuplicates = 0.0;
        timeFsCopy = 0.0;
        timeFiltration = 0.0;
        timeVerification = 0.0;
        maqMapping = false;
    }

};

template <typename TSpec = RazerSSpec<> >
struct RazerSOptions : RazerSCoreOptions<TSpec>
{
    typedef RazerSCoreOptions<TSpec> TCoreOptions;
    SeqFileIn readFile;           // left read's SeqFile (we have to keep it open and store it here to stream it only once)
};

template <typename TSpec>
String<unsigned char> RazerSCoreOptions<TSpec>::errorCutOff;

template <typename TSpec>
typename RazerSCoreOptions<TSpec>::TPreprocessing RazerSCoreOptions<TSpec>::forwardPatterns;

//////////////////////////////////////////////////////////////////////////////
// Typedefs

enum RAZERS_ERROR
{
    RAZERS_INVALID_OPTIONS = -1,
    RAZERS_READS_FAILED    = -2,
    RAZERS_GENOME_FAILED   = -3,
    RAZERS_INVALID_SHAPE   = -4
};

//////////////////////////////////////////////////////////////////////////////
// Definitions

template <typename TReadSet, typename TShape, typename TSpec>
struct Cargo<Index<TReadSet, IndexQGram<TShape, TSpec> > >
{
    typedef struct
    {
        double      abundanceCut;
        int         _debugLevel;
    } Type;
};

//////////////////////////////////////////////////////////////////////////////
// Memory tuning

#ifdef RAZERS_MEMOPT

template <typename TReadSet, typename TShape, typename TSpec>
struct SAValue<Index<TReadSet, IndexQGram<TShape, TSpec> > >
{
    typedef Pair<
        unsigned,
        unsigned,
        BitCompressed<24, 8>    // max. 16M reads of length < 256
        > Type;
};

#else

template <typename TReadSet, typename TShape, typename TSpec>
struct SAValue<Index<TReadSet, IndexQGram<TShape, TSpec> > >
{
    typedef Pair<
        unsigned,               // many reads
        unsigned                // of arbitrary length
        > Type;
};

#endif

template <>
struct Size<Dna5String>
{
    typedef unsigned Type;
};

template <typename TReadSet, typename TShape>
struct Size<Index<TReadSet, IndexQGram<TShape> > >
{
    typedef unsigned Type;
};

template <typename TReadSet, typename TShape>
struct Size<Index<TReadSet, IndexQGram<TShape, OpenAddressing> > >
{
    typedef unsigned Type;
};


#ifdef RAZERS_PRUNE_QGRAM_INDEX

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TReadSet, typename TShape, typename TSpec>
inline bool _qgramDisableBuckets(Index<TReadSet, IndexQGram<TShape, TSpec> > & index)
{
    typedef Index<TReadSet, IndexQGram<TShape, TSpec> > TReadIndex;
    typedef typename Fibre<TReadIndex, QGramDir>::Type      TDir;
    typedef typename Iterator<TDir, Standard>::Type         TDirIterator;
    typedef typename Value<TDir>::Type                      TSize;

    TDir & dir    = indexDir(index);
    bool result  = false;
    unsigned counter = 0;
    TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
    if (thresh < 100)
        thresh = 100;

    TDirIterator it = begin(dir, Standard());
    TDirIterator itEnd = end(dir, Standard());
    for (; it != itEnd; ++it)
        if (*it > thresh)
        {
            *it = (TSize) - 1;
            result = true;
            ++counter;
        }

    if (counter > 0 && cargo(index)._debugLevel >= 1)
        std::cerr << "Removed " << counter << " k-mers" << std::endl;

    return result;
}

#endif


template <
    typename TFragmentStore_,
    typename TMatches_,
    typename TRazerSOptions_,
    typename TRazerSMode_,
    typename TFilterPattern_,
    typename TCounts_
    >
struct MatchVerifier
{
    typedef TFragmentStore_                                 TFragmentStore;
    typedef TMatches_                                       TMatches;
    typedef TRazerSOptions_                                 TOptions;
    typedef TRazerSMode_                                    TRazerSMode;
    typedef TFilterPattern_                                 TFilterPattern;
    typedef TCounts_                                        TCounts;

    typedef typename TRazerSMode::TMatchNPolicy             TMatchNPolicy;

    typedef typename TFragmentStore::TReadSeqStore          TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type const       TRead;
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;
    typedef typename TFragmentStore::TContigSeq             TContigSeq;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedRead;
    typedef typename Value<TAlignQualityStore>::Type        TAlignQuality;
    typedef typename Size<TContigSeq>::Type                 TSize;
    typedef typename MakeSigned_<TSize>::Type               TContigPos;
    typedef ModifiedString<TRead, ModReverse>               TRevRead;

    typedef typename Value<TMatches>::Type TMatchRecord;

#ifdef RAZERS_BANDED_MYERS
    typedef PatternState_<TRead, Myers<AlignTextBanded<FindInfix, TMatchNPolicy, TMatchNPolicy>, True, void> > TPatternState;
    typedef PatternState_<TRevRead, Myers<AlignTextBanded<FindPrefix, TMatchNPolicy, TMatchNPolicy>, True, void> > TRPatternState;
#else  // #ifdef RAZERS_BANDED_MYERS
    typedef Pattern<TRead, Myers<FindInfix, False, void> >      TMyersPattern;
    typedef Pattern<TRevRead, Myers<FindInfix, False, void> >   TRevMyersPattern;
    typedef typename PatternState<TMyersPattern>::Type          TPatternState;
    typedef typename PatternState<TRevMyersPattern>::Type       TRPatternState;
#endif  // #ifdef RAZERS_BANDED_MYERS

    TMatches        * matches;
    TOptions        * options;              // RazerS options
    TFilterPattern  * filterPattern;
    TCounts         * cnts;

    TMatchRecord    m;
    TPatternState   patternState;
    TRPatternState  revPatternState;
    TSize           genomeLength;
    TSize           rightClip;
    TContigPos      sinkPos;
    bool            onReverseComplement;
    bool            oneMatchPerBucket;

    double compactionTime;

    MatchVerifier() :
        genomeLength(0), rightClip(0), sinkPos(std::numeric_limits<TContigPos>::max()), onReverseComplement(false), oneMatchPerBucket(false), compactionTime(0) {}

    MatchVerifier(TMatches_ & _matches, TOptions & _options, TFilterPattern & _filterPattern, TCounts & _cnts) :
        matches(&_matches),
        options(&_options),
        filterPattern(&_filterPattern),
        cnts(&_cnts)
    {
        genomeLength = 0;
        rightClip = 0;
        sinkPos = std::numeric_limits<TContigPos>::max() >> 1;
        onReverseComplement = false;
        oneMatchPerBucket = false;
        compactionTime = 0;
    }

    inline void push()
    {
        if (onReverseComplement)
        {
            // transform coordinates to the forward strand
            m.beginPos = genomeLength - m.beginPos;
            m.endPos = genomeLength - m.endPos;
            std::swap(m.beginPos, m.endPos);
            m.orientation = 'R';
        }
        else
        {
            m.orientation = 'F';
        }

//SEQAN_OMP_PRAGMA(critical)
// begin of critical section
        {
            if (!options->spec.DONT_DUMP_RESULTS)
            {
//                    std::cout << "begin: "<<m.beginPos <<"\tendPos: "<<m.endPos << "\terrors: "<<m.score <<std::endl;
                appendValue(*matches, m, Generous());

                if ((int64_t)length(*matches) > options->compactThresh)
                {
                    double beginTime = sysTime();
                    typename Size<TMatches>::Type oldSize = length(*matches);

                    if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE || options->threshold == 0)
                        maskDuplicates(*matches, *options, TRazerSMode());      // overlapping parallelograms cause duplicates
                    // SEQAN_ASSERT_MSG((back(*matches).endPos - back(*matches).beginPos == 100), "len == %d", int(m.endPos - m.beginPos));

                    compactMatches(*matches, *cnts, *options, TRazerSMode(), *filterPattern, COMPACT);
                    // SEQAN_ASSERT_MSG((back(*matches).endPos - back(*matches).beginPos == 100), "len == %d", int(m.endPos - m.beginPos));

                    if (length(*matches) * 4 > oldSize)                 // the threshold should not be raised
                    {       // fprintf(stderr, "[raising threshold]");
                            // options->compactThresh += (options->compactThresh >> 1);	// if too many matches were removed
                        options->compactThresh = (int64_t)(options->compactThresh * options->compactMult);
                    }

//						if (options._debugLevel >= 2)
//							std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
                    double endTime = sysTime();
                    compactionTime += (endTime - beginTime);
                }
            }
            ++options->countVerification;
        }
// end of critical section
    }

};



//////////////////////////////////////////////////////////////////////////////
// Read a list of genome file names
template<typename TSpec>
int getGenomeFileNameList(CharString filename, StringSet<CharString> & genomeFileNames, RazerSCoreOptions<TSpec> &options)
{
	std::ifstream file;
	file.open(toCString(filename), std::ios_base::in | std::ios_base::binary);
	if (!file.is_open())
		return RAZERS_GENOME_FAILED;

    DirectionIterator<std::fstream, Input>::Type reader(file);
    if (!atEnd(reader))
        return 0;

    clear(genomeFileNames);
	if (*reader == '>' && *reader != '@')	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		if(options._debugLevel >=1)
			std::cout << std::endl << "Reading multiple genome files:" << std::endl;
		
		unsigned i = 1;
        CharString line;
		while (!atEnd(reader))
		{
            readLine(line, reader);
            cropOuter(line, IsWhitespace());
			appendValue(genomeFileNames, line);
			if(options._debugLevel >=2)
				std::cout <<"Genome file #"<< i <<": " << back(genomeFileNames) << std::endl;
			++i;
		}
		if(options._debugLevel >=1)
			std::cout << i-1 << " genome files total." << std::endl;
	}
	else		//if file starts with a fasta header --> regular one-genome-file input
		appendValue(genomeFileNames, filename, Exact());
	file.close();
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences with or w/o quality values
template <typename TFSSpec, typename TFSConfig, typename TRazerSOptions>
bool loadReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
	SeqFileIn &seqFile,
    TRazerSOptions & options)
{
    bool countN = !(options.matchN || options.outputFormat == 1);

    String<uint64_t> qualSum;
    String<Dna5Q>    seq;
    CharString       qual;
    CharString       seqId;

    unsigned seqCount = 0;
    unsigned kickoutcount = 0;

    while (!atEnd(seqFile))
    {
        readRecord(seqId, seq, qual, seqFile);
        ++seqCount;

        if (options.readNaming == 0 || options.readNaming == 3)
        {
            if (!options.fullFastaId)
                cropAfterFirst(seqId, IsWhitespace());  // read Fasta id up to the first whitespace
        }
        else
        {
            clear(seqId);
        }

        if (countN)
        {
            int count = 0;
            int cutoffCount = (int)(options.errorRate * length(seq));
            for (unsigned j = 0; j < length(seq); ++j)
                if (getValue(seq, j) == 'N')
                    if (++count > cutoffCount)
                    {
                        clear(seq);
                        clear(seqId);
                        clear(qual);  // So no qualities are assigned below.
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

        // append read to fragment store
        appendRead(store, seq, seqId);

        unsigned len = length(seq);
        if (length(qualSum) <= len)
        {
            resize(qualSum, len, 0u);
            resize(options.readLengths, len + 1, 0u);
        }
        ++options.readLengths[len];

        for (unsigned i = 0; i < len; ++i)
            qualSum[i] += getQualityValue(seq[i]);
    }

    // memory optimization
    // we store reads in a concat-direct stringset and can shrink its size
    shrinkToFit(store.readSeqStore.concat);
    shrinkToFit(store.readSeqStore.limits);

    // compute error probabilities
    resize(options.avrgQuality, length(qualSum));
    unsigned coverage = 0;
    for (unsigned i = length(qualSum); i != 0; )
    {
        coverage += options.readLengths[i];
        --i;
        options.avrgQuality[i] = (double)qualSum[i] / (double)coverage;
    }
    estimateErrorDistributionFromQualities(options);

    typedef Shape<Dna, SimpleShape> TShape;
    typedef typename SAValue<Index<StringSet<Dna5String>, IndexQGram<TShape, OpenAddressing> > >::Type TSAValue;
    TSAValue sa(0, 0);
    sa.i1 = ~sa.i1;
    sa.i2 = ~sa.i2;

    if ((unsigned)sa.i1 < length(store.readSeqStore) - 1)
    {
        std::cerr << "Maximal read number of " << (unsigned)sa.i1 + 1 << " exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << std::endl;
        seqCount = 0;
    }
    if ((unsigned)sa.i2 < length(options.readLengths) - 2)
    {
        std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << std::endl;
        seqCount = 0;
    }

    if (options._debugLevel > 1 && kickoutcount > 0)
        std::cerr << "Ignoring " << kickoutcount << " low quality reads.\n";

    if (options._debugLevel > 1)
    {
        std::cerr << std::endl;
        std::cerr << "Average quality profile:" << std::endl;
        unsigned cols = length(options.avrgQuality);
        if (cols > 40)
            cols = 40;

        for (int j = 20; j >= 0; --j)
        {
            std::cout.width(3);
            if (j % 5 == 0)
                std::cout << 2 * j;
            else
                std::cout << ' ';
            std::cout << " | ";
            for (unsigned i = 0; i < cols; ++i)
            {
                unsigned c = i * (length(options.avrgQuality) - 1) / (cols - 1);
                double x = options.avrgQuality[c];
                std::cout << ((2 * j + 1 <= x) ? '*' : ' ');
            }
            std::cout << std::endl;
        }
        std::cout << "    +-";
        for (unsigned i = 0; i < cols; ++i)
            std::cout << ((i % 5 == 0) ? '+' : '-');
        std::cout << std::endl << "  ";
        for (unsigned i = 0; i < cols; ++i)
        {
            unsigned c = i * (length(options.avrgQuality) - 1) / (cols - 1);
            std::cout.width(5);
            if (i % 5 == 0)
                std::cout << c + 1;
        }
        std::cout << std::endl;
    }

    return seqCount > 0;
}

//////////////////////////////////////////////////////////////////////////////
// Read the first sequence of a multi-sequence file
// and return its length
inline int estimateReadLength(SeqFileIn &seqFile)
{
	if (atEnd(seqFile))
		return 0;

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

// Comparators for RazerS1-style matches.

// TODO(holtgrew): Slightly different comparators than in previous RazerS 3 version, add back the additional checks?

template <typename TReadMatch>
struct LessBeginPos :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    inline bool operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // genome position and orientation
        if (a.contigId < b.contigId) return true;
        if (a.contigId > b.contigId) return false;

        if (a.beginPos < b.beginPos) return true;

        return false;
    }

};

template <typename TReadMatch>
struct LessRNoBeginPos :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    inline bool operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (a.readId < b.readId) return true;
        if (a.readId > b.readId) return false;

        // genome position and orientation
        if (a.contigId < b.contigId) return true;
        if (a.contigId > b.contigId) return false;

        if (a.beginPos < b.beginPos) return true;
        if (a.beginPos > b.beginPos) return false;

        if (a.orientation == '-') return false;
        if (b.orientation == '-') return true;

        if (a.orientation < b.orientation) return true;
        if (a.orientation > b.orientation) return false;

        // quality
        if (a.score > b.score) return true;
        if (b.score > a.score) return false;

        if (a.endPos > b.endPos) return true;

        return false;
    }

};

template <typename TReadMatch>
struct LessRNoBeginPosMP :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    inline bool operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (a.readId < b.readId) return true;
        if (a.readId > b.readId) return false;

        // genome position and orientation
        if (a.contigId < b.contigId) return true;
        if (a.contigId > b.contigId) return false;

        if (a.beginPos < b.beginPos) return true;
        if (a.beginPos > b.beginPos) return false;

        if (a.orientation == '-') return false;
        if (b.orientation == '-') return true;

        if (a.orientation < b.orientation) return true;
        if (a.orientation > b.orientation) return false;

        // quality
        if (a.pairScore > b.pairScore) return true;
        if (a.pairScore < b.pairScore) return false;

        if (a.libDiff < b.libDiff) return true;
        if (a.libDiff > b.libDiff) return false;

        if (a.endPos > b.endPos) return true;

        return false;
    }

};

// ... to sort matches and remove duplicates with equal gEnd
template <typename TReadMatch>
struct LessRNoEndPos :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    inline bool operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (a.readId < b.readId) return true;
        if (a.readId > b.readId) return false;

        // genome position and orientation
        if (a.contigId < b.contigId) return true;
        if (a.contigId > b.contigId) return false;

        if (a.endPos   < b.endPos) return true;
        if (a.endPos   > b.endPos) return false;

        if (a.orientation == '-') return false;
        if (b.orientation == '-') return true;

        if (a.orientation < b.orientation) return true;
        if (a.orientation > b.orientation) return false;

        // quality
        if (a.score > b.score) return true;
        if (b.score > a.score) return false;

        if (a.beginPos < b.beginPos) return true;

        return false;
    }

};

template <typename TReadMatch>
struct LessRNoEndPosMP :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    int libSize;
    LessRNoEndPosMP(int _libSize) :
        libSize(_libSize) {}

    inline bool operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (a.readId < b.readId) return true;
        if (a.readId > b.readId) return false;

        // genome position and orientation
        if (a.contigId < b.contigId) return true;
        if (a.contigId > b.contigId) return false;

        if (a.endPos   < b.endPos) return true;
        if (a.endPos   > b.endPos) return false;

        if (a.orientation == '-') return false;
        if (b.orientation == '-') return true;

        if (a.orientation < b.orientation) return true;
        if (a.orientation > b.orientation) return false;

        // quality
        if (a.pairScore > b.pairScore) return true;
        if (a.pairScore < b.pairScore) return false;

        if (a.libDiff < b.libDiff) return true;
        if (a.libDiff > b.libDiff) return false;

        if (a.beginPos < b.beginPos) return true;

        return false;
    }

};

template <typename TReadMatch>
struct LessScoreBackport :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    inline bool operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (a.readId < b.readId) return true;
        if (a.readId > b.readId) return false;

        // quality
        if (a.orientation == '-') return false;
        if (b.orientation == '-') return true;

        if (a.score > b.score) return true;
        if (b.score > a.score) return false;

        // Sort by leftmost begin pos, longest end pos on ties.
        if (a.contigId < b.contigId) return true;
        if (a.contigId > b.contigId) return false;

        if (a.orientation < b.orientation) return true;
        if (a.orientation > b.orientation) return false;

        if (a.beginPos < b.beginPos) return true;
        if (a.beginPos > b.beginPos) return false;

        if (a.endPos < b.endPos) return false;
        if (a.endPos > b.endPos) return true;

        return false;
    }

};

// TODO(holtgrew): Merge with above.

template <typename TReadMatch>
struct LessScoreBackport3Way :
    public std::binary_function<TReadMatch, TReadMatch, int>
{
    inline int operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (a.readId < b.readId) return -1;
        if (a.readId > b.readId) return 1;

        // quality
        if (a.orientation != '-' || b.orientation != '-')
        {
            if (a.orientation == '-') return -1;
            if (b.orientation == '-') return 1;
        }

        if (a.score > b.score) return -1;
        if (b.score > a.score) return 1;

        // Sort by leftmost begin pos, longest end pos on ties.
        if (a.contigId < b.contigId) return -1;
        if (a.contigId > b.contigId) return 1;

        if (a.orientation < b.orientation) return -1;
        if (a.orientation > b.orientation) return 1;

        if (a.beginPos < b.beginPos) return -1;
        if (a.beginPos > b.beginPos) return 1;

        if (a.endPos < b.endPos) return 1;
        if (a.endPos > b.endPos) return -1;

        return 0;
    }

};


// Comparators for Fragment Store

template <typename TAlignedReadStore, typename TLessScore>
struct LessRNoGPos :
    public std::binary_function<typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool>
{
    typedef typename Value<TAlignedReadStore>::Type TAlignedRead;
    TLessScore lessScore;

    LessRNoGPos(TLessScore const & _lessScore) :
        lessScore(_lessScore) {}

    inline bool operator()(TAlignedRead const & a, TAlignedRead const & b) const
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

        int result = lessScore.compare(a, b);
        if (result == 0)
        {
            // prefer reads that support more of the reference
            return _max(a.beginPos, a.endPos) > _max(b.beginPos, b.endPos);
        }
        return result == -1;
    }

};

// ... to sort matches and remove duplicates with equal gEnd
template <typename TAlignedReadStore, typename TLessScore>
struct LessRNoGEndPos :
    public std::binary_function<typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool>
{
    typedef typename Value<TAlignedReadStore>::Type TAlignedRead;
    TLessScore lessScore;

    LessRNoGEndPos(TLessScore const & _lessScore) :
        lessScore(_lessScore) {}

    inline bool operator()(
        typename Value<TAlignedReadStore>::Type const & a,
        typename Value<TAlignedReadStore>::Type const & b) const
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

        int result = lessScore.compare(a, b);
        if (result == 0)
        {
            // prefer reads that support more of the reference
            return _min(a.beginPos, a.endPos) < _min(b.beginPos, b.endPos);
        }
        return result == -1;
    }

};

template <typename TAlignedReadStore, typename TAlignedReadQualityStore, typename TRazerSMode>
struct LessScore :
    public std::binary_function<typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool>
{
    TAlignedReadQualityStore & qualStore;

    LessScore(TAlignedReadQualityStore & _qualStore) :
        qualStore(_qualStore) {}

    inline int compare(
        typename Value<TAlignedReadStore>::Type const & a,
        typename Value<TAlignedReadStore>::Type const & b) const
    {
        typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

        // read number
        if (a.readId < b.readId) return -1;
        if (a.readId > b.readId) return 1;

        // quality
        if (a.id == TAlignedRead::INVALID_ID) return 1;
        if (b.id == TAlignedRead::INVALID_ID) return -1;

        typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
        typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
        if (qa.pairScore > qb.pairScore) return -1;
        if (qa.pairScore < qb.pairScore) return 1;

        if (qa.score > qb.score) return -1;
        if (qb.score > qa.score) return 1;

        return 0;
    }

    inline bool operator()(
        typename Value<TAlignedReadStore>::Type const & a,
        typename Value<TAlignedReadStore>::Type const & b) const
    {
        return compare(a, b) == -1;
    }

};

// longest prefix mapping
template <typename TAlignedReadStore, typename TAlignedReadQualityStore, typename TGapMode, typename TScoreMode, typename TMatchNPolicy>
struct LessScore<TAlignedReadStore, TAlignedReadQualityStore, RazerSMode<RazerSPrefix, TGapMode, TScoreMode, TMatchNPolicy> >:
    public std::binary_function<typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool>
{
    TAlignedReadQualityStore & qualStore;

    LessScore(TAlignedReadQualityStore & _qualStore) :
        qualStore(_qualStore) {}

    inline int compare(
        typename Value<TAlignedReadStore>::Type const & a,
        typename Value<TAlignedReadStore>::Type const & b) const
    {
        typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

        // read number
        if (a.readId < b.readId) return -1;
        if (a.readId > b.readId) return 1;

        // quality
        if (a.id == TAlignedRead::INVALID_ID) return 1;
        if (b.id == TAlignedRead::INVALID_ID) return -1;

        typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
        typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
        if (qa.errors < qb.errors) return -1;
        if (qa.errors > qb.errors) return 1;

        if (qa.score > qb.score) return -1;
        if (qb.score > qa.score) return 1;

        return 0;
    }

    inline bool operator()(
        typename Value<TAlignedReadStore>::Type const & a,
        typename Value<TAlignedReadStore>::Type const & b) const
    {
        return compare(a, b) == -1;
    }

};

//////////////////////////////////////////////////////////////////////////////

template <typename TAlignedReadQualityStore, typename TRazerSMode>
struct BinFunctorDefault
{
    TAlignedReadQualityStore & qualStore;

    BinFunctorDefault(TAlignedReadQualityStore & _qualStore) :
        qualStore(_qualStore) {}

    template <typename TAlignedRead>
    inline int operator()(TAlignedRead & alignedRead) const
    {
        return qualStore[alignedRead.id].errors;
    }

};


//////////////////////////////////////////////////////////////////////////////
// Mark duplicate matches for deletion
template <typename TMatches, typename TIterator, typename TOptions, typename TRazerSMode>
void maskDuplicates(TMatches &, TIterator const itBegin, TIterator const itEnd, TOptions & options, TRazerSMode)
{
    typedef typename Value<TMatches>::Type  TMatch;
    typedef typename TMatch::TContigPos     TContigPos;

    TContigPos  beginPos, endPos;
    unsigned    contigId, readId;
    char        orientation;
    unsigned    masked;
    TIterator   it;
    double      beginTime = sysTime();

    //////////////////////////////////////////////////////////////////////////////
    // remove matches with equal ends

    // we can skip one sort step in no-gap mode and with pigeonhole filter
    if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE || options.threshold != 0)
    {
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_SORT);
#endif
        if (options.libraryLength >= 0)
            std::sort(itBegin, itEnd, LessRNoEndPosMP<TMatch>(options.libraryLength));
        else
            std::sort(itBegin, itEnd, LessRNoEndPos<TMatch>());
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_SORT);
#endif

        beginPos = -1;
        endPos = -1;
        contigId = TMatch::INVALID_ID;
        readId = TMatch::INVALID_ID;
        orientation = '-';
        masked = 0;
        it = itBegin;

        for (; it != itEnd; ++it)
        {
            if ((*it).pairMatchId != TMatch::INVALID_ID && (it->readId & 1) != 0)
                continue;                                                                   // remove only single reads or left mates

            TContigPos itEndPos = _max((*it).beginPos, (*it).endPos);
            if (endPos == itEndPos && orientation == (*it).orientation &&
                contigId == (*it).contigId && readId == (*it).readId)
            {
                (*it).orientation = '-';
                masked += 1;
                continue;
            }
            readId = (*it).readId;
            contigId = (*it).contigId;
            endPos = itEndPos;
            orientation = (*it).orientation;
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    // remove matches with equal begins

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
    if (options.libraryLength >= 0)
        std::sort(itBegin, itEnd, LessRNoBeginPosMP<TMatch>());
    else
        std::sort(itBegin, itEnd, LessRNoBeginPos<TMatch>());
    // std::cerr << "(SORTING " << itEnd-itBegin << " MATCHES)";
    // sortAlignedReads(store.alignedReadStore, TLessBeginPos(TLessScore(store.alignQualityStore)));
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

    beginPos = -1;
    endPos = -1;
    contigId = TMatch::INVALID_ID;
    readId = TMatch::INVALID_ID;
    orientation = '-';
    masked = 0;
    it = itBegin;

    for (; it != itEnd; ++it)
    {
        if ((*it).orientation == '-')
            continue;
        if ((*it).pairMatchId != TMatch::INVALID_ID && (it->readId & 1) != 0)
            continue;                                                                   // remove only single reads or left mates

        TContigPos itBeginPos = _min((*it).beginPos, (*it).endPos);
        if (beginPos == itBeginPos && readId == (*it).readId &&
            contigId == (*it).contigId && orientation == ((*it).beginPos < (*it).endPos))
        {
            (*it).orientation = '-';
            masked += 1;
            continue;
        }
        readId = (*it).readId;
        contigId = (*it).contigId;
        beginPos = itBeginPos;
        orientation = (*it).beginPos < (*it).endPos;
    }

#ifdef RAZERS_DEFER_COMPACTION
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
    //////////////////////////////////////////////////////////////////////////////
    // sort matches by begin position when using defered compaction
    std::sort(itBegin, itEnd, LessBeginPos<TMatch>());
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
#endif  // #ifdef RAZERS_DEFER_COMPACTION

    options.timeMaskDuplicates += sysTime() - beginTime;
    if (options._debugLevel >= 2)
        fprintf(stderr, " [%u matches masked]", masked);
}

template <typename TMatches, typename TOptions, typename TRazerSMode>
void maskDuplicates(TMatches & matches, TOptions & options, TRazerSMode const & mode)
{
    maskDuplicates(matches, begin(matches, Standard()), end(matches, Standard()), options, mode);
}

/*
//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template <typename TFragmentStore, typename TCounts, typename TBinFunctor, typename TRazerSMode>
void countMatches(TFragmentStore &store, TCounts &cnt, TBinFunctor &binF, TRazerSMode)
{
    typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
    typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;

    typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
    typedef typename Value<TCounts>::Type							TRow;
    typedef typename Value<TRow>::Type								TValue;

    sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));

    TIterator it = begin(store.alignedReadStore, Standard());
    TIterator itEnd = end(store.alignedReadStore, Standard());

    unsigned readId = TAlignedRead::INVALID_ID;
    int lastBin = -1;
    int64_t count = 0;

    String<TValue> row, empty;
    for (; it != itEnd; ++it)
    {
        if ((*it).id == TAlignedRead::INVALID_ID) continue;
        int bin = binF((*it).id);

        if (readId == (*it).readId)
        {
            if (lastBin == bin)
                ++count;
            else
            {
                appendValue(row, TValue(bin, count), Generous());
                lastBin = bin;
                count = 1;
            }
        }
        else
        {
            while (length(cnt) < readId)
                appendValue(cnt, empty, Generous());
            appendValue(cnt, row, Generous());
            clear(row);
            readId = (*it).readId;
            lastBin = bin;
            count = 1;
        }
    }
    while (length(cnt) < readId)
        appendValue(cnt, empty, Generous());
    appendValue(cnt, row, Generous());
}
*/
//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template <typename TFragmentStore, typename TCounts, typename TRazerSMode>
void countMatches(TFragmentStore & store, TCounts & cnt, TRazerSMode const &)
{
    typedef typename TFragmentStore::TAlignedReadStore              TAlignedReadStore;
    typedef typename TFragmentStore::TAlignQualityStore             TAlignQualityStore;

    typedef typename Value<TAlignedReadStore>::Type                 TAlignedRead;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TIterator;
    typedef typename Value<TCounts>::Type                           TRow;
    typedef typename Value<TRow>::Type                              TValue;

    TIterator it = begin(store.alignedReadStore, Standard());
    TIterator itEnd = end(store.alignedReadStore, Standard());

    unsigned readId = TAlignedRead::INVALID_ID;
    short errors = -1;
    int64_t count = 0;
    int64_t maxVal = std::numeric_limits<TValue>::max();

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
    std::sort(begin(store.alignedReadStore, Standard()), end(store.alignedReadStore, Standard()), LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
    //sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

    for (; it != itEnd; ++it)
    {
        if ((*it).id == TAlignedRead::INVALID_ID)
            continue;
        if (readId == (*it).readId && errors == store.alignQualityStore[(*it).id].errors)
            ++count;
        else
        {
            if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < length(cnt))
                cnt[errors][readId] = (maxVal < count) ? (TValue)maxVal : (TValue)count;
            readId = (*it).readId;
            errors = store.alignQualityStore[(*it).id].errors;
            count = 1;
        }
    }
    if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < length(cnt))
        cnt[errors][readId] = (TValue)count;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFilterPattern, typename TReadNo, typename TMaxErrors>
inline void
setMaxErrors(TFilterPattern &, TReadNo, TMaxErrors)
{}

template <typename TIndex, typename TPigeonholeSpec, typename TReadNo, typename TMaxErrors>
inline void
setMaxErrors(Pattern<TIndex, Pigeonhole<TPigeonholeSpec> > & filterPattern, TReadNo readNo, TMaxErrors maxErrors)
{
    if (maxErrors < 0)
        maskPatternSequence(filterPattern, readNo, maxErrors >= 0);
}

template <typename TIndex, typename TSwiftSpec, typename TReadNo, typename TMaxErrors>
inline void
setMaxErrors(Pattern<TIndex, Swift<TSwiftSpec> > & filterPattern, TReadNo readNo, TMaxErrors maxErrors)
{
    // if (readNo==643)
    //  std::cout<<"dman"<<std::endl;
    int minT = _qgramLemma(filterPattern, readNo, maxErrors);
    if (minT > 1)
    {
//		std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
        if (maxErrors < 0)
            minT = std::numeric_limits<int>::max();
        setMinThreshold(filterPattern, readNo, (unsigned)minT);
    }
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template <
    typename TMatches,
    typename TCounts,
    typename TSpec,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TSwift,
    typename TMatchNPolicy
    >
void compactMatches(
    TMatches & matches,
    TCounts &,
    RazerSCoreOptions<TSpec> & options,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const &,
    TSwift & swift,
    CompactMatchesMode compactMode)
{
    // fprintf(stderr, "[compact]");
    double beginTime = sysTime();
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename Iterator<TMatches, Standard>::Type TIterator;
    //typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;

    unsigned readNo = -1;
    unsigned hitCount = 0;
    unsigned hitCountCutOff = options.maxHits;
    int scoreCutOff = std::numeric_limits<int>::min();
    int scoreRangeBest = (IsSameType<TAlignMode, RazerSGlobal>::VALUE && !IsSameType<TScoreMode, RazerSScore>::VALUE) ? -(int)options.scoreDistanceRange : std::numeric_limits<int>::max();
    ignoreUnusedVariableWarning(scoreRangeBest);
    ignoreUnusedVariableWarning(compactMode);

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
#ifdef RAZERS_EXTERNAL_MATCHES
    if (compactMode == COMPACT_FINAL_EXTERNAL)
    {
        typedef Pipe<TMatches, Source<> > TSource;
        typedef LessScoreBackport3Way<TMatch> TLess;
        typedef Pool<TMatch, SorterSpec<SorterConfigSize<TLess, typename Size<TSource>::Type> > > TSorterPool;

        TSource source(matches);
        TSorterPool sorter;
        sorter << source;
        matches << sorter;

        for (unsigned i = 1; i < length(matches); ++i)
            SEQAN_ASSERT_LEQ(TLess() (matches[i - 1], matches[i]), 0);
    }
    else
    {
#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
    std::sort(begin(matches, Standard()), end(matches, Standard()), LessScoreBackport<TMatch>());
    // sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
#ifdef RAZERS_EXTERNAL_MATCHES
}

#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

    TIterator it = begin(matches, Standard());
    TIterator itEnd = end(matches, Standard());
    TIterator dit = it;
    TIterator ditBeg = it;
    // fprintf(stderr, "[%u matches to compact]", unsigned(itEnd - it));
    unsigned disabled = 0;

    for (; it != itEnd; ++it)
    {
        if ((*it).orientation == '-')
            continue;
        int score = (*it).score;
        int errors = -(*it).score;
        ignoreUnusedVariableWarning(errors);

        //if (readNo == 643) std::cerr <<"["<<score<<","<<errors<<"] "<<std::flush;
        if (readNo == (*it).readId && (*it).pairMatchId == TMatch::INVALID_ID)  // Only compact unpaired matches.
        {
            if (score <= scoreCutOff)
            {
//    std::cout<<"decreased errCutOff["<<readNo<<"] from "<<(unsigned)options.errorCutOff[readNo];
                options.errorCutOff[readNo] = -scoreCutOff;
//    std::cout<<"to "<<(unsigned)options.errorCutOff[readNo]<<std::endl;
                setMaxErrors(swift, readNo, -scoreCutOff - 1);
                while (it != itEnd && readNo == (*it).readId)
                    ++it;
                --it;
                continue;
            }
            if (++hitCount >= hitCountCutOff)
            {
#ifdef RAZERS_MASK_READS
                if (hitCount == hitCountCutOff)
                {
                    // we have enough, now look for better matches
                    if (options.purgeAmbiguous && (options.scoreDistanceRange == 0 || score > scoreRangeBest))
                    {
                        // std::cerr << "PURGED " << readNo << std::endl;
//    std::cout<<"decreased errCutOff["<<readNo<<"] from "<<(unsigned)options.errorCutOff[readNo];
                        options.errorCutOff[readNo] = 0;
//    std::cout<<"to "<<(unsigned)options.errorCutOff[readNo]<<std::endl;
                        setMaxErrors(swift, readNo, -1);
                        ++disabled;
                        // if (options._debugLevel >= 2)
                        //  std::cerr << "(read #" << readNo << " disabled)";
                    }
                    else
                    // we only need better matches
                    if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
                    {
                        // std::cerr << "LIMITED " << readNo << std::endl;
//    std::cout<<"decreased errCutOff["<<readNo<<"] from "<<(unsigned)options.errorCutOff[readNo];
                        options.errorCutOff[readNo] = errors;
//    std::cout<<"to "<<(unsigned)options.errorCutOff[readNo]<<std::endl;
                        setMaxErrors(swift, readNo, errors - 1);
                        if (errors == 0)
                            ++disabled;
                        // if (errors == 0 && options._debugLevel >= 2)
                        //  std::cerr << "(read #" << readNo << " disabled)";
                    }

                    if (options.purgeAmbiguous && (compactMode == COMPACT_FINAL || compactMode == COMPACT_FINAL_EXTERNAL))
                    {
                        if (options.scoreDistanceRange == 0 || score > scoreRangeBest || compactMode == COMPACT_FINAL || compactMode == COMPACT_FINAL_EXTERNAL)
                        {
                            dit = ditBeg;
                        }
                        else
                        {
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
            if (options.scoreDistanceRange > 0)
                scoreCutOff = score - options.scoreDistanceRange;
            ditBeg = dit;
        }
        *dit = *it;
        ++dit;
    }
    unsigned origSize = length(matches);
    resize(matches, dit - begin(matches, Standard()));
    // compactAlignedReads(store);
    options.timeCompactMatches += sysTime() - beginTime;
    // fprintf(stderr, "[compacted in %f s]", endTime - beginTime);
    unsigned newSize = length(matches);
    if (options._debugLevel >= 2)
    {
        fprintf(stderr, " [%u matches removed]", unsigned(origSize - newSize));
        if (disabled > 0)
            fprintf(stderr, " [%u reads disabled]", disabled);
    }
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template <
    typename TMatches,
    typename TCounts,
    typename TSpec,
    typename TGapMode,
    typename TSwift,
    typename TMatchNPolicy
    >
void compactMatches(
    TMatches & matches,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> &,
    RazerSMode<RazerSGlobal, TGapMode, RazerSQuality<RazerSMAQ>, TMatchNPolicy> const &,
    TSwift &,
    CompactMatchesMode compactMode)
{
    typedef typename Value<TMatches>::Type                      TMatch;
    typedef typename Iterator<TMatches, Standard>::Type         TIterator;
    //typedef RazerSMode<RazerSGlobal, TGapMode, RazerSQuality<RazerSMAQ>, TMatchNPolicy> TRazerSMode;

    unsigned readNo = -1;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
    std::sort(
        begin(matches, Standard()),
        end(matches, Standard()),
        LessScoreBackport<TMatch>());
    // sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

    TIterator it = begin(matches, Standard());
    TIterator itEnd = end(matches, Standard());
    TIterator dit = it;

    //number of errors may not exceed 31!
    bool second = true;
    for (; it != itEnd; ++it)
    {
        if ((*it).orientation == '-')
            continue;
        if (readNo == (*it).readId)
        {
            //second best match
            if (second)
            {
                second = false;
                if ((cnts[1][(*it).readId] & 31)  > (*it).editDist)
                {
                    //this second best match is better than any second best match before
                    cnts[1][(*it).readId] = (*it).editDist; // second best dist is this editDist
                    // count is 0 (will be updated if countFirstTwo)
                }
                if (compactMode == COMPACT_FINAL)
                    if ((cnts[1][(*it).readId] >> 5) != 2047)
                        cnts[1][(*it).readId] += 32;
            }
            else
            {
                if ((*it).editDist <= (cnts[0][(*it).readId] & 31))
                    if (cnts[0][(*it).readId] >> 5 != 2047)
                        cnts[0][(*it).readId] += 32;
                if ((*it).editDist <= (cnts[1][(*it).readId] & 31))
                    if ((cnts[1][(*it).readId] >> 5) != 2047)
                        cnts[1][(*it).readId] += 32;
                continue;
            }
        }
        else //best match
        {
            second = true;
            readNo = (*it).readId;
            //cnts has 16bits, 11:5 for count:dist
            if ((cnts[0][(*it).readId] & 31)  > (*it).editDist)
            {
                //this match is better than any match before
                cnts[1][(*it).readId] = cnts[0][(*it).readId]; // best before is now second best
                // (count will be updated when match is removed)
                cnts[0][(*it).readId] = (*it).editDist; // best dist is this editDist
                // count is 0 (will be updated if countFirstTwo)
            }
            if (compactMode == COMPACT_FINAL)
                if ((cnts[0][(*it).readId] >> 5) != 2047)
                    cnts[0][(*it).readId] += 32;
            // shift 5 to the right, add 1, shift 5 to the left, and keep count
        }
        *dit = *it;
        ++dit;
    }

    resize(matches, dit - begin(matches, Standard()));
}

/* // fallback
template <
    typename TFragmentStore,
    typename TCounts,
    typename TSpec,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TSwift
>
void compactMatches(
    TFragmentStore &,
    TCounts &,
    RazerSCoreOptions<TSpec> &,
    RazerSMode<TAlignMode, TGapMode, TScoreMode> const,
    TSwift &,
    CompactMatchesMode)
{
}
*/

//////////////////////////////////////////////////////////////////////////////
// Best Hamming prefix verification
template <
    typename TMatchVerifier,
    typename TGenome,
    typename TRead,
    typename TMatchNPolicy>
inline bool
matchVerify(
    TMatchVerifier & verifier,
    Segment<TGenome, InfixSegment> inf,                                 // potential match genome region
    unsigned readId,                                                    // read number
    TRead const & read,                                                 // read
    RazerSMode<RazerSPrefix, RazerSUngapped, RazerSErrors, TMatchNPolicy> const &)  // Hamming only
{
    typedef Segment<TGenome, InfixSegment>                  TGenomeInfix;
    typedef typename Iterator<TGenomeInfix, Standard>::Type TGenomeIterator;
    typedef typename Iterator<TRead const, Standard>::Type  TReadIterator;

//	unsigned maxErrors = (unsigned)(verifier.options->prefixSeedLength * verifier.options->errorRate);
    unsigned maxErrors = verifier.options->errorCutOff[readId];
    if (maxErrors == 0)
        return false;

    --maxErrors;

    // verify
    TReadIterator ritBeg    = begin(read, Standard());
    TReadIterator ritEnd    = end(read, Standard());
    unsigned ndlLength      = ritEnd - ritBeg;

    if (length(inf) < ndlLength)
        return false;

    TGenomeIterator git     = begin(inf, Standard());
    TGenomeIterator gitEnd  = end(inf, Standard()) - (ndlLength - 1);

    unsigned errorThresh = (verifier.oneMatchPerBucket) ? std::numeric_limits<unsigned>::max() : maxErrors;
    unsigned minErrors = maxErrors + 2;
    int bestHitLength = 0;

    for (; git < gitEnd; ++git)
    {
        unsigned errors = 0;
        TReadIterator r = ritBeg;
        TGenomeIterator g = git;
        for (; r != ritEnd; ++r, ++g)
            if ((verifier.options->compMask[ordValue(*g)] & verifier.options->compMask[ordValue(*r)]) == 0)
            {
                if (r - ritBeg < (int)verifier.options->prefixSeedLength)   // seed
                {
                    if (++errors > maxErrors)               // doesn't work for islands with errorThresh > maxErrors
                        break;
                }
                else
                    break;
            }

        if (errors < minErrors)
        {
            minErrors = errors;
            bestHitLength = r - ritBeg;
            verifier.m.beginPos = git - begin(host(inf), Standard());
        }
        else if (errors == minErrors && bestHitLength < r - ritBeg)
        {
            bestHitLength = r - ritBeg;
            verifier.m.beginPos = git - begin(host(inf), Standard());
        }
        else if (errorThresh < errors)
        {
            if (minErrors <= maxErrors)
            {
                verifier.m.endPos = verifier.m.beginPos + bestHitLength;
                verifier.q.pairScore = verifier.q.score = bestHitLength;
                verifier.q.errors = minErrors;

                if (maxErrors > minErrors + verifier.options->dRange)
                    maxErrors = minErrors + verifier.options->dRange;
                minErrors = maxErrors + 2;

                verifier.push();
            }
        }
    }

    if (minErrors <= maxErrors)
    {
        verifier.m.endPos = verifier.m.beginPos + bestHitLength;
        verifier.q.pairScore = verifier.q.score = bestHitLength;
        verifier.q.errors = minErrors;
        if (!verifier.oneMatchPerBucket)
        {
            if (maxErrors > minErrors + verifier.options->dRange)
                maxErrors = minErrors + verifier.options->dRange;
            verifier.push();
            // update maximal errors per read
            unsigned cutOff = maxErrors + 1;
            if ((unsigned)verifier.options->errorCutOff[readId] > cutOff)
                verifier.options->errorCutOff[readId] = cutOff;
        }
        return true;
    }
    return false;
}

template <typename TRazerSMode>
struct UseQualityValues__
{
    enum
    {
        VALUE = false
    };
};
template <typename TAlignMode, typename TGapMode, typename TSpec, typename TMatchNPolicy>
struct UseQualityValues__<RazerSMode<TAlignMode, TGapMode, RazerSQuality<TSpec>, TMatchNPolicy> >
{
    enum
    {
        VALUE = true
    };
};

//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
    typename TMatchVerifier,
    typename TGenome,
    typename TRead,
    typename TScoreMode,
    typename TMatchNPolicy>
inline bool
matchVerify(
    TMatchVerifier & verifier,
    Segment<TGenome, InfixSegment> inf,                             // potential match genome region
    unsigned readId,                                            // read number
    TRead const & read,                                             // reads
    RazerSMode<RazerSGlobal, RazerSUngapped, TScoreMode, TMatchNPolicy> const &) // Semi-global, no gaps
{
    typedef Segment<TGenome, InfixSegment>                          TGenomeInfix;
    typedef typename Iterator<TGenomeInfix, Standard>::Type         TGenomeIterator;
    typedef typename Iterator<TRead const, Standard>::Type          TReadIterator;
    typedef RazerSMode<RazerSGlobal, RazerSUngapped, TScoreMode, TMatchNPolicy> TRazerSMode;

#ifdef RAZERS_DEBUG
    std::cout << "Verify: " << std::endl;
    std::cout << "Genome: " << inf << "\t" << beginPosition(inf) << "," << endPosition(inf) << std::endl;
    std::cout << "Read:   " << read << std::endl;
#endif

    int mismatchDelta, scoreInit;
    int minScore;
    if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
    {
        minScore = verifier.options->errorCutOff[readId];
        if (minScore == 0)
            return false;

        minScore = -minScore + 1;
    }
    else if (UseQualityValues__<TRazerSMode>::VALUE)
        minScore = -verifier.options->absMaxQualSumErrors;
    else if (IsSameType<TScoreMode, RazerSScore>::VALUE)
    {
        minScore = verifier.options->minScore;
        mismatchDelta = scoreMatch(verifier.options->scoringScheme) - scoreMismatch(verifier.options->scoringScheme);
        scoreInit = scoreMatch(verifier.options->scoringScheme) * length(read);
    }

    // verify
    TReadIterator ritBeg    = begin(read, Standard());
    TReadIterator ritEnd    = end(read, Standard());
    unsigned ndlLength      = ritEnd - ritBeg;

    if (length(inf) < ndlLength)
        return false;

    TGenomeIterator git     = begin(inf, Standard());
    TGenomeIterator gitEnd  = end(inf, Standard()) - (ndlLength - 1);

    int maxScore = minScore - 1;
    int scoreThresh = (verifier.oneMatchPerBucket) ? std::numeric_limits<int>::max() : minScore;
    int score, errors;

    for (; git < gitEnd; ++git)
    {
        if (!IsSameType<TScoreMode, RazerSScore>::VALUE)
            score = 0;
        else
            score = scoreInit;

        if (!IsSameType<TScoreMode, RazerSErrors>::VALUE)
            errors = 0;

        TGenomeIterator g = git;
        for (TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
            if ((verifier.options->compMask[ordValue(*g)] & verifier.options->compMask[ordValue(*r)]) == 0)
            {
                if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
                {
                    // A. Count mismatches only
                    --score;
                }
                else
                {
                    ++errors;
                    if (UseQualityValues__<TRazerSMode>::VALUE)
                        // B. Count mismatches and mismatch qualities
                        score -= getQualityValue(*g);
                    else if (IsSameType<TScoreMode, RazerSScore>::VALUE)
                        // C. Count mismatches and alignment score
                        score -= mismatchDelta;
                    else
                        SEQAN_FAIL("Unsupported score mode!");
                }
                if (score < minScore)   // doesn't work for islands with errorThresh > maxErrors
                    break;
            }

        if (score > maxScore)
        {
            maxScore = score;
            if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
                verifier.m.score = score;
            else
                verifier.m.score = errors;
            verifier.m.beginPos = git - begin(host(inf), Standard());
        }
#ifdef RAZERS_ISLAND_CRITERION
        else if (scoreThresh > score)
        {
            if (maxScore >= minScore)
            {
                // for RazerSErrors bestErrors == -maxScore
                verifier.m.endPos = verifier.m.beginPos + ndlLength;
                verifier.m.pairScore = verifier.m.score = maxScore;
                if (!verifier.oneMatchPerBucket)
                {
                    if (minScore < maxScore - verifier.options->dRange)
                        minScore = maxScore - verifier.options->dRange;
                    verifier.push();
                }
                maxScore = minScore - 1;
            }
        }
#else
        (void)scoreThresh;
#endif
    }

    if (maxScore >= minScore)
    {
        verifier.m.endPos = verifier.m.beginPos + ndlLength;
        verifier.m.pairScore = verifier.m.score = maxScore;
        if (!verifier.oneMatchPerBucket)
        {
            if (minScore < maxScore - verifier.options->dRange)
                minScore = maxScore - verifier.options->dRange;
            verifier.push();

            // update maximal errors per read
            unsigned cutOff = -minScore + 1;
            if ((unsigned)verifier.options->errorCutOff[readId] > cutOff)
                verifier.options->errorCutOff[readId] = cutOff;
        }
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
    typename TMatchVerifier,
    typename TGenome,
    typename TRead,
    typename TMatchNPolicy>
inline bool
matchVerify(
    TMatchVerifier & verifier,
    Segment<TGenome, InfixSegment> inf,                                 // potential match genome region
    unsigned readId,                                                    // read number
    TRead const & read,                                                 // reads
    RazerSMode<RazerSGlobal, RazerSGapped, RazerSErrors, TMatchNPolicy> const &) // Mismatches and Indels
{
    if (empty(inf))
        return false;

    typedef Segment<TGenome, InfixSegment>                  TGenomeInfix;
    typedef typename Prefix<TRead const>::Type              TReadPrefix SEQAN_UNUSED_TYPEDEF;
    typedef typename Position<TGenomeInfix>::Type           TPosition;
    typedef typename MakeSigned_<TPosition>::Type           TDistance;

    // find read match end
    typedef Finder<TGenomeInfix>                            TMyersFinder;
    typedef typename TMatchVerifier::TOptions::TMyersPattern TMyersPattern SEQAN_UNUSED_TYPEDEF;
    typedef typename TMatchVerifier::TPatternState          TPatternState;

    // find read match begin
    // TODO(holtgrew): Use reverse-search here, as well!
    typedef ModifiedString<TGenomeInfix, ModReverse>        TGenomeInfixRev;
    typedef Finder<TGenomeInfixRev>                         TMyersFinderRev;

#ifdef RAZERS_NOOUTERREADGAPS
    typedef ModifiedString<TReadPrefix, ModReverse>         TReadRev;
#else
    typedef ModifiedString<TRead const, ModReverse>         TReadRev;
#endif
    typedef Pattern<TReadRev, MyersUkkonenGlobal>           TMyersPatternRev;

    unsigned ndlLength = length(read);
    int maxScore = std::numeric_limits<int>::min();
    int minScore = verifier.options->errorCutOff[readId];
    if (minScore == 0)
        return false;

    minScore = -minScore + 1;

    TDistance minSinkDistance = std::numeric_limits<TDistance>::max();
    TPosition maxPos = 0;
    TPosition lastPos = length(inf);
#ifdef RAZERS_ISLAND_CRITERION
    unsigned minDistance = (verifier.oneMatchPerBucket) ? lastPos : 1;
#endif

#ifdef RAZERS_NOOUTERREADGAPS
    TGenomeInfix origInf(inf);
    setEndPosition(inf, endPosition(inf) - 1);
    --ndlLength;
    TReadPrefix readPrefix = prefix(read, ndlLength);
#else
    TRead readPrefix(read);  // here only infixes (no sequence) is copied
#endif

    TMyersFinder myersFinder(inf);
#ifndef RAZERS_BANDED_MYERS
    TMyersPattern & myersPattern = verifier.options->forwardPatterns[readId];
#endif  // #ifdef RAZERS_BANDED_MYERS
    TPatternState & state = verifier.patternState;

#ifdef RAZERS_DEBUG
    std::cout << "Verify: " << std::endl;
    std::cout << "Genome: " << inf << "\t" << beginPosition(inf) << "," << endPosition(inf) << std::endl;
    std::cout << "Read:   " << read << "(id: " << readId << ")" << std::endl;
#endif

    // find end of best semi-global alignment
#ifdef RAZERS_BANDED_MYERS
    while (find(myersFinder, readPrefix, state, minScore))
#else  // #ifdef RAZERS_BANDED_MYERS
    while (find(myersFinder, myersPattern, state, minScore))
#endif  // #ifdef RAZERS_BANDED_MYERS
    {
        TPosition const pos = position(hostIterator(myersFinder));
        int score = getScore(state);

#ifdef RAZERS_NOOUTERREADGAPS
        // Manually align the last base of the read
        //
        // In this case myersPattern contains the whole read without the
        // last base. We compare the bases and adjust the score.
        // We also have to adjust inf and remove the last base of the
        // genomic region that has to be verified.
        SEQAN_ASSERT_LT(pos + 1, length(origInf));
        if ((verifier.options->compMask[ordValue(origInf[pos + 1])] & verifier.options->compMask[ordValue(back(read))]) == 0)
            if (--score < minScore)
                continue;
#endif
#ifdef RAZERS_ISLAND_CRITERION
        if (lastPos + minDistance < pos)
        {
            if (minScore <= maxScore)
            {
                verifier.m.endPos = beginPosition(inf) + maxPos + 1;
                verifier.m.pairScore = verifier.m.score = maxScore;
//              verifier.m.errors = -maxScore;

                if (maxScore == 0)
                    verifier.m.beginPos = verifier.m.endPos - ndlLength;
                else
                {
                    // find beginning of best semi-global alignment
                    TPosition infBeginPos = beginPosition(inf);
                    TPosition infEndPos = endPosition(inf);
                    TPosition newInfEndPos = verifier.m.endPos;

#ifdef RAZERS_BANDED_MYERS
                    verifier.revPatternState.leftClip = infEndPos - newInfEndPos + verifier.rightClip;
#endif
                    setEndPosition(inf, newInfEndPos);

//					// limit the beginning to needle length plus errors (== -maxScore)
//					if (length(inf) > ndlLength - maxScore)
//						setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);

                    // we eventually have to search before the beginning of our parallelogram
                    // otherwise alignments of an island in the previous parallelogram
                    // could be cut and prevent that an island in this parallelgram is found
                    if (endPosition(inf) > (unsigned)(ndlLength - maxScore))
                        setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
                    else
                        setBeginPosition(inf, 0);

#ifdef RAZERS_NOOUTERREADGAPS
                    // The best score must be corrected to hold the score of the prefix w/o the last read base
                    if ((verifier.options->compMask[ordValue(origInf[maxPos + 1])] & verifier.options->compMask[ordValue(back(read))]) == 0)
                        ++maxScore;
#endif

                    TReadRev            readRev(readPrefix);
                    TGenomeInfixRev     infRev(inf);
                    TMyersPatternRev    myersPatternRev(readRev);
                    TMyersFinderRev     myersFinderRev(infRev);

                    verifier.m.beginPos = verifier.m.endPos;
#ifdef RAZERS_BANDED_MYERS
                    while (find(myersFinderRev, readRev, verifier.revPatternState, maxScore))
                    {
                        if (maxScore <= getScore(verifier.revPatternState))
                        {
                            maxScore = getScore(verifier.revPatternState);
                            verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);
                        }
                    }
#else
                    _patternMatchNOfPattern(myersPatternRev, verifier.options->matchN);
                    _patternMatchNOfFinder(myersPatternRev, verifier.options->matchN);
                    while (find(myersFinderRev, myersPatternRev, maxScore))
                    {
                        if (maxScore <= getScore(myersPatternRev))
                        {
                            maxScore = getScore(myersPatternRev);
                            verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);
                        }
                    }
#endif
                    setBeginPosition(inf, infBeginPos);
                    setEndPosition(inf, infEndPos);
#ifdef RAZERS_BANDED_MYERS
                    if (verifier.m.beginPos == verifier.m.endPos)
                        continue;
#endif
                }
                // minDistance implicitly forbids to get here with verifier.oneMatchPerBucket == true
                SEQAN_ASSERT_NOT(verifier.oneMatchPerBucket);
                SEQAN_ASSERT_LT(verifier.m.beginPos, verifier.m.endPos);

#ifdef RAZERS_NOOUTERREADGAPS
                // The match end position must be increased by the omitted base.
                ++verifier.m.endPos;
#endif
                if (minScore < verifier.m.score - verifier.options->dRange)
                    minScore = verifier.m.score - verifier.options->dRange;

                verifier.push();
                maxScore = minScore - 1;
                minSinkDistance = std::numeric_limits<TDistance>::max();
            }
        }
#endif  // #ifdef RAZERS_ISLAND_CRITERION

        // minimize distance between sink and estimated match begin
        TDistance sinkDistance = verifier.sinkPos - ((TDistance)(beginPosition(inf) + pos) - (TDistance)ndlLength);
//        sinkDistance = maxValue(sinkDistance);
        if (sinkDistance < (TDistance)0)
            sinkDistance = -sinkDistance;

        if (score > maxScore || (score == maxScore && sinkDistance <= minSinkDistance))
        {
            maxScore = score;
            minSinkDistance = sinkDistance;
            maxPos = pos;
        }
        lastPos = pos;
    }

    if (minScore <= maxScore)
    {
        verifier.m.endPos = beginPosition(inf) + maxPos + 1;
        verifier.m.pairScore = verifier.m.score = maxScore;
//		verifier.m.errors = -maxScore;

        if (maxScore == 0)
            verifier.m.beginPos = verifier.m.endPos - ndlLength;
        else
        {
            // find beginning of best semi-global alignment
            TPosition newInfEndPos = verifier.m.endPos;

#ifdef RAZERS_BANDED_MYERS
            verifier.revPatternState.leftClip = endPosition(inf) - newInfEndPos + verifier.rightClip;
#endif
            setEndPosition(inf, newInfEndPos);

//					// limit the beginning to needle length plus errors (== -maxScore)
//					if (length(inf) > ndlLength - maxScore)
//						setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);

            // we eventually have to search before the beginning of our parallelogram
            // otherwise alignments of an island in the previous parallelogram
            // could be cut and prevent that an island in this parallelgram is found
            if (endPosition(inf) > (unsigned)(ndlLength - maxScore))
                setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
            else
                setBeginPosition(inf, 0);

#ifdef RAZERS_NOOUTERREADGAPS
            // The best score must be corrected to hold the score of the prefix w/o the last read base
            if ((verifier.options->compMask[ordValue(origInf[maxPos + 1])] & verifier.options->compMask[ordValue(back(read))]) == 0)
                ++maxScore;
#endif

            TReadRev            readRev(readPrefix);
            TGenomeInfixRev     infRev(inf);
            TMyersPatternRev    myersPatternRev(readRev);
            TMyersFinderRev     myersFinderRev(infRev);

            verifier.m.beginPos = verifier.m.endPos;
#ifdef RAZERS_BANDED_MYERS
            while (find(myersFinderRev, readRev, verifier.revPatternState, maxScore))
            {
                if (maxScore <= getScore(verifier.revPatternState))
                {
                    maxScore = getScore(verifier.revPatternState);
                    verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);
                }
            }
#else
            _patternMatchNOfPattern(myersPatternRev, verifier.options->matchN);
            _patternMatchNOfFinder(myersPatternRev, verifier.options->matchN);
            while (find(myersFinderRev, myersPatternRev, maxScore))
            {
                if (maxScore <= getScore(myersPatternRev))
                {
                    maxScore = getScore(myersPatternRev);
                    verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);
                }
            }
#endif
#ifdef RAZERS_BANDED_MYERS
            if (verifier.m.beginPos == verifier.m.endPos)
            {
#ifdef RAZERS_DEBUG
                std::cout << "FAILED2" << std::endl;
#endif
                return false;
            }
#endif // RAZERS_BANDED_MYERS
        }
        SEQAN_ASSERT_LT(verifier.m.beginPos, verifier.m.endPos);

#ifdef RAZERS_NOOUTERREADGAPS
        // The match end position must be increased by the omitted base.
        ++verifier.m.endPos;
#endif

        if (!verifier.oneMatchPerBucket)
        {
            if (minScore < verifier.m.score - verifier.options->dRange)
                minScore = verifier.m.score - verifier.options->dRange;
            verifier.push();

            // update maximal errors per read
            unsigned cutOff = -minScore + 1;
            if ((unsigned)verifier.options->errorCutOff[readId] > cutOff)
            {
//    std::cout<<"maxScore="<<verifier.m.score  << " minScore=" << minScore<<std::endl;
//    std::cout<<"decreased errCutOff["<<readId<<"] from "<<(unsigned)verifier.options->errorCutOff[readId];
                verifier.options->errorCutOff[readId] = cutOff;
//    std::cout<<"to "<<(unsigned)verifier.options->errorCutOff[readId]<<std::endl;
            }
        }

#ifdef RAZERS_DEBUG
        std::cout << "OK" << std::endl;
#endif
        return true;
    }
#ifdef RAZERS_DEBUG
    std::cout << "FAILED3" << std::endl;
#endif
    return false;
}

template <
    typename TMatchVerifier,
    typename TGenome,
    typename TRead,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TMatchNPolicy>
inline bool
matchVerify(
    TMatchVerifier &,
    Segment<TGenome, InfixSegment>,                             // potential match genome region
    unsigned,                                                   // read number
    TRead const &,                                              // read
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const &)
{
    SEQAN_FAIL("Verification not implemented!");
    return false;
}

template <typename TOptions>
inline void
estimateErrorDistributionFromQualities(TOptions & options)
{
    typedef typename TOptions::TProb TFloat;

    resize(options.errorProb, length(options.avrgQuality));
//	std::cout<< "ERROR PROBS:"<<std::endl;
    for (unsigned i = 0; i < length(options.avrgQuality); ++i)
    {
        //    qual = -10 log_10 p
        //         = -10 log p / log 10
        //       p = exp^(qual * log 10 / -10)
        //   log p = qual * log 10 / -10;
        double e = options.avrgQuality[i] * log(10.0) / -10.0;
        TFloat sequencingError;
        sequencingError.data_value = e;
//		sequencingError = exp(e);
        options.errorProb[i] = (TFloat)1.0 - ((TFloat)1.0 - sequencingError) * ((TFloat)1.0 - options.mutationRate);
//		std::cout<<e<<':'<<options.errorProb[i]<<'\t';
    }
//	std::cout<<std::endl<<std::endl;
}

//////////////////////////////////////////////////////////////////////////////
// Customize filter
template <typename TDelta, typename TOptions>
void computeQGramLengths(TDelta & minDelta, TOptions const & options)
{
    const unsigned maxLength1 = length(options.readLengths);
//    const unsigned maxLength = maxLength1 - 1;
//    const unsigned maxErrors = (unsigned) floor(options.errorRate * maxLength);
//    const unsigned maxErrors1 = maxErrors + 1;

    unsigned seqCount = 0;
    String<unsigned> maxDelta;
    resize(minDelta, options.maxOverlap + 1, std::numeric_limits<unsigned>::max());
    resize(maxDelta, options.maxOverlap + 1, 3);

    // compute delta (==stepSize) for different overlaps
    for (unsigned len = 0; len < maxLength1; ++len)
    {
        if (options.readLengths[len] == 0)
            continue;

        seqCount += options.readLengths[len];
        for (unsigned ol = 0; ol <= options.maxOverlap; ++ol)
        {
            // sequence must have sufficient length
            if (len <= ol)
                continue;

            // cut overlap many characters from the end
            unsigned errors = (unsigned) floor(options.errorRate * len);
            unsigned delta = (len - ol) / (errors + 1);

            // ignore too short delta-grams
            if (delta < 3)
                continue;
            if (minDelta[ol] > delta)
                minDelta[ol] = delta;
            if (maxDelta[ol] < delta)
                maxDelta[ol] = delta;
        }
    }

//	std::cout<< "Deltas:"<<std::endl;
//    for (unsigned ol = 0; ol < length(minDelta); ++ol)
//	{
//        if (minDelta[ol] < 3) minDelta[ol] = maxDelta[ol];
//		std::cout<<minDelta[ol]+ol<<'\t';
//	}
//	std::cout<<std::endl<<std::endl;
}

template <typename TEstLosses, typename TDelta, typename TOptions>
unsigned estimatePigeonholeLosses(TEstLosses & estLosses, TDelta const & delta, TOptions const & options)
{
    typedef typename Value<TEstLosses>::Type TFloat;

    const unsigned maxLength1 = length(options.readLengths);
    const unsigned maxLength = maxLength1 - 1;
    const unsigned maxErrors = (unsigned) floor(options.errorRate * maxLength);
    const unsigned maxErrors1 = maxErrors + 1;

    // ------------------------------------------------------------------------------------------------
    // compute probs to have 0,1,...,maxErrors errors in the prefix of length 1,...,maxLength
    // ------------------------------------------------------------------------------------------------
    String<TFloat> p_prefix;
    resize(p_prefix, maxLength * maxErrors1);
    p_prefix[0] = (TFloat)1.0 - options.errorProb[0];   // 0 error
    p_prefix[1] = options.errorProb[0];                 // 1 error
    for (unsigned e = 2; e <= maxErrors; ++e)           // 2 errors are not possible in a sequence of length 1
        p_prefix[e] = 0.0;
    for (unsigned i = 1, idx = maxErrors1; i < maxLength; ++i)
    {
        p_prefix[idx] = p_prefix[idx - maxErrors1] * ((TFloat)1.0 - options.errorProb[i]);
        ++idx;
        for (unsigned e = 1; e <= maxErrors; ++e, ++idx)
            p_prefix[idx] = p_prefix[idx - maxErrors1] * ((TFloat)1.0 - options.errorProb[i]) + p_prefix[idx - maxErrors1 - 1] * options.errorProb[i];
    }

    // ------------------------------------------------------------------------------------------------
    // compute expected numbers of reads with 0,1,...,maxErrors errors
    // ------------------------------------------------------------------------------------------------

    clear(estLosses);
    resize(estLosses, maxErrors1 + 2, 0.0);
    double maxLoss = 0.0;
//		std::cout << "len:"<<0<<"  "<<options.readLengths[0]<<std::endl;
    for (unsigned len = 1; len <= maxLength; ++len)
    {
        if (options.readLengths[len] == 0)
            continue;
//		std::cout << "len:"<<len<<"  "<<options.readLengths[len]<<'\t'<<(TFloat)options.readLengths[len]<<'\t'<<(double)(TFloat)options.readLengths[len]<<'\t'<<std::endl;

        unsigned errors = (unsigned) floor(options.errorRate * len);
        TFloat sum = 0.0;
        for (unsigned e = 0; e <= errors; ++e)
        {
            sum += p_prefix[(len - 1) * maxErrors1 + e];
            estLosses[2 + e] += p_prefix[(len - 1) * maxErrors1 + e] * (TFloat)options.readLengths[len];
        }
        maxLoss += (double)sum * options.readLengths[len];
    }
    maxLoss *= options.lossRate;

    unsigned overlap = 0;

    // ------------------------------------------------------------------------------------------------
    // compute loss for each overlap and select best the shape
    // ------------------------------------------------------------------------------------------------
    for (unsigned ol = 0; ol < length(delta); ++ol)
    {
        unsigned stepSize = delta[ol];
        unsigned q = stepSize + ol;

        // if we had the same q-gram before
        if (ol > 0 && delta[ol - 1] + ol - 1 == q)
            continue;

        // don't allow more than two overlapping q-grams
        if (2 * ol > stepSize)
            break;

        // ------------------------------------------------------------------------------------------------
        // compute probs to have 0,1,...,maxErrors errors in a segment
        // ------------------------------------------------------------------------------------------------

        String<TFloat> p;
        resize(p, maxLength * maxErrors1);
        for (unsigned i = 0, idx = 0; i < maxLength; ++i)
        {
            if (i % stepSize == 0 || (i >= stepSize && i % stepSize == ol))
            {
                // initialization
                p[idx++] = (TFloat)1.0 - options.errorProb[i];      // 0 error
                p[idx++] = options.errorProb[i];                    // 1 error
                for (unsigned e = 2; e <= maxErrors; ++e, ++idx)    // 2 errors are not possible in a sequence of length 1
                    p[idx] = 0.0;
            }
            else
            {
                p[idx] = p[idx - maxErrors1] * ((TFloat)1.0 - options.errorProb[i]);
                ++idx;
                for (unsigned e = 1; e <= maxErrors; ++e, ++idx)
                {
                    // a sequence of length i with e errors has either:
                    //   - e errors in prefix i-1 and no error at position i
                    //   - e-1 errors in prefix i-1 and an error at position i
                    p[idx] = p[idx - maxErrors1] * ((TFloat)1.0 - options.errorProb[i]) + p[idx - maxErrors1 - 1] * options.errorProb[i];
                }
            }
        }

        String<TFloat> p_last;
        resize(p_last, maxLength * maxErrors1);
        for (unsigned i = 0, idx = 0; i < maxLength; ++i)
        {
            if (i == 0 || (i >= stepSize && i % stepSize == ol))
            {
                // initialization
                p_last[idx++] = (TFloat)1.0 - options.errorProb[i]; // 0 error
                p_last[idx++] = options.errorProb[i];               // 1 error
                for (unsigned e = 2; e <= maxErrors; ++e, ++idx)    // 2 errors are not possible in a sequence of length 1
                    p_last[idx] = 0.0;
            }
            else
            {
                p_last[idx] = p_last[idx - maxErrors1] * ((TFloat)1.0 - options.errorProb[i]);
                ++idx;
                for (unsigned e = 1; e <= maxErrors; ++e, ++idx)
                    p_last[idx] = p_last[idx - maxErrors1] * ((TFloat)1.0 - options.errorProb[i]) + p_last[idx - maxErrors1 - 1] * options.errorProb[i];
            }
        }

        // ------------------------------------------------------------------------------------------------
        // compute probs to lose every prefix of pigeonhole segments (q-grams) with 0,...,maxErrors errors
        // ------------------------------------------------------------------------------------------------
        unsigned maxXErrors1 = _min(maxErrors, ol) + 1;
        unsigned maxSegments = (maxLength - ol) / stepSize;
        String<TFloat> M;
        resize(M, maxSegments * maxErrors1 * maxXErrors1, 0.0);

        unsigned idx = maxXErrors1;
        for (unsigned e0 = 1; e0 <= maxErrors; ++e0, idx += maxXErrors1)
        {
            if (q != stepSize)
            {
                unsigned maxX = _min(e0, ol);
                for (unsigned x0 = 0; x0 <= maxX; ++x0)
                {
                    M[idx + x0] = p[(stepSize - 1) * maxErrors1 + (e0 - x0)] * p[(q - 1) * maxErrors1 + x0];
//                  std::cout<<"M[0,"<<e0<<','<<x0<<"]="<<M[idx + x0]<<std::endl;
                }
            }
            else
                M[idx] = p[(stepSize - 1) * maxErrors1 + e0];
//		std::cout<<"p(C_"<<0<<'='<<e0<<")="<<p[( (stepSize - 1)) * maxErrors1 + e0]<<std::endl;
//		std::cout<<"p(X_"<<0<<'='<<e0<<")="<<p[( (q - 1)) * maxErrors1 + e0]<<std::endl;
        }
/*
    SEQAN_OMP_PRAGMA(critical)
    {
        for (unsigned e0 = 1; e0 <= maxErrors; ++e0)
            for (unsigned i = 0; i < maxLength; ++i)
                std::cerr<<"p("<<i<<'='<<e0<<")="<<p[i * maxErrors1 + e0]<<std::endl;
        for (unsigned e0 = 1; e0 <= maxErrors; ++e0)
            for (unsigned i = 0; i < maxLength; ++i)
                std::cerr<<"p_last("<<i<<'='<<e0<<")="<<p_last[i * maxErrors1 + e0]<<std::endl;
    }
*/
//		std::cout<<"=============================="<<std::endl;

        // M[s,e,x] at position ((s * maxErrors1) + e) * maxErrors1 + x
        for (unsigned s = 1; s < maxSegments; ++s)
            for (unsigned e = 0; e <= maxErrors; ++e, idx += maxXErrors1)
            {
                if (ol != 0)
                {
                    unsigned maxX = _min(e, ol);
                    for (unsigned x = 0; x <= maxX; ++x)
                    {
                        TFloat sum = 0.0;
                        for (unsigned _e = 0; _e <= e - x; ++_e)
                        {
                            unsigned max_X = _min(_e, ol);
                            for (unsigned _x = 0; _x <= max_X; ++_x)
                                if (_e - _x < e)    // at least one error in seed X_{s-1} C_s X_s
                                    sum +=   M[((s - 1) * maxErrors1 + _e) * maxXErrors1 + _x]
                                           * p[(stepSize * s + (stepSize - 1)) * maxErrors1 + (e - _e - x)]
                                           * p[(stepSize * s + (q - 1)) * maxErrors1 + x];
                        }

                        M[idx + x] = sum;
//                        std::cerr<<"M["<<s<<','<<e<<','<<x<<"]="<<M[idx + x]<<std::endl;
                    }
//                  std::cerr<<"p(C_"<<s<<'='<<e<<")="<<p[(stepSize * s + (stepSize - 1)) * maxErrors1 + e]<<std::endl;
//                  std::cerr<<"p(X_"<<s<<'='<<e<<")="<<p[(stepSize * s + (q - 1)) * maxErrors1 + e]<<std::endl;
                }
                else
                {
                    TFloat sum = 0.0;
                    for (unsigned _e = 0; _e < e; ++_e)
                        sum +=   M[((s - 1) * maxErrors1 + _e) * maxXErrors1]
                               * p[(stepSize * s + (stepSize - 1)) * maxErrors1 + (e - _e)];

                    M[idx] = sum;
//                    std::cerr<<"M["<<s<<','<<e<<",0]="<<M[idx]<<std::endl;
                }
            }

//		// no overlap => no loss => no estimation required
//		if (ol == 0) continue;

        // ------------------------------------------------------------------------------------------------
        // sum up the expected numbers of lost reads for every read length 0,...,maxLength
        // ------------------------------------------------------------------------------------------------
        double loss = 0.0;
//		std::cout<< "LOSSES FOR OVERLAP " << ol <<":"<<std::endl;

        appendValue(estLosses, (double)ol);
        appendValue(estLosses, (double)q);
        resize(estLosses, length(estLosses) + maxErrors1, 0.0);
        for (unsigned len = 1; len <= maxLength; ++len)
        {
            if (options.readLengths[len] == 0)
                continue;
            if (len < q)
                continue;

            unsigned errors = (unsigned) floor(options.errorRate * len);
            unsigned segments = (len - ol) / stepSize;

//			TFloat divider = p_prefix[(len - 1) * maxErrors1];
//			std::cout<<"DIVIDER"<<0<< ':'<<p_prefix[(len - 1) * maxErrors1]<<std::endl;
//			for (unsigned e = 1; e <= errors; ++e)
//			{
//				divider += p_prefix[(len - 1) * maxErrors1 + e];
//				std::cout<<"DIVIDER"<<e<< ':'<<p_prefix[(len - 1) * maxErrors1 + e]<<std::endl;
//			}
//			std::cout<<"DIVIDERSUM:"<<divider<<std::endl;

            TFloat lossPerLength = 0.0;
            unsigned segmentedLen = segments * stepSize + ol;
            for (unsigned e1 = 0; e1 <= errors; ++e1)
            {
                unsigned maxX = _min(e1, ol);
                for (unsigned x = 0; x <= maxX; ++x)
                {
                    TFloat sum = 0.0;
                    if (segmentedLen < len)
                    {
                        unsigned maxE2 = _min(errors - e1, len - segmentedLen);
                        for (unsigned e2 = 0; e2 <= maxE2; ++e2)
                        {
                            sum += p_last[(len - 1) * maxErrors1 + e2];
                            estLosses[length(estLosses) - maxErrors1 + e1 + e2] += (M[((segments - 1) * maxErrors1 + e1) * maxXErrors1 + x] * p_last[(len - 1) * maxErrors1 + e2] /* / divider */) * options.readLengths[len];
                        }
                    }
                    else
                    {
                        sum = 1.0;
                        estLosses[length(estLosses) - maxErrors1 + e1] += (M[((segments - 1) * maxErrors1 + e1) * maxXErrors1 + x] /* / divider */) * options.readLengths[len];
                    }
                    lossPerLength += M[((segments - 1) * maxErrors1 + e1) * maxXErrors1 + x] * sum;
                }
            }

//			lossPerLength /= divider;
            loss += options.readLengths[len] * (double)lossPerLength;
//			std::cout<<len<<':'<<options.readLengths[len] * (double)lossPerLength<<'\t';

//			if (len >= delta[ol] * (errors + 2) + ol) continue;
//			loss += options.readLengths[len] * (q0[len * options.maxOverlap + ol - 1] / p[len * maxErrors1 + errors]);
        }
        if (loss > maxLoss)
            break;
//		std::cout<<std::endl<<std::endl;

        overlap = ol;
    }
    return overlap;
}

template <typename TIndex, typename TPigeonholeSpec, typename TOptions>
void _applyFilterOptions(Pattern<TIndex, Pigeonhole<TPigeonholeSpec> > & filterPattern, TOptions const & options)
{
    if (options.lossRate == 0.0)
    {
        filterPattern.params.delta = options.delta;
        filterPattern.params.overlap = options.overlap;
        _patternInit(filterPattern, options.errorRate);

        if (options._debugLevel >= 2)
        {
            SEQAN_OMP_PRAGMA(critical)
            {
                CharString str;
                shapeToString(str, filterPattern.shape);
                std::cout << std::endl << "Pigeonhole settings:" << std::endl;
                std::cout << "  shape:    " << length(str) << '\t' << str << std::endl;
                std::cout << "  stepsize: " << getStepSize(host(filterPattern)) << std::endl;
            }
            return;
        }
    }

    typedef typename TOptions::TProb TFloat;

    String<TFloat> estLosses;
    const unsigned maxErrors = (unsigned) floor(options.errorRate * length(options.readLengths));
    const unsigned maxErrors1 = maxErrors + 1;

    String<unsigned> delta;
    computeQGramLengths(delta, options);
    unsigned maxWeight = _pigeonholeMaxShapeWeight(indexShape(host(filterPattern)));

    if (delta[0] < maxWeight && options.lossRate != 0.0)
    {
        // lossy filtration

        // we can stop if the proposed weight exceeds or reaches the maximum
        for (unsigned ol = 0; ol < length(delta); ++ol)
            if (delta[ol] + ol > maxWeight)
            {
                resize(delta, ol);
                break;
            }

        filterPattern.params.overlap = estimatePigeonholeLosses(estLosses, delta, options);
        _patternInit(filterPattern, options.errorRate);

        if (options._debugLevel >= 2)
        {
            SEQAN_OMP_PRAGMA(critical)
            {
                std::cout << "     e | e error reads | loss ol =" << std::setw(2) << estLosses[maxErrors1 + 2];
                for (unsigned i = 2 * (maxErrors1 + 2); i < length(estLosses); i += maxErrors1 + 2)
                    std::cout << " |      ol =" << std::setw(2) << estLosses[i];
                std::cout << std::endl;
                std::cout << "       |               |       q =" << std::setw(2) << estLosses[maxErrors1 + 3];
                for (unsigned i = 2 * (maxErrors1 + 2); i < length(estLosses); i += maxErrors1 + 2)
                    std::cout << " |       q =" << std::setw(2) << estLosses[i + 1];
                std::cout << std::endl;
                std::cout << " ------+---------------";
                for (unsigned i = maxErrors1 + 2; i < length(estLosses); i += maxErrors1 + 2)
                    std::cout << "+-------------";
                std::cout << std::endl;
                std::cout.setf(std::ios::fixed);
                TFloat sumReads;
                for (unsigned e = 0; e <= maxErrors; ++e)
                {
                    std::cout.fill(' ');
                    std::cout << std::setprecision(0);
                    std::cout << std::setw(6) << e;
                    std::cout << " |";
                    std::cout << std::setw(14) << (unsigned)estLosses[e + 2];
                    std::cout << std::setprecision(2);
                    for (unsigned i = maxErrors1 + 2 + e + 2; i < length(estLosses); i += maxErrors1 + 2)
                        std::cout << " |" << std::setw(12) << estLosses[i];
                    std::cout << std::endl;
                    sumReads += estLosses[e + 2];
                }
                std::cout << " ------+---------------";
                for (unsigned i = maxErrors1 + 2; i < length(estLosses); i += maxErrors1 + 2)
                    std::cout << "+-------------";
                std::cout << std::endl;

                std::cout << std::setprecision(0);
                std::cout << " total |" << std::setw(14) << (unsigned)sumReads << ' ';
                std::cout << std::setprecision(2);
                for (unsigned i = maxErrors1 + 4; i < length(estLosses); i += maxErrors1 + 2)
                {
                    sumReads = 0.0;
                    for (unsigned e = 0; e <= maxErrors; ++e)
                        sumReads += estLosses[i + e];
                    std::cout << '|' << std::setw(12) << sumReads;
                    std::cout << ((filterPattern.params.overlap == (unsigned)estLosses[i - 2]) ? '*' : ' ');
                }
                std::cout << std::endl;

                CharString str;
                shapeToString(str, filterPattern.shape);
                std::cout << std::endl << "Pigeonhole settings:" << std::endl;
                std::cout << "  shape:    " << length(str) << '\t' << str << std::endl;
                std::cout << "  stepsize: " << getStepSize(host(filterPattern)) << std::endl;
            }
        }
    }
    else
    {
        // lossless filtration
        filterPattern.params.overlap = 0;
    }

}

template <typename TIndex, typename TSwiftSpec, typename TOptions>
void _applyFilterOptions(Pattern<TIndex, Swift<TSwiftSpec> > & filterPattern, TOptions const & options)
{
    filterPattern.params.minThreshold = options.threshold;
    filterPattern.params.tabooLength = options.tabooLength;
    _patternInit(filterPattern, options.errorRate, 0);
}

//////////////////////////////////////////////////////////////////////////////
// Find read matches in a single genome sequence
template <
    typename TMatches,
    typename TFragmentStore,
    typename TReadIndex,
    typename TFilterSpec,
    typename TCounts,
    typename TRazerSOptions,
    typename TRazerSMode>
void _mapSingleReadsToContig(
    TMatches & matches,
    TFragmentStore & store,
    unsigned                                  contigId,             // ... and its sequence number
    Pattern<TReadIndex, TFilterSpec> & filterPattern,
    TCounts & cnts,
    char                                      orientation,              // q-gram index of reads
    TRazerSOptions & options,
    TRazerSMode                       const & mode)
{
    // FILTRATION
    typedef typename TFragmentStore::TContigSeq             TContigSeq;
    typedef Finder<TContigSeq, TFilterSpec>                 TFilterFinder;
    typedef Pattern<TReadIndex, TFilterSpec>                TFilterPattern;

    // VERIFICATION
    typedef MatchVerifier<
        TFragmentStore,
        TMatches,
        TRazerSOptions,
        TRazerSMode,
        TFilterPattern,
        TCounts>                                            TVerifier;
    typedef typename Fibre<TReadIndex, FibreText>::Type     TReadSet;

    // iterate all genomic sequences
    if (options._debugLevel >= 1)
    {
        std::cerr << std::endl << "Process genome seq #" << contigId;
        if (orientation == 'F')
            std::cerr << "[fwd]";
        else
            std::cerr << "[rev]";
    }
    lockContig(store, contigId);
    TContigSeq & contigSeq = store.contigStore[contigId].seq;
    if (orientation == 'R')
        reverseComplement(contigSeq);

    TReadSet & readSet = indexText(host(filterPattern));
    TFilterFinder   filterFinder(contigSeq, options.repeatLength, 1);
    TVerifier       verifier(matches, options, filterPattern, cnts);

#ifdef RAZERS_BANDED_MYERS
    typename Size<TContigSeq>::Type contigLength = length(contigSeq);
#endif

    // initialize verifier
    verifier.onReverseComplement = (orientation == 'R');
    verifier.genomeLength = length(contigSeq);
    verifier.m.contigId = contigId;

    double beginTime = sysTime();
    // Build q-gram index separately, so we can better compute the time for it.
    indexRequire(host(filterPattern), QGramSADir());
    options.timeBuildQGramIndex += sysTime() - beginTime;

    // iterate all verification regions returned by SWIFT
    while (find(filterFinder, filterPattern, options.errorRate))
    {
        // std::cout << "read id = " << (*filterFinder.curHit).ndlSeqNo << ", " << beginPosition(filterFinder) << std::endl;

//        if (length(infix(filterFinder)) < length(readSet[(*filterFinder.curHit).ndlSeqNo]))
//            continue;  // Skip if hit length < read length.  TODO(holtgrew): David has to fix something in banded myers to make this work.
#ifdef RAZERS_BANDED_MYERS
        verifier.patternState.leftClip = (beginPosition(filterFinder) >= 0) ? 0 : -beginPosition(filterFinder);           // left clip if match begins left of the genome
        verifier.rightClip = (endPosition(filterFinder) <= contigLength) ? 0 : endPosition(filterFinder) - contigLength;  // right clip if match end right of the genome
#endif
        verifier.m.readId = (*filterFinder.curHit).ndlSeqNo;
        if (!options.spec.DONT_VERIFY)
            matchVerify(verifier, infix(filterFinder), verifier.m.readId, readSet[verifier.m.readId], mode);
        ++options.countFiltration;
    }
    if (!unlockAndFreeContig(store, contigId))                          // if the contig is still used
        if (orientation == 'R')
            reverseComplement(contigSeq);
    // we have to restore original orientation
}

//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TCounts,
    typename TSpec,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TReadIndex,
    typename TMatchNPolicy,
    typename TFilterSpec>
int _mapSingleReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy>  const & mode,
    TReadIndex & readIndex,
    TFilterSpec)
{
    typedef FragmentStore<TFSSpec, TFSConfig>           TFragmentStore;
    //typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;
    typedef Pattern<TReadIndex, TFilterSpec>            TFilterPattern; // filter

    //typedef typename Value<TReadSeqStore>::Type const   TRead;
    //typedef Pattern<TRead, MyersUkkonen>                TMyersPattern;  // verifier
    // typedef Pattern<TRead, Myers<FindInfix, False, void> >	TMyersPattern;	// verifier

    typedef typename TFragmentStore::TContigSeq TContigSeq;
    typedef typename Position<TContigSeq>::Type TContigPos;
    typedef MatchRecord<TContigPos> TMatchRecord;

    // configure Swift pattern
    TFilterPattern filterPattern(readIndex);

    // configure filter pattern
    // (if this is a pigeonhole filter, all sequences must be appended first)
    _applyFilterOptions(filterPattern, options);
    filterPattern.params.printDots = options._debugLevel > 0;

    // clear stats
    options.countFiltration = 0;
    options.countVerification = 0;
    options.timeMapReads = 0;
    options.timeDumpResults = 0;
    options.timeBuildQGramIndex = 0;
    options.timeCompactMatches = 0;
    options.timeMaskDuplicates = 0;
    options.timeFsCopy = 0;
    options.timeFiltration = 0;
    options.timeVerification = 0;
    SEQAN_PROTIMESTART(find_time);

    // We collect the matches in a more compact data structure than the
    // AlignedReadStoreElement from FragmentStore.
    String<TMatchRecord> matches;

    // iterate over genome sequences
    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId)
    {
        // lock to prevent releasing and loading the same contig twice
        // (once per _mapSingleReadsToContig call)
        lockContig(store, contigId);
#ifndef RAZERS_WINDOW
        if (options.forward)
            _mapSingleReadsToContig(matches, store, contigId, filterPattern, cnts, 'F', options, mode);
        if (options.reverse)
            _mapSingleReadsToContig(matches, store, contigId, filterPattern, cnts, 'R', options, mode);
#else
        if (options.forward)
            _mapSingleReadsToContigWindow(store, contigId, filterPattern, cnts, 'F', options, mode);
        if (options.reverse)
            _mapSingleReadsToContigWindow(store, contigId, filterPattern, cnts, 'R', options, mode);
#endif
        unlockAndFreeContig(store, contigId);
    }

    options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
    if (options._debugLevel >= 1)
        std::cerr << std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << std::endl;

    double beginCopyTime = sysTime();
    // Final mask duplicates and compact matches.
    Nothing nothing;
    if (IsSameType<TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
        maskDuplicates(matches, options, mode);
    compactMatches(matches, cnts, options, mode, nothing, COMPACT_FINAL);
    // Write back to store.
    reserve(store.alignedReadStore, length(matches), Exact());
    reserve(store.alignQualityStore, length(matches), Exact());
    typedef typename Iterator<String<TMatchRecord>, Standard>::Type TIterator;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
    typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
    for (TIterator it = begin(matches), itEnd = end(matches); it != itEnd; ++it)
    {
        SEQAN_ASSERT_NEQ(it->orientation, '-');
        SEQAN_ASSERT_LEQ(it->beginPos, it->endPos);
        if (it->orientation == 'R')
            std::swap(it->beginPos, it->endPos);
        appendValue(store.alignedReadStore, TAlignedReadStoreElem(length(store.alignQualityStore), it->readId, it->contigId, it->beginPos, it->endPos));
        appendValue(store.alignQualityStore, TAlignedQualStoreElem(it->pairScore, it->score, -it->score));
    }
    options.timeFsCopy = sysTime() - beginCopyTime;

    if (options._debugLevel >= 1)
    {
        std::cerr << "Masking duplicates took          \t" << options.timeMaskDuplicates << " seconds" << std::endl;
        std::cerr << "Compacting matches took          \t" << options.timeCompactMatches << " seconds" << std::endl;
        std::cerr << "Building q-gram index took       \t" << options.timeBuildQGramIndex << " seconds" << std::endl;
        std::cerr << "Time for copying back            \t" << options.timeFsCopy << " seconds" << std::endl;
        std::cerr << "Time for filtration              \t" << options.timeFiltration << " seconds" << std::endl;
        std::cerr << "Time for verifications           \t" << options.timeVerification << " seconds" << std::endl;
        std::cerr << std::endl;
        std::cerr << "___FILTRATION_STATS____" << std::endl;
        std::cerr << "Filtration counter:      " << options.countFiltration << std::endl;
        std::cerr << "Successful verifications: " << options.countVerification << std::endl;
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for SWIFT (default)
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TCounts,
    typename TSpec,
    typename TShape,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TMatchNPolicy,
    typename TFilterSpec>
int _mapSingleReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & mode,
    TFilterSpec)
{
    typedef FragmentStore<TFSSpec, TFSConfig>           TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;
#ifndef RAZERS_OPENADDRESSING
    typedef Index<TReadSeqStore, IndexQGram<TShape> >   TIndex;         // q-gram index
#else
    typedef Index<TReadSeqStore, IndexQGram<TShape, OpenAddressing> >   TIndex;
#endif

    // configure q-gram index
    TIndex swiftIndex(store.readSeqStore, shape);
#ifdef RAZERS_OPENADDRESSING
    swiftIndex.alpha = options.loadFactor;
#endif
    cargo(swiftIndex).abundanceCut = options.abundanceCut;
    cargo(swiftIndex)._debugLevel = options._debugLevel;

    return _mapSingleReads(store, cnts, options, mode, swiftIndex, TFilterSpec());
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for SWIFT with Micro RNA
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TCounts,
    typename TSpec,
    typename TShape,
    typename TGapMode,
    typename TScoreMode,
    typename TMatchNPolicy,
    typename TFilterSpec>
int _mapSingleReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<RazerSPrefix, TGapMode, TScoreMode, TMatchNPolicy> const & mode,
    TFilterSpec)
{
    typedef FragmentStore<TFSSpec, TFSConfig>               TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore          TReadSeqStore;

//	typedef typename Value<TReadSeqStore>::Type				TRead;
//	typedef typename Infix<TRead>::Type						TReadInfix;
//	typedef StringSet<TReadInfix>							TReadSet;
    typedef TReadSeqStore                                   TReadSet;
    typedef Index<TReadSet, IndexQGram<TShape> >            TIndex;         // q-gram index

//	TReadSet readSet;
//	unsigned readCount = length(store.readSeqStore);
//	resize(readSet, readCount, Exact());
//
//	for (unsigned i = 0; i < readCount; ++i)
//		assign(readSet[i], prefix(store.readSeqStore[i], _min(length(store.readSeqStore[i]), options.prefixSeedLength)));
//
//	// configure q-gram index
//	TIndex swiftIndex(readSet, shape);
    TIndex swiftIndex(store.readSeqStore, shape);
    cargo(swiftIndex).abundanceCut = options.abundanceCut;
    cargo(swiftIndex)._debugLevel = options._debugLevel;

    return _mapSingleReads(store, cnts, options, mode, swiftIndex, TFilterSpec());
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different filters specs
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TCounts,
    typename TSpec,
    typename TShape,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TMatchNPolicy>
int _mapSingleReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & mode)
{
    if (options.threshold > 0)
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
        return _mapSingleReads(store, cnts, options, shape, mode, Swift<TSwiftSpec>());
    }
    else
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, void, Hamming_>::Type TPigeonholeSpec;
        return _mapSingleReads(store, cnts, options, Shape<Dna, OneGappedShape>(), mode, Pigeonhole<TPigeonholeSpec>());
    }
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different filters specs
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TCounts,
    typename TSpec,
    typename TShape,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TMatchNPolicy>
int _mapMatePairReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & mode)
{
    if (options.threshold > 0)
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
        return _mapMatePairReads(store, cnts, options, shape, mode, Swift<TSwiftSpec>());
    }
    else
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, void, Hamming_>::Type TPigeonholeSpec;
        return _mapMatePairReads(store, cnts, options, Shape<Dna, OneGappedShape>(), mode, Pigeonhole<TPigeonholeSpec>());
    }
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for single/mate-pair mapping
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TCounts,
    typename TSpec,
    typename TShape,
    typename TRazerSMode>
int _mapReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    TRazerSMode                       const & mode)
{
    //typedef FragmentStore<TFSSpec, TFSConfig>           TFragmentStore;
    //typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;
    //typedef typename Value<TReadSeqStore>::Type         TRead;
    //typedef Pattern<TRead const, MyersUkkonen>          TMyersPattern;  // verifier


    options.dRange = options.scoreDistanceRange;
    if (options.dRange == 0)
        options.dRange = 1 << 30;
    if (options.maxHits == 1 && !options.purgeAmbiguous)
        options.dRange = -1;
    else
        --options.dRange;


    int readCount = length(store.readSeqStore);
    resize(options.errorCutOff, readCount, Exact());

#ifndef RAZERS_BANDED_MYERS
    if (options.gapMode == RAZERS_GAPPED)
        resize(options.forwardPatterns, readCount, Exact());
#endif

    // set absolute error cut-offs
    Dna5String tmp;
    
    SEQAN_OMP_PRAGMA(parallel for private(tmp))
    for (int i = 0; i < readCount; ++i)
    {
        unsigned err = (unsigned)(length(store.readSeqStore[i]) * options.errorRate);
        options.errorCutOff[i] = (err < 255) ? err + 1 : 255;

        // optionally compute preprocessing information
#ifndef RAZERS_BANDED_MYERS
        if (!empty(store.readSeqStore[i]))
            continue;
        _patternMatchNOfPattern(options.forwardPatterns[i], options.matchN);
        _patternMatchNOfFinder(options.forwardPatterns[i], options.matchN);

#ifdef RAZERS_NOOUTERREADGAPS
        if (options.libraryLength >= 0 && (i & 1) == 1)
        {
            tmp = store.readSeqStore[i];
            reverseComplement(tmp);
            setHost(options.forwardPatterns[i], prefix(tmp, length(tmp) - 1));
        }
        else
            setHost(options.forwardPatterns[i], prefix(store.readSeqStore[i], length(store.readSeqStore[i]) - 1));
#else
        if (options.libraryLength >= 0 && (i & 1) == 1)
        {
            tmp = store.readSeqStore[i];
            reverseComplement(tmp);
            setHost(options.forwardPatterns[i], tmp);
        }
        else
            setHost(options.forwardPatterns[i], store.readSeqStore[i]);
#endif
#endif  // #ifdef RAZERS_BANDED_MYERS
    }

#ifdef _OPENMP
    if (options.threadCount == 0 || length(store.readNameStore) < MIN_PARALLEL_WORK)
#endif
    {
        // Sequential RazerS
        #ifdef RAZERS_MATEPAIRS
        if (options.libraryLength >= 0)
            return _mapMatePairReads(store, cnts, options, shape, mode);
        else
        #endif  // #ifndef RAZERS_MATEPAIRS
        return _mapSingleReads(store, cnts, options, shape, mode);
    }
#ifdef _OPENMP
    else
    {
        // Parallel RazerS
        #ifdef RAZERS_MATEPAIRS
        if (options.libraryLength >= 0)
            return _mapMatePairReadsParallel(store, cnts, options, shape, mode);
        else
        #endif  // #ifndef RAZERS_MATEPAIRS
        return _mapSingleReadsParallel(store, cnts, options, shape, mode);
    }
#endif
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different shapes
template <typename TFSSpec, typename TFSConfig, typename TCounts, typename TSpec, typename TRazersMode>
int _mapReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TRazersMode                       const & mode)
{
    Shape<Dna, SimpleShape>     ungapped;
    Shape<Dna, OneGappedShape>  onegapped;
    Shape<Dna, GenericShape>    gapped;

    // 2x3 SPECIALIZATION

    // select best-fitting shape
    if (stringToShape(ungapped, options.shape))
        return _mapReads(store, cnts, options, ungapped, mode);

//	if (stringToShape(onegapped, options.shape))
//		return _mapReads(store, cnts, options, onegapped, mode);
    if (stringToShape(gapped, options.shape))
        return _mapReads(store, cnts, options, gapped, mode);

    return RAZERS_INVALID_SHAPE;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different score modes
template <typename TFSSpec, typename TFSConfig, typename TCounts, typename TSpec, typename TAlignMode, typename TGapMode, typename TMatchNPolicy>
int _mapReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    RazerSMode<TAlignMode, TGapMode, Nothing, TMatchNPolicy> const)
{
    if (options.scoreMode == RAZERS_ERRORS)
        return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSErrors, TMatchNPolicy>());

/*	if (options.scoreMode == RAZERS_SCORE)
        return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSScore, TMatchNPolicy>());
    if (options.scoreMode == RAZERS_QUALITY)
        return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSQuality<>, TMatchNPolicy>());
*/  return RAZERS_INVALID_OPTIONS;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different gap and align modes
template <typename TFSSpec, typename TFSConfig, typename TCounts, typename TSpec>
int _mapReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options)
{
    if (options.matchN)
    {
        if (options.gapMode == RAZERS_GAPPED)
        {
//             if (options.alignMode == RAZERS_LOCAL)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSGapped, Nothing, NMatchesAll_>());
//             if (options.alignMode == RAZERS_PREFIX)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSGapped, Nothing, NMatchesAll_>());
            if (options.alignMode == RAZERS_GLOBAL)
                return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSGapped, Nothing, NMatchesAll_>());
        }
        else
        {
//             if (options.alignMode == RAZERS_LOCAL)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSUngapped, Nothing, NMatchesAll_>());
//             if (options.alignMode == RAZERS_PREFIX)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSUngapped, Nothing, NMatchesAll_>());
            if (options.alignMode == RAZERS_GLOBAL)
                return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSUngapped, Nothing, NMatchesAll_>());
        }
    }
    else
    {
        if (options.gapMode == RAZERS_GAPPED)
        {
            //     if (options.alignMode == RAZERS_LOCAL)
            //         return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSGapped, Nothing, NMatchesNone_>());
            //     if (options.alignMode == RAZERS_PREFIX)
            //         return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSGapped, Nothing, NMatchesNone_>());
            if (options.alignMode == RAZERS_GLOBAL)
                return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSGapped, Nothing, NMatchesNone_>());
        }
        else
        {
            // if (options.alignMode == RAZERS_LOCAL)
            //     return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSUngapped, Nothing, NMatchesNone_>());
            // if (options.alignMode == RAZERS_PREFIX)
            //     return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSUngapped, Nothing, NMatchesNone_>());
            if (options.alignMode == RAZERS_GLOBAL)
                return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSUngapped, Nothing, NMatchesNone_>());
        }
    }
    return RAZERS_INVALID_OPTIONS;
}

}

#endif
