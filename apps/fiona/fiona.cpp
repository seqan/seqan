// We are using the following parameters in the paper:
//
// -nt 8 -g 100286070 -id 1 -e 0.01 -dsr 3 -krr 0.01 -os 0.3 -or 1 -l 0 0 -i 1000
//
// Setting options triggers quite complex behaviours that we need to simplify later.

//#define FIONA_NOERROROPTIMIZATION   //enable mode to emulate error correction by random encounter
#define SEQAN_PROFILE		// enable time measuring
//#define FIONA_MEMOPT		// small suffix array values (<16mio reads of length <256)
#define FIONA_USE_SA        // use binary search in a suffix array for traversal
#define FIONA_REDUCE_MEMORY
//#define FIONA_OVERLAP_WITH_EDIT_DISTANCE  // allow indels in the overlap (instead of only mismatches)
#define FIONA_CONSENSUS_REDUCE_MEMORY

//#define FIONA_MAX_CORRECTIONS_PER_BASE 3  // record and limit the number of found corrections per base

#define FIONA_MATCH_N
#define FIONA_MAXIMIZE_OVERLAPSUM   // instead of maximizing the SUM of left and right, simply use the MAXIMUM of left and right
#define FIONA_NO_SEPARATE_OVERLAPSUM //in this mode just save the max over left and right in the linked list of corrections and only keep one correction per position
#define FIONA_INTERNAL_MEMORY

// debugging
//#define SEQAN_DEBUG_INDEX
//#define SEQAN_DEBUG
//#define SEQAN_VERBOSE
//#define SEQAN_VVERBOSE

// iodebug
//#define SEQAN_DEBUG_OR_TEST_
//#define SEQAN_HEADER_PIPE_DEBUG



#define FIONA_ALLOWINDELS	// allow for indels (chooses a less compact FragmentStore)

//    // currently, consensus works only without indels
//    #ifndef FIONA_OVERLAP_WITH_EDIT_DISTANCE
//        #define FIONA_CONSENSUS
//    #endif
//
//    //#define FIONA_FIXED_OVERLAP_ERRORS  // use fixed (ISMB) instead of error rate dependent threshold for overlap errors
//
//    // Dave's proposal to locally chose the operation with maximal support
//    #define FIONA_MAXIMIZE_SUPPORT
//
//    //#define FIONA_DISTANCE_BASED_ERROR_OPTIMIZATION // in this mode errors are corrected independent of type (mismatch or indels) but constraining the corrections to be apart at least min_k distance on a read


#ifdef FIONA_ILLUMINA

//  Illumina settings (currently set above)

  #define FIONA_FIXED_OVERLAP_ERRORS
  #define FIONA_CONSENSUS
  #undef FIONA_OVERLAP_WITH_EDIT_DISTANCE
  #undef FIONA_MAXIMIZE_SUPPORT
  #undef FIONA_DISTANCE_BASED_ERROR_OPTIMIZATION
  #define FIONA_BINARY_NAME "fiona_illumina"

#else

//  Indel settings

  #undef FIONA_FIXED_OVERLAP_ERRORS
  #undef FIONA_CONSENSUS
  #define FIONA_OVERLAP_WITH_EDIT_DISTANCE
  #define FIONA_MAXIMIZE_SUPPORT
  #define FIONA_DISTANCE_BASED_ERROR_OPTIMIZATION
  #define FIONA_BINARY_NAME "fiona"

#endif

#include <seqan/platform.h>

#if defined(_OPENMP)
    #include <omp.h>
    #define SEQAN_PARALLEL      // Only enable parallelism in fiona if OpenMP is enabled.
    #define FIONA_PARALLEL		// divide suffix tree into subtrees for each possible 3-gram

    #if defined(STDLIB_GNU)
        #include <parallel/algorithm>
    #endif
#else
    #pragma message("Please enable OpenMP.")
#endif  // #ifdef _OPENMP

// The q-gram length used for the q-gram index.  This has to be hard-coded as a precompiler definition since it is part
// of the template parameters for the indices.
#ifndef QGRAM_LENGTH
#define QGRAM_LENGTH 10                    // must be less or equal to fromLevel
#endif

// The hardcoded maximal indel length.
const unsigned int MAX_INDEL_LENGTH = 4;
// The hardcoded maximal number of rounds when using auto round detection.
const unsigned int MAX_NUM_ROUND = 6;

// The program's version and release date are kept here at the top of the file for easier maintenance.
const char * PROGRAM_VERSION = "0.2";

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
//#include <sys/resource.h>

// TODO (hugues) 
// 1_ Update all formulas of mixed poisson and binomial to use
// boost library
// 2_ Use the boost functions for all computations involving binomial/poisson

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/version.h>

#if SEQAN_VERSION_MAJOR == 1 && SEQAN_VERSION_MINOR == 4
// in SeqAn 1.4.x the SA tree and the View classes were in extras
#include <../../include/seqan/index/index_sa_stree.h>
#include <../../include/seqan/basic/basic_view.h>
#include <../../include/seqan/sequence/iterator_range.h>
#endif

#include "index_qgram_parallel.h"

//Boost Math headers
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/binomial.hpp>

using boost::math::binomial;

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;



// TODO(holtgrew): This raises a warning with Boost 1.42. Deactivate warnings, activate again afterwards. The correct #pragma has to be used for each supported compiler.
//#include <boost/math/distributions/normal.hpp>

//#define MEDIAN

using namespace seqan;

#ifdef FIONA_ALLOWINDELS

	// NOTE:
	// Currently we have to change the StringSet spec of the readSeqStore
	// to Owner as ConcatDirect<> (default) is not able to notice if a read
	// changes its size (after correction).
	struct FionaStoreConfig:
		public FragmentStoreConfig<>
	{
        //typedef String<Dna5, Packed<> > TReadSeq;
        typedef String<Dna5> TReadSeq;
		typedef Owner<>	TReadSeqStoreSpec;
        typedef Owner<> TReadNameStoreSpec;
	};

	typedef FragmentStore<void, FionaStoreConfig> TFionaFragStore;
    typedef Value<TFionaFragStore::TReadSeqStore>::Type TRead;
    typedef Infix<TRead>::Type TReadPrefix;

#else

	struct FionaStoreConfig:
		public FragmentStoreConfig<>
	{
        //typedef String<Dna5, Packed<> > TReadSeq;
        typedef String<Dna5> TReadSeq;
        typedef Owner<> TReadNameStoreSpec;
	};

	typedef FragmentStore<void, FionaStoreConfig> TFionaFragStore;
    typedef Value<TFionaFragStore::TReadSeqStore>::Type TReadPrefix;

#endif

typedef StringSet<TReadPrefix> TReadPrefixes;
#ifdef FIONA_USE_SA
typedef Index<TFionaFragStore::TReadSeqStore, IndexSa<> > TFionaIndex;
#else
typedef Index<TFionaFragStore::TReadSeqStore, IndexWotd<> > TFionaIndex;
#endif
typedef Index<TReadPrefixes, IndexQGram< Shape<Dna5, UngappedShape<QGRAM_LENGTH> > > > TFionaQgramIndex;
//typedef Index<TFionaFragStore::TReadSeqStore, IndexQGram< Shape<Dna5, UngappedShape<5> > > > TFionaQgramIndex;


struct FionaPoisson_;
struct FionaExpected_;
struct FionaCount_;
struct FionaPoissonSens_;
struct FionaPoissonClassif_;

typedef Tag<FionaPoisson_> const FionaPoisson;
typedef Tag<FionaPoissonSens_> const FionaPoissonSens;
typedef Tag<FionaExpected_> const FionaExpected;
typedef Tag<FionaCount_> const FionaCount;
typedef Tag<FionaPoissonClassif_> const FionaPoissonClassif;


struct FionaCorrectedError
{
	unsigned int   correctReadId;
//	unsigned int   occurrences;
	unsigned short errorPos;
	unsigned short correctPos;
	unsigned short overlap;
	signed   char  indelLength;		// 0..mismatch, <0..deletion, >0..insertion
//	unsigned char  mismatches;
};

// Enum for representing the Fiona method.

enum FionaMethod
{
    CONTROL_FP,
    EXPECTED,
    CONTROL_FN,
    COUNT,
    CLASSIFIER
};

// Return a string with the name for the given Fiona method.

char const * methodName(FionaMethod m)
{
    switch (m)
    {
        case CONTROL_FP:
            return "CONTROL FALSE POSITIVES";
        case EXPECTED:
            return "EXPECTED";
        case CONTROL_FN:
            return "CONTROL FALSE NEGATIVES";
        case COUNT:
            return "COUNT";
        case CLASSIFIER:
            return "CLASSIFIER";
    }
    return "INVALID";
}

// Convert a valid method name into a FionaMethod.

FionaMethod methodForName(seqan::CharString const & str)
{
    if (str == "control_fp")
        return CONTROL_FP;
    else if (str == "expected")
        return EXPECTED;
    else if (str == "control_fn")
        return CONTROL_FN;
    else if (str == "count")
        return COUNT;
    else
        return CLASSIFIER;
}

struct FionaOptions
{
    // Verbosity:  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

	int64_t genomeLength;
	double strictness;
	unsigned acceptedMismatches;
	int maxIndelLength;
	bool autolevel;
	int fromLevel;
	int toLevel;
	unsigned cycles;
    unsigned cycle;
	double errorrate;
	double overlap_errorrate;
	double oddserrorreads;
	double wovsum;
    int debugRead, corrRead;
    unsigned packagesPerThread;
    int loopLevel;
    // Whether or not to append correction information into the FASTA output headers.
    bool appendCorrectionInfo;
    double kmerAbundanceCutoff;
    double kmerStdDevCutOff;
    int depthSampleRate;
    double timeComputeOverlapSum;

    // The number of errors to correct per read relative to the read length.
    double relativeErrorsToCorrect;

    CharString inputFilename;
    CharString outputFilename;

    FionaMethod method;
    int numThreads;

// internal parameters

	String<double> expectedTheoretical;
    String<int> errorCutoffs;
    String<unsigned> repeatCutoffs;
	matrix<double> overlapSumCutoffs;
	//create new String with allowed errors per Read that is used to skip further corrections on reads
	String<unsigned char> allowedCorrectionsPerRead;

	bool limitCorrPerRound;
	bool trimNsOnOutput;
    unsigned numSuperPackages;

	FionaOptions()
	{
        verbosity = 0;
        method = CLASSIFIER;
        numThreads = 1;
        limitCorrPerRound = true;
        trimNsOnOutput = true;
		genomeLength = 0;
		strictness = 0.0001;
		acceptedMismatches = 1;
		maxIndelLength = 1;
		cycles = 6;
        cycle = 1;
		autolevel = false;
		fromLevel = 0;
		toLevel = 0;
#ifdef FIONA_ILLUMINA
		errorrate = 0.01;
        kmerAbundanceCutoff = 0.01;
#else
		errorrate = 0.05;
        kmerAbundanceCutoff = 0.05;
#endif
		overlap_errorrate = 0;
		oddserrorreads = 0;
		wovsum = 0.3;
        debugRead = -1;
        corrRead = -1;
        packagesPerThread = 100;
        kmerStdDevCutOff = 2.0;
        depthSampleRate = 3;
        relativeErrorsToCorrect = 0.02;
        timeComputeOverlapSum = 0;
        loopLevel = -1;
        appendCorrectionInfo = false;
        numSuperPackages = 10;
	}
};

// Return C-style string "YES"/"NO" depending on the value of b.  Useful for printing options.

char const * yesNo(bool b)
{
    return b ? "YES" : "NO";
}

// Print options to out.

void printOptions(std::ostream & out, FionaOptions const & options)
{
    out << "__OPTIONS________________________________________________________\n"
        << "\n"
        << "MODE INFORMATION\n"
        << "\n";
#if defined(FIONA_PARALLEL)
    out << "  PARALLEL MODE            YES\n";
#else  // #if defined(FIONA_PARALLEL)
    out << "  PARALLEL MODE            NO\n";
#endif  // #if defined(FIONA_PARALLEL)
#if defined(FIONA_NOERROROPTIMIZATION)
    out << "  RANDOM ENCOUNTER MODE    YES\n";
#else  // #if defined(FIONA_NOERROROPTIMIZATION)
    out << "  RANDOM ENCOUNTER MODE    NO\n";
#endif  // #if defined(FIONA_NOERROROPTIMIZATION)
#if defined(FIONA_ALLOWINDELS)
    out << "  ALLOW INDELS MODE        YES\n";
#else  // #if defined(FIONA_ALLOWINDELS)
    out << "  ALLOW INDELS MODE        NO\n";
#endif  // #if defined(FIONA_ALLOWINDELS)
#if defined(FIONA_USE_SA)
    out << "  SUFFIX ARRAY MODE        YES\n";
#else  // #if defined(FIONA_USE_SA)
    out << "  SUFFIX ARRAY MODE        NO\n";
#endif  // #if defined(FIONA_USE_SA)
    out << "\n"
        << "CONSTANTS\n"
        << "\n"
        << "  K-MER LENGTH             " << QGRAM_LENGTH << "\n"
        << "  MAX INDEL LENGTH         " << MAX_INDEL_LENGTH << "\n"
        << "  MAX NUM ROUNDS           " << MAX_NUM_ROUND << "\n"
        << "\n"
        << "OPTIONS\n"
        << "\n"
        << "  METHOD                   " << methodName(options.method) << "\n"
        << "  GENOME LENGTH            " << options.genomeLength << "\n"
        << "  STRICTNESS               " << options.strictness << "\n"
        << "  ACCEPTED MISMATCHES      " << options.acceptedMismatches << "\n"
        << "  MAX INDEL LENGTH         " << options.maxIndelLength << "\n"
        << "  CYCLES                   " << options.cycles << "\n"
        << "  CYCLE                    " << options.cycle << "\n"
        << "  AUTOLEVEL                " << yesNo(options.autolevel) << "\n"
        << "  FROM LEVEL               " << options.fromLevel << "\n"
        << "  TO LEVEL                 " << options.toLevel << "\n"
        << "  ERROR RATE               " << options.errorrate << "\n"
        << "  ODDS ERROR READS         " << options.oddserrorreads << "\n"
        << "  WOV SUM                  " << options.wovsum << "\n"
        << "  RELATIVE ERRORS TO CORR. " << options.relativeErrorsToCorrect << "\n";
    if (options.debugRead != -1)
        out << "  DEBUG READ               " << options.debugRead << "\n";
    if (options.corrRead != -1)
        out << "  CORR READ                " << options.corrRead << "\n";
    out << "  PACKAGES PER THREAD      " << options.packagesPerThread << "\n"
        << "  KMER ABOUNDANCE CUTOFF   " << options.kmerAbundanceCutoff << "\n"
        << "  DEPTH SAMPLE RATE        " << options.depthSampleRate << "\n"
        << "  TIME COMPUTE OVERLAP SUM " << options.timeComputeOverlapSum << "\n";
    if (options.loopLevel != -1)
        out << "  LOOP LEVEL               " << options.loopLevel << "\n";
    out << "  APPEND CORRECTION INFO   " << yesNo(options.appendCorrectionInfo) << "\n"
        << "  NUM THREADS              " << options.numThreads << "\n"
        << "\n";
}

// used for profiling
struct FionaResources
{
    unsigned long   bucketBegin;
    unsigned long   bucketEnd;
    double          cpuTime;
    unsigned        investigatedNodes;
    unsigned        putCorrections;

    FionaResources():
        bucketBegin(0),
        bucketEnd(0),
        cpuTime(0),
        investigatedNodes(0),
        putCorrections(0) {}

    bool operator < (FionaResources const &right) const
    {
        return cpuTime > right.cpuTime;
    }
};

//Temp global variable num of families
unsigned int nfamilies = 0;


/*
new struct for Fiona to record several errors and their corrections per read
with support for different indels and both strands.
The struct is saved in a String<Correction> that essentially is a save efficient 
linked list of these structs. Therefore the first variable nextCorrection points to the 
position in the String where the next correction for that read position can be found.
The last linked list item for a read position holds the maxValue for int.
The array overlap records the maximum overlap sum observed for the forward and the reverse strand.
Note that for indel events special care must be taken if both strands are considered
as the indel position may differ relative to orientation. 
Note that the maximum recordable overlapsum is bounded by 65.535.

*/
struct CorrectionIndelPos
{
    unsigned int    nextCorrection;    // -1u..last correction
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
    unsigned int    correctReadId;
    unsigned short  correctPos;
#endif
    unsigned short  errorPos;
    unsigned char   onReverse:1;        // 0 if forward strand, 1 if reverse strand similar to boolean strand variable
    unsigned char   foundCorrections:7;
    unsigned short  overlap[2];         // 0..forward, 1..reverse
    Dna5            correctSeq[MAX_INDEL_LENGTH];
    signed char     indelLength;        // -x..insertion
                                        //  0..mismatch
                                        //  x..deletion
};

namespace seqan
{

	struct CargoQgramIndex
	{ 
        FionaOptions* optionsPtr;
	};

/*restriction for certain levels - between max and min, also table with frequency may be to use eventually TODO*/
	struct FionaNodeConstraints
	{ 
        unsigned replen_max; 
        unsigned replen_min;
//        String<unsigned> *repeatCutoffs;  
//        std::map<unsigned,double> frequency;
	};

	template <>
	struct Cargo<TFionaQgramIndex>
    {
		typedef CargoQgramIndex Type;
	}; 

	template <>
	struct Cargo<TFionaIndex> { 
		typedef FionaNodeConstraints Type; 
	}; 

#ifdef FIONA_MEMOPT

	typedef Pair<
		unsigned,
		unsigned,
		BitCompressed<24, 8>	// max. 16M reads of length < 256
	> TSAValue;

#else

	typedef Pair<
		unsigned,				// many reads
		unsigned short,			// of arbitrary length
		Pack
	> TSAValue;

#endif

	template <>
	struct SAValue< TFionaIndex > 
	{
		typedef TSAValue Type;
	};

	template <>
	struct SAValue< TFionaQgramIndex > 
	{
		typedef TSAValue Type;
	};

#ifndef FIONA_INTERNAL_MEMORY
	// use a mmap string for storing the q-grams
	template <>
	struct Fibre< TFionaQgramIndex, FibreSA > 
	{
#ifdef FIONA_REDUCE_MEMORY
		typedef String<TSAValue, External<ExternalConfigLarge<> > > Type;
#else
		typedef String<TSAValue, MMap<> > Type;
#endif
	};
#endif

#ifdef FIONA_PARALLEL

	template <>
	struct Fibre< TFionaIndex, FibreSA >
	{
		typedef Fibre< TFionaQgramIndex, QGramSA >::Type TSA;
#ifdef FIONA_REDUCE_MEMORY
        typedef Value<TSA>::Type TValue;
        typedef Range<TValue*> Type;
#else
		typedef Infix<TSA>::Type Type;
#endif
	};

#endif

    // Copy from store_io.h, could/should go into library.
    template <typename TFSSpec, typename TFSConfig, typename TFileName>
    bool loadReadsNoNames(FragmentStore<TFSSpec, TFSConfig> &store, TFileName &fileName, FionaOptions const & options)
    {
        //StringSet<String<Dna5, Packed<> >, Owner<ConcatDirect<> > > reads;

        SeqFileIn seqFile;
        if (!open(seqFile, toCString(fileName)))
            return false;

        // read sequences
    //    String<Dna5Q> seq;
    //    CharString qual;

        CharString _id;
        String<Dna5> seq;
        unsigned i = 0;

        while (!atEnd(seqFile))
        {
            readRecord(_id, seq, seqFile);
//            appendRead(store, seq/*, _id*/);
            appendValue(store.readSeqStore, seq);
            if (options.verbosity >= 1 && ++i % 100000 == 0)
                std::cout<<'.'<<std::flush;
        }
//        store.readSeqStore = reads;
        return true;
    }

    template <typename THash>
    inline bool
    hashContainsN(THash h)
    {
        // we assume that the hash value was computed for a Dna5 q-gram (sigma=5)
        while (h != 0)
        {
            if (h % 5 == 4)
                return true;
            h /= 5;
        }
        return false;
    }

    template <typename TDir>
    struct GreaterBucketSize
    {
        TDir const &dir;
        
        GreaterBucketSize(TDir const &dir_): dir(dir_) {}
        
        inline bool
        operator () (unsigned a, unsigned b)
        {
            return dir[a] > dir[b];
        }
    };

    template <typename TDir>
    inline void
    maskRepeatBuckets(TDir &dir, FionaOptions const & options)
    {
        typedef typename Value<TDir>::Type                              TDirValue;
        typedef typename MakeSigned<typename Size<TDir>::Type>::Type    TDirSize;

        const TDirValue PURGE_BUCKET = (TDirValue)-1;
        TDirSize dirLen = length(dir);

        if (options.verbosity >= 1)
            std::cerr << "Purge repetitive k-mers ........... " << std::flush;

        // extract bucket numbers
        uint64_t suffixCount = 0;
        String<unsigned> bktIdx;

//        #pragma omp parallel for reduction(+:suffixCount)
        for (TDirSize i = 0; i < dirLen - 1; ++i)
            if (dir[i] != PURGE_BUCKET && dir[i] > 0)
            {
                suffixCount += dir[i];
//                SEQAN_OMP_PRAGMA(critical)
                appendValue(bktIdx, i);
            }
        
        // sort them descendingly by bucket size
        sort(bktIdx, GreaterBucketSize<TDir>(dir), Parallel());
        
        // mask for removal of the largest buckets that 
        // contain overall at most 2% of all suffixes
        if (options.verbosity >= 1)
            std::cerr << " suffixes: " << suffixCount << std::endl;
        uint64_t threshN1 = (uint64_t)(suffixCount * options.kmerAbundanceCutoff);
        Dna5String kmer;
        suffixCount = 0;
        for (unsigned i = 0; i < length(bktIdx); ++i)
        {
            unsigned bkt = bktIdx[i];
            suffixCount += dir[bkt];
            if (suffixCount < threshN1)
            {
                dir[bkt] = PURGE_BUCKET;
                unhash(kmer, bkt, QGRAM_LENGTH);
                if (options.verbosity >= 2)
                    std::cerr << kmer << ' ' << std::flush;
            } else
                break;
        }
        if (options.verbosity >= 2)
            std::cerr << std::endl;
    }

    template <typename TDir>
    inline void
    maskRepeatBuckets2(TDir &dir, FionaOptions const & options)
    {
        typedef typename Value<TDir>::Type                              TDirValue;
        typedef typename MakeSigned<typename Size<TDir>::Type>::Type    TDirSize;
        typedef typename Iterator<String<unsigned> >::Type              TBktIter;

        const TDirValue PURGE_BUCKET = (TDirValue)-1;
        TDirSize dirLen = length(dir);

        if (options.verbosity >= 1)
            std::cerr << "Purge repetitive k-mers ........... " << std::flush;

        // extract bucket numbers
        uint64_t suffixCount = 0;
        String<unsigned> bktIdx;

//        #pragma omp parallel for reduction(+:suffixCount)
        for (TDirSize i = 0; i < dirLen - 1; ++i)
            if (dir[i] != PURGE_BUCKET && dir[i] > 0)
            {
                suffixCount += dir[i];
//                SEQAN_OMP_PRAGMA(critical)
                appendValue(bktIdx, i);
            }
        
        // sort them descendingly by bucket size
        sort(bktIdx, GreaterBucketSize<TDir>(dir), Parallel());

        TBktIter itFirst = begin(bktIdx, Standard());
        TBktIter itLast = end(bktIdx, Standard());

        // omit the disabled buckets in front
        while (dir[*itFirst] == PURGE_BUCKET && itFirst != itLast)
            ++itFirst;

        // omit the empty buckets in the back
        do {
            --itLast;
        } while (dir[*itLast] == 0 && itLast != itFirst);
        ++itLast;

        // compute k-mer median
        unsigned n = itLast - itFirst;
        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<dir[*itFirst]<<'\t'<<dir[*(itFirst+1)]<<'\t'<<dir[*(itFirst+2)]<<std::endl;
        std::cout<<dir[*(itLast-1)]<<'\t'<<dir[*(itLast-2)]<<'\t'<<dir[*(itLast-3)]<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::endl;
        TDirValue median = dir[*(itFirst + n / 2)];
        double stdDev = 0;
        for (TBktIter it = itFirst; it != itLast; ++it)
        {
            double diff = (double)dir[*it] - (double)median;
            stdDev += diff * diff;
        }
        
        // compute k-mer standard deviation
        stdDev = std::sqrt(stdDev / (n - 1));
        TDirValue cutOff = median + stdDev * options.kmerStdDevCutOff;

        if (options.verbosity >= 1)
        {
            std::cerr << " k-mer median: " << median << std::endl;
            std::cerr << " k-mer stddev: " << stdDev << std::endl;
            std::cerr << " k-mer cutoff: " << cutOff << std::endl;
        }

        // remove bucket above the cut-off
        Dna5String kmer;
        for (TBktIter it = itFirst; it != itLast; ++it)
        {
            if (dir[*it] > cutOff)
            {
                dir[*it] = PURGE_BUCKET;
                unhash(kmer, *it, QGRAM_LENGTH);
                if (options.verbosity >= 2)
                    std::cerr << kmer << ' ' << std::flush;
            }
        }

        if (options.verbosity >= 2)
            std::cerr << std::endl;
    }

    inline bool
    _qgramDisableBuckets(TFionaQgramIndex &index)
    {
        typedef TFionaQgramIndex                    TIndex;
        typedef Fibre<TIndex, QGramDir>::Type       TDir;
        typedef Fibre<TIndex, QGramShape>::Type     TShape;
        typedef Value<TDir>::Type                   TDirValue;
        typedef MakeSigned<Size<TDir>::Type>::Type  TDirSize;
        typedef Host<TShape>::Type                  TValue;
        typedef ValueSize<TValue>::Type             TValueSize;

        const TDirValue PURGE_BUCKET = (TDirValue)-1;

        TDir &dir = indexDir(index);
        TDirSize dirLen = length(dir);
        TShape &shape = indexShape(index);
        String<TValue> kmer;

        // 1. manually remove all homopolymer repeats
        for (TValueSize x = 0; x < ValueSize<TValue>::VALUE; ++x)
        {
            clear(kmer);
            resize(kmer, length(shape), (TValue)x);
            dir[hash(shape, begin(kmer, Standard()))] = PURGE_BUCKET;
        }

        // 2. mask k-mers with Ns
        SEQAN_OMP_PRAGMA(parallel for)
        for (TDirSize i = 0; i < dirLen - 1; ++i)
            if (hashContainsN(i))
                dir[i] = PURGE_BUCKET;

        // 3. mask k-mers from repeat regions
        maskRepeatBuckets(dir, *(cargo(index).optionsPtr));
//        maskRepeatBuckets2(dir, *(cargo(index).optionsPtr));

        return true;
    }

    // check whether a string can be exactly be overlapped with itself
    // try different overlap offsets 1,...,6
    template <typename TString>
    inline bool
    isSelfRepetitive(TString const &str)
    {
        typedef typename Iterator<TString const, Standard>::Type TIterator;
        TIterator itBegin = begin(str, Standard());
        TIterator itEnd = end(str, Standard());
        
        unsigned maxOverlap = _min(6, (itEnd - itBegin) / 2);
        for (unsigned ofs = 1; ofs <= maxOverlap; ++ofs)
        {
            TIterator it1 = itBegin;
            TIterator it2 = itBegin + ofs;
            for (; it2 != itEnd; ++it1, ++it2)
                if (*it1 != *it2) break;
            if (it2 == itEnd)
            {
//                std::cerr << ofs << '\t' << str << " skipped." << std::endl;
                return true;
            }
        }
        return false;
    }


	/*TODO THIS FONCTION CAN BE CHANGED FOR THE FREQUENCY - here just one experience*/
	/*by the use also the frequency for A,T,G,C*/
	/*higher frequency - high level as min in which we will begin the searching*/

	/*hide the node between certain level*/
	template <typename TSpec>
	inline bool nodeHullPredicate(Iter<TFionaIndex, VSTree<TopDown<TSpec> > > &it)
	{	
		//Hugues: to parse all the node levels, I would use nodeDepth, e.g.
		//return nodeDepth(it) < cargo(container(it)).replen_max;
		return parentRepLength(it) <= cargo(container(it)).replen_max;
    }

	template <typename TSpec>
	inline bool nodePredicate(Iter<TFionaIndex, VSTree<TopDown<TSpec> > > &it)
	{
//		return true;
		FionaNodeConstraints &cons = cargo(container(it));
		unsigned repLen = parentRepLength(it);
		//the same here why isn't it nodeDepth ?
		//unsigned repLen = nodeDepth(it);
		/*TODO may utilise >=*/
		return cons.replen_min <= repLen && repLen <= cons.replen_max;
	}

}  // namespace seqan


// fill an array of type alphabet with the correction string
// of the correct read, in case the reverse sequence is seeked the 
// readID is changed accordingly.
template <typename TAlphabet, typename TValue,typename TValue2,typename TFragmentStore>
inline void getCorrectionString(TAlphabet array[], //the array has to be at least of length abs(indelLength)
                                signed char indelLength,
                                TValue correctReadId,
                                TValue2 correctPos,
                                bool strand,
                                TFragmentStore &store)
{
    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type         TRead;
    typedef typename Iterator<TRead, Standard>::Type    TReadIterator;

    SEQAN_ASSERT_LEQ(indelLength, 0);
    if (strand)
    {
        // change the pos in the correctRead and the readID
        unsigned readCount = length(store.readSeqStore) / 2;

        // switch to reverse-complements
        if (correctReadId < readCount)
            correctReadId += readCount;
        else
            correctReadId -= readCount;

        // mirror positions
        correctPos = length(store.readSeqStore[correctReadId]) - correctPos;

		// for insertions the position is correct already
        if (indelLength == 0)
            --correctPos;
    }

    // from here on the correctPos and the readID should refer
    // to the forward strand respective to the error position
    TReadIterator it = begin(store.readSeqStore[correctReadId], Standard()) + correctPos;
    array[0] = *it;    // always copy first char (also for mismatches when indelLength==0)
    for (int i = 1; i < -indelLength; ++i)
        array[i] = *(++it); //store.readSeqStore[correctReadId][correctPos+i];
}

/*
update a Correction entry with a higher overlap sum if necessary
*/
	
template <typename TCorrection,typename TValue,typename TValueShort,typename TValue2,typename TAlphabet>
inline void updateCorrectionEntry(String<TCorrection> &correctionList,
                                  TValue posInCorrectionList,
                                  //			TValue2 correctReadId,
                                  //			TValue3 correctPos,
                                  TValueShort overlap,
                                  bool strand,
                                  TValue correctReadId,
				  TValue2 indelLength,
				  TAlphabet &correctSeq,
                                  TValue2 previousIndelLength)
{
    TValue pos = (strand)? 1: 0;
#ifdef FIONA_NO_SEPARATE_OVERLAPSUM
     pos=0;
     bool change = false;
#else
	(void) previousIndelLength;
        (void) correctSeq;
        (void) indelLength;
#endif

    //because we are querying this position for a potential update we had found a correction
    correctionList[posInCorrectionList].foundCorrections++;
    if (correctionList[posInCorrectionList].overlap[pos] < overlap)
    {
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
        correctionList[posInCorrectionList].correctReadId = correctReadId;  // for debugging only
#else
        (void)correctReadId;
#endif
        //			correctionList[posInCorrectionList].correctPos= correctPos;
#ifdef FIONA_NO_SEPARATE_OVERLAPSUM     
       change = true;
#endif
        correctionList[posInCorrectionList].overlap[pos] = overlap;
    }
#ifdef FIONA_NO_SEPARATE_OVERLAPSUM
         //check if the overlap sum is the same put the type of correction preferred
        else if(correctionList[posInCorrectionList].overlap[pos] == overlap){
	   if(previousIndelLength ==0)
                  return;
           if(indelLength ==0)
		change =true;
	   else if(indelLength > previousIndelLength)
		change= true;
	}
      if(change){
        //update the indeLength and string as well
	correctionList[posInCorrectionList].indelLength = indelLength;
         // add seq to array in the struct
     correctionList[posInCorrectionList].correctSeq[0] = correctSeq[0];    // we always copy the first char (even for deletions, which are rare)
     for(int i = 1; i < -indelLength; ++i)           // copy insertion string (indelLength < 0)
        correctionList[posInCorrectionList].correctSeq[i] = correctSeq[i];
        }
#endif 

    return;
}

/*
Basic function to fill a new CorrectionIndelPos with values.
The nextCorrection field is initialized with highest integer value.
*/

template <
    typename TCorrection,
    typename TPos1,
    typename TPos2,
    typename TAlphabet,
    typename TOverlap,
    typename TIndel,
    typename TSize,
    typename TId >
inline void fillCorrection(TCorrection &newCorrection,
                           //	TValue correctReadId,
                           TPos1 correctPos,
                           TPos2 errorPos,
                           TAlphabet correctSeq[],  //this should always be in forward direction
                           TOverlap overlap,
                           bool strand,
                           TIndel indelLength,
                           TSize readLength,
                           TId correctReadId)
{
    //TValue empty=maxValue(readLength);
    //fill Correction struct
    newCorrection.nextCorrection = std::numeric_limits<unsigned>::max();  // it will be the last correction in the linked list
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
    newCorrection.correctReadId = correctReadId;        // only for debugging purposes
    newCorrection.correctPos = correctPos;
#else
    (void)correctReadId;
    (void)correctPos;
#endif
    newCorrection.indelLength = indelLength;
    newCorrection.foundCorrections=1;
    if (!strand)
    {
        // forward
        newCorrection.errorPos = errorPos;
        newCorrection.overlap[0] = overlap;      // 0..forward, 1..reverse
        newCorrection.overlap[1] = 0;
	newCorrection.onReverse = 0;
    }
    else
    {
        // mirror positions
        newCorrection.errorPos = readLength - errorPos;
	newCorrection.onReverse = 1;
        // for insertions the position is correct already
        if (indelLength == 0)
            --newCorrection.errorPos;
        else if (indelLength > 0)  // deletion
            newCorrection.errorPos -= indelLength;

        newCorrection.overlap[1] = overlap;
        newCorrection.overlap[0] = 0;      // 0..forward, 1..reverse
    }
#ifdef FIONA_NO_SEPARATE_OVERLAPSUM
        //put is always in overlap[0]
        newCorrection.overlap[0] = overlap;      // 0..forward, 1..reverse
        newCorrection.overlap[1] = 0;
#endif

    // add seq to array in the struct
    newCorrection.correctSeq[0] = correctSeq[0];    // we always copy the first char (even for deletions, which are rare)
    for(int i = 1; i < -indelLength; ++i)           // copy insertion string (indelLength < 0)
        newCorrection.correctSeq[i] = correctSeq[i];
}


/*
	existsCorrectionAtPos computes given the error type and length
	if it corresponds to the same position.  
	The position has to be converted
	depending if on a different strand. 
	Example (forward sequence):
	A C G T A C G 
	0 1 2 3 4 5 6
	
	mismatch at forward position 2 
	Pos(reverse) = readLength - Pos(forward) -1 = 7-2-1=4

	insertion at forward position 2
	Pos(reverse) = readlength - Pos(forward) = 7-2 = 5 
	
	deletion of size del at position 2
	Pos(reverse) = readlength - Pos(forward) - del = 7-2-2 = 3
	
	We assume that the strand of the previousError that is given to the function
	is forward and as such that previousErrorPos is on the forward strand of the read
	If the position is the same the character(s) for the correction are compared with the entry 

*/
template <typename TValueShort, typename TValueShort2, typename TValueLength,typename TAlphabet, typename TAlphabet2>
inline bool existsCorrectionAtPos(TValueShort previousErrorPos,
                                  //	TValue2 previousCorrectReadId,
                                  signed char previousIndelLength,
                                  TAlphabet2 previousCorrectSeq[],
                                  TValueShort2 errorPos,
                                  //	TValue3 correctReadId,
                                  bool strand,
                                  signed char indelLength,
                                  TValueLength readLength,
                                  TAlphabet correctSeq[])
{
    (void)previousCorrectSeq;
    (void)correctSeq;
    //std::cerr << "in existsCorrectionAtPos, previousindelLength, indelLength"<<(int)previousIndelLength<<","<<(int)indelLength<<std::endl;

#ifdef FIONA_NO_SEPARATE_OVERLAPSUM
        //make extra part for new mode to avoid clash with FIONA_NOERROROPTIMIZATION  
	//here change strand pos already
        if (strand)
        {
            // mirror positions
            errorPos = readLength - errorPos;

            //for insertions the position is correct already
            if (indelLength == 0)
                --errorPos;
            else if (indelLength > 0)  // deletion
                errorPos -= indelLength;
        }

        if (errorPos == previousErrorPos)
            return true;
	else
            return false;

#endif

    if (previousIndelLength == indelLength)
    {

        //here change strand pos already
        if (strand)
        {
            // mirror positions
            errorPos = readLength - errorPos;

            //for insertions the position is correct already
            if (indelLength == 0)
                --errorPos;
            else if (indelLength > 0)  // deletion
                errorPos -= indelLength;
        }

#ifdef FIONA_NOERROROPTIMIZATION
        if (errorPos == previousErrorPos)
            return true;

#else    //this is the original part
        if (indelLength == 0)
        {
            //check here if position and alphabet character for correction is the same
            return errorPos == previousErrorPos && correctSeq[0] == previousCorrectSeq[0];
        }
        else if (indelLength > 0) //deletion in read no string to compare
        {
            return errorPos == previousErrorPos;
        }
        else
        {
            SEQAN_ASSERT_LT(indelLength, 0);
            //insertion in read
            if (previousErrorPos == errorPos)
            {
                int lcp = 0;
                for (; lcp < -indelLength; ++lcp)
                    if (correctSeq[lcp] != previousCorrectSeq[lcp])
                        break;
                return lcp == -indelLength;
            }
        }
#endif
    }
    return false;
}


// this lock augments a string class by thread-safety as follows:
//  - supports multiple concurrent readers (possibly waiting for writer to finish)
//  - supports only a single writer at a time (possibly waiting for readers or other writers to finish)
//  - the writer has higher priority than all readers
struct StringLock
{
    volatile unsigned readers;
    volatile unsigned writers;

    StringLock():
        readers(0),
        writers(0)
    {
    }
};

static StringLock correctionListLock;


inline void
lockReading(StringLock &lock)
{
    do
    {
        // wait for the end of a write access
        while (lock.writers != 0)
        {}

//        SEQAN_OMP_PRAGMA(atomic)
//        ++lock.readers;
        atomicInc(lock.readers);

        if (lock.writers == 0)
            break;

        // writer hasn't noticed us -> retry
//        SEQAN_OMP_PRAGMA(atomic)
//        --lock.readers;
        atomicDec(lock.readers);
    } while (true);
}

inline void
unlockReading(StringLock &lock)
{
//    SEQAN_OMP_PRAGMA(atomic)
//    --lock.readers;
    atomicDec(lock.readers);
}

inline void
lockWriting(StringLock &lock)
{
    // wait until we are the only writer
    while (atomicCas(lock.writers, 0u, 1u) != 0)
    {}

    // wait until all readers are done
    while (lock.readers != 0)
    {}
}

inline void
unlockWriting(StringLock &lock)
{
    lock.writers = 0;
}




/*
new function that returns the number of corrections already entered for a read position
	
*/
template <typename TCorrection,typename TValue, typename TId1, typename TPos,typename TReadStore>
inline unsigned
getFoundCorrections(
    String<TCorrection> const &correctionList,
    String<TValue> &firstCorrectionForRead,
    TId1 erroneousReadId,
    TPos errorPos,
    TReadStore & store)
{
    unsigned numReads = length(store.readSeqStore) / 2;

    // project errorPos on forward strand
    if (erroneousReadId >= numReads)
	{
        erroneousReadId -= numReads;
        // mirror positions
        errorPos = length(store.readSeqStore[erroneousReadId]) - errorPos;
	}

    unsigned numCorrections = 0;

    TValue currentPos = firstCorrectionForRead[erroneousReadId];
    while (currentPos != std::numeric_limits<TValue>::max())
    {
        if (correctionList[currentPos].errorPos == errorPos)
        {
            numCorrections = correctionList[currentPos].foundCorrections;
            break;
        }
        currentPos = correctionList[currentPos].nextCorrection;
    }
    return numCorrections;
}

/*
the add new linked list item append a new entry to the correctionList and 
takes care that the corresponding entries in firstCorrectionForRead are updated	
*/
template <typename TCorrection,typename TValue, typename TId1, typename TId2, typename TPos1,typename TPos2,typename TOverlap, typename TReadStore, typename TAlphabet>
inline void addCorrectionEntry(String<TCorrection> &correctionList,
                               String<TValue> &firstCorrectionForRead,
                               TId1 erroneousReadId,
                               TId2 correctReadId,
                               TPos1 correctPos,
                               TPos2 errorPos,
                               TOverlap overlap,
                               bool strand,
                               signed char indelLength,
                               //	TValueLength readLength,
                               TReadStore & store,
                               TAlphabet & correctSeq)
{
    //get Alphabet type from String here
    //Dna5 correctSeq[MAX_INDEL_LENGTH];
    //std::cerr <<abs( (int)indelLength) << " lll " << MAX_INDEL_LENGTH<<std::endl;
    // get the correction sequence
    //if(indelLength <= 0) //only get string if insertion in read or mismatch
    //	getCorrectionString(correctSeq,indelLength,correctReadId,correctPos,strand,store);
    /*	//DEBUG
     unsigned length2=abs((int)indelLength);
     if(indelLength ==0){
     length2=1;
     }
     std::cerr << "In add CorrectionEntry retrieved CorrectionSequence:";
     for(unsigned i =0;i<length2;++i)
     std::cerr << correctSeq[i];
     std::cerr << "errorPos,strand" <<errorPos<<strand<< std::endl;
     //DEBUG */
    //change the erroneousReadId if necessary
    unsigned numReads = length(store.readSeqStore) / 2;

    // switch to reverse-complements for strand
    if (erroneousReadId >= numReads)
        erroneousReadId -= numReads;
    //std::cerr << " erroneousReadID and dummy  " <<erroneousReadId<< ","<< dummyErroneousReadId<< std::endl;
    //TValue empty=maxValue(erroneousReadId);
    //first check if a Correction for erroneousReadId exists already

    TValue insertLinkAt = std::numeric_limits<TValue>::max();
    TValue currentPos = firstCorrectionForRead[erroneousReadId];
    while (currentPos != std::numeric_limits<TValue>::max())
    {
        TCorrection &corr = correctionList[currentPos];
        if (existsCorrectionAtPos(
                corr.errorPos,
                corr.indelLength,
                corr.correctSeq,
                errorPos,
                strand,
                indelLength,
                length(store.readSeqStore[erroneousReadId]),
                correctSeq))
        {
#ifndef FIONA_NOERROROPTIMIZATION
            updateCorrectionEntry(
                correctionList,
                currentPos,
                overlap,
                strand,
                correctReadId,
                indelLength,
                correctSeq,
                correctionList[currentPos].indelLength);
#endif
            return;
        }
        insertLinkAt = currentPos;
        currentPos = correctionList[currentPos].nextCorrection;
    }
    
    TCorrection newCorrection;
    fillCorrection(
        newCorrection,
        correctPos,
        errorPos,
        correctSeq,
        overlap,
        strand,
        indelLength,
        length(store.readSeqStore[erroneousReadId]),
        correctReadId);

    if (length(correctionList) == capacity(correctionList))
    {
        lockWriting(correctionListLock);
        reserve(correctionList, capacity(correctionList) + 1, Generous());
        unlockWriting(correctionListLock);
    }
    appendValue(correctionList, newCorrection);

    if (insertLinkAt == std::numeric_limits<TValue>::max())
        firstCorrectionForRead[erroneousReadId] = length(correctionList) - 1;
    else
        correctionList[insertLinkAt].nextCorrection = length(correctionList) - 1;

/*
    if(firstCorrectionForRead[dummyErroneousReadId] == maxValue<TValue>())
    {	//we are going to add the first Correction element for the read at the end later
        firstCorrectionForRead[dummyErroneousReadId] = (TValue) length(correctionList);
    }
    else
    {
        TValue currentPos = firstCorrectionForRead[dummyErroneousReadId];
        bool nextElem =true;
        while(nextElem)
        {
			if(existsCorrectionAtPos(correctionList[currentPos].errorPos,correctionList[currentPos].indelLength,correctionList[currentPos].correctSeq, errorPos, strand, indelLength, length(store.readSeqStore[erroneousReadId]),correctSeq) )
            {
				// add the new information if the overlap sum is bigger decided by addCorrectionInfo




#ifndef FIONA_NOERROROPTIMIZATION
				updateCorrectionEntry(correctionList,currentPos,overlap,strand,correctReadId,indelLength,correctSeq,correctionList[currentPos].indelLength);
#endif
				return;
			}
			if (correctionList[currentPos].nextCorrection != maxValue<TValue>())
            {
				currentPos = correctionList[currentPos].nextCorrection;
			}
            else
            {
				nextElem =false;
				//we have not found any similar correction item so we extend the last elem
				// in the linked list
                correctionList[currentPos].nextCorrection = (TValue) length(correctionList);
			}
        }
    }


// if we land here than either there was never an entry or no Correction entry
// of the same pos or error type. Therefore we
// create new Correction struct and add that to the list.

    //std::cerr<< "Put new correction:"<<correctReadId<<","<<correctPos<<","<<errorPos<<","<<overlap<<","<<strand<<","<<(int)indelLength<<std::endl;
    TCorrection newCorrection;
    fillCorrection(
        newCorrection,
        correctPos,
        errorPos,
        correctSeq,
        overlap,
        strand,
        indelLength,
        length(store.readSeqStore[erroneousReadId]),
        correctReadId);
    
    appendValue(correctionList,newCorrection);
*/
    return;
}

   template <typename TCorrection,typename TValue,typename TReadStore>
        inline void  _testCorrectionStruct(String<TCorrection> &correctionList, String<TValue> &firstCorrectionForRead,TReadStore &store)
	{
		// we assume we work with three reads here
		TValue empty = std::numeric_limits<TValue>::max();
		appendValue(firstCorrectionForRead,empty);
		appendValue(firstCorrectionForRead,empty);
		appendValue(firstCorrectionForRead,empty);
		signed char indelLength = 0;
		unsigned short  correctPos=12;
		unsigned short errorPos =24;
		unsigned short overlap =5;
		bool strand = false;
		unsigned int correctReadId=4;
		unsigned int erroneousReadId=0;
		unsigned int readLength =36; 
		TCorrection tester;
		Dna5 correctSeq[MAX_INDEL_LENGTH];
		correctSeq[0]='A';
		correctSeq[1]='T';
		String <int> Ovsumcutoffs;
		clear(Ovsumcutoffs);
		resize(Ovsumcutoffs, readLength+1, 0);
		// test to fill a Corrections we will work with all the time
		fillCorrection(tester,correctPos,errorPos,correctSeq,overlap,strand,indelLength,readLength,correctReadId); //mismatch
		// add correction to the overall list
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[0],strand,tester.indelLength,store,correctSeq);
		
		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		// create same entry with higher overlap sum
		tester.overlap[0]=12;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[0],strand,tester.indelLength,store,correctSeq);
		
		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		//test adding same correction on different strand
		strand = true;
		erroneousReadId = 6; //because of reverse strand
		tester.errorPos=readLength-24-1;
		correctPos = readLength-12-1;
		tester.overlap[1]=8;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);
		
		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		//test adding new correction with indel at same position
		strand = false;
		erroneousReadId = 0; //because of forward strand
		correctPos=12;
		tester.errorPos=24;
		tester.overlap[0]=11;
		tester.indelLength=-2;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[0],strand,tester.indelLength,store,correctSeq);
		
		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		//test adding new correction with indel at same position but on opposite strand
		strand = true;
		erroneousReadId = 6; //because of reverse strand
		correctPos= readLength - 12 -2;
		tester.errorPos = readLength - tester.errorPos - 2;
		tester.overlap[1]=13;
		tester.indelLength=-2;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);

		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		//test adding new correction with insertion at same position
		strand = false;
		erroneousReadId = 0; //because of forward strand
		correctPos = 12;
		tester.errorPos = 24;
		tester.overlap[0]=3;
		tester.indelLength=2;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[0],strand,tester.indelLength,store,correctSeq);

		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		//test adding new correction with insertion at different position but on opposite strand
		strand = true;
		erroneousReadId = 6; //because of reverse strand
		correctPos = readLength - correctPos;
		tester.errorPos = readLength - tester.errorPos;
		tester.overlap[1]=7;
		tester.indelLength=2;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);

		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		//test adding new correction with insertion at same position but different correction sequence
		strand = false;
		erroneousReadId = 0; //because of forward strand
		correctPos = 2;
		tester.errorPos = 24;
		tester.overlap[1]=10;
		tester.indelLength=2;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);

		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	

		//add insertion Correction on reverse strand for new readID
		erroneousReadId=2;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);

		//_dumpCorrectionList(correctionList, firstCorrectionForRead);	
		//add deletion Correction on reverse strand for new readID
		tester.errorPos = readLength - 24 - 2;
		tester.indelLength=-2;
		tester.overlap[1]=15;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);


		//adding another mismatch for non-conficting indel correction
		erroneousReadId =0;
		strand = false;
		correctPos = 13;
		tester.errorPos = 13;
		tester.overlap[1]=97;
		tester.indelLength=0;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);

		//adding yet another mismatch to test if several corrections are done for non-conficting indel correction
		erroneousReadId =0;
		strand = false;
		correctPos = 5;
		tester.errorPos = 7;
		tester.overlap[1]=90;
		tester.indelLength=0;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);
		
		//adding another deletion for testing of output
		erroneousReadId =1;
		strand = false;
		correctPos = 12;
		tester.errorPos = 11;
		tester.overlap[1]=20;
		tester.indelLength=2;
		addCorrectionEntry(correctionList,firstCorrectionForRead,erroneousReadId,correctReadId,correctPos,tester.errorPos,tester.overlap[1],strand,tester.indelLength,store,correctSeq);
        _dumpCorrectionList(correctionList, firstCorrectionForRead,store);	

		//test non Conflicitng Correction selection
		//add several other corrections 
		FionaOptions options;

//unsigned numberCorrections = applyReadErrorCorrections(correctionList,firstCorrectionForRead,store, Ovsumcutoffs,options);
//		std::cout << "did "<< numberCorrections << " many corrections."<<std::endl;
		return;
	}
/*
	function that goes through all reads and prints there Corrections listed 
*/
template <typename  TCorrection,typename TValue, typename TStore>
inline void _dumpCorrectionList(
	String<TCorrection> const &correctionList,
	String<TValue> const &firstCorrectionForRead,
	TStore & store)
{
	// go through all reads and show which Corrections are listed
    for (unsigned int i =0;i< length(firstCorrectionForRead);i++)
    {
        std::cerr << "Found "<<length(correctionList)<<" corrections. Look at readID: " <<i<<std::endl;
        if (firstCorrectionForRead[i] != std::numeric_limits<TValue>::max())
        {
            std::cerr << "found Correction for read "<<i<< " at pos (in String<corrections>) "<< firstCorrectionForRead[i] << std::endl;
            _dumpCorrectionIndelPos(correctionList[firstCorrectionForRead[i]],i,store);
            TValue next = correctionList[firstCorrectionForRead[i]].nextCorrection;
            while (next != std::numeric_limits<TValue>::max())
            {
                std::cerr << "found Correction for read "<<i<< " at pos (in String<corrections>) "<< next << std::endl;
                _dumpCorrectionIndelPos(correctionList[next],i,store);
                next = correctionList[next].nextCorrection;
            }
        }	
    }
}
template <typename TCorrection, typename TStore>
inline void _dumpCorrectionIndelPos(
	TCorrection const &correction,
	unsigned errorReadId,
	TStore& store)
{
	std::cerr << "error___read_id\t" << errorReadId << std::endl;
	std::cerr << "error_pos      \t" << correction.errorPos << std::endl;
	std::cerr << "next correction\t" << correction.nextCorrection << std::endl;
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
	std::cerr << "correct_read_id\t" << correction.correctReadId << std::endl;
	std::cerr << "correct_pos    \t" << correction.correctPos << std::endl;
#endif
	std::cerr << "overlap        \t F:" << correction.overlap[0] << " R: " << correction.overlap[1] << std::endl;
	std::cerr << "indel_len      \t" << (int)correction.indelLength << std::endl;
	unsigned length=abs((int)correction.indelLength);
	if(correction.indelLength ==0){
		length=1;
	}
	std::cerr << "CorrectionSequence:";
	for(unsigned i =0;i<length;++i)
		std::cerr << correction.correctSeq[i];
	std::cerr << std::endl;
	std::cerr << "error___read   \t" << store.readSeqStore[errorReadId][correction.errorPos] << '\t';
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
         for (unsigned i = 0; i < correction.correctPos; ++i)
                 std::cerr << ' ';
#endif
     std::cerr << store.readSeqStore[errorReadId] << std::endl;
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
         std::cerr << "correct_read   \t" << store.readSeqStore[correction.correctReadId][correction.correctPos] << '\t';
         for (unsigned i = 0; i < correction.errorPos; ++i)
                 std::cerr << ' ';
     std::cerr << store.readSeqStore[correction.correctReadId] << std::endl;
#endif
}

/*matching string*/
inline bool strContains(std::string const & inputStr, std::string const & searchStr)
{
	return inputStr.find(searchStr) != std::string::npos;
}


template <typename TCorrection>
struct LessOverlap : public std::binary_function<TCorrection, TCorrection, bool >
{
    inline bool operator() (TCorrection const &a, TCorrection const &b) const
    {
        // sort by overlap
        if (a.overlap[0] != b.overlap[0])
            return a.overlap[0] > b.overlap[0];

        // if overlap equal prefer mismatches
        // if no mismatch prefer deletion
        if (a.indelLength == 0 || b.indelLength == 0)
            return a.indelLength == 0 && b.indelLength != 0;    // return (abs((int)a.indelLength) < abs((int)b.indelLength));
        else
            return a.indelLength > b.indelLength;
    }
};

template <typename TCorrection>
struct LessPositionOverlap : public std::binary_function<TCorrection, TCorrection, bool >
{
    inline bool operator() (TCorrection const &a, TCorrection const &b) const
    {
        // check difference in position
        if (a.errorPos != b.errorPos)
            return a.errorPos > b.errorPos;

        // if position is the same sort by overlap
        if (a.overlap[0] != b.overlap[0])
            return a.overlap[0] > b.overlap[0];

        // if overlap equal prefer mismatches
        // if no mismatch prefer deletion
        if (a.indelLength == 0 || b.indelLength == 0)
            return a.indelLength == 0 && b.indelLength != 0;    // return (abs((int)a.indelLength) < abs((int)b.indelLength));
        else
            return a.indelLength > b.indelLength;
    }
};




template<typename TCorrection,typename TValue,typename TReadStore>
inline unsigned applyReadErrorCorrections(String<TCorrection> const &correctionList,
                                          String<TValue> &firstCorrectionForRead,
                                          TReadStore & store,
                                          FionaOptions & options)
{
    if (options.verbosity >= 1)
    {
        std::cerr << "Find non conflicting Corrections\n"
                  << "Length of the linked list is " << length(correctionList) << "\n"
                  << "  => length " << length(correctionList) * sizeof(TCorrection) << " bytes\n"
                  << "Capacity of the linked list is " << capacity(correctionList) << "\n"
                  << "  => capacity " << capacity(correctionList) * sizeof(TCorrection) << " bytes\n";
    }

    unsigned numberCorrections = 0;
    String<TCorrection> possibleCorrections;
    String<unsigned int> correctionsToSave;

    SEQAN_OMP_PRAGMA(parallel for schedule(guided) private(correctionsToSave, possibleCorrections) reduction(+ : numberCorrections))
    for (int readId = 0; readId < (int)length(firstCorrectionForRead); ++readId)
    {
        // descend if correction exists
        TValue corrId = firstCorrectionForRead[readId];
        if (corrId == std::numeric_limits<TValue>::max())
            continue;

        // step through linked list and collect corrections
        clear(possibleCorrections);
        do
        {
            appendValue(possibleCorrections, correctionList[corrId]);
            TCorrection & correction = back(possibleCorrections);

            //add the overlap values for both strand for the new object in the string
#ifdef FIONA_MAXIMIZE_OVERLAPSUM
            correction.overlap[0] = _max(correction.overlap[0], correction.overlap[1]);
#else
            correction.overlap[0] += correction.overlap[1];
#endif
            //_dumpCorrectionIndelPos(correction, readId);
            corrId = correction.nextCorrection;
        } while (corrId != std::numeric_limits<TValue>::max());

        //sorting by Position first to get the best correction per Position
        //sorting is done arbitrarily from large to small(right to left)

        sort(possibleCorrections, LessPositionOverlap<CorrectionIndelPos>(), Parallel());

#ifndef FIONA_NO_SEPARATE_OVERLAPSUM
	//only remove if several corrections per position are saved

        if (readId == options.debugRead)
        {
            std::cerr << "output after sorting for positions and removal of positions"<<std::endl;
            for (unsigned j=0;j<length(possibleCorrections);++j)
            {
                std::cerr << "Output for read "<<readId<<"entry "<<j<<std::endl;
                _dumpCorrectionIndelPos(possibleCorrections[j],readId,store);
            }
        }  //output for read

        //go through all Correction struct and keep only the first per Position
        if (length(possibleCorrections) > 1)
        {
            unsigned lastErrorPos = possibleCorrections[0].errorPos;
            unsigned differentPos = 1;
            //std::cerr <<"length if corectios: "<<length(possibleCorrections)<<std::endl;
            for (unsigned j = 1; j < length(possibleCorrections); ++j)
            {
                // if position decreases means we found a new error position
                if (possibleCorrections[j].errorPos < lastErrorPos)
                {
                    lastErrorPos = possibleCorrections[j].errorPos;
                    //if we skipped Correction entries that had been on the same position
                    //move the Correction with new position forward to differentPos
                    if (j != differentPos)
                        possibleCorrections[differentPos] = possibleCorrections[j];
                    ++differentPos;
                }
            }
            //finally resize to the length of different Positions found
            resize(possibleCorrections, differentPos);
        }//finished the cases with more than one Correction

#endif //FIONA_NO_SEPARATE_OVERLAPSUM

        /*/output for debug
         std::cerr << "output after sorting for positions and removal of positions"<<std::endl;
         for(unsigned j=0;j<length(possibleCorrections);++j){
         std::cerr << "entry "<<j<<std::endl;
         _dumpCorrectionIndelPos(possibleCorrections[j],readId);

         }*/

        //sorting by overlap sum now
#ifndef FIONA_NOERROROPTIMIZATION    //dont sort in random encounter (local) mode
        sort(possibleCorrections, LessOverlap<CorrectionIndelPos>(), Parallel());
#endif
        //go through all Correction struct and keep the ones with highest overlapsum
        // and without conflict in terms of error type
        // execute all correctins that are mismatches and have a higher overlapsum than the best indel correction
        //create a new array that saves which correction should be used to accomodate letter N corrections
        clear(correctionsToSave);
        appendValue(correctionsToSave, 0);

#ifndef FIONA_NOERROROPTIMIZATION

        unsigned errorReadLength = length(store.readSeqStore[readId]);

#ifndef FIONA_DISTANCE_BASED_ERROR_OPTIMIZATION
        // if indel use just the first correction except there are Ns
        bool notFoundCorrectionLimit = (possibleCorrections[0].indelLength == 0);
	
	if(notFoundCorrectionLimit){  //if indel dont do anything more
        for (unsigned j = 1; j < length(possibleCorrections); ++j)
        {
            if (possibleCorrections[j].indelLength != 0)
                continue;

            // Todo (Hugues) integrate varying k when available.
            if (notFoundCorrectionLimit
                && possibleCorrections[j].overlap[0] > options.overlapSumCutoffs(errorReadLength, possibleCorrections[j].errorPos))
            {
                // if not too many corrections are found do the correction and proceed
                if (length(correctionsToSave) < (unsigned) options.allowedCorrectionsPerRead[readId])
                {
                    appendValue(correctionsToSave, j);
                    continue;
                }

                // dont pick any further corrections
                notFoundCorrectionLimit = false;
            }

            // special test to check for letter N corrections which where not detected above
            if (ordValue(store.readSeqStore[readId][possibleCorrections[j].errorPos]) == 4)
                appendValue(correctionsToSave, j);
        }
        }//endif

#else //new mode that uses distances between corrections 
    unsigned int next = 1;
	while (next < length(possibleCorrections))
	{
        if (possibleCorrections[next].overlap[0] <= options.overlapSumCutoffs(errorReadLength, possibleCorrections[next].errorPos))
        {
            ++next;  // skip this correction
            continue;
        }

        bool add = true;
		for(unsigned int i = 0; i < length(correctionsToSave); i++){
			//check if the distance between correction positions is apart min_k to each correction accepted
			if(abs(possibleCorrections[correctionsToSave[i]].errorPos-possibleCorrections[next].errorPos) <= options.fromLevel){
				add=false;
				break;
			}
		}
		if (add)
			appendValue(correctionsToSave, next);

        next++;  //check the next possible correction in the list
	}
#endif


#else 
        //do all corrections if no global mode
        resize(correctionsToSave, length(possibleCorrections));
        for (unsigned j = 0; j < length(possibleCorrections); ++j)
            correctionsToSave[j] = j;
#endif

        //decrease the number of allowed Corrections
        options.allowedCorrectionsPerRead[readId] -= length(correctionsToSave);

        // apply the corrections this part is very similar to the original end of the function applyReadErrorCorrections
        // in the first Fiona version but that we dont need to copy the original reads as we have saved the correction string in the
        // CorrectionIndelPos struct.
        numberCorrections += length(correctionsToSave);

        unsigned correction;
        //for (unsigned int correction =0;correction < numberOfValidCorrections; ++correction)
        for (unsigned count = 0; count < length(correctionsToSave); ++count)
        {
            correction = correctionsToSave[count]; //set the position in the CorrectionString we are executing now
//	    SEQAN_OMP_PRAGMA(critical) {	
//	    std::cout << "found " << (unsigned ) getFoundCorrections(correctionList,firstCorrectionForRead,
//                                 (unsigned) readId, possibleCorrections[correction].errorPos, false, store)<< "corrections read:" << readId<< " pos:" << possibleCorrections[correction].errorPos << std::endl;
//	    }

            if ((unsigned)readId == (unsigned) options.debugRead)
                //if(readId == 18473)
            {
                std::cerr << "BEFORE:" << std::endl;
                _dumpCorrectionIndelPos(possibleCorrections[correction], readId,store);
            }
            //SEQAN_OMP_PRAGMA(critical) {
            //    _dumpCorrectionsIndelPos(possibleCorrections[]
            //}

            std::ostringstream m;
            if (options.appendCorrectionInfo)
            {
                // Start the read correction information with the common prefix " corrected:\t" if it is empty.
                if (empty(store.readNameStore[readId]))
                    m << " corrected:\t";
                else
                    m << "\t";
            }

            if (possibleCorrections[correction].indelLength == 0)
            {
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
                m << possibleCorrections[correction].errorPos << "(" << options.cycle <<',' <<possibleCorrections[correction].overlap[0] << ","<<possibleCorrections[correction].correctReadId<<"):" << store.readSeqStore[readId][possibleCorrections[correction].errorPos] << "->" << possibleCorrections[correction].correctSeq[0];
#else
                m << possibleCorrections[correction].errorPos << "(" << options.cycle <<',' <<possibleCorrections[correction].overlap[0] << ",???"<<"):" << store.readSeqStore[readId][possibleCorrections[correction].errorPos] << "->" << possibleCorrections[correction].correctSeq[0];
#endif
                store.readSeqStore[readId][possibleCorrections[correction].errorPos] = possibleCorrections[correction].correctSeq[0];
            }
#ifdef FIONA_ALLOWINDELS
            else if (possibleCorrections[correction].indelLength > 0)
            {
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
                m << possibleCorrections[correction].errorPos << "(" << options.cycle <<',' <<possibleCorrections[correction].overlap[0] << ","<<possibleCorrections[correction].correctReadId<<"):-" << infix(store.readSeqStore[readId], possibleCorrections[correction].errorPos, possibleCorrections[correction].errorPos + possibleCorrections[correction].indelLength);
#else
                m << possibleCorrections[correction].errorPos << "(" << options.cycle <<',' <<possibleCorrections[correction].overlap[0] << ",???"<<"):-" << infix(store.readSeqStore[readId], possibleCorrections[correction].errorPos, possibleCorrections[correction].errorPos + possibleCorrections[correction].indelLength);
#endif
                erase(store.readSeqStore[readId], possibleCorrections[correction].errorPos, possibleCorrections[correction].errorPos + possibleCorrections[correction].indelLength);
            } else {
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
                m << possibleCorrections[correction].errorPos << "(" << options.cycle <<',' <<possibleCorrections[correction].overlap[0] <<","<<possibleCorrections[correction].correctReadId<<"):+";
#else
                m << possibleCorrections[correction].errorPos << "(" << options.cycle <<',' <<possibleCorrections[correction].overlap[0] <<",??"<<"):+";
#endif
                String<Dna5> proxy;  //TODO: construct is ugly
                resize(proxy,abs(possibleCorrections[correction].indelLength),Exact());
                proxy=infix(possibleCorrections[correction].correctSeq,0,(unsigned) abs(possibleCorrections[correction].indelLength));
                m << proxy;
                insert(store.readSeqStore[readId], possibleCorrections[correction].errorPos,proxy);
            }
#endif
#ifdef FIONA_DISTANCE_BASED_ERROR_OPTIMIZATION
            //correct offset of subsequent error positions downstream of the current pos after indel corrections
            if(possibleCorrections[correction].indelLength != 0){
		//SEQAN_OMP_PRAGMA(critical) {     
		//std::cerr << "Did (read/length) " << readId << "/"<<length(store.readSeqStore[readId])<< " ) e-pos/indel " << possibleCorrections[correction].errorPos << "/" << (int)possibleCorrections[correction].indelLength << std::endl;
		//}
                  for(unsigned nextcount = count+1; nextcount < length(correctionsToSave); ++nextcount)
                  {
             	      unsigned nextcorrection = correctionsToSave[nextcount];
                      if(possibleCorrections[nextcorrection].errorPos > possibleCorrections[correction].errorPos){
                          possibleCorrections[nextcorrection].errorPos -= (int) possibleCorrections[correction].indelLength;  
		//SEQAN_OMP_PRAGMA(critical) {     std::cerr << readId << " changed occurrence " << " ) e-pos/indel " << possibleCorrections[nextcorrection].errorPos << "/" << (int)possibleCorrections[nextcorrection].indelLength << std::endl; }
                        }
                  }
            }
#endif  //DISTANCE_BASED_ERROR_OPTIMIZATION

            // Extend the read correction information if configured to do so.
            if (options.appendCorrectionInfo)
                append(store.readNameStore[readId], m.str());
        }
    }
    //give her the number of corrections
    return numberCorrections;
}




///Functions interpretation of results
//Proba of a read with no error
template <typename TErrorRate, typename TPrefixLen>
inline double probabilityNoError(TErrorRate perrorrate, TPrefixLen k)
{
	return pow(1.0 - perrorrate, (double) k);
}

// factorial
template <typename TValue>
inline double factorial(TValue n)
{
	double fact = 1.0;
	for (TValue i = 2; i <= n; ++i)
		fact *= i;
	return fact;
}

#ifdef CAN_BE_REMOVED
// Expected value - general, use if all reads have the same length
template <typename TExpectedValues, typename TReadSet, typename TGenomeSize>
void expectedValueEqualReadsLength(TExpectedValues & expected, TReadSet const & readSet, TGenomeSize const genomeLength)
{
	//
	//	E(m) = (read_length - suffix_length + 1) * numberReads / genomeLength
	//

	// without reverse complement
	unsigned readCount = length(readSet) / 2;
	unsigned readLength = length(readSet[0]);

	clear(expected);
	for (unsigned suffixLength = 0; suffixLength <= readLength; ++suffixLength)
		append(expected, (readLength - suffixLength + 1) * readCount / (double)genomeLength);
}
#endif

// Calculate a read length histogram
// Excluding reverse complements
template <typename TReadLengthHist, typename TReadSet>
void computeReadLengthHistogram(TReadLengthHist &readLenHist, TReadSet const &readSet)
{
    clear(readLenHist);
    unsigned numReads = length(readSet) / 2;
	for (unsigned i = 0; i < numReads; ++i)
	{
		unsigned readLength = length(readSet[i]);
		if (readLength >= length(readLenHist))
			resize(readLenHist, readLength + 1, 0);
		++readLenHist[readLength];
	}    
}

// Expected value for set of reads with different length
// Also gets the expected number of erroneous reads
template <typename TExpectedValues, typename TNumRead, typename TGenomeSize, typename TErrorrate>
double expectedValueTheoretical(TExpectedValues & expected, String<TNumRead> const &readLenHist, TGenomeSize const genomeLength, TErrorrate const errorrate)
{
	//
	//	E(m) = (8 - suffix_length + 1) * numberReads / genomeLength
	//
    	
	double nerrreads = 0.0;

	// a = read_length - suffix_length + 1
	clear(expected);
	resize(expected, length(readLenHist), 0.0);
	for (unsigned readLen = 1; readLen < length(readLenHist); ++readLen)
	{
        TNumRead numReads = readLenHist[readLen];
        if (numReads == 0)
            continue;

		if (errorrate != 0.0)
			nerrreads += numReads * (1 - probabilityNoError(errorrate, readLen));

		for (unsigned suffixLength = 0; suffixLength <= readLen; ++suffixLength)
		{
			double a = readLen - suffixLength + 1;
			expected[suffixLength] += a * (double)numReads / (double)genomeLength;
		}
	}
	return nerrreads;
}

/* Standard Deviation */
template <typename TDeviationValues, typename TReadSet, typename TGenomeSize>
void standardDeviation(TDeviationValues & deviation, TReadSet const & readSet, TGenomeSize const genomeLength)
{
	// without reverse complement
	unsigned readCount = length(readSet) / 2;
	unsigned readLength = length(readSet[0]);

	//
	//	SD(m)= numberReads*((read_length - suffix_length + 1)/genomeLength 
	//			- (read_length - suffix_length + 1)^2/genomeLength ^2)
	//

	double valueFirst;
	double valueSecond;
	resize(deviation, readLength + 1);
	for (unsigned suffixLength = 0; suffixLength <= readLength; ++suffixLength)
	{
		valueFirst  = (readLength - suffixLength + 1) / (double)genomeLength;
		valueSecond = valueFirst * valueFirst;
		deviation[suffixLength] = sqrt((valueFirst - valueSecond) * readCount);	
	}
}


//P(X = k) for X ~ Poisson
////  WARNING These computation of poisson values are not accurate 
////  for low probabilities !!! 
////  around 50% relative error at 1e-5 proba
template <typename TValue, typename TMean>
inline double dpois(TValue k, TMean mean)
{
	return pow(mean, (double)k) * exp(-mean) / factorial(k);
}

// cumulative poisson distribution
// P(X <= k) with X ~ Poisson(mean)
template <typename TValue, typename TMean>
inline double ppois(TValue k, TMean mean)
{
	// return gsl_cdf_poisson_P(k,mean);
	double negExp = exp(-mean);
	double pValue = 0.0;
	double pow = 1.0;
	double fact = 1.0;
	for (TValue i = 0; i <= k; ++i, fact *= i){
		pValue += pow * negExp / fact;
		pow *= mean;
	}
	return pValue;
}

// give the highest k such that P(X < c) <= alpha with X ~ Poisson(lambda)
// carefull it's exclusive (for compatibility with the cutoff mode)
template <typename TPValue, typename TMean>
inline unsigned qpois(TPValue alpha, TMean lambda){
	double negExp = exp(-lambda);
	double pValue = 0.0;
	double pow = 1.0;
	double fact = 1.0;
	unsigned k = 0;
	while (pValue <= alpha){
		pValue += pow * negExp / fact;
		k++;
		pow *= lambda;
		fact *= k;
	}
	return k;
}

//Estimate the two first terms of the error part into the mixture
//this is only valid for error rate sufficiently low (ie  < 5%..)
//
template <typename TValue, typename TMean, typename TErrorrates, typename TPrefixLen>
inline double dpoismixerror(TValue k, TMean lambda, TErrorrates errorrate, TPrefixLen prefixlen)
{
	double noerrlm2 = pow(1.0 - errorrate, prefixlen - 2.0);
	double noerrlm1 = noerrlm2 * (1-errorrate);
	double errexp1err = lambda * noerrlm1 * (errorrate/3);
	double errexp2err = lambda * noerrlm2 * (errorrate/3) * (errorrate/3);
//	double perr1 = ((double) prefixlen) * noerrlm1* errorrate ;
//	double perr2 = ((double) prefixlen * (prefixlen -1.0) / 2 ) * noerrlm2 * (errorrate) *(errorrate);
	double perr1 =  noerrlm1 ;
	double perr2 =  ((prefixlen -1.0) / 2 )* noerrlm2 * errorrate;
	double sc = perr1 + perr2;
	perr1 /= sc; 
	perr2 /= sc; 
	double proba;
	proba = perr1 * dpois(k, errexp1err) + perr2 * dpois(k, errexp2err);
	return(proba);
}
			 
// cdf for a poisson mixture of errors, given the 
// expected value on the real reads before sequencing
//// WARNING THIS IS extremely sensitive  !! 
///
template <typename TValue, typename TMean, typename TErrorrates, typename TPrefixLen>
inline double ppoismixerror(TValue k, TMean lambda, TErrorrates errorrate, TPrefixLen prefixlen)
{
	double noerrlm2 = pow(1.0 - errorrate, prefixlen - 2.0);
	double noerrlm1 = noerrlm2 * (1-errorrate);
	double errexp1err = lambda * noerrlm1 * (errorrate/3);
	double errexp2err = lambda * noerrlm2 * (errorrate/3) * (errorrate/3);
	double perr1 = prefixlen * noerrlm1* (errorrate) ;
	double perr2 = (prefixlen * (prefixlen -1) / 2 ) * noerrlm2 * (errorrate) *(errorrate);
	double sc = perr1 + perr2;
	perr1 /= (sc); 
	perr2 /= (sc); 
	double negExp1err = exp(-errexp1err); //
	double negExp2err = exp(-errexp2err);
	double pow1err = 1.0;
	double pow2err = 1.0;
	double fact = 1.0;
	double pValue = 0.0;
//	std::cerr << "Level " << prefixlen << " and params for correc\n";
//	std::cerr << "Expected and comp: " << lambda << " - " << errexp1err << " - " << errexp2err << "\n";
//	std::cerr << "the observed count: " << k << ", perr is " << perr1 << " - " << perr2  << "\n";
	TValue i = 0;
    for (i = 0; i <= k; ++i, fact *= i){
        pValue += perr1 * pow1err * negExp1err / fact;
		pValue += perr2 * pow2err * negExp2err / fact;
		pow1err *= errexp1err; 
		pow2err *= errexp2err;
//		std::cerr << "alpha =" << pValue << std::endl;
	}
	return pValue;
}

///Similar to qpois for the mixture of poisson distributions with error proportions
template <typename TPValue, typename TMean, typename TErrorrates, typename TPrefixLen>
inline unsigned qpoismixerror(TPValue alpha, TMean lambda, TErrorrates errorrate, TPrefixLen prefixlen)
{
	double noerrlm2 = pow(1.0 - errorrate, prefixlen - 2.0);
	double noerrlm1 = noerrlm2 * (1-errorrate);
	double errexp1err = lambda * noerrlm1 * (errorrate/3);
	double errexp2err = lambda * noerrlm2 * (errorrate/3) * (errorrate/3);
	double perr1 = prefixlen * noerrlm1* (errorrate) ;
	double perr2 = (prefixlen * (prefixlen -1) / 2 ) * noerrlm2 * (errorrate) *(errorrate);
	double sc = perr1 + perr2;
	perr1 /= (sc); 
	perr2 /= (sc); 
	double negExp1err = exp(-errexp1err); //
	double negExp2err = exp(-errexp2err);
	double pow1err = 1.0;
	double pow2err = 1.0;
	double fact = 1.0;
	double pValue = 0.0;
	unsigned k = 0;
	//std::cerr << "Level " << prefixlen << " and params for correc\n";
	//std::cerr << "Expected and comp: " << lambda << " - " << errexp1err << " - " << errexp2err << "\n";
	//std::cerr << "the observed count: " << k << ", perr is " << perr1 << " - " << perr2  << "\n";
	while (pValue <= alpha){
        pValue += perr1 * pow1err * negExp1err / fact;
		pValue += perr2 * pow2err * negExp2err / fact;
		k++;
		pow1err *= errexp1err; 
		pow2err *= errexp2err;
		//std::cerr << "alpha =" << pValue << std::endl;
		fact *= k;
	}
	return k;
}

//
//Cutoff value under a classification scheme according to sign( log {P(X = c| error) / P(X = c |no error)} + log(prior)) 
//prior = 0 does automatic adjustment of priors according to the number of expected errors
template <typename TOddsRatio, typename TMean, typename TErrorrates, typename TPrefixLen>
inline unsigned PoisClassifCutoff(TOddsRatio prior, TMean lambda, TErrorrates errorrate, TPrefixLen prefixlen){
	double noerr = probabilityNoError(errorrate, prefixlen);	
	double noerrlm2 = pow(1.0 - errorrate, prefixlen - 2.0);
	double noerrlm1 = noerrlm2 * (1-errorrate);
	double errexpnoerr = lambda * noerr;
	double errexp1err = lambda * noerrlm1 * (errorrate/3);
	double errexp2err = lambda * noerrlm2 * (errorrate/3) * (errorrate/3);
	double perr1 = prefixlen * noerrlm1* (errorrate) ;
	double perr2 = (prefixlen * (prefixlen -1) / 2 ) * noerrlm2 * (errorrate) *(errorrate);
	double sc = perr1 + perr2;
	perr1 /= (sc); 
	perr2 /= (sc); 
	double negExpnoerr = exp(-errexpnoerr);
	double negExp1err = exp(-errexp1err); //
	double negExp2err = exp(-errexp2err);
	double pow1err = 1.0;
	double pow2err = 1.0;
	double pownoerr = 1.0;
	double fact = 1.0;
	double Ppostnoerr = 1.0;
	double Pposterr = 0.0;
	if (prior == 0){
		prior = (1-noerr) / noerr;
	}else {
		prior = prior * (1-noerr) / noerr; 
	}

	unsigned k = 0;
	//std::cerr << "prefix of len " << prefixlen << " expected " << lambda << " prior " << prior << " error-rate " <<errorrate << std::endl;
	//we impose that k shall not be higher than the expected count for correct reads.
	unsigned kquart = (unsigned)round(errexpnoerr);
	
	//while ((log((Pposterr / (Ppostnoerr)) * prior) > 0) || k == 0){
	while ((k < kquart && log(Pposterr / Ppostnoerr * prior) > 0) || k == 0){
		//not as efficient
		//Ppostnoerr = dpois(k, lambda*noerr);
		//Pposterr = dpoismixerror( k, lambda, errorrate, prefixlen);
		//std::cerr << "(" << k << " , "  << Ppostnoerr << " , " << Pposterr << ")" << std::endl;
		Ppostnoerr = pownoerr * negExpnoerr / fact;
        Pposterr = (perr1 * pow1err * negExp1err + perr2 * pow2err * negExp2err) / fact;
		pownoerr *= errexpnoerr;
		pow1err *= errexp1err; 
		pow2err *= errexp2err;
		k++;
		fact *= k;
	}
	//std::cerr << "(" << k << " , "  << Ppostnoerr << " , " << Pposterr << ")" << std::endl;
	return k;
	
}

// The probability distribution of a repeat, given the probability of the word
// !!! We assume simply non overlapping word and do a poisson approximation
inline double drepeat(int nrep, double pword, long genomelength){
	return dpois(nrep, pword * genomelength);
}


//Compute the odds of having a repeat for all counts from kmin to kmax 
//given an (iid M00) genome, a prefix length, the errorrate and the expected coverage
template <typename TPosterior, typename TCount, typename TMean, typename TErrorrates , typename TPrefixLen, typename TGenomeLen>
void OddsRepeat(TPosterior & post1occ, TCount cmin, TCount cmax, TMean lambda, TErrorrates errorrate, TPrefixLen prefixlen, TGenomeLen genomelength){
	double pword = 1.0 / pow(4.0, (double)prefixlen);
	int nrmax = 10;
	//prior on having a repeat for a random genome, not used now.
	String <double> prepeats ;
	resize(prepeats, nrmax+1, 0.0);
	double pnoerr = probabilityNoError(errorrate, prefixlen);
	for (int i = 0; i <= nrmax; i++) {
		prepeats[i] = drepeat(i, pword, genomelength);
	}
	String <double> posteriors;
	clear(post1occ);
	resize(post1occ, cmax, 0.0);
	for (int c = cmin; c < cmax; c++) {
		clear(posteriors);
		resize( posteriors, nrmax+1, 0.0);
		double sum = 0.0;
		posteriors[0] = 0.0;
		for (int nr=1; nr <= nrmax; nr++) {
			posteriors[nr] = ((pnoerr) * ppois(c, lambda*nr*pnoerr) + (1-pnoerr)*ppoismixerror(c, lambda*nr, errorrate, prefixlen)) * prepeats[nr];
			sum += posteriors[nr]; 
		}
		//rescale
		for (int nr=1; nr <= nrmax; nr++){
			posteriors[nr] /= sum;
		}
		post1occ[c] = (1-posteriors[1])/ posteriors[1] ;
	}
	
}

//The same methods, but computes the cutoff k where the odds of repeat over non repeat is > odds
template <typename TOdds, typename TMean, typename TErrorrate , typename TPrefixLen, typename TGenomeLen>
inline int OddsRepeatCutoff(
		TOdds  odds, 
		TMean lambda, 
		TErrorrate errorrate,
		TPrefixLen prefixlen, 
		TGenomeLen genomelength){
	double pword = 1.0 / pow(4.0, (double)prefixlen);
	int nrmax = 10;
	String <double> posteriors;
	clear(posteriors);
	resize( posteriors, nrmax+1, 0.0);
	//prior on having a repeat for a random genome, not used for now.
	String <double> prepeats ;
	resize(prepeats, nrmax+1, 0.0);
	double pnoerr = probabilityNoError(errorrate, prefixlen);
	for (int i = 0; i <= nrmax; i++) {
		prepeats[i] = drepeat(i, pword, genomelength);
	}
	//
	double post1occ = 0.0;
	int c = (odds < 1) ? 0 : (int)(lambda * pnoerr);
	int cmax = (int)(100.0 * lambda * pnoerr);
	while(post1occ < odds && c <cmax){
		clear(posteriors);
		resize( posteriors, nrmax+1, 0.0);
		double sum = 0.0;
		posteriors[0] = 0.0;
		for (int nr=1; nr <= nrmax; nr++) {
			posteriors[nr] = ((pnoerr) * ppois(c, lambda*nr*pnoerr) + (1-pnoerr)*ppoismixerror(c, lambda*nr, errorrate, prefixlen)) * prepeats[nr];
			sum += posteriors[nr]; 
		}
	//rescale
		for (int nr=1; nr <= nrmax; nr++){
			posteriors[nr] /= sum;
		}
		post1occ = (1-posteriors[1])/ posteriors[1] ;
		c++;
	}
	if (c == cmax) return 0;
	return c;
}
	
template <typename TCuttoffs, typename TOdds, typename TExpectedValues, typename TErrorrate , typename TPrefixLen, typename TGenomeLen>
void computeCutoffRepeats(
	TCuttoffs & thresholds,
	TOdds const odds,
	TExpectedValues const & expected,
	TErrorrate errorrate,
	TPrefixLen kmin,
	TPrefixLen kmax,
	TGenomeLen genomelength)
{
	clear(thresholds);
	resize(thresholds, kmax + 1, std::numeric_limits<typename Value<TCuttoffs>::Type>::max());
	for (int i = kmin; i <= kmax; i++)
		thresholds[i] = OddsRepeatCutoff(odds, expected[i], errorrate, i, genomelength); // Dave: I added "/ 3", otherwise this cutoff seems to have no effect
}



//Compute the number of ways of placing at most e errors in a read of length l such that
//any interval of length k has at least one error
inline void
precomputeOverlapCombinatorics(
    String<double> &expectedCorrectOverlapSum,
    String<double> &expectedFPOverlapSum,
    int maxNonSeedOverlap,
    int k,
	FionaOptions & options)
{
    clear(expectedCorrectOverlapSum);
    clear(expectedFPOverlapSum);
	resize(expectedCorrectOverlapSum, maxNonSeedOverlap + 1);
	resize(expectedFPOverlapSum, maxNonSeedOverlap + 1);

	SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1))
    for (int nonSeedOverlap = maxNonSeedOverlap; nonSeedOverlap >= 0; --nonSeedOverlap)
    {
        binomial Zerr(nonSeedOverlap, options.errorrate);
        binomial Zotherov(nonSeedOverlap, (1 + options.errorrate) / 4);

		unsigned maxErrors = (unsigned)(options.overlap_errorrate * (nonSeedOverlap + k + 1));

        double expovsumT1 = 0;
        double expovsumT2 = 0;
        
        for (unsigned numError = 0; numError < maxErrors; ++numError)
        {
            expovsumT1 += (k + 1 + nonSeedOverlap - numError) * pdf(Zerr, numError);
            expovsumT2 += (k + 1 + numError)                  * pdf(Zotherov, numError);
        }
        
        expectedCorrectOverlapSum[nonSeedOverlap] = expovsumT1;
        expectedFPOverlapSum[nonSeedOverlap] = expovsumT2;
    }
}

//
// Model for putting a threshold on the OverlapSum 
// We compute the expected overlapsum count for the real and the erroneous reads
//
template <typename TMean, typename TErrorrate, typename TWeight, typename TReadLengthHist>
inline int OddsOverlapSumCutoff(
	int i,          // position of the error
	TMean lambda,   // expected coverage per position
	int lenError,   // error read length
	int k,          // seed length
	TErrorrate errorrate,
	TWeight w,
    String<double> & expectedCorrectOverlapSum,
    String<double> & expectedFPOverlapSum,
    TReadLengthHist & readLenHist,
	FionaOptions & /*options*/)
{
	double pnoerr = probabilityNoError(errorrate, k+1); 
	double potherpos = probabilityNoError(errorrate, k) * (3./4); 
	double pword = 1.0 / pow(4.0, (double)k);
    
	double totalExpectedCorrectOverlapSum = 0.0;    // expected value for bona fide reads
	double totalExpectedFPOverlapSum = 0.0;         // expected value for reads with errors
	//double varovsum = 0; //we would like also to compute the variance (later).
	//count the forward strand overlaps

    for (int lenCorrect = 1; lenCorrect < (int)length(readLenHist); ++lenCorrect)
    {
        unsigned numReads = readLenHist[lenCorrect];
        if (numReads == 0)
            continue;

        double localExpectedCorrectOverlapSum = 0;
        double localExpectedFPOverlapSum = 0;
        int stepSize = 1 + (lenCorrect / 150);

        // seed is left of the error
        if (i >= k)
        {
            // j is the start position of the seed
            for (int j = 0; j < lenCorrect - k; j += stepSize)
            {
                //                   |left overlap|   |========== right overlap =======================|
                int nonSeedOverlap = _min(i - k, j) + _min(lenError - (i + 1), lenCorrect - (j + k + 1));

                localExpectedCorrectOverlapSum += stepSize * expectedCorrectOverlapSum[nonSeedOverlap];
                localExpectedFPOverlapSum += stepSize * expectedFPOverlapSum[nonSeedOverlap];
            }
        }

        // seed is right of the error
        if (i < lenError - k)
        {
            // j is the position of the "error" in the correct read
            for (int j = 0; j < (int)lenCorrect - k; j += stepSize)
            {
                //                   |left overlap|   |========== right overlap ===========================|
                int nonSeedOverlap = _min(i, j)     + _min(lenError - (i + k + 1), lenCorrect - (j + k + 1));
                
                localExpectedCorrectOverlapSum += stepSize * expectedCorrectOverlapSum[nonSeedOverlap];
                localExpectedFPOverlapSum += stepSize * expectedFPOverlapSum[nonSeedOverlap];
            }
        }

        totalExpectedCorrectOverlapSum += localExpectedCorrectOverlapSum * numReads;
        totalExpectedFPOverlapSum      += localExpectedFPOverlapSum      * numReads;
    }

    totalExpectedCorrectOverlapSum *= lambda * pnoerr;
    totalExpectedFPOverlapSum      *= lambda * potherpos * pword;

//    if (options.verbosity >= 2)
//    {
//        SEQAN_OMP_PRAGMA(critical(outputOS))
//        {
//            std::cerr << "(pos:" << i << ")(rl:" << lenError << "-" << k << ")"
//                      << " -- lambda - OvSum noerror/errors: (" << lambda << " , " << totalExpectedCorrectOverlapSum
//                      << " -- " <<  totalExpectedFPOverlapSum << ")" << std::endl;
//        }
//    }

	int cutoff = (int)((1 - w) * totalExpectedCorrectOverlapSum + w * totalExpectedFPOverlapSum);
	return _max(cutoff, 5);


//	SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) reduction(+:expovsumT1) reduction(+:expovsumT2))
//	for (int j = 1; j <= i - k; ++j)
//    {
//		int mnov = l - j - k ; //max n. overlapping bases with the read
//		for (int nerr = 0; nerr < mnov; nerr++)
//        {
//			expovsumT1 += (k + 1 + mnov - nerr) * matZerr(nerr, mnov);
//			expovsumT2 += (k + 1 + nerr) * matZotherov(nerr, mnov);
//		}
//	}
//	//Now the reverse strand
//	SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) reduction(+:expovsumT1) reduction(+:expovsumT2))
//	for (int j = l; j >= i + k; j--)
//    {
//		int mnov = j - k - 1;
//		for (int nerr = 0; nerr < mnov; nerr++)
//        {
//			expovsumT1 += (k + 1 + mnov - nerr) * matZerr(nerr, mnov);
//			expovsumT2 += (k + 1 + nerr) * matZotherov(nerr, mnov);
//		}
//	}
//	expovsumT1 *= lambda * pnoerr;
//	expovsumT2 *= lambda * potherpos * pword;
//    if (options.verbosity >= 2)
//        std::cerr << "(pos:" << i << ")(rl:" << l << "-" << k << ")"
//                  << " -- lambda - OvSum noerror/errors: (" << lambda << " , " << expovsumT1
//                  << " --" <<  expovsumT2 << ")" << std::endl;
	//Todo (Hugues) find a more rational way for this parameter, this one could be too stringent
	//one could hypothesize that the variance is equal to expectation and use a binomial approx
	//on the T2 values
//	int cutoff = ceil((1-w)*expovsumT1 + w*expovsumT2);
//	cutoff = (cutoff < 5) ? 5 : cutoff;
//	return (cutoff);
}
//
//
////
//// Compute the overlapsum cutoff for each position in a read
template <typename TCutOffMatrix, typename TPrefixLen, typename TReadLengthHist, typename TGenomeLen>
void ComputeCutoffOverlapSum(
	TCutOffMatrix & thresholds,
	TPrefixLen k,
    TReadLengthHist &readLenHist,
	TGenomeLen genomeLength,
    FionaOptions & options)
{
    unsigned maxReadLength = length(readLenHist) - 1;

    String<double> expectedCorrectOverlapSum;
    String<double> expectedFPOverlapSum;

    precomputeOverlapCombinatorics(
        expectedCorrectOverlapSum,
        expectedFPOverlapSum,
        maxReadLength - k - 1,  // the maximal value of mnov is readlen-k
        k,
        options);

	thresholds.resize(maxReadLength + 1, maxReadLength);  // readLength, position
	thresholds.clear();

	for (unsigned lenError = 1; lenError <= maxReadLength; ++lenError)
    {
        if (readLenHist[lenError] == 0)
            continue;
        //std::cout << '.' << std::flush;

        int stepSize = 1 + (lenError / 150);

        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1))
        for (int posError = 0; posError < (int)lenError / 2; ++posError)
        {
            double thresh = -1;

            if (posError <= k || posError % stepSize == 0)
                thresh = OddsOverlapSumCutoff(
                            posError,
                            1.0 / genomeLength,
                            lenError,
                            k,
                            options.errorrate,
                            options.wovsum,
                            expectedCorrectOverlapSum,
                            expectedFPOverlapSum,
                            readLenHist,
                            options);

            thresholds(lenError, posError) = thresh;
            thresholds(lenError, lenError - posError - 1) = thresh;
        }
        for (int posError = 1; posError < (int)lenError; ++posError)
        {
            if (thresholds(lenError, posError) == -1)
                thresholds(lenError, posError) = thresholds(lenError, posError - 1);
        }
    }

//	//for (int k = kmin; k <= kmax; k++)
//	double expected = (double)readcount / genomelength;
//	unsigned pad = readLen / 150 + 1;
//	for (unsigned len = 1; len <= readLen; ++len)
//    {
//		//std::cerr << "Computing for length " << i << "expected is " << expected << std::endl;
//		if (len >= readLen / 2 + 1)
//			thresholds[len] = thresholds[readLen - len + 1];
//		else {
//			if (len == 1 || len % pad == 0u)
//				thresholds[len] = OddsOverlapSumCutoff(len, expected, readLen, k, errorrate, w, matZerr, matZotherov, verbosity);
//			else
//				thresholds[len] = thresholds[len - 1];
//		}
//	}
}




/*estimated a median value for a given level*/
template < typename TIndex, class TSpec >
double medianLevel(Iter<TIndex, VSTree<TSpec> > iter){

	double totalOccs = 0.0;
	double sumMedian = 0.0;
	double median = 0.0;
	double mediumTotalOccs = 0.0;

  std::map<unsigned, unsigned> vectorOccurrences;

	goBegin(iter);
	for (; !atEnd(iter); ++iter)
	{
		unsigned numOccs = countOccurrences(iter);
		++vectorOccurrences[numOccs];
		totalOccs += numOccs;
	}

	mediumTotalOccs = totalOccs / 2.0;

  std::map<unsigned,unsigned>::iterator iterMap;
	for (iterMap = vectorOccurrences.begin (); iterMap != vectorOccurrences.end (); ++iterMap)
	{
		sumMedian += iterMap->second*iterMap->first;
		if (sumMedian >= mediumTotalOccs)
		{
			median = iterMap->first;
			break;
		}
	}
	return median;
}

template <typename TPercentage, typename TSize>
inline double probabilityOneError(TPercentage percentageErr, TSize repLen)
{
	return 1.0 - pow(1.0 - percentageErr, (double)repLen);
}

//Compute the number of ways of placing at most e errors in a read of length l such that
//any interval of length k has at least one error
template <typename TMatrix>
void CombinatoricsNoSeed(TMatrix &m, int lread, int kmax, int nerrmax)
{
	int m1 =  2*kmax-1 < lread ? (2*kmax-1) : lread ;
	int i,k;
	//Verify the size
	//std::cerr << "resizing\n";
	m.resize(nerrmax + 1, lread + 1);
	m.clear();
	for (i = 0; i <= lread; ++i)
		for (k=0; k<= kmax; ++k)
			m(k,i) = 0;
	//std::cerr << "all values to 0\n";	
	for (i = 0; i <= lread; ++i)
		m(0,i) = m(1,i) = 0;
	for (i = 0; i < kmax; ++i)
		m(1,i) = i;
	for (i = kmax; i <= m1; ++i)
		m(1, i) = 2 * kmax - i;
	for (int ne = 2; ne <= nerrmax; ++ne)
    {
		for (k = 0; k < ne; ++k)
			m(ne, k) = 0;
		for (k = ne; k <= kmax; ++k)
			m(ne, k) = (int)boost::math::binomial_coefficient<double>(k, ne);
		int mpos = (ne + 1) * kmax;
		for (k = mpos; k <= lread; ++k)
			m(ne, k) = 0;
		int mposlread = mpos < (lread +1) ? mpos : (lread+1);
		for (k = kmax + 1; k < mposlread; ++k)
        {
			unsigned csum = 0;
			for (int j = k-kmax; j<k; j++)
				csum += m(ne-1, j);
			m(ne,k) = csum;
		}
	}
	//should print here for debug
//	std::cerr << "Computed the Combinatorial matrix\nk = " << kmax << ", l = " << lread << std::endl;
//	for (unsigned j = 1; j <= nerrmax; ++j){
//		for (i = 1; i <= lread; ++i)
//			std::cerr << "\t" << m(j,i);
//		std::cerr << std::endl;
//	}		
}


//Expected number of correctable reads


//Expected number of uncorrectable reads
template <typename TUExpCounts, typename TPrefixLen, typename TErrorRate, typename TReadLen, typename TNReads>
void UncorrectableExpected(
		TUExpCounts & UncorrExp, 
		TPrefixLen const kmin, 
		TPrefixLen const kmax,
		TReadLen const lread, 
		TNReads const nreads,
		TErrorRate const perr){
	clear(UncorrExp);
	resize(UncorrExp, kmax+1, 0.0);
	//Proba place exactly i errors
	String <double> pkerrs;
	clear(pkerrs);
	resize( pkerrs, kmax+1, 0.0);
	for (int i =0; i<=kmax; ++i)
		pkerrs[i] = pow(perr, (double)i) * pow(1.0 - perr, (double)(lread - i));
	matrix <unsigned> MatNoSeed (kmax+1, lread+1);
	for (int k = kmin; k <= kmax; k++){
		CombinatoricsNoSeed(MatNoSeed, lread, k, k);
		for (int nerr = 1; nerr < k; ++nerr){
			UncorrExp[k] += MatNoSeed(nerr,lread) * pkerrs[nerr] * nreads;
		}
	}
}

template <typename TUExpCounts, typename TPrefixLen, typename TErrorRate, typename TNumRead>
void UncorrectableExpectedBases(
		TUExpCounts & uncorrExp, 
		TPrefixLen const kmin, 
		TPrefixLen const kmax,
		String<TNumRead> readLenHist,
		TErrorRate const perr)
{
	clear(uncorrExp);
	resize(uncorrExp, kmax + 1, 0.0);

	//Proba place exactly i errors
	matrix<double> pkerrs(kmax + 1, length(readLenHist));

    // precompute the probabilities to have k errors in a read of length readLen
	SEQAN_OMP_PRAGMA(parallel for)
	for (int readLen = 1; readLen < (int)length(readLenHist); ++readLen)
        for (TPrefixLen k = 1; k <= kmax; ++k)
            pkerrs(k, readLen) = pow(perr, (double)k) * pow(1.0 - perr, (double)(readLen - k));

    // distribute the work (interval [kmin..kmax+1)) over different threads
    Splitter<int> splitter(kmin, kmax + 1);
	SEQAN_OMP_PRAGMA(parallel for num_threads(length(splitter)) schedule(static))
	for (int i = 0; i < (int)length(splitter); ++i)
    {
        // add for each anchor k the number of uncorrectable reads
        matrix<unsigned> matNoSeed;
        for (int k = splitter[i]; k < splitter[i + 1]; ++k)
        {
            CombinatoricsNoSeed(matNoSeed, length(readLenHist) - 1, k, k);
            for (unsigned readLen = 1; readLen < length(readLenHist); ++readLen)
            {
                TNumRead numReads = readLenHist[readLen];
                if (numReads == 0)
                    continue;

                for (int numErrors = 1; numErrors < k; ++numErrors)
                    uncorrExp[k] += matNoSeed(numErrors, readLen) * pkerrs(numErrors, readLen) * (double) numReads * (double) readLen;
            }
        }
	}
}


//Expected number of destructible reads
template <typename TDestCounts, typename TPrefixLen, typename TErrorRate, typename TReadLen, typename TNReads, typename TGenomeLen>
void DestructibleExpected(
		TDestCounts & DestrExp, 
		TPrefixLen const kmin, 
		TPrefixLen const kmax,
		TReadLen const lread, 
		TNReads const nreads,
		TErrorRate const perr, 
		TGenomeLen const genomelen){
	clear(DestrExp);
	resize(DestrExp, kmax+1, 0.0);
	double muw = pow(4.0, kmin - 1.0);
	for (int k = kmin; k <= kmax; ++k){
		muw *= 4;
		double qw = (1- pow(1-perr, (double)k))*(1-perr)*(1-pow(1-1.0/muw, (double)genomelen))* (3.0/4);
		DestrExp[k] = (1- pow(1-qw, (double)(lread-k)))*pow(1-perr, (double)lread) * nreads;
	}
}

//Expected number of destructible reads
template <typename TDestCounts, typename TPrefixLen, typename TErrorRate, typename TNumRead, typename TGenomeLen>
void DestructibleExpectedBases(
		TDestCounts & destrExp, 
		TPrefixLen const kmin, 
		TPrefixLen const kmax,
		String<TNumRead> readLenHist,
		TErrorRate const perr, 
		TGenomeLen const genomelen)
{
	clear(destrExp);
	resize(destrExp, kmax + 1, 0.0);

	SEQAN_OMP_PRAGMA(parallel for)
    for (int k = kmin; k <= kmax; ++k)
    {
        double muw = pow(4.0, (double)k);
        double qw = (1 - pow(1 - perr, (double)k)) * (1 - perr) * (1 - pow(1 - 1.0 / muw, (double)genomelen)) * 0.75;

        for (unsigned readLen = 1; readLen < length(readLenHist); ++readLen)
        {
            uint64_t numReads = readLenHist[readLen];
            if (numReads == 0)
                continue;

            destrExp[k] += (1 - pow(1 - qw, (double)(readLen - k))) * pow(1 - perr, (double)readLen) * (double) numReads * (double) readLen;
        }
    }
}

template <typename TDestCounts, typename TPrefixLen, typename TErrorRate, typename TReadLen, typename TNReads, typename TGenomeLen>
void DestructibleExpectedFiona(
						  TDestCounts & DestrExp, 
						  TPrefixLen const kmin, 
						  TPrefixLen const kmax,
						  TReadLen const /*lread*/, 
						  TNReads const /*nreads*/,
						  TErrorRate const /*perr*/, 
						  TGenomeLen const /*genomelen*/){
	clear(DestrExp);
	resize(DestrExp, kmax+1, 0.0);
//	double muw = pow(4.0, kmin -1);
	for (unsigned k = kmin; k <= kmax; ++k){
		//This version should account for the fact that Fiona gets the best correction over the range of Ks
	}
}



/*precomputation of thresholds for varying k*/
template <typename TThresholds, typename TExpectedValues, typename TStrictness, typename TErrorRate, typename TOddsError, typename TPrefixLen>
void ComputeCutoffErroneous(
		TThresholds & thresholds,
		TExpectedValues & ,
		TStrictness const cutoff,
		TErrorRate  const,
		TOddsError  const,
		TPrefixLen  const kmin,
		TPrefixLen  const kmax,
		FionaCount  const)
{
	clear(thresholds);
	resize(thresholds, kmax+1, 0);
	for (int k = kmin; k <= kmax; k++){
		thresholds[k] = (int)cutoff;
	}
}


//Fixed count mode
template <typename TThresholds, typename TExpectedValues, typename TStrictness, typename TErrorRate, typename TOddsError, typename TPrefixLen>
void ComputeCutoffErroneous(
	TThresholds & thresholds,
	TExpectedValues & expected,
	TStrictness const,
	TErrorRate const,
	TOddsError const,
	TPrefixLen const kmin,
	TPrefixLen const kmax,
	FionaExpected const)
{
		clear(thresholds);
		resize(thresholds, kmax+1, 0);
	for (int k = kmin; k <= kmax; k++){
		thresholds[k] = (int)expected[k];
	}
}

//Poisson pvalue mode
template <typename TThresholds, typename TExpectedValues, typename TStrictness, typename TErrorRate, typename TOddsError, typename TPrefixLen>
void ComputeCutoffErroneous(
	TThresholds & thresholds,
	TExpectedValues & expected,
	TStrictness const strictness,
	TErrorRate  const,
	TOddsError  const,
	TPrefixLen  const kmin,
	TPrefixLen  const kmax, 
	FionaPoisson const)
	{
		clear(thresholds);
		resize(thresholds, kmax+1, 0);
		for (int k = kmin; k <= kmax; k++){
			thresholds[k] = qpois(strictness, expected[k]);
		}
	}


//Poisson sensitivity mode
template <typename TThresholds, typename TExpectedValues, typename TStrictness, typename TErrorRate, typename TOddsError, typename TPrefixLen>
void ComputeCutoffErroneous(
	TThresholds & thresholds,
	TExpectedValues & expected,
	TStrictness const falsenegrate,
	TErrorRate const errorrate,
	TOddsError const,
	TPrefixLen const kmin,
	TPrefixLen const kmax,
	FionaPoissonSens const)
{
	clear(thresholds);
	resize(thresholds, kmax+1, 0);
	for (int k = kmin; k<= kmax; k++){
		thresholds[k] = 1 + qpoismixerror(1-falsenegrate, expected[k], errorrate, k); 
	}
}

// Poisson Classification mode
// priorerror is the parameter for more stringent (>1) more loose (<1) error detection. 
// oddserrorreads is the computed ratio of erroneous over correct reads
// value of 0 turns automatic prior computation of the amount of erroneous reads.
template <typename TThresholds, typename TExpectedValues, typename TStrictness, typename TErrorRate, typename TOddsError, typename TPrefixLen>
void ComputeCutoffErroneous(
	TThresholds & thresholds,
	TExpectedValues & expected,
	TStrictness const priorerror, 
	TErrorRate const errorrate,
	TOddsError const oddserrorreads, 
	TPrefixLen const kmin,
	TPrefixLen const kmax,
	FionaPoissonClassif const)
{
	clear(thresholds);
	resize(thresholds, kmax+1, 0);
	for (int k = kmin; k<= kmax; k++){
		thresholds[k] = PoisClassifCutoff(priorerror*oddserrorreads, expected[k], errorrate, k); 
		//(hugues) should we always add 1 just in case ?
		//thresholds[k]++;
	}
}


/*
* Linear Model fitting used for determination of number of Rounds by fitting
* log data to capture deviation from a exponential distribution
*/
struct LinearModel
{
	double		intercept;
	double 		slope;
	unsigned int	numberObservations;
	unsigned int	numberPredictors;
};


/* compute a fitted value under a linear regression model */
template <typename LinearModel, typename TValue>
inline TValue fittedValue(
        LinearModel const &linearModel,
	TValue  x) 
{
	return((TValue)(linearModel.intercept + (linearModel.slope * (double)x)));
}

	/* use the standard maximum likelihood estimator for linear regression 
	 * for datapoints x=(x_1, ..., x_n) and y=(y_1, .., y_n)  */
template <typename LinearModel, typename TValue>
inline void linearRegression(
        LinearModel &linearModel,
        String<TValue> const  &x,
        String<TValue> const &y)
{
	/* get the means first  */
	TValue meanX=0;
	TValue meanY=0;
	for(unsigned int i=0;i<length(x);i++)
	{
		meanX += x[i];
		meanY += y[i];
	}
	meanX = meanX/(TValue)length(x);
	meanY = meanY/(TValue)length(y);

	/* use the standard maximum likelihood estimator for linear regression  
	 * and compute first the slope than the intercept of the linear function */
	TValue covarianceXY = 0;
	TValue varianceX = 0;
    for(unsigned int i=0;i<length(x);i++)
    {
        covarianceXY += (x[i] - meanX) * (y[i] - meanY);
        varianceX    += (x[i] - meanX) * (x[i] - meanX);
    }

	/* save the parameters in the model */
	linearModel.slope     = (double) covarianceXY/varianceX;
	linearModel.intercept = (double) (meanY - ((TValue)linearModel.slope * meanX) );
	linearModel.numberObservations = (unsigned int) length(x);
    linearModel.numberPredictors  = (unsigned int) 1;
}

       /* compute R-Square (or Coefficient of Determination) for a set of values that
	* have been fit using a linear regression model */ 
template <typename LinearModel, typename TValue>
inline TValue RSquare(
        LinearModel const &linearModel,
        String<TValue> const &x,
        String<TValue> const &y)
{
        /* get the mean of y  */
        TValue meanY=0;
        for(unsigned int i=0;i<length(y);i++)
        {
                meanY += y[i];
        }
        meanY = meanY/(TValue)length(y);
	/* R-Square is defined as  1 - SSerror/SStotal, where SSerror is the sum of residual errors
	 * and SStotal is the variance of the y values	*/
	TValue SSerror = 0;
	TValue SStotal =0;
    for(unsigned int i=0;i<length(x);i++)
    {
        SStotal  += (y[i] - meanY) * (y[i] - meanY);
        SSerror  += pow((y[i] - fittedValue(linearModel,x[i])), 2.0);
    }
	return((TValue) ( 1- (SSerror/SStotal)));
}

       /* compute the adjusted R-Square for a set of values that
	* have been fit using a linear regression model */ 
template <typename LinearModel, typename TValue>
inline TValue adjustedRSquare(
        LinearModel const &linearModel,
        String<TValue> const &x,
        String<TValue> const &y)
{
	/* The adjusted R-Square value is nothing but the R-Square value corrected by the number
	 * of observed variables n and the number of predictors k
	 * AdjRSquare = 1 - (1-RSquare)*(n-1)/(n-k-1) */

	TValue R_2 = RSquare(linearModel,x,y);
	return( (TValue) (1 - (TValue)(1-R_2) * (TValue)(linearModel.numberObservations -1)/(TValue)(linearModel.numberObservations - linearModel.numberPredictors -1)));
}

template <typename TFragmentStore, typename TCorrection>
inline void _dumpCorrection(
	TFragmentStore &store,
	TCorrection const &correction,
	unsigned errorReadId)
{
	std::cerr << std::endl;
	std::cerr << "error___read_id\t" << errorReadId << std::endl;
    std::cerr << "error_pos      \t" << correction.errorPos << std::endl;
    std::cerr << "error_read length \t"<<length(store.readSeqStore[errorReadId])<<std::endl;
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
    std::cerr << "correct_read_id\t" << correction.correctReadId << std::endl;
	std::cerr << "correct_pos    \t" << correction.correctPos << std::endl;
    std::cerr << "correct_read length \t"<<length(store.readSeqStore[correction.correctReadId])<<std::endl;
#endif
	std::cerr << "overlap        \t" << correction.overlap << std::endl;
	std::cerr << "indel_len      \t" << (int)correction.indelLength << std::endl;
	std::cerr << "error___read   \t" << store.readSeqStore[errorReadId][correction.errorPos] << '\t';
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
	for (unsigned i = 0; i < correction.correctPos; ++i)
		std::cerr << ' ';
#endif
    std::cerr << store.readSeqStore[errorReadId] << std::endl;
#ifndef FIONA_CONSENSUS_REDUCE_MEMORY
	std::cerr << "correct_read   \t" << store.readSeqStore[correction.correctReadId][correction.correctPos] << '\t';
	for (unsigned i = 0; i < correction.errorPos; ++i)
		std::cerr << ' ';
    std::cerr << store.readSeqStore[correction.correctReadId] << std::endl;
#endif
}	

/*change the erroneous nucleotide in all reads identify with errors*/
template <typename TFragmentStore, typename TCorrections>
void applyReadErrorCorrections(
	TFragmentStore &store,
	TCorrections const &corrections,
    FionaOptions const & options)
{
	typedef typename Value<TCorrections>::Type TCorrection;
	int readCount = length(corrections);

	// we have to make a temp-copy in order to use original (not corrected) reads for correction
	StringSet<typename TFragmentStore::TReadSeq> originalReads;
	resize(originalReads, length(store.readSeqStore), Exact());

	SEQAN_OMP_PRAGMA(parallel for schedule(guided))
	for (int readId = 0; readId < readCount; ++readId)
	{
		TCorrection const &corr = corrections[readId];
		if (corr.overlap != 0 && corr.indelLength <= 0)
			originalReads[corr.correctReadId] = store.readSeqStore[corr.correctReadId];
	}

	SEQAN_OMP_PRAGMA(parallel for)
	for (int readId = 0; readId < readCount; ++readId)
	{
#ifdef SEQAN_VERBOSE
		std::cerr << "at readID: "<<readId<<std::endl;
#endif
		TCorrection const &corr = corrections[readId];
		if (corr.overlap == 0) continue;

        if (readId == options.debugRead)
        {
            std::cerr << "BEFORE:" << std::endl;
            _dumpCorrection(store, corr, readId);
        }

    std::ostringstream m;
		if (strContains(toCString(store.readNameStore[readId]), "corrected"))
			m << "\t";
		else
			m << " corrected:\t";

#ifdef SEQAN_VERBOSE
		_dumpCorrection(store, corr, readId);
#endif

		if (corr.indelLength == 0)
		{	

			m << corr.errorPos <<  "(" << options.cycle <<',' <<corr.overlap << ","<< corr.correctReadId<<"):" << store.readSeqStore[readId][corr.errorPos] << "->" << originalReads[corr.correctReadId][corr.correctPos];
			store.readSeqStore[readId][corr.errorPos] = originalReads[corr.correctReadId][corr.correctPos];
		}
#ifdef FIONA_ALLOWINDELS
		else if (corr.indelLength > 0)
		{
			m << corr.errorPos << "(" << options.cycle <<',' <<corr.overlap << ","<< corr.correctReadId<<"):-" << infix(store.readSeqStore[readId], corr.errorPos, corr.errorPos + corr.indelLength);
            erase(store.readSeqStore[readId], corr.errorPos, corr.errorPos + corr.indelLength);
		} else {
			m << corr.errorPos << "(" << options.cycle <<',' <<corr.overlap << ","<< corr.correctReadId<<"):+" << infix(originalReads[corr.correctReadId], corr.correctPos, corr.correctPos + -corr.indelLength);
			insert(store.readSeqStore[readId], corr.errorPos, infix(originalReads[corr.correctReadId], corr.correctPos, corr.correctPos + -corr.indelLength));
		}
#endif
		append(store.readNameStore[readId], m.str());
        if (readId == options.debugRead)
        {
            std::cerr << "AFTER:" << std::endl;
            _dumpCorrection(store, corr, readId);
        }
#ifdef SEQAN_VERBOSE
		std::cerr << "corrected:";
		for (unsigned i = 0; i < corr.correctPos; ++i)
			std::cerr << ' ';
		std::cerr << store.readSeqStore[readId] << std::endl;
#endif
	}
}

template <typename TObservedCount, typename TCutoffCount>
inline bool
potentiallyErroneousNode(TObservedCount observed, TCutoffCount cutoff)
{
	return (long)observed < cutoff;
}


//template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
//inline bool potentiallyErroneousNode(
//	TObserved observed,
//	TExpected expected,
//	TStrictness strictness,
//	TErrorRate,
//	TPrefixLen,
//	FionaPoisson const)
//{
//	// compare the cumulative poisson distribution with the p-value (strictness)
//    double negExp = exp(-expected);
//    double pValue = 0.0;
//	double pow = 1.0;
//	double fact = 1.0;
//
//    for (TObserved i = 0; i <= observed && pValue <= strictness; ++i, fact *= i){
//        pValue += pow * negExp / fact;
//		pow *= expected;
//	}
//	return pValue <= strictness;
//}
//
//template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
//inline bool potentiallyErroneousNode(
//	TObserved observed,
//	TExpected expected,
//	TStrictness,
//	TErrorRate,
//	TPrefixLen,
//	FionaExpected const)
//{
//	// compare the weight for a node with a cutoff given by the strictness param.
//	return observed < expected;
//}
//
//template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
//inline bool potentiallyErroneousNode(
//	TObserved observed,
//	TExpected,
//	TStrictness cutoff,
//	TErrorRate,
//	TPrefixLen,
//	FionaCount const)
//{
//	// compare the weight for a node with a fixed cutoff
//	return observed < cutoff;
//}
//
//
//
//template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
//inline bool potentiallyErroneousNode(
//	TObserved observed,
//	TExpected expected,
//	TStrictness falsenegrate,
//	TErrorRate errorrate,
//    TPrefixLen prefixlen,
//	FionaPoissonSens const)
//{
//	// Poisson based threshold, given a fixed percentage of missed errors (1 - min sensitivity)
//	// consider only the cases with one and two errors
//	// the average error rate and the expected value allow to compute the expected count for an error.
//	
//	//special case when current node count ==1, we always consider that as an error
//	if(observed == (TObserved) 1) return true;
//
//	double sensitivity = 1 - falsenegrate;
//	double noerrlm2 = pow(1-errorrate, prefixlen - 2.0);
//	double noerrlm1 = noerrlm2 * (1-errorrate);
//	double errexp1err = expected * noerrlm1 * errorrate/3;
//	double errexp2err = expected * noerrlm2 * (errorrate/3) *(errorrate/3);
//	double perr1 = prefixlen * noerrlm1* errorrate ;
//	double perr2 = (prefixlen * (prefixlen -1) / 2 ) * noerrlm2 * errorrate * errorrate;
//	double sc = perr1 + perr2;
//	perr1 /= (sc); 
//	perr2 /= (sc); 
//	double negExp1err = exp(-errexp1err); //
//	double negExp2err = exp(-errexp2err);
//    double probaerror = 0.0;
//	double pow1err = 1.0;
//	double pow2err = 1.0;
//	double fact = 1.0;
//	//std::cerr << "Level " << prefixlen << " and params for correc\n";
//	//std::cerr << "Expected and comp: " << expected << " - " << errexp1err << " - " << errexp2err << "\n";
//	//std::cerr << "the basic observed count: " << observed << ", perr is " << perr1 << " - " << perr2  << "\n";
//	TObserved i = 0;
//    for (i = 0; i <= observed && probaerror <= sensitivity; ++i, fact *= i){
//        probaerror += perr1 * pow1err * negExp1err / fact;
//		probaerror += perr2 * pow2err * negExp2err / fact;
//		pow1err *= errexp1err; 
//		pow2err *= errexp2err;
//	}
//	//std::cerr << "Stopped at observed value **" << i-1 << "** for a sens. of " << probaerror << "and a thr at " << sensitivity << "\n";
//	return probaerror <= sensitivity;
//}
//
//template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
//inline bool potentiallyErroneousNode(
//	TObserved observed,
//	TExpected expected,
//	TStrictness priorerror,//the a priori odds of errors pi_err/(1-pi_err) (default should be 1)
//	TErrorRate errorrate,
//	TPrefixLen prefixlen,
//	FionaPoissonClassif const)
//{
//	// Poisson based threshold, we compute the logodds of being an error vs a genuine read
//	// consider only the cases with one and two errors
//	// the average error rate and the expected value allow to compute the expected count for an error.
//	
//	//special case when current node count ==1, we always consider that as an error
//	if(observed == (TObserved) 1) return true;
//	
//	double noerr = pow(1-errorrate, prefixlen);
//	double noerrlm2 = pow(1-errorrate, prefixlen - 2.0);
//	double noerrlm1 = noerrlm2 * (1-errorrate);
////	double errexpnoerr = expected * noerr;
////	double errexp1err = expected * noerrlm1 * (errorrate/3);
////	double errexp2err = expected * noerrlm2 * (errorrate/3) *(errorrate/3);
//	double perr1 = prefixlen * noerrlm1* (errorrate) ;
//	double perr2 = (prefixlen * (prefixlen -1) / 2 ) * noerrlm2 * (errorrate) *(errorrate);
//	double sc = perr1 + perr2;
//	perr1 /= (sc); 
//	perr2 /= (sc); 
////	double negExpnoerr = exp(-errexpnoerr);
////	double negExp1err = exp(-errexp1err); //
////	double negExp2err = exp(-errexp2err);
////	double pow1err = 1.0;
////	double pow2err = 1.0;
////	double pownoerr = 1.0;
////	double fact = 1.0;
//    double probaerror = 0.0;
//	double probanoerror = 0.0;
//	//std::cerr << "Level " << prefixlen << " and params for correc\n";
//	//std::cerr << "Expected and comp: " << errexpnoerr << " - " << errexp1err << " - " << errexp2err << "\n";
//	//std::cerr << "the basic observed count: " << observed << ", perr is " << perr1 << " - " << perr2  << "\n";
//	//TObserved i = 0;
//	probaerror = dpoismixerror( observed, expected, errorrate, prefixlen);
//	probanoerror = dpois(observed, expected*noerr);
////    for (i = 0; i <= observed; ++i, fact *= i){
////        probaerror += perr1 * pow1err * negExp1err / fact;
////		probaerror += perr2 * pow2err * negExp2err / fact;
////		pow1err *= errexp1err; 
////		pow2err *= errexp2err;
////		probanoerror += pownoerr * negExpnoerr / fact;
////	}
//	//std::cerr << "Stopped at observed value **" << i-1 << "with proba of err " << probaerror << "and probanoerr " << probanoerr << "\n";
//	return (log((probaerror)/ (probanoerror) * priorerror) > 0);
//}


struct Overlap  // for each operation (substitution/deletion/insert) there is one Overlap entry
{
    unsigned overlapSumLeft;
    unsigned overlapSumRight;
    unsigned readId;                // the readId of one correct candidate
    unsigned short correctPos;      // the position of the correct base in this candidate (reads shouldn't be longer than 65536bp)
#ifdef FIONA_MAXIMIZE_SUPPORT
    unsigned short _errorsRight;    // temporary
    unsigned support;               // absolute number of reads supporting this operation (with minimal errors right of the anchor)
#endif
#ifdef FIONA_CONSENSUS
    typedef ProfileChar<Dna5, unsigned short>   TProfileValue;
    typedef String<TProfileValue>               TConsensus;
    
    TConsensus consensus;
#endif
};

// this function is used to extend the read suffixes after a common seed with at most #maxErrors errors
template <typename TIter>
inline unsigned
_comparePrefixesWithEditDistance(
	Range<TIter> seq1,
	Range<TIter> seq2,
	unsigned maxErrors)
{
    typedef Range<TIter> TSeq;
    typename Size<TSeq>::Type len1 = length(seq1);
    typename Size<TSeq>::Type len2 = length(seq2);

    PatternState_<TSeq, Myers<AlignTextBanded<FindPrefix, NMatchesAll_, NMatchesAll_>, True, void> > state;
    state.maxErrors = maxErrors;
//    state.leftClip = maxErrors / 2;
    state.leftClip = maxErrors;

    typename Iterator<TSeq, Rooted>::Type iter;
    unsigned minErrors = maxErrors + 1;

    if (len1 <= len2)
    {
        if (len1 == 0)
            return 0;
        
        seq2.end = seq2.begin + _min(len2, len1 + state.leftClip);
        iter = begin(seq2, Rooted());
        if (_patternInitSmallStateBanded(iter, seq1, state))
            for (; !atEnd(iter) && _findMyersSmallPatternsBanded(iter, seq1, state, True()); goNext(iter))
                if (minErrors > state.errors)
                    minErrors = state.errors;
    }
    else
    {
        seq1.end = seq1.begin + _min(len1, len2 + state.leftClip);
        iter = begin(seq1, Rooted());
        if (_patternInitSmallStateBanded(iter, seq2, state))
            for (; !atEnd(iter) && _findMyersSmallPatternsBanded(iter, seq2, state, True()); goNext(iter))
                if (minErrors > state.errors)
                    minErrors = state.errors;
    }
    return minErrors;
}

template <typename TIter>
inline unsigned
_comparePrefixesWithEditDistanceReverse(
	Range<TIter> seq1,
	Range<TIter> seq2,
	unsigned maxErrors)
{
    typedef Range<TIter> TSeq;
    typename Size<TSeq>::Type len1 = length(seq1);
    typename Size<TSeq>::Type len2 = length(seq2);
    typedef ModifiedString<TSeq, ModReverse> TRev;

    PatternState_<TRev, Myers<AlignTextBanded<FindPrefix, NMatchesAll_, NMatchesAll_>, True, void> > state;
    state.maxErrors = maxErrors;
//    state.leftClip = maxErrors / 2;
    state.leftClip = maxErrors;

    typename Iterator<TRev, Rooted>::Type iter;
    unsigned minErrors = maxErrors + 1;

    if (len1 <= len2)
    {
        if (len1 == 0)
            return 0;
        
        seq2.begin = seq2.end - _min(len2, len1 + state.leftClip);

        TRev rseq1(seq1);
        TRev rseq2(seq2);
        iter = begin(rseq2, Rooted());
        if (_patternInitSmallStateBanded(iter, rseq1, state))
            for (; !atEnd(iter) && _findMyersSmallPatternsBanded(iter, rseq1, state, True()); goNext(iter))
                if (minErrors > state.errors)
                    minErrors = state.errors;
    }
    else
    {
        seq1.begin = seq1.end - _min(len1, len2 + state.leftClip);

        TRev rseq1(seq1);
        TRev rseq2(seq2);
        iter = begin(rseq1, Rooted());
        if (_patternInitSmallStateBanded(iter, rseq2, state))
            for (; !atEnd(iter) && _findMyersSmallPatternsBanded(iter, rseq2, state, True()); goNext(iter))
                if (minErrors > state.errors)
                    minErrors = state.errors;
    }
    return minErrors;
}

//newfunction for traversing and filling the linked list with/*detect and repair the reads with errors*/
template <
    int LOOP_LEVEL,
    typename TTreeIterator,
    typename TFragmentStore,
    typename TCorrections,
    typename TValueId,
    typename TAlgorithm >
void traverseAndSearchCorrections(
	TTreeIterator iter,
	TFragmentStore &store,
	String<TCorrections> & correctionList,
	String<TValueId> &firstCorrectionForRead,
	FionaOptions & options,
	Tag<TAlgorithm> const,
	unsigned readLength,
    FionaResources &resources)
{
	typedef typename Container<TTreeIterator>::Type TFionaIndex;
	typedef typename Fibre<TFionaIndex, FibreText>::Type TReadSet;
	typedef typename Fibre<TFionaIndex, FibreSA>::Type TSA;
	typedef typename Infix<TSA const>::Type TOccs;
	typedef typename Iterator<TOccs, Standard>::Type TOccsIterator;
	typedef typename Value<TReadSet>::Type TRead;
	typedef typename Value<TRead>::Type TValue;
	typedef typename Iterator<TRead, Standard>::Type TReadIterator;

    double start = omp_get_wtime();
    TFionaIndex &index = container(static_cast<TTreeIterator&>(iter));
	unsigned readCount = length(store.readSeqStore) / 2;
	//for debugging of read data
	String<TOccs, Array<4> > correctCandidates;     // there are at most 4 correcting branches
    Overlap bestCorrection[1+2*MAX_INDEL_LENGTH];   // for a given branch store for every indel-size the best correction
    Dna5 correctSeq[MAX_INDEL_LENGTH];
	const TValue unknownChar = unknownValue<TValue>();

	//compute maximum allowed mismatches per read
	//We allow up to maxAcceptedMismatches between two reads. The threshold is decided 
	//by controlling the proportion pNeig of the reads that are expected with that many mismatches (given the error rate)
    unsigned replen_min = cargo(index).replen_min;
    unsigned cycle = options.cycle;
	float pNeig = 0.95f;
	binomial Nmismatch(readLength, options.errorrate);
	unsigned maxAcceptedMismatches = _max((unsigned) ceil(quantile(Nmismatch, pNeig)), 2u);   //unlikely that 2 reads share the same error         (weese:) don't understand the comment
//	double oldAcceptedMismatches = (options.errorrate * readLength);

//	std::cout << std::endl << "acceptedMismatches: "  << maxAcceptedMismatches << " readLength " << readLength << "  error " << options.errorrate << std::endl;
//	std::cout << std::endl << "Old acceptedMismatches value: "  << oldAcceptedMismatches << std::endl;

//	for (goBegin(iter); !atEnd(iter); )                     // do a DFS walk
	while (!atEnd(iter))                                        // do a DFS walk
    {
//        if (value(iter).range.i1==1798 && value(iter).range.i2==1811 && value(iter).repLen==23)
//        {
//            std::cerr<<std::endl<<"START";
//            std::cerr << representative(iter)<<std::endl;
//            for(int i=0; i<countOccurrences(iter); ++i)
//                std::cerr << getOccurrences(iter)[i] << '\t';
//            std::cerr<<"END"<<std::endl;
//        }
		unsigned commonPrefix = parentRepLength(iter);			// length of parent label
//        if (seqan::range(iter).i1==127 && prefix(representative(iter),commonPrefix) == "AAAAACAAAAACA")
//        std::cout<<"HERE"<<std::endl;

        if ((int)commonPrefix < options.fromLevel)
        {
            goNext(iter);
            continue;
        }

#ifdef FIONA_USE_SA
        if ((int)commonPrefix > options.toLevel)
            SEQAN_FAIL("how can that be?");
#endif

		SEQAN_ASSERT_LT(commonPrefix + 1, length(options.expectedTheoretical));
//		SEQAN_ASSERT_LEQ(options.fromLevel, (int)commonPrefix); // doesn't hold for the first node (=root node)

        // only examine branches where the string depth is a multiple of depthSampleRate
		TValue firstEdgeChar = parentEdgeFirstChar(iter);
        bool skipNode = ((commonPrefix + cycle - replen_min) % options.depthSampleRate != 0);
        bool leaveNodeRight = (int)commonPrefix >= options.toLevel || countOccurrences(iter) < 3 || firstEdgeChar == unknownChar;

#ifdef FIONA_USE_SA
        if (!skipNode && isFirstChild(iter) && isLastChild(iter))
            skipNode = true;
#endif

		if (!skipNode && firstEdgeChar != unknownChar &&   // N is always an error
			!potentiallyErroneousNode(countOccurrences(iter), options.errorCutoffs[commonPrefix+1])) //New test is easier
			//!potentiallyErroneousNode(countOccurrences(iter), options.expectedTheoretical[commonPrefix+1], options.strictness, options.errorrate, commonPrefix + 1 ,alg))
		{
            skipNode = true;
		}

        // don't descent over repeats
        if (!skipNode)
        {
            if (/*countOccurrences(iter) >= options.repeatCutoffs[repLength(iter)] &&*/
//                isSelfRepetitive(representative(iter)))
                isSelfRepetitive(infixWithLength(indexText(index), getOccurrence(iter), commonPrefix + 1)))
            {
                // anchor plus first edge char is repetitive -> skip the whole parent
                skipNode = true;
                leaveNodeRight = true;
            }
        }

        if (skipNode)
        {
            if (leaveNodeRight)
                goNextRight(iter);
            else
                goNext(iter);
            continue;
        }


        ++resources.investigatedNodes;
		//
		//	get the id and position (suffix begin) for suspected nodes
		//	for which we can find a more optimal correction
		//
		TOccs errorCandidates = getOccurrences(iter);

		/*copy the iterator for iterate over the siblings*/
//		typename Iterator<TFionaIndex, TopDown<> >::Type iterSibling(index, nodeUp(iter));
		TTreeIterator iterSibling(iter);
		goUp(iterSibling);

		//
		//	potential reads for make the correction,
		//	because at the same level, with the same prefix
		//

		clear(correctCandidates);
		if (!goDown(iterSibling))
			SEQAN_ASSERT_FAIL("going up and down failed!?");


if (LOOP_LEVEL != 0)
{
        
		// pick potentially correct reads
        unsigned long thickestBranchCount = length(errorCandidates);
        TOccs thickestBranchOccs;
		do
		{
            TValue siblingFirstEdgeChar = parentEdgeFirstChar(iterSibling);
			if (siblingFirstEdgeChar != firstEdgeChar && siblingFirstEdgeChar != unknownChar)
            {
                unsigned long siblingOccCount = countOccurrences(iterSibling);
                // record the thickest branch
                if (thickestBranchCount < siblingOccCount)// && siblingOccCount < options.repeatCutoffs[commonPrefix+1])
                {
                    thickestBranchCount = siblingOccCount;
                    thickestBranchOccs = getOccurrences(iterSibling);
                }
                
                if (!potentiallyErroneousNode(siblingOccCount, options.errorCutoffs[commonPrefix + 1]) &&
                    siblingOccCount < options.repeatCutoffs[commonPrefix+1])
                    //!potentiallyErroneousNode(countOccurrences(iterSibling), options.expectedTheoretical[commonPrefix + 1], options.strictness, options.errorrate,commonPrefix +1, alg))
                {
                    // save the id and position(where the suffix begin in the reads) in the table of IDs correct
                    // also the number of occurrences
                    appendValue(correctCandidates, getOccurrences(iterSibling));
                }
            }
            if (ordValue(siblingFirstEdgeChar) == 3) break; // we ignore the N-node right of the T
		} while (goRight(iterSibling));

        // if there aren't any, try the thickest branch
        if (empty(correctCandidates)  &&  thickestBranchCount > length(errorCandidates))
            appendValue(correctCandidates, thickestBranchOccs);
}

		// continue only if we have found a correct read
		if (!empty(correctCandidates))
        {
//std::cout << seqan::range(iter) << '\t' << parentEdgeFirstChar(iter) << '\t' << prefix(representative(iter),commonPrefix) << '\t';
//for (unsigned j = 0; j < length(correctCandidates); ++j)
//    std::cout << beginPosition(correctCandidates[j]) << ',' << endPosition(correctCandidates[j]) << '\t';
//std::cout << '\n';

            double computeOverlapSumStart = omp_get_wtime();

            // make the comarison between the substrings(the suffix after the position of error
            TOccsIterator errorRead = begin(errorCandidates, Standard());
            TOccsIterator errorReadEnd = end(errorCandidates, Standard());
            for (; errorRead != errorReadEnd; ++errorRead)
            {
                unsigned errorReadId = (*errorRead).i1; //debug?
                unsigned fwdReadId = errorReadId;
                // swap to forward read ID because the allowed errors are only saved for the forward read IDs
                if (fwdReadId >= readCount)
                      fwdReadId -= readCount;

                //check here if the max number of corrections was made already 
                if (options.allowedCorrectionsPerRead[fwdReadId] == 0)
                    continue;

                unsigned positionError = (*errorRead).i2 + commonPrefix;

#ifdef FIONA_MAX_CORRECTIONS_PER_BASE
                lockReading(correctionListLock);
                unsigned foundCorr = getFoundCorrections(correctionList, firstCorrectionForRead, errorReadId, positionError, store);
                unlockReading(correctionListLock);
                if (foundCorr >= FIONA_MAX_CORRECTIONS_PER_BASE)
                    continue;
#endif

                TReadIterator itEBegin = begin(store.readSeqStore[errorReadId], Standard());
                TReadIterator itEPrefixBegin = itEBegin + (*errorRead).i2;
                TReadIterator itEEnd = end(store.readSeqStore[errorReadId], Standard());

    if (LOOP_LEVEL == 1)
        continue;

                for (unsigned c = 0; c < length(correctCandidates); ++c)
                {
                    TOccsIterator corrRead = begin(correctCandidates[c], Standard());
                    TOccsIterator corrReadEnd = end(correctCandidates[c], Standard());
                    
                /*
                Here should go the new Branch and Bound algorithm that
                looks if the current error pos has already a higher overlaps sum
                than what is possible to achieve with an SA infix from correctCandidates
                by considering the current strand. Although the Overlapsum might be
                computed before the for loop with c
                */

                    /////////////////////////////////////////////////////////////////////////////////
                    // reset overlap sums
                    for (int i = 0; i < options.maxIndelLength * 2 + 1; ++i)
                    {
                        Overlap &ov = bestCorrection[i];
                        ov.overlapSumLeft = 0;
                        ov.overlapSumRight = 0;
                    #ifdef FIONA_MAXIMIZE_SUPPORT
                        ov.support = 0;
                    #endif
                    #ifdef FIONA_CONSENSUS
                        clear(ov.consensus);
                        resize(ov.consensus, (itEEnd - itEPrefixBegin) - commonPrefix + options.maxIndelLength + 1);
                    #endif
                    }
                    
    if (LOOP_LEVEL == 2)
        continue;

//               bool debug = (errorReadId == (unsigned int)options.debugRead || errorReadId-length(store.readSeqStore)/2 == (unsigned int)options.debugRead);
//               if (debug)
//               {
//               std::cout<<"found"<<std::endl;
//               }

                    for (; corrRead != corrReadEnd; ++corrRead)
                    {					
                        /////////////////////////////////////////////////////////////////////////////////
                        // compare overlap left of the common prefix
                        //this part can be done without considering the type of indel
                        //if the left part has already too many errors no further investigation necessary

                        TReadIterator itE = itEBegin;
                        TReadIterator itCBegin = begin(store.readSeqStore[(*corrRead).i1], Standard());
                        TReadIterator itCLeft = itCBegin;
                        TReadIterator itCEnd = end(store.readSeqStore[(*corrRead).i1], Standard());

                        int delta = (*errorRead).i2 - (*corrRead).i2;
                        unsigned overlapLeft = positionError;

                        if (delta > 0)
                        {
                            // erroneous reads starts left of correct read
                            // EEEEEEEEEEEEEEE
                            //          CCCCCCCCCCCCCCC
                            overlapLeft -= delta;
                            itE += delta;
                        }
                        else
                        {
                            // erroneous reads starts right or with correct read
                            //         EEEEEEEEEEEEEEE
                            // CCCCCCCCCCCCCCC
                            itCLeft += -delta;
                        }
#ifdef FIONA_FIXED_OVERLAP_ERRORS
                        unsigned acceptedMismatchesLeft = maxAcceptedMismatches;
#else
                        // TOTAL NUMBER OF ERRORS IN THE OVERLAP
                        unsigned acceptedMismatchesLeft = std::max(2u, (unsigned)(options.overlap_errorrate * _min(itEEnd - itE, itCEnd - itCLeft)));
#endif
                        // COMPARE READ OVERLAPS WITH ERRORS
                    #ifdef FIONA_OVERLAP_WITH_EDIT_DISTANCE
                        // ALLOW INDELS

                        unsigned maxLen = _min((*errorRead).i2, (*corrRead).i2);
                        itE = itEPrefixBegin;
                        TReadIterator itC = itCBegin + (*corrRead).i2;

//                        if (debug)
//                        {
//                        #pragma omp critical
//                        {
//                        std::cout<< std::endl<< Range<TReadIterator>(itEBegin, itE) << std::endl;
//                        std::cout<< Range<TReadIterator>(itCBegin, itC) << std::endl<< std::endl;
//                        }
//                        }

                        for (; maxLen != 0; --maxLen)
                        {
                            --itC;
                            --itE;
                        #ifdef FIONA_MATCH_N
                            if (ordValue(*itE) == 4) continue;
                        #endif
                            if (*itE != *itC)
                                break;
                        }

                        unsigned errorsLeft =_comparePrefixesWithEditDistanceReverse(
                            Range<TReadIterator>(itEBegin, itE),
                            Range<TReadIterator>(itCBegin, itC),
                            acceptedMismatchesLeft);

                        // too many mismatches right of the common prefix?
                        if (acceptedMismatchesLeft < errorsLeft)
                            continue;
                        acceptedMismatchesLeft -= errorsLeft;
                    #else

                        for (; itE < itEPrefixBegin; ++itE, ++itCLeft)
                        {
                        #ifdef FIONA_MATCH_N
                            if (ordValue(*itE) == 4) continue;
                        #endif
                            if (*itE != *itCLeft)
                                if (--acceptedMismatchesLeft == std::numeric_limits<unsigned>::max()) break;
                        }
                        
                        // too many mismatches left of the common prefix?
                        if (acceptedMismatchesLeft == std::numeric_limits<unsigned>::max())
                            continue;
                    #endif

                        if (overlapLeft > (maxAcceptedMismatches - acceptedMismatchesLeft)) // kann weg (immer true)
                            overlapLeft -= (maxAcceptedMismatches - acceptedMismatchesLeft);

                        // the position in the read until which there is the same prefix
                        unsigned positionCorrect = (*corrRead).i2 + commonPrefix;


    if (LOOP_LEVEL == 3)
        continue;
                    #ifdef FIONA_MAXIMIZE_SUPPORT
                        unsigned overallMinErrorsRight = acceptedMismatchesLeft + 1;
                    #endif

                        /////////////////////////////////////////////////////////////////////////////////
                        // compare overlap right of the common prefix by considering all types of allowed errors
                        for (int indel = -options.maxIndelLength; indel <= options.maxIndelLength; ++indel)
                        {
                            Overlap &overlap = bestCorrection[((indel >= 0)? indel * 2: -indel * 2 - 1)];
                        #ifdef FIONA_MAXIMIZE_SUPPORT
                            overlap._errorsRight = acceptedMismatchesLeft + 2;  // initialize with "infinity"
                        #endif
                            itE = itEBegin + positionError;
                            TReadIterator itC = itCBegin + positionCorrect;

                            if (indel == 0)
                            {
                                // mismatch
                                ++itE;
                                ++itC;
                            }
                            else if (indel > 0)
                            {
                                // gap in correct read (must be deleted in err. read)
                                itE += indel;
                                if (itE + indel >= itEEnd || itC + indel >= itCEnd) continue;
                            }
                            else
                            {
                                // gap in erroneous read (must be filled with an insertion in err. read)
                                itC += -indel;
                                if (itC + -indel >= itCEnd || itC + -indel >= itCEnd) continue;
                            }
                            
    if (LOOP_LEVEL == 4)
        continue;

                            // COMPARE READ OVERLAPS WITH ERRORS
                            unsigned rightOverlapLen = _min(itEEnd - itE, itCEnd - itC);
                            unsigned acceptedMismatches = acceptedMismatchesLeft;
                            unsigned errorsRight;
                        #ifdef FIONA_OVERLAP_WITH_EDIT_DISTANCE
                            // ALLOW INDELS

                            // first extend seed without errors
                            for (unsigned maxLen = rightOverlapLen; maxLen != 0; --maxLen, ++itE, ++itC)
                            {
                            #ifdef FIONA_MATCH_N
                                if (ordValue(*itE) == 4) continue;
                            #endif
                                if (*itE != *itC)
                                    break;
                            }

                            // now continue with banded Myers
                            errorsRight =_comparePrefixesWithEditDistance(
                                Range<TReadIterator>(itE, itEEnd),
                                Range<TReadIterator>(itC, itCEnd),
                                acceptedMismatches);

                            // too many mismatches right of the common prefix?
                            if (acceptedMismatches < errorsRight)
                                continue;
                            acceptedMismatches -= errorsRight;
                        #else
                            // ALLOW ONLY MISMATCHES

                            // manually count mismatches
                            for (unsigned maxLen = rightOverlapLen; maxLen != 0; --maxLen, ++itE, ++itC)
                            {
                            #ifdef FIONA_MATCH_N
                                if (ordValue(*itE) == 4) continue;
                            #endif
                                if (*itE != *itC)
                                    if (--acceptedMismatches == std::numeric_limits<unsigned>::max())
                                        break;
                            }

                            // too many mismatches right of the common prefix?
                            if (acceptedMismatches == std::numeric_limits<unsigned>::max())
                                continue;

                            errorsRight = acceptedMismatchesLeft - acceptedMismatches;
                        #endif

                        #ifdef FIONA_MAXIMIZE_SUPPORT
                            overlap._errorsRight = errorsRight;         // store number errors for this operation
                            if (overallMinErrorsRight > errorsRight)    // update minimum over all operations
                                overallMinErrorsRight = errorsRight;
                        #else
                            ignoreUnusedVariableWarning(errorsRight);
                        #endif

//               bool debug = (errorReadId == (unsigned int)options.debugRead || errorReadId-length(store.readSeqStore)/2 == (unsigned int)options.debugRead);
//                            if(debug)
//                            {
//                                 std::cerr << "error_read     \t" << indel << '\t' << suffix(store.readSeqStore[errorReadId],positionError) << '\t' << errorReadId << '\t' << (acceptedMismatchesLeft-acceptedMismatches) << '\n';
//                                 std::cerr << "correct_read   \t" << indel << '\t' << suffix(store.readSeqStore[(*corrRead).i1],positionCorrect) << '\t' << (*corrRead).i1 << '\n';
//                                 //std::cerr << store.readSeqStore[(*corrRead).i1] << std::endl;
//                            }

                            // correction candidate
                            overlap.overlapSumLeft += overlapLeft;
                            int overlapRight = rightOverlapLen;
                            if (overlapRight + _min(indel,0) >= 0)
                                overlapRight += _min(indel,0); // hier immer indel laenge abziehen damit Mismatches nicht benachteiligt werden
                            if (indel == 0)
                                ++overlapRight;
                            if (overlapRight >= (int)(acceptedMismatchesLeft - acceptedMismatches))
                                overlapRight -= (int)(acceptedMismatchesLeft - acceptedMismatches);
                            overlap.overlapSumRight += overlapRight;

                            overlap.readId = (*corrRead).i1;
                            overlap.correctPos = positionCorrect;

//                            if (debug)
//                            {
//                                std::cerr << "overlap      \t" << overlapLeft << '+' << overlapRight << '\t' << indel << std::endl;
//                            }

                        #ifdef FIONA_CONSENSUS
                            // (only for mismatches/insertions) compute consensus of correct reads
                            // that fulfill the acceptedMismatch-critereon
                            if (indel <= 0)
                            {
                                // reset itC and repeat comparison to increase consensus counters
                                typename Iterator<Overlap::TConsensus, Standard>::Type itCons = begin(overlap.consensus, Standard());

                                // there are rightOverlapLen bases right of the error
                                itC = itCLeft + commonPrefix;
                                itCEnd = itC + rightOverlapLen;

                                // take the error (mismatch/gap in error read) into account for consensus computation
                                if (indel == 0)
                                    ++itCEnd;
                                else if (indel < 0)
                                    itCEnd += -indel;
                                
                                for (; itC < itCEnd; ++itC, ++itCons)
                                    ++(*itCons).count[ordValue(*itC)];
                                SEQAN_ASSERT_LEQ(itCons, end(overlap.consensus, Standard()));
                            }
                        #endif
                        }

                    #ifdef FIONA_MAXIMIZE_SUPPORT
                        // update support for all operations
                        for (int i = 0; i < options.maxIndelLength * 2 + 1; ++i)
                        {
                            Overlap &overlap = bestCorrection[i];
                            if (overlap._errorsRight == overallMinErrorsRight)  // increase support for operation with minimal errors
                                ++overlap.support;
                        }
                    #endif
                    }


                #ifdef FIONA_MAXIMIZE_SUPPORT
                    // determine operation with maximal support
                    unsigned bestOperation = 0;
                    unsigned maxSupport = bestCorrection[0].support;
                    for (int i = 1; i < options.maxIndelLength * 2 + 1; ++i)
                    {
                        Overlap &overlap = bestCorrection[i];
                        if (maxSupport < overlap.support)
                        {
                            maxSupport = overlap.support;
                            bestOperation = i;
                        }
                    }
                    if (maxSupport == 0)
                        continue;
                #endif

                    // try finding the best correction with the highest overlap sum
                    // here is where we start to differ as we record the Corrections in the new
                    // new linked list (correctionList) with CorrectionIndelPos structs
                    // but instead of choosing the Correction with the highest Overlapsum from each error type
                    // we treat each of theses cases by using addCorrectionEntry
                    bool strand = (errorReadId >= readCount);

                #ifdef FIONA_MAXIMIZE_SUPPORT
                    int i = bestOperation;
                #else
                    for (int i = 0; i < options.maxIndelLength * 2 + 1; ++i)
                #endif
                    {
                        unsigned short overlapSum = bestCorrection[i].overlapSumLeft + bestCorrection[i].overlapSumRight;

                        if (overlapSum == 0)
                            continue;

                        //create variables to use addCorrectionEntry this is critical for OMP
                        signed char indel = ((i & 1) == 0)? i / 2: -((i + 1) / 2);

#ifdef FIONA_CONSENSUS
                        TReadIterator itE = itEPrefixBegin + commonPrefix;
                        typename Iterator<Overlap::TConsensus, Standard>::Type itCons = begin(bestCorrection[i].consensus, Standard());

                        if (indel == 0)     // mismatch
                        {
                            // extract first consensus base
                            correctSeq[0] = *itCons;
                            if (strand)
                                correctSeq[0] = FunctorComplement<Dna5>()(correctSeq[0]);
                            ++itE;
                            ++itCons;
                        }
                        else if (indel > 0) // gap in correct read (must be deleted in err. read)
                        {
                            itE += indel;
                        }
                        else                // gap in erroneous read (must be filled with an insertion in err. read)
                        {
                            SEQAN_ASSERT_LT(indel, 0);
                            // extract consensus bases to determine insert
                            if (strand)
                            {
                                for (int l = -indel; l > 0; ++itCons)
                                    correctSeq[--l] = FunctorComplement<Dna5>()((Dna5)*itCons);
                            }
                            else
                            {
                                for (int l = 0; l < -indel; ++l, ++itCons)
                                    correctSeq[l] = *itCons;
                            }
                        }
                        SEQAN_ASSERT_LEQ(itCons, end(bestCorrection[i].consensus, Standard()));

                 /*      SEQAN_OMP_PRAGMA(critical(TestConsensusOverlapsum))
                                    {
                                           std::cout << "normCorrect: " << errorReadId<< " "<< positionError<< " " << overlapSum<< " " <<  options.overlapSumCutoffs[positionError+1]<< " " << cycle<< std::endl;
                                    } */

                        // 1. add major mismatch/indel correction
                        SEQAN_OMP_PRAGMA(critical(addCorrection))
                        {
                            ++resources.putCorrections;
                            addCorrectionEntry(
                                correctionList,
                                firstCorrectionForRead,
                                errorReadId,
                                bestCorrection[i].readId,
                                bestCorrection[i].correctPos,
                                positionError,
                                overlapSum,
                                strand,
                                indel,
                                store,
                                correctSeq);
                        }

                        if (indel == 0)
                        {
                            // 2. add minor consensus mismatch corrections
                            //bool debugConsensus = false;
                            for (; itE < itEEnd; ++itE, ++itCons)
                            {
                                unsigned maxBase = getMaxIndex(*itCons);
                                if (maxBase == 4) continue;     // skip if N is the consensus base
                                unsigned frequency = (*itCons).count[maxBase];
                                if (frequency < 2) break;       // at least 2 suffixes need to vote for that base
                                if ((*itCons).count[ordValue(*itE)] < frequency)
                                {
                                    // TODO:
                                    // (weese:) I'm not sure if we better scale the overlap-sum relative to the coverage at the anchor (front(consensus))
                                    //          or relative to the coverage at the current base (*itCons)
                                    unsigned totalCounts = totalCount(*itCons);
                                    if(totalCounts == 0)
                                         continue;
                                  //totalCount(*itCons /*front(consensus)*/));
                                    unsigned consOverlapSum = (unsigned)((frequency * overlapSum) / totalCounts);
				   if(consOverlapSum >1) //penalize consensus correction by one to always to major correction first
					--consOverlapSum; 
                                    //unsigned testOverlapSum = (unsigned)((frequency * overlapSum) / totalCount(front(itCons)));
                                  /*  SEQAN_OMP_PRAGMA(critical(TestConsensusOverlapsum))
                                    {
                                           std::cout << "overCorrect: " << errorReadId<< " "<< maxBase << " " <<itE-itEBegin<< " " << overlapSum << " " << consOverlapSum << " "   << frequency <<" overcutoff: " << options.overlapSumCutoffs[itE-itEBegin+1]<< " " << cycle << std::endl;
                                    } */
                                    correctSeq[0] = (Dna5)maxBase;
                                    if (strand)
                                        correctSeq[0] = FunctorComplement<Dna5>()(correctSeq[0]);
                                    SEQAN_OMP_PRAGMA(critical(addCorrection))
                                    {
                                        ++resources.putCorrections;
                                        addCorrectionEntry(correctionList,firstCorrectionForRead,errorReadId,bestCorrection[i].readId,bestCorrection[i].correctPos,itE-itEBegin,consOverlapSum,strand,indel,store,correctSeq);
                                    }

                                    //std::cout << "replace " << *itE << " at position " << (itE - itEBegin) << " in read " << errorReadId;
                                    //std::cout << " by " << (Dna5)maxBase << " (support=" << (*itCons).count[maxBase] << ")" << std::endl;
                                    //debugConsensus = true;
                                }
                            }//go over all positions in consensus
                        }//no indel
/*
                        if (debugConsensus)
                        {
                            for(int divider = 10; divider != 0; divider /= 10)
                            {
                                for(int xx=0; xx<positionError+1+indel; ++xx) std::cout << ' ';
                                for(unsigned xa=0; xa<length(bestCorrection[i].consensus); ++xa)
                                {
                                    unsigned maxBase = getMaxIndex(bestCorrection[i].consensus[xa]);
                                    if (bestCorrection[i].consensus[xa].count[maxBase] == 0) break;
                                    std::cout << (bestCorrection[i].consensus[xa].count[maxBase] / divider) % 10;
                                }
                                std::cout << std::endl;
                            }


                            for(int xx=0; xx<indel; ++xx) std::cout << ' ';
                            std::cout << prefix(store.readSeqStore[errorReadId], positionError) << ' ';
                            for(unsigned xa=0; xa<length(bestCorrection[i].consensus); ++xa)
                            {
                                unsigned maxBase = getMaxIndex(bestCorrection[i].consensus[xa]);
                                if (bestCorrection[i].consensus[xa].count[maxBase] == 0) break;
                                std::cout << (Dna5)maxBase;
                            }
                            std::cout << std::endl;

                            std::cout << prefix(store.readSeqStore[errorReadId], positionError) << ' ';
                            for(int xy=0; xy<-indel; ++xy) std::cout << ' ';
                            std::cout << suffix(store.readSeqStore[errorReadId], positionError + _max(0, -indel)) << std::endl;
                            std::cout << std::endl;
                        }
*/
                       
#else // FIONA_CONSENSUS
                        if (indel <= 0) //only get string if insertion in read or mismatch
                            getCorrectionString(correctSeq,indel,bestCorrection[i].readId,bestCorrection[i].correctPos,strand,store);
                        SEQAN_OMP_PRAGMA(critical(addCorrection))
                        {
             /*std::ofstream myfile;
myfile.open ("EntryOutput.txt",std::ios::app);
myfile << "beforeEntry " << " "<< errorReadId<< " "<<bestCorrection[i].readId<< " "<<bestCorrection[i].correctPos<< " "<<positionError<< " "<<(int)indel << " " << overlapSum << " " <<representative(iter) <<"endEntry" << std::endl;
myfile.close();*/
                            ++resources.putCorrections;
                            addCorrectionEntry(correctionList,firstCorrectionForRead,errorReadId,bestCorrection[i].readId,bestCorrection[i].correctPos,positionError,overlapSum,strand,indel,store,correctSeq);
                        }
#endif // FIONA_CONSENSUS
                    }

                }//finished errorCandidate

            }//finished analysis error read

            options.timeComputeOverlapSum += (omp_get_wtime() - computeOverlapSumStart);
        } // if (!empty(correctCandidates))

        if (leaveNodeRight)
            goNextRight(iter);
        else
            goNext(iter);
	}
    resources.cpuTime = omp_get_wtime() - start;
}

// CorrectionIndelPos  


/*GC-content*/
/*fonction which allow to determine the frequency for each nucleotide*/
template < typename TFionaIndex, typename TSpec>
String<double, Array<5> >
determineFrequency(Iter< TFionaIndex, VSTree<TSpec> > iter)
{
    /*calculate the frequency for each nucleotide*/
    /*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
    String<double, Array<5> > frequency;
    resize(frequency, 5, 0);

    goBegin(iter);

//    int position = 0;
    /*nombre total nucleotides*/
    double total = countOccurrences(iter);
    /*the first is A (alphabetical ordre)*/
    goDown(iter);
    do {
        frequency[ordValue(parentEdgeFirstChar(iter))] = countOccurrences(iter) / total;
    } while (goRight(iter));

    /*table of frequency for each nucleotide*/
    return frequency;
}

/*construction Suffix Array */
template <typename TFragmentStore, typename TAlgorithm>
unsigned correctReads(
	TFragmentStore & store,
	FionaOptions & options,
	Tag<TAlgorithm> const alg)
{
	/*iterator with restrictions*/
	typedef Iterator<TFionaIndex, TopDown<ParentLinks<Preorder> > >::Type TConstrainedIterator;

	// append their reverse complements
	unsigned readCount = length(store.readSeqStore);
    if (options.verbosity >= 2)
        std::cerr << "Add reverse complements: " << std::flush;

    Dna5String tmp;
    //if (options.genomeLength != 1)
    //(Hugues) we need to get the reads and their reverse complements 
    //when estimating the error rates (but the letter substitutions should be reversed).
    SEQAN_PROTIMESTART(timeRevComp);
    for (unsigned i = 0; i < readCount; ++i)
    {
        tmp = store.readSeqStore[i];
        reverseComplement(tmp, Serial());
        appendValue(store.readSeqStore, tmp);
    }

    if (options.verbosity >= 2)
        std::cerr << SEQAN_PROTIMEDIFF(timeRevComp) << " seconds." << std::endl;

    /*table with the theoretical values*/

    String<uint64_t> readLengthHist;
    computeReadLengthHistogram(readLengthHist, store.readSeqStore);
    unsigned maxReadLength = length(readLengthHist) - 1;
    if (options.verbosity >= 2)
        std::cerr << "Maximal read length:" << maxReadLength << std::endl;

    double experrreads = expectedValueTheoretical(options.expectedTheoretical, readLengthHist, options.genomeLength, options.errorrate);

    if (IsSameType<TAlgorithm, FionaExpected_>::VALUE)
    {
        String<double> sd;
        standardDeviation(sd, store.readSeqStore, options.genomeLength);

        /*The strictness value allows one to estimate the confidence intervall*/
        for (unsigned i = 0; i < length(options.expectedTheoretical); ++i)
        {
            double expectedTemporary = options.expectedTheoretical[i] - options.strictness * sd[i];

            /*If the connfidential intervall take value less than 1 ??? not sure for that*/
            /*if(expectedTemporary < 1){
                options.expectedTheoretical[i] = 1.1;
            }else{*/
                options.expectedTheoretical[i] = expectedTemporary;
            //}
        }
    }

    if (IsSameType<TAlgorithm, FionaCount_>::VALUE)
        for (unsigned i = 0; i < length(options.expectedTheoretical); ++i)
            options.expectedTheoretical[i] = options.strictness;

    // std::cerr << " run the multiple-Correction-per-Round Fiona method " <<std::endl;
    // the linked list for read corrections
    String<CorrectionIndelPos> correctionList;
    // the first correction occurrence (if any) for a read
    // if no occurrence exists set the enrty to maxINt Value
    String<unsigned int> firstCorrectionForRead; // should this really be always created anew for different cycles, could be created outside correctReads?
    // we assume we work with three reads here
    resize(firstCorrectionForRead, readCount, std::numeric_limits<unsigned>::max(), Exact());

    // Determine the number of allowed corrections per round per read, depending on the read length and the
    // configuration in options.relativeErrorsToCorrect.  We set a hard lower limit of 2.
    unsigned const MIN_ALLOWED_CORRECTIONS = 2;
    if (options.limitCorrPerRound || length(options.allowedCorrectionsPerRead) == 0u)
    {
        resize(options.allowedCorrectionsPerRead, readCount, Exact());
        for (unsigned i = 0; i < readCount; ++i)
        {
            unsigned readLength = length(store.readSeqStore[i]);
            options.allowedCorrectionsPerRead[i] = _max((unsigned) ceil(options.relativeErrorsToCorrect * readLength), MIN_ALLOWED_CORRECTIONS);
        }
    }

    /*
    //copy Hugues code to get the allowed errors
	binomial Nmismatch(maxReadLength, options.errorrate);
    unsigned maxAcceptedMismatches = (unsigned) ceil(quantile(Nmismatch, 0.95));
    if (maxAcceptedMismatches < 2)
        maxAcceptedMismatches = 2; //unlikely that 2 reads share the same error
    if (options.verbosity >= 3)
        std::cerr << "Max number of corrections in a read: " << maxAcceptedMismatches << std::endl;
    if (options.limitCorrPerRound)
        clear(options.allowedCorrectionsPerRead);
    resize(options.allowedCorrectionsPerRead, readCount, maxAcceptedMismatches, Exact());
    */
	
    if (options.verbosity >= 2)
    {
        if (IsSameType<TAlgorithm, FionaExpected_>::VALUE)
            std::cerr << std::endl << "Method with expected value for each level" << std::endl;
        if (IsSameType<TAlgorithm, FionaPoisson_>::VALUE) 
            std::cerr << std::endl << "Method with p-value and Poisson distribution" << std::endl;
        if (IsSameType<TAlgorithm, FionaPoissonSens_>::VALUE)
            std::cerr << std::endl << "Method with sensitivity and Poisson distribution" << std::endl;
        if (IsSameType<TAlgorithm, FionaCount_>::VALUE)
            std::cerr << std::endl << "Method with fixed count for each level" << std::endl;
        if (IsSameType<TAlgorithm, FionaPoissonClassif_>::VALUE)
            std::cerr << std::endl << "Log-odds method assuming a Poisson coverage distribution" << std::endl;
    }
	
	///Give the set of cut-off for the various methods
	double oddserrreads = (experrreads /  (readCount - experrreads));
	options.oddserrorreads = oddserrreads;
	if (options.verbosity >= 2 && experrreads != 0)
    {
		std::cerr << "Error rate provided: " << options.errorrate << std::endl;
		std::cerr << "Expected number of erroneous reads (percent): " << floor(experrreads) << " ("  << floor(experrreads/readCount *100) << ")" << std::endl;
		std::cerr << "Odds of erroneous/correct reads: " << oddserrreads << std::endl;
	}
	
	if (options.autolevel){
        if (options.verbosity >= 2)
            std::cerr << "Setting automatic from/to level with number of Correctables/Uncorrectables reads" << std::endl;
		//std::cerr << "Computing the expected number of uncorrectable (Uk) and Destructible (Dk) reads" << std::endl;
		//
        // TODO(holtgrew): This logic is broken.
//		float errrate = options.errorrate;
		if (options.verbosity >= 2 && options.errorrate == 0)
        {
			std::cerr << "Error rate not given, set to default value of 1%" << std::endl;
//			errrate = 0.01;
		}
		float minexpcov = 5.; //ask for mincov of 5, we could just do + 10
		int toplevel = (int) (maxReadLength - minexpcov * options.genomeLength / readCount) - 1;
		int mink = 5;
		int maxk = std::min((int)maxReadLength - 2, 50);
		toplevel = (toplevel > (int)maxReadLength) ? (int)maxReadLength : toplevel;

        String<double> uncorrectables;
        String<double> destructibles;

		UncorrectableExpectedBases(uncorrectables, mink, maxk, readLengthHist, options.errorrate);
		DestructibleExpectedBases(destructibles, mink, maxk, readLengthHist, options.errorrate, options.genomeLength);

        if (options.verbosity >= 2)
        {
            std::cerr << "max. read length: " << maxReadLength << std::endl;
            std::cerr << "k\tUncorrectable\tDestructible\tUc+Dc" << std::endl;
        }
		float TmpExp = uncorrectables[mink] + destructibles[mink];
		unsigned HiTEC_mink = mink;
		for (int i = mink + 1; i <= maxk; i++)
        {
			double tmp = uncorrectables[i] + destructibles[i];
            
            if (options.verbosity >= 2)
                std::cerr << i << "\t" << uncorrectables[i] << "\t" << destructibles[i] << "\t" << tmp << std::endl;
            
			if (tmp < TmpExp)
            {
				HiTEC_mink = i;
				TmpExp = tmp; 
			}
		}
        if (options.verbosity >= 1)
            std::cerr << "Determination of minimum tree level according to HiTEC strategy." << std::endl;
        double mink_genome = log(200.0 * options.genomeLength) / log(4.0);
		options.fromLevel = (HiTEC_mink < mink_genome) ? HiTEC_mink : (int)mink_genome;
        
		int upl = (options.fromLevel + 10) < (int)maxReadLength ? (options.fromLevel + 10) : (int)maxReadLength;
		options.toLevel = upl; // OLD let to high toplevels: toplevel < upl ? upl : toplevel;
        if (options.verbosity >= 1)
            std::cerr << "The estimated top level is " << options.fromLevel << " and the down level is " << options.toLevel << std::endl;
	}
	
    if (options.verbosity >= 2)
    {
        std::cerr << "Expected coverage of k-mers before sequencing (k, coverage):" << std::endl;
        for (int i = options.fromLevel; i<= options.toLevel; i++){
            std::cerr << "(" <<  i << " , " << options.expectedTheoretical[i] << ")  ";	
        }
        std::cerr << std::endl;
    }
	
	ComputeCutoffErroneous(options.errorCutoffs, options.expectedTheoretical, options.strictness, options.errorrate, options.oddserrorreads, options.fromLevel, options.toLevel + 1, alg);
    if (options.verbosity >= 2)
    {
        std::cerr << "Computed cutoffs for errors (level, cutoff):" << std::endl;
        for (int i = options.fromLevel; i<= options.toLevel; i++)
            std::cerr << "(" <<  i << " , " << options.errorCutoffs[i] << ")  ";
        std::cerr << std::endl;
    }
	//(Hugues) shouldn't we increment all automatics cutoffs by 1 for errors of one value (sensitivity + oddsratio ?) 
	
	if (options.errorrate != 0)
    {
		computeCutoffRepeats(options.repeatCutoffs, 1, options.expectedTheoretical, options.errorrate, options.fromLevel, options.toLevel + 1, options.genomeLength);
        if (options.verbosity >= 2)
        {
            std::cerr << "Computed cutoffs for repeats (level, cutoff):" << std::endl;
            for (int i = options.fromLevel; i<= options.toLevel; i++){
                std::cerr << "(" <<  i << " , " << options.repeatCutoffs[i] << ")  ";
            }
            std::cerr << std::endl;
        }
	}
    else
    {
        resize(options.repeatCutoffs, options.toLevel + 2, std::numeric_limits<unsigned>::max());
    }

    if (options.verbosity >= 2)
        std::cerr << std::endl;

    if (options.errorrate != 0 && options.wovsum != 0)
    {
        SEQAN_PROTIMESTART(timeCutoffComp);
        ComputeCutoffOverlapSum(options.overlapSumCutoffs, options.fromLevel, readLengthHist, options.genomeLength, options);
        /*int hack[100] = {  876, 876, 873, 870, 866, 863, 859, 856, 852, 848, 843, 839, 834, 830, 842, 855, 867, 878, 890, 901, 911, 922, 932, 941, 950, 959, 968, 976, 984, 991, 998, 1005, 1011, 1017, 1023, 1028, 1033, 1038, 1042, 1046, 1049, 1053, 1055, 1058, 1060, 1062, 1063, 1064, 1065, 1065, 1065, 1065, 1064, 1063, 1062, 1060, 1058, 1055, 1053, 1049, 1046, 1042, 1038, 1033, 1028, 1023, 1017, 1011, 1005, 998, 991, 984, 976, 968, 959, 950, 941, 932, 922, 911, 901, 890, 878, 867, 855, 842, 830, 834, 839, 843, 848, 852, 856, 859, 863, 866, 870, 873, 876, 8763 };
        for (unsigned i=0;i<100;++i)
            options.overlapSumCutoffs(100, i) = hack[i];
        */
        if (options.verbosity >= 2)
        {
            unsigned len = _min(100u, maxReadLength);
            std::cerr << "Computed cutoffs for nb overlap bp (pos, cutoff)" << std::endl;
            for (unsigned i = 0; i < len; ++i)
                std::cerr << "(" << i + 1 << " , " << options.overlapSumCutoffs(len, i) << ")  ";
            std::cerr << std::endl;
            std::cerr << "Time required for cutoffs computation: " << SEQAN_PROTIMEDIFF(timeCutoffComp) << " seconds." << std::endl;
        }
	}
    else
    {
        options.overlapSumCutoffs.resize(length(readLengthHist), length(readLengthHist) - 1);  // readLength, position
        for (unsigned readLen = 0; readLen < length(readLengthHist); ++readLen)
            for (unsigned pos = 0; pos < readLen; ++pos)
                options.overlapSumCutoffs(readLen, pos) = 3; //default parameter is 3, that is 3 bp are needed to correct.
    }
	
//	std::ofstream pdist("proba_dist.txt");
	
//	pdist << "count";
//	for (int i = options.fromLevel ; i<= options.toLevel; i++) 
//			pdist << "\t" << i << ".expected.cov\t" <<  i << ".dpois\t" <<  i << ".dpoismix\t" << i << ".ppoismix\t";
//	pdist << std::endl;
//	for (int k = 0; k <= 50; k++){
//		pdist << k ;
//		for (int i = options.fromLevel; i<= options.toLevel; i++){
//			double noerr = probabilityNoError(options.errorrate, i);
//			pdist << "\t" << options.expectedTheoretical[i] << "\t" << dpois(k, options.expectedTheoretical[i] * noerr) << "\t" << dpoismixerror(k, options.expectedTheoretical[i], options.errorrate, i) << "\t" << 1-ppoismixerror(k, options.expectedTheoretical[i], options.errorrate, i) ;
//		}
//		pdist << std::endl;
//	}
//	pdist.close();

	//exit(0);

    options.timeComputeOverlapSum = 0;
    if (options.verbosity >= 1)
        std::cerr << "Searching..." << std::endl;
	SEQAN_PROTIMESTART(search);

#ifndef FIONA_PARALLEL
	// FIONA NON-PARALLEL SEARCH
		
	// construct suffix array of the set of reads
    if (options.verbosity >= 1)
        std::cerr << "Construct suffix array" << std::endl;
	SEQAN_PROTIMESTART(construct);
	TFionaIndex myIndex(store.readSeqStore);
	TConstrainedIterator myConstrainedIterator(myIndex);

	/*calculate the frequency for each nucleotide, didn't use for the moment*/
	/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
//	String<double, Array<5> > frequency = determineFrequency(myConstrainedIterator);

    if (options.verbosity >= 1)
        std::cerr << "Time required for suffix array construction : " << SEQAN_PROTIMEDIFF(construct) << " seconds." << std::endl;

	/*restrictions just for estimating the genome length if there is no data*/

#ifdef MEDIAN
	unsigned level = fromLevel;
	ofstream out("medianLevels.txt"); 
	ofstream median("medianForEachLevel.txt");
#endif // #ifdef MEDIAN

	if (options.genomeLength == 1)
	{
        if (options.verbosity >= 1)
        {
            std::cerr << "Generating Hugues' stats file." << std::endl;
            std::cerr << "Between levels " << options.fromLevel << " and " << options.toLevel << std::endl;
        }
        std::ofstream stats("stats.txt");
		Iterator<TFionaIndex, TopDown<ParentLinks<Preorder> > >::Type it(myIndex);
		//change the iterator definition here
		cargo(myIndex).replen_min = options.fromLevel;
		cargo(myIndex).replen_max = options.toLevel; 	
		//cargo(myIndex).frequency = frequency;
		TConstrainedIterator myItStat(myIndex);
		goBegin(myItStat);
		CharString tmp;
		//get the length of the reads.
		unsigned readMaxLength = 0;
		for (unsigned i = 0; i < readCount; ++i)
		{
			unsigned readLength = length(store.readSeqStore[i]);
			readMaxLength = (readMaxLength < readLength) ? readLength : readMaxLength;
		}
        if (options.verbosity >= 1)
            std::cerr << "Max observed read length : " << readMaxLength << std::endl;
		//String<int> freq;
		std::map < unsigned, std::vector <int> > freqpos; //
		for (unsigned i = 0; i < 6; i++) freqpos[i].assign(readMaxLength, 0);
		std::vector <int> freqmarginal(5, 0 );
		
		
		//stats << "branch\tlength\ttreeDep\tletter\treadPos\tfreq" << std::endl;
		stats << "prefix\ttreeDepth\treadPos\tna\tnc\tng\tnt\tnn\tntot\tnfather" << std::endl;
		
		while (!atEnd(myItStat))
		{
			//unsigned ofs = parentRepLength(myItStat);
			//tmp = parentEdgeLabel(myItStat);
			CharString tmp2 = representative(myItStat);
			unsigned cfather = countOccurrences(myItStat);
			if (parentRepLength(myItStat) == 0) { 
				goNext(myItStat);
				continue;} 
			//std::cerr << "******* Now in prefix : " << toCString(tmp) << " of the string: " << tmp2 << std::endl; 
			if (isLeaf(myItStat)){ 
				//	std::cerr << "which is a leaf, next one" << std::endl;
				goNext(myItStat);
				continue;
			}
			//copy the iterator to get all siblings
			Iterator<TFionaIndex, TopDown<ParentLinks<Preorder> > >::Type it2(myItStat);
			goDown(it2);
			//goDown(it2);
			unsigned let = 0;
			//count each position and read.
			do {
				//Get count per letter per pos
				/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3 and 'N'= 4*/
				//WE NEED TO GET THE VALUE OF LET HERE
				for (unsigned j = 0; j < countOccurrences(it2); ++j)
				{
					unsigned readID = getOccurrences(it2)[j].i1 ; 
					unsigned rl = length(store.readSeqStore[readID]);
					unsigned thei2 = getOccurrences(it2)[j].i2; //check i2
					unsigned posInRead = (readID < readCount) ? thei2 : rl - thei2 - 1; 
					//unsigned posInRead =  0;  //  ofs + length(tmp); //changer ici, que vaut i2 ?
					//std::cerr << "Looking at sequence repres: " << representative(it2) << std::endl ; 
					//std::cerr << "we have length=" << rl << " for read number " << readID << std::endl;
					//std::cerr << " i2=" << thei2 << " or with dir access " << getOccurrences(it2)[j].i2 << std::endl;
					//std::cerr <<"letter:" << let << " and with ofs "	<< ofs << " and length(tmp): " << length(tmp) << std::endl;
					//if (length(freqpos[let]) <= posInRead)
					//	freqpos[let][posInRead + 1] = 0;
					++freqpos[let][posInRead];
					
					//cl = cl < posInRead ? posInRead : cl;
				}
				//we should check the letter here  
				++let;
			}  while( goRight(it2));
			
			freqmarginal.assign(6,0);
			for (unsigned i= 0; i < readMaxLength; ++i)
				for (int cl=0; cl < 5; ++cl){
					freqpos[5][i] += freqpos[cl][i]; 
					freqmarginal[cl] += freqpos[cl][i];
				}
			for (int cl=0; cl < 5; ++cl)
				freqmarginal[5] += freqmarginal[cl]; 
			//also log the marginal, easier
			stats << tmp2 << '\t' << nodeDepth(myItStat) << '\t' << -1 << '\t';
			stats << freqmarginal[0] << '\t' << freqmarginal[1] << '\t'; 
			stats << freqmarginal[2] << '\t' << freqmarginal[3] << '\t' << freqmarginal[4] << '\t' ;
			stats << freqmarginal[5] << '\t' << cfather <<  std::endl;	
			
			
			for (unsigned cpos = 0; cpos < readMaxLength; ++cpos){
				
				if (freqpos[5][cpos] > 0) {
					stats << tmp2 << '\t' << nodeDepth(myItStat) << '\t' << cpos << '\t';
					stats << freqpos[0][cpos] << '\t' << freqpos[1][cpos] << '\t'; 
					stats << freqpos[2][cpos] << '\t' << freqpos[3][cpos] << '\t' << freqpos[4][cpos] << '\t' ;
					stats << freqpos[5][cpos] << '\t' << cfather <<  std::endl;	
				}
			}
			
			for (unsigned i = 0; i < 6; i++) freqpos[i].assign(readMaxLength, 0);
			
			goNext(myItStat);		
		}
		stats.close();
        if (options.verbosity >= 1)
            std::cerr << "  Done." << std::endl;
		exit(0);
	}
			
	if (options.genomeLength == 0)
	{
#ifdef MEDIAN		
		for (; level < toLevel; ++level)
		{
			//int logRation = (log10(static_cast<double>(length(setReads)/2))/(log10(4.0)));
			//int l = logRation + 1;
			//std::cerr << l << std::endl;
			
			cargo(myIndex).replen_min = level;
			cargo(myIndex).replen_max = level+2; 	
//			cargo(myIndex).frequency = frequency;
			double numOccs = 0.0;
			
			median << level << " " << medianLevel(myConstrainedIterator) << std::endl;
			goBegin(myConstrainedIterator);
			while (!atEnd(myConstrainedIterator))
			{
				if (parentRepLength(myConstrainedIterator) > level)
				{
					numOccs = countOccurrences(myConstrainedIterator);
					out << level << " " << numOccs << std::endl; 
				}
				++myConstrainedIterator;
			}
		}
#else // #ifdef MEDIAN
		//int logRation = log10(static_cast<double>(readCount)) / log10(4.0);
		//int l = logRation + 1;
		//std::cerr << l << std::endl;
		cargo(myIndex).replen_min = options.fromLevel;
		cargo(myIndex).replen_max = options.fromLevel + 2; 	
//		cargo(myIndex).frequency = frequency;

		double expectedValueGivenLevel = medianLevel(myConstrainedIterator);

        // Compute mean read length as an estimate.  This will be the read length for Illumina data and for now
        // a good enough value for 454 data.
        // TODO: Think of something more clever in the future.
		uint64_t readLengthSum = 0;
		for (unsigned i = 0; i < readCount; ++i)
		    readLengthSum += length(store.readSeqStore[i]);
		unsigned readLength = readLengthSum / readCount;
        if (options.verbosity >= 1)
            std::cerr << "Average read length " << readLength << "\n";

		/* a = readLength - path_label + 1 */
		/*here plus 1 also because the level is between fromLevel and toLevel*/
		double a = readLength - options.fromLevel + 2;
		options.genomeLength = static_cast<int64_t>(readCount * a / expectedValueGivenLevel);
        if (options.verbosity >= 1)
        {
            std::cerr << "Expected median coverage :" << expectedValueGivenLevel << " for k-mer of length:" << options.fromLevel << std::endl; 
            std::cerr << "The estimated genome length is " << options.genomeLength << std::endl;
        }
#endif // #ifdef MEDIAN
	}

	/*restrictions for the searching levels*/
	cargo(myIndex).replen_min = options.fromLevel;
	cargo(myIndex).replen_max = options.toLevel;
//    cargo(myIndex).frequency = frequency;
//    cargo(myIndex).repeatCutoffs = options.repeatCutoffs;

	/*the core of the correction method*/
    FionaResources resources;
    traverseAndSearchCorrections<-1>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg,maxReadLength, resources);

    if (options.verbosity >= 2)
        std::cerr << "Time for searching between given levels: "<< SEQAN_PROTIMEDIFF(search) << " seconds." << std::endl
                  << "Time for compute overlap sums:           "<< options.timeComputeOverlapSum << " seconds." << std::endl;

#else // #ifndef FIONA_PARALLEL
	// FIONA PARALLEL SEARCH

	// construct q-gram index
//	TFionaQgramIndex qgramIndex(store.readSeqStore);
    TReadPrefixes prefixes;
    unsigned cutLength = options.fromLevel - 5; // we only need suffixes of length >= fromLevel (and their q-gram anchor)
    for (unsigned i = 0; i < length(store.readSeqStore); ++i)
        if (length(store.readSeqStore[i]) >= (unsigned)options.fromLevel)
            appendValue(prefixes, infix(store.readSeqStore[i], 0, length(store.readSeqStore[i]) - cutLength));
        else
            appendValue(prefixes, infix(store.readSeqStore[i], 0, 1));  // to short suffixes don't need to appear in the q-gram index
    TFionaQgramIndex qgramIndex(prefixes);
    cargo(qgramIndex).optionsPtr = &options;

	String<uint64_t> packages;
    SEQAN_PROTIMESTART(constructQgramExt);
    if (options.verbosity >= 1)
        std::cerr << "Construct external q-gram index ... " << std::flush;

#ifdef FIONA_INTERNAL_MEMORY

    SEQAN_PROTIMESTART(countQgramInt);
    if (options.verbosity >= 1)
        std::cerr << std::endl << "Counting phase of internal q-gram index ... " << std::flush;

    // 1. count q-grams
    resize(indexDir(qgramIndex), _fullDirLength(qgramIndex), Exact());
    _qgramClearDir(indexDir(qgramIndex), qgramIndex.bucketMap, Parallel());
    _qgramCountQGrams(indexDir(qgramIndex), qgramIndex.bucketMap, indexText(qgramIndex), indexShape(qgramIndex), getStepSize(qgramIndex), Parallel());
    _qgramDisableBuckets(qgramIndex);

    typedef typename Fibre<TFionaQgramIndex, QGramDir>::Type TQGramDir;
    typedef typename Value<TQGramDir>::Type TQGramDirValue;
    typedef typename Size<TQGramDir>::Type TQGramDirSize;
    TQGramDir origDir = indexDir(qgramIndex);

    // 2. create super packages for multiple q-gram index creations
    unsigned dirLen = length(origDir);
    uint64_t numSuffixes = 0;
    SEQAN_OMP_PRAGMA(parallel for reduction(+ : numSuffixes))
    for (int i = 0; i < (int)dirLen; ++i)
        if (origDir[i] != (TQGramDirValue)-1)
            numSuffixes += origDir[i];

	String<uint64_t> superPackages;
    uint64_t sumSuffixes = 0;
    uint64_t nextThresh = 0;
    for (unsigned i = 0; i < dirLen; ++i)
    {
        if (nextThresh <= sumSuffixes)
        {
            appendValue(superPackages, i);
            nextThresh = (length(superPackages) * numSuffixes) / (uint64_t)options.numSuperPackages;
        }
        if (origDir[i] != (TQGramDirValue)-1)
            sumSuffixes += origDir[i];
    }

    if (options.verbosity >= 1)
        std::cerr << "done. (" << SEQAN_PROTIMEDIFF(countQgramInt) << " seconds, " << sumSuffixes << " kmers)" << std::endl;

    // 3. create partial q-gram index and iterate over its packages
    for (unsigned superPackage = 0; superPackage < options.numSuperPackages; ++superPackage)
    {
        SEQAN_PROTIMESTART(fillQgramSAInt);
        if (options.verbosity >= 1)
        {
            std::cerr << "Create partial internal q-gram index (" << superPackage + 1 << " of " << options.numSuperPackages << ", ";
            std::cerr << "buckets " << superPackages[superPackage] << '-' << superPackages[superPackage+1] <<") ... " << std::flush;
        }

        // 4. keep only a portion of buckets with ids in [beginBucket..endBucket)
        TQGramDir &dir = indexDir(qgramIndex);
        dir = origDir;
        TQGramDirSize beginBucket = superPackages[superPackage];
        TQGramDirSize endBucket = superPackages[superPackage + 1];

        for (TQGramDirSize i = 0; i < beginBucket; ++i)
            dir[i] = (TQGramDirValue)-1;
        for (TQGramDirSize i = endBucket; i < dirLen - 1; ++i)
            dir[i] = (TQGramDirValue)-1;

        resize(indexSA(qgramIndex), _qgramCummulativeSum(indexDir(qgramIndex), True(), True(), Unsigned<1>(), Parallel()), Exact());
        _qgramFillSuffixArray(indexSA(qgramIndex), indexText(qgramIndex), indexShape(qgramIndex), indexDir(qgramIndex), qgramIndex.bucketMap, getStepSize(qgramIndex), True(), Parallel());
        _qgramPostprocessBuckets(indexDir(qgramIndex), Parallel());
        
        if (options.verbosity >= 1)
            std::cerr << "done. (" << SEQAN_PROTIMEDIFF(fillQgramSAInt) << " seconds, " << back(indexDir(qgramIndex)) << " kmers, ";

#else
    resize(indexSA(qgramIndex), _qgramQGramCount(qgramIndex), Exact());
    resize(indexDir(qgramIndex), _fullDirLength(qgramIndex), Exact());
    createQGramIndexExt(qgramIndex);

//	createQGramIndexExtSA(qgramIndex);  // alternative but (unfortunately) slower variant using a Mapper
	resize(indexSA(qgramIndex), back(indexDir(qgramIndex)), Exact());
#endif

#if defined(FIONA_REDUCE_MEMORY) && !defined(FIONA_INTERNAL_MEMORY)
    flush(indexSA(qgramIndex));
#endif

    if (options.verbosity >= 1)
        std::cerr << "done. (" << SEQAN_PROTIMEDIFF(constructQgramExt) << " seconds, " << back(indexDir(qgramIndex)) << " kmers, ";
#ifndef FIONA_INTERNAL_MEMORY
	clear(indexText(qgramIndex));
#endif

//    unsigned dirLen = length(indexDir(qgramIndex));
//	SEQAN_PROTIMESTART(purgeNBuckets);
//	typedef typename Fibre<TFionaIndex, FibreSA>::Type TSA;
//    typedef typename Iterator<TSA, Standard>::Type TSAIter;
//    TSAIter beginSA = begin(indexSA(qgramIndex), Standard());
//    TSAIter srcIt = beginSA;
//    TSAIter dstIt = beginSA;
//    unsigned i;
//    resize(maskedBuckets, dirLen - 1, false);
//
//    // 1. mask k-mers with Ns
//    for (i = 0; i < dirLen - 1; ++i)
//        maskedBuckets[i] = hashContainsN(i);
//
//    // 2. mask k-mers that are trivial repeats, e.g. X^k
//    maskTrivialRepeats(maskedBuckets, indexShape(qgramIndex));
//
//    // 3. mask k-mers from repeat regions
//    maskRepeatBuckets(indexDir(qgramIndex));
//
//    // 3. remove all masked k-mers from k-mer index
//	std::cerr << "Purge repetitive k-mers ........... " << std::flush;
//    for (i = 0; i < dirLen - 1; ++i)
//    {
//        uint64_t bucketLen = indexDir(qgramIndex)[i + 1] - indexDir(qgramIndex)[i];
//        indexDir(qgramIndex)[i] = dstIt - beginSA;
//        // copy bucket unless it is marked for removal
//        if (!maskedBuckets[i])
//        {            
//            if (dstIt != srcIt)
//                std::copy(srcIt, srcIt + bucketLen, dstIt);
//            dstIt += bucketLen;
//        }
//        srcIt += bucketLen;
//    }
//    indexDir(qgramIndex)[i] = dstIt - beginSA;
//    resize(indexSA(qgramIndex), dstIt - beginSA);
//	std::cerr << "done. (" << SEQAN_PROTIMEDIFF(purgeNBuckets) << " seconds)" << std::endl;

    // distribute q-gram buckets over work packages
#ifdef FIONA_USE_SA
    packages = indexDir(qgramIndex);
#else
    unsigned dirLen = length(indexDir(qgramIndex));
    uint64_t numPacks = options.packagesPerThread * omp_get_max_threads();
    uint64_t numSuffixes = back(indexDir(qgramIndex));
    uint64_t nextThresh = numSuffixes / numPacks;
	appendValue(packages, 0);
    for (unsigned i = 1; i < dirLen; ++i)
    {
        if (nextThresh <= indexDir(qgramIndex)[i])
        {
            appendValue(packages, indexDir(qgramIndex)[i]);
            nextThresh = (length(packages) * numSuffixes) / numPacks;
        }
    }
#endif
    if (options.verbosity >= 1)
        std::cerr << length(packages) << " packages)" << std::endl;

    // don't need the q-gram dir any more, from now we use packages (multiple q-gram buckets)
    clear(indexDir(qgramIndex));
    shrinkToFit(indexDir(qgramIndex));
    unsigned finished = 0;
    bool inTerm = isatty(fileno(stdout));

    if (options.verbosity >= 1)
        std::cerr << "Parallel suffix tree traversal .... ";
    if (inTerm && options.verbosity >= 2)
        std::cerr << "  0%";
    if (options.verbosity >= 2)
        std::cerr << std::flush;

    //investigatedNodes=0;
    //putCorrections=0;

    String<FionaResources> resourcesPerPackage;
    resize(resourcesPerPackage, length(packages) - 1, Exact());

#if defined(FIONA_REDUCE_MEMORY) && !defined(FIONA_INTERNAL_MEMORY)
    flush(indexSA(qgramIndex));
    FileMapping<File<> > mapping;
    open(mapping, indexSA(qgramIndex).file);
#endif

    double startTime = sysTime();
    std::vector<double> done(omp_get_max_threads(), 0);

    // this must be done before and here (out of the parallel section)
    _refreshStringSetLimits(store.readSeqStore);

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic,1))
	for (int i = 1; i < (int)length(packages); ++i)
	{
        typedef uint64_t                            TFileSize;

        SEQAN_OMP_PRAGMA(atomic)
        ++finished;

        TFileSize bktBegin = packages[i-1];
        TFileSize bktEnd = packages[i];
		if (bktBegin + 3 >= bktEnd)     // we need at least 3 suffixes to distinguish correct from incorrect bases
            continue;

		TFionaIndex myIndex(store.readSeqStore);

#if defined(FIONA_REDUCE_MEMORY) && !defined(FIONA_INTERNAL_MEMORY)
        typedef Fibre<TFionaIndex, FibreSA>::Type   TSA;
        typedef typename Value<TSA>::Type           TSAValue;
        typedef typename Size<TSA>::Type            TSASize;

        TFileSize mapOfs = bktBegin & ~(TFileSize)0xfff;
        TFileSize mapSize = (TFileSize)sizeof(TSAValue) * (bktEnd - mapOfs);

        TSAValue *mapPtr = (TSAValue*)mapFileSegment(mapping, (TFileSize)sizeof(TSAValue) * mapOfs, mapSize, MAP_COPYONWRITE | MAP_RDWR);
        indexSA(myIndex).begin = mapPtr + (bktBegin - mapOfs);
        indexSA(myIndex).end   = mapPtr + (bktEnd   - mapOfs);
#else
		indexSA(myIndex) = toRange(infix(indexSA(qgramIndex), bktBegin, bktEnd));
#endif

		cargo(myIndex).replen_min = options.fromLevel;
		cargo(myIndex).replen_max = options.toLevel;
//        cargo(myIndex).repeatCutoffs = options.repeatCutoffs;

		if (inTerm && options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical(progressOutput))
            {
                for (int u = 0; u <= omp_get_thread_num(); ++u)
                    std::cerr << '\n';
                std::cerr << "thread " << omp_get_thread_num() << "\t: ";
                std::cerr << prefix(suffix(store.readSeqStore, front(indexSA(myIndex))), QGRAM_LENGTH) << " - ";
                std::cerr << prefix(suffix(store.readSeqStore, back(indexSA(myIndex))), QGRAM_LENGTH);
                std::cerr << '\t' << (bktEnd - bktBegin) << "          ";
                std::cerr << (char)27 << '[' << (omp_get_thread_num() + 1) << 'F' << std::flush;
            }
        }

		TConstrainedIterator myConstrainedIterator(myIndex);

#ifdef FIONA_USE_SA
        // extend prefix sorting to the k-max
        _refineQGramIndexBucket(
            indexSA(myIndex),
            indexText(myIndex),
            QGRAM_LENGTH,
            options.toLevel + 1);   // sort by one more character, as the toLevel restriction is applied to the parentRepLength

        value(myConstrainedIterator).range.i1 = 0;
        _setSizeInval(value(myConstrainedIterator).range.i2);
        value(myConstrainedIterator).parentRight = value(myConstrainedIterator).range.i2;
        value(myConstrainedIterator).repLen = QGRAM_LENGTH;
        value(myConstrainedIterator).lastChar = suffix(store.readSeqStore, front(indexSA(myIndex)))[QGRAM_LENGTH - 1];
#endif
        FionaResources &resources = resourcesPerPackage[i - 1];
        resources.bucketBegin = bktBegin;
        resources.bucketEnd = bktEnd;
        traverseAndSearchCorrections<-1>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg, maxReadLength, resources);
/*    if (options.loopLevel == 0)
                traverseAndSearchCorrections<0>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg, maxReadLength, resources);
    if (options.loopLevel == 1)
                traverseAndSearchCorrections<1>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg, maxReadLength, resources);
    if (options.loopLevel == 2)
                traverseAndSearchCorrections<2>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg, maxReadLength, resources);
    if (options.loopLevel == 3)
                traverseAndSearchCorrections<3>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg, maxReadLength, resources);
    if (options.loopLevel == 4)
                traverseAndSearchCorrections<4>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg, maxReadLength, resources);
    if (options.loopLevel == 5)
                traverseAndSearchCorrections<5>(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg, maxReadLength, resources);
*/
//		traverseAndSearchCorrections(myConstrainedIterator, store, correctionList, firstCorrectionForRead, options, alg,readLength, resources);
//		mmapAdvise(indexSA(qgramIndex), MAP_DONTNEED, bktBegin, bktEnd);
#if defined(FIONA_REDUCE_MEMORY) && !defined(FIONA_INTERNAL_MEMORY)
        unmapFileSegment(mapping, mapPtr, mapSize);
#endif
        if (inTerm && options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical(progressOutput))
            {
				std::cerr << "Parallel suffix tree traversal .... " << std::setw(3) << (100 * finished) / (length(packages) - 1) << '%' << std::flush;

                for (int u = 0; u <= omp_get_thread_num(); ++u)
                    std::cerr << '\n';
                std::cerr << "thread " << omp_get_thread_num() << "\t: done                                 ";
                std::cerr << (char)27 << '[' << (omp_get_thread_num() + 1) << 'F' << std::flush;
            }
		}
		
		done[omp_get_thread_num()] = sysTime() - startTime;
	}
	if (options.verbosity >= 1)
	{
	    for (unsigned i = 0; i < length(done); ++i)
    	    std::cerr << "Thread " << i << " took " << done[i] << " wallclock seconds." << std::endl;
	}
	
	if (inTerm && options.verbosity >= 2)
    {
		std::cerr << "\b\b\b\b";
        for (int u = 0; u <= omp_get_max_threads() + 1; ++u)
            std::cerr << std::endl;
    }
	if (inTerm && options.verbosity >= 1)
        std::cerr << "done. (" << SEQAN_PROTIMEDIFF(search) << " seconds)" << std::endl;

    sort(resourcesPerPackage, Parallel());
	if (inTerm && options.verbosity >= 2)
    {
        std::cerr << std::endl;
        std::cerr << "kmerFirst\tkmerLast\ttime\t\tnodes\tcorrections\tsuffixes" << std::endl;
        for (unsigned i = 0; i < 30; ++i)
        {
            FionaResources &resources = resourcesPerPackage[i];
            std::cerr << prefix(suffix(store.readSeqStore, indexSA(qgramIndex)[resources.bucketBegin]), QGRAM_LENGTH) << '\t';
            std::cerr << prefix(suffix(store.readSeqStore, indexSA(qgramIndex)[resources.bucketEnd - 1]), QGRAM_LENGTH) << '\t';
            std::cerr << resources.cpuTime << '\t' << resources.investigatedNodes << '\t' << resources.putCorrections << '\t' << (resources.bucketEnd - resources.bucketBegin) << std::endl;
//          std::cerr << investigatedNodes << " many error nodes have been investigated" <<std::endl;
//          std::cerr << putCorrections << " many Corrections where found and potentially saved "<< std::endl;
        }
        std::cerr << std::endl;
    }
#ifdef FIONA_INTERNAL_MEMORY
    }   // for-loop over all superpackages
#endif

#endif // #ifndef FIONA_PARALLEL


	unsigned totalCorrections = 0;
	//get the number of corrections from the next function as more than one correction
	//per read might occur
	totalCorrections = applyReadErrorCorrections(correctionList,firstCorrectionForRead,store,options);
	unsigned readCorrections=0;
    for (unsigned a = 0; a < length(firstCorrectionForRead); ++a)
    {
        if (firstCorrectionForRead[a] == std::numeric_limits<unsigned>::max()) continue;
        ++readCorrections;
    }

    if (options.verbosity >= 1)
        std::cerr << "Total corrected reads number is "<< readCorrections << std::endl
                  << "Number of total corrections is "<< totalCorrections << std::endl;

    // remove reverse complements
    resize(store.readSeqStore, readCount);
    return totalCorrections;
}

// Write output and return 0 on success, a different value on errors.  Update numCorrected to reflect the number of
// corrected reads.

template <typename TFragmentStore>
int writeOutput(unsigned & numCorrected, TFragmentStore const & store, FionaOptions const & options)
{
    // Write out the corrected reads and stream through input file for getting the read ids.
    SeqFileIn inFile;
    std::cerr << "Opening input " << options.inputFilename << "\n";
    if (!open(inFile, toCString(options.inputFilename)))
    {
        std::cerr << "ERROR: Could not open " << options.inputFilename << " for reading.\n";
        return 1;
    }
    bool success;
    SeqFileOut outFile;
    if (options.outputFilename != "-")
        success = open(outFile, toCString(options.outputFilename));
    else
        success = open(outFile, std::cout, Fasta());
    std::cerr << "Opening output " << options.outputFilename << "\n";

    if (!success)
    {
        std::cerr << "ERROR: Could not open " << options.outputFilename << " for writing.\n";
        return 1;
    }

    // Buffer variables for FASTA id/sequence.
    CharString id;
    Dna5String seq2, seq;

    numCorrected = 0;
    // disable linebreak in output file
    context(outFile).options.lineLength = 0;

    for (unsigned i = 0; i < length(store.readSeqStore); ++i)
    {
        // Read record from input (for the id only).
        readRecord(id, seq, inFile);

        // Overwrite the sequence from the input file with the corrected read sequence.
        seq = store.readSeqStore[i];
        // Append correction information to the read id if we collected any.
        if (options.appendCorrectionInfo)
        {
            append(id, store.readNameStore[i]);
            // TODO(holtgrew): Value of numCorrected is not increased if options.appendCorrectionInfo is false.
            numCorrected += !empty(store.readNameStore[i]);
        }

        // Trim leading and trailing Ns that could not be substituted.
        if (options.trimNsOnOutput)
        {
            int beginPos = 0, endPos = length(seq);
            for (; beginPos < (int)length(seq); ++beginPos)
                if (seq[beginPos] != 'N')
                    break;
            for (; endPos > 0; --endPos)
                if (seq[endPos - 1] != 'N')
                    break;
            if (beginPos > endPos)
                endPos = beginPos;

            if (options.appendCorrectionInfo && (beginPos != 0 || endPos != (int)length(seq)))
            {
                std::stringstream ss;
                ss << " trimmed to [" << beginPos << ", " << endPos << ")";
                append(id, ss.str().c_str());
            }

            seq2 = infix(seq, beginPos, endPos);
            seq = seq2;
        }

        // Write out the FASTA record to the output file.
        writeRecord(outFile, id, seq);
    }

    std::cerr << "Wrote " << length(store.readSeqStore) << " sequences\n";

    return 0;
}

// Parse the command line and return the status of the parsing.
seqan::ArgumentParser::ParseResult
parseCommandLine(FionaOptions & options, int argc, char const ** argv)
{
    // Setup command line parser.
    seqan::ArgumentParser parser(FIONA_BINARY_NAME);

    // Set short description, version, and date.
    setShortDescription(parser, "Parallel and automatic read error correction");
    setCategory(parser, "Error Correction");

    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-g\\fP \\fIGENOME_LEN\\fP \\fIIN.{fq,fa}\\fP \\fIOUT.fa\\fP");
    addDescription(parser,
                   "Fiona is a tool for the correction of NGS read data sets.  It uses a novel "
                   "statistical approach for high quality and state-of-the art data structures "
                   "for low resource consumptions and features a good parallelization.");
    addDescription(parser,
                   "You have to specify the approximate genome length of the donor in \\fIGENOME_LEN\\fP. "
                   "The reads are read from the file \\fIIN.{fq,fa}\\fP and are written to \\fIOUT.fa\\fP.");

    // Fiona gets two parameters:  The paths to the input and the output files.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, 0, "fa fasta fq fastq");
    setHelpText(parser, 0, "An input file with reads to be corrected.");
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setValidValues(parser, 1, "fa fasta fq fastq");
    setHelpText(parser, 1, "An output file to store the corrected reads.");

    // General Options.
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "More verbose output."));
    addOption(parser, seqan::ArgParseOption("", "correction-infos",
                                            "Enable embedding of correction information in the output file."));

    // Internal Options
    addSection(parser, "Internal");
    addOption(parser, seqan::ArgParseOption("", "global-corr-limit", "Limit corrections globally and not per round."));
    addOption(parser, seqan::ArgParseOption("", "no-final-trim-ns", "Disable trimming of Ns at the end."));

    // Dataset Properties.
    addSection(parser, "Dataset Properties");

    addOption(parser, seqan::ArgParseOption("g", "genome-length", "Approximate length of the underlying genome.",
                                            seqan::ArgParseOption::INT64, "LEN"));
    setMinValue(parser, "genome-length", "1");
    setRequired(parser, "genome-length");

    addOption(parser, seqan::ArgParseOption("e", "error-rate",
                                            "Approximate per-base error rate in the read set. A slight "
                                            "overestimation gives better results.",
                                            seqan::ArgParseOption::DOUBLE, "ERATE"));
    setMinValue(parser, "error-rate", "0");
    setMaxValue(parser, "error-rate", "1");
    setDefaultValue(parser, "error-rate", options.errorrate);

    addOption(parser, seqan::ArgParseOption("oe", "overlap-error-scale",
                                            "The \\fIerror-rate\\fP is multiplied by this scale to define the error rate cutoff in the pairwise read overlap.",
                                            seqan::ArgParseOption::DOUBLE, "ERATE"));
    setMinValue(parser, "overlap-error-scale", "0");
    setDefaultValue(parser, "overlap-error-scale", "2");

    // Tree Iteration Options.
    addSection(parser, "Tree Iteration Options");

    // TODO(holtgrew): The old help for --levels was "set to 0 <x> for auto HiTEC, 1 <x> for auto fiona".. Fix this here?
    addOption(parser, seqan::ArgParseOption("fl", "from-level",
                                            "Set the lower bound on the level for suffix tree DFS.  Use "
                                            "\\fI0\\fP for both \\fIfrom-level\\fP and \\fIto-level\\fP "
                                            "to get automatic level detection.",
                                            seqan::ArgParseOption::INTEGER, "LEVEL"));
    setMinValue(parser, "from-level", "0");
    setDefaultValue(parser, "from-level", options.fromLevel);

    // TODO(holtgrew): The old help for --levels was "set to 0 <x> for auto HiTEC, 1 <x> for auto fiona".. Fix this here?
    addOption(parser, seqan::ArgParseOption("tl", "to-level",
                                            "Set the upper bound on the level for suffix tree DFS.  Use "
                                            "\\fI0\\fP for both \\fIto-level\\fP and \\fIto-level\\fP "
                                            "to get automatic level detection.",
                                            seqan::ArgParseOption::INTEGER, "LEVEL"));
    setMinValue(parser, "to-level", "0");
    setDefaultValue(parser, "to-level", options.toLevel);

    addOption(parser, seqan::ArgParseOption("dsr", "depth-sample-rate", "The depth sampling rate factor.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "depth-sample-rate", "1");
    setDefaultValue(parser, "depth-sample-rate", options.depthSampleRate);

    // Repeat Masking Options
    addSection(parser, "Repeat Masking Options");

    addOption(parser, seqan::ArgParseOption("", "no-mask-repeats", "Turn off automatic repeat masking."));
    //hideOption(parser, "no-mask-repeats");

    addOption(parser, seqan::ArgParseOption("krr", "kmer-repeat-ratio",
                                            "The fraction of k-mers that are considered as repeats.",
                                            seqan::ArgParseOption::DOUBLE, "RATIO"));
    setMinValue(parser, "kmer-repeat-ratio", "0");
    setMaxValue(parser, "kmer-repeat-ratio", "1");
    setDefaultValue(parser, "kmer-repeat-ratio", options.kmerAbundanceCutoff);

    addOption(parser, seqan::ArgParseOption("krsd", "kmer-repeat-std-dev",
                                            "Multiples of standard deviation (for k-mer repeat cut-off).",
                                            seqan::ArgParseOption::DOUBLE, "SCALE"));
    setMinValue(parser, "kmer-repeat-std-dev", "0");
    setDefaultValue(parser, "kmer-repeat-std-dev", options.kmerStdDevCutOff);

    // Correction Algorithm Options.
    addSection(parser, "Correction Algorithm Options");

    addOption(parser, seqan::ArgParseOption("", "method", "Selects the correction method to use.",
                                            seqan::ArgParseOption::STRING, "NAME"));
    setValidValues(parser, "method", "classifier control_fp control_fn expected count");
    setDefaultValue(parser, "method", "classifier");

    addOption(parser, seqan::ArgParseOption("i", "iterations",
                                            "Number of iterations.  Use \\fI0\\fP for auto-detection.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "iterations", "0");
    setDefaultValue(parser, "iterations", options.cycles);

    addOption(parser, seqan::ArgParseOption("f", "expected",
                                            "Use expected value correction with the given strictness cutoff "
                                            "for the \\fIexpected\\fP method.",
                                            seqan::ArgParseOption::DOUBLE, "CUTOFF"));
    setMinValue(parser, "expected", "0");
    setDefaultValue(parser, "expected", "1");

    addOption(parser, seqan::ArgParseOption("c", "count", "Use fixed count correction cutoff.",
                                            seqan::ArgParseOption::DOUBLE, "CUTOFF"));
    setMinValue(parser, "count", "0");
    setDefaultValue(parser, "count", "7");

    addOption(parser, seqan::ArgParseOption("or", "odds-ratio", "Odds-ratio for the \\fIclassifier\\fP method.",
                                            seqan::ArgParseOption::DOUBLE, "RATIO"));
    setMinValue(parser, "odds-ratio", "0");
    setDefaultValue(parser, "odds-ratio", "1");

    addOption(parser, seqan::ArgParseOption("p", "p-value", "The p value for the \\fIexpected\\fP mode. In "
                                            "sensitivity mode, this is the false discovery rate.",
                                            seqan::ArgParseOption::DOUBLE, "P-VALUE"));
    setMinValue(parser, "p-value", "0");
    setDefaultValue(parser, "p-value", "1");

    addOption(parser, seqan::ArgParseOption("m", "mismatches", "The number of accepted mismatches per read.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "mismatches", "0");
    setDefaultValue(parser, "mismatches", options.acceptedMismatches);

    addOption(parser, seqan::ArgParseOption("os", "overlap-sum",
                                            "Filter on the number of overlapping bp needed to correct an "
                                            "erroneous bp.  A smaller value leads to lower sensitivity, a "
                                            "higher value leads to higher sensitivity.",
                                            seqan::ArgParseOption::DOUBLE, "P-VALUE"));
    setMinValue(parser, "overlap-sum", "0");
    setMaxValue(parser, "overlap-sum", "1");
    setDefaultValue(parser, "overlap-sum", options.wovsum);

#ifdef FIONA_ALLOWINDELS
    addOption(parser, seqan::ArgParseOption("id", "indel-length", "Maximal indel length.  Use \\fI0\\fP for "
                                            "correcting only substitutions and \\fI1\\fP for edit distance "
                                            "corrections on Illumina reads.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "indel-length", "0");
    {
        std::stringstream tmp;
        tmp << MAX_INDEL_LENGTH;
        setMaxValue(parser, "indel-length", tmp.str().c_str());
    }
    setDefaultValue(parser, "indel-length", options.maxIndelLength);
#endif

    // DEBUG Options

    addOption(parser, seqan::ArgParseOption("", "loop-level", "For time measurements.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "loop-level", options.loopLevel);
    hideOption(parser, "loop-level");

    addOption(parser, seqan::ArgParseOption("", "debug-read", "Dump information for a read given by its id.",
                                            seqan::ArgParseOption::INTEGER, "ID"));
    hideOption(parser, "debug-read");

    addOption(parser, seqan::ArgParseOption("", "corr-read", "Dump information for a correcting read.",
                                            seqan::ArgParseOption::INTEGER, "ID"));
    hideOption(parser, "corr-read");

    // Parallelization Options.
    addSection(parser, "Parallelization Options");

    addOption(parser, seqan::ArgParseOption("nt", "num-threads", "Number of threads to use (default 1).",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "num-threads", "1");
    setDefaultValue(parser, "num-threads", options.numThreads);

#ifdef FIONA_INTERNAL_MEMORY
    addOption(parser, seqan::ArgParseOption("", "super-packages", "Number of internal q-gram index creation runs.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "super-packages", "1");
    setDefaultValue(parser, "super-packages", options.numSuperPackages);
#endif


    addOption(parser, seqan::ArgParseOption("ppt", "packages-per-thread",
                                            "Set the number of work packages per thread.  More packages result "
                                            "lower memory consumption but possibly a longer running time.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "packages-per-thread", "1");
    setDefaultValue(parser, "packages-per-thread", options.packagesPerThread);

    // Documentation on the correction methods.
    addTextSection(parser, "Method Description");
    addText(parser, "The \\fB--method\\fP parameter can be used to select one of the following methods to decide on the coverage cutoff");
    addListItem(parser, "\\fIclassifier\\fP (default)", "Default method, all k-mers are modeled as a mixture of Poisson distributions.  Detect k-mers with errors using a naive bayes classifier (parameters: \\fB--error-rate\\fP, \\fB--odds-ratio\\fP)");
    addListItem(parser, "\\fIcontrol_fp\\fP (type I mode)", "The cutoff is set to control the probability of a false positive detection for each k, using the Poisson distribution for the reads without errors. (parameter: FP probability \\fB--p-value\\fP)");
    addListItem(parser, "\\fIcontrol_fn\\fP (type II mode)", "The cutoff is set to control the probability of a false negative detection (parameters: FN proba or 1-FN proba through \\fB--error-rate\\fP)");
    addListItem(parser, "\\fIexpected\\fP", "The cutoff is set to lambda-alpha*lambda^2, where lambda is the expected coverage of reads before sequencing.");
    addListItem(parser, "\\fIcount\\fP", "The cutoff is set manually with \\fB--count\\fP and the same for all values of k.");

    // Usage Examples.
    addTextSection(parser, "Examples");
    addText(parser,
            "Most users will only have to specify the genome length using the mandatory \\fB-g\\fP parameter and "
            "enable multi-threading using the \\fB-nt\\fP option.  For best performance, use as many threads "
            "as you have (virtual) cores in your machine.");

    std::string toolName = "\\fB" FIONA_BINARY_NAME "\\fP";
    addListItem(parser, toolName + " \\fB-g\\fP 4639675 IN.fq OUT.fq",
                "Correct reads in \\fIIN.fq\\fP with one thread and write the results to \\fIOUT.fq\\fP. "
                "The estimated genome length fits for E.coli.");

    addListItem(parser, toolName + " \\fB-nt\\fP 16 \\fB-g\\fP 4639675 IN.fq OUT.fq",
                "Same as above, but use \\fI16\\fP threads.");

    addListItem(parser, toolName + " \\fB-id\\fP 0 \\fB-g\\fP 4639675 IN.fq OUT.fq",
                "Sequential correction that corrects only mismatches no indels.");

    addListItem(parser, toolName + " \\fB-e\\fP 0.02 \\fB-g\\fP 4639675 IN.fq OUT.fq",
                "Sequential correction using an expected base error rate of 2%.");

    // Environment Variables.
    addTextSection(parser, "Environment Variables.");
    addText(parser,
            "Fiona uses the \\fBTMPDIR\\fP environment variable for creating temporary files.  If not set then "
            "\\fI/tmp\\fP is used which is fine for desktop settings.  In a compute server/data center setup, "
            "ask your administrator for the appropriate value.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine().
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract Argument and Option Value

    getArgumentValue(options.inputFilename, parser, 0);
    getArgumentValue(options.outputFilename, parser, 1);

    if (isSet(parser, "verbose"))
        options.verbosity = 1;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 2;

    options.trimNsOnOutput = !isSet(parser, "no-final-trim-ns");
    options.limitCorrPerRound = !isSet(parser, "global-corr-limit");

    seqan::CharString tmp;
    getOptionValue(tmp, parser, "method");
    options.method = methodForName(tmp);

    // The option member strictness is double-typed and stores the cutoff-like values.  We branch here between the
    // different methods and get the appropriate option into strictness.
    switch (options.method)
    {
        case EXPECTED:
            getOptionValue(options.strictness, parser, "expected");
            break;
        case COUNT:
            getOptionValue(options.strictness, parser, "count");
            break;
        case CONTROL_FP:
            getOptionValue(options.strictness, parser, "p-value");
            break;
        case CONTROL_FN:
            // Empty on purpose, --error-rate is always interpreted.
            break;
        case CLASSIFIER:
            getOptionValue(options.strictness, parser, "odds-ratio");
            break;
    }

    double overlapErrorScale = 2;
    getOptionValue(options.errorrate, parser, "error-rate");
    getOptionValue(overlapErrorScale, parser, "overlap-error-scale");
    options.overlap_errorrate = overlapErrorScale * options.errorrate;

    options.appendCorrectionInfo = isSet(parser, "correction-infos");
    getOptionValue(options.fromLevel, parser, "from-level");
    getOptionValue(options.toLevel, parser, "to-level");
    getOptionValue(options.packagesPerThread, parser, "packages-per-thread");
	getOptionValue(options.wovsum, parser, "overlap-sum");
    getOptionValue(options.genomeLength, parser, "genome-length");
	getOptionValue(options.acceptedMismatches, parser, "mismatches");
	getOptionValue(options.cycles, parser, "iterations");
    getOptionValue(options.loopLevel, parser, "loop-level");
	getOptionValue(options.kmerAbundanceCutoff, parser, "kmer-repeat-ratio");
	getOptionValue(options.kmerStdDevCutOff, parser, "kmer-repeat-std-dev");
	getOptionValue(options.depthSampleRate, parser, "depth-sample-rate");

#ifdef FIONA_ALLOWINDELS
	getOptionValue(options.maxIndelLength, parser, "indel-length");
#endif
    getOptionValue(options.numThreads, parser, "num-threads");
#ifdef FIONA_INTERNAL_MEMORY
    getOptionValue(options.numSuperPackages, parser, "super-packages");
#endif
    getOptionValue(options.debugRead, parser, "debug-read");
    getOptionValue(options.corrRead, parser, "corr-read");

    // Check Arguments.

	if (options.packagesPerThread >= _intPow((unsigned)ValueSize<Dna5>::VALUE, QGRAM_LENGTH))
    {
		std::cerr << "warning: packages-per-thread parameter is decreased to " << _intPow((unsigned)ValueSize<Dna5>::VALUE, QGRAM_LENGTH) - 1 << std::endl;
        // nothing more needs to be done, the parameter is implicitly decreased
    }

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, const char* argv[]) 
{
    // Declare options variable and parse command line.
    FionaOptions options;
    seqan::ArgumentParser::ParseResult parseRes = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0 (e.g. help
    // was printed).
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;

    // Set number of threads to use from the command line.
#if defined(_OPENMP)
    omp_set_num_threads(options.numThreads);
#endif

    if (options.verbosity >= 1)
        std::cerr << "FIONA - Read Correction\n"
                  << "=======================\n\n";

    bool autoCycles = (options.cycles == 0 || options.cycles == 1000);
    bool bestExpFit = (options.cycles == 0);
    if (autoCycles) options.cycles = MAX_NUM_ROUND;
    if (options.verbosity >= 1)
        printOptions(std::cerr, options);

    SEQAN_PROTIMESTART(correction);

    // Load original set of reads without read names.  When collecting correction information for debugging, we will
    // collect the correction string in the readNameStore and allocate space for this below.  When writing out, we will
    // stream through the input file again and get the read ids from there.
    TFionaFragStore store;
    if (options.verbosity >= 1)
        std::cerr << "Loading reads from " << options.inputFilename << "\n";
    if (!loadReadsNoNames(store, options.inputFilename, options))
    {
        std::cerr << "Failed to open reads file " << options.inputFilename << "\n"
                    << "Exiting ...\n";
        return 1;
    }
    else
    {
        if (options.verbosity >= 1)
            std::cerr << "Loaded " << length(store.readSeqStore) << " reads (" << SEQAN_PROTIMEDIFF(correction) << ".\n";
    }
    // If we collect correction information in the sequence ids/headers then we need to allocate space for them in the
    // readNameStore.
    if (options.appendCorrectionInfo)
        resize(store.readNameStore, length(store.readSeqStore), Exact());
    if (options.verbosity >= 1)
        std::cerr << "Done loading " << length(store.readSeqStore) << " sequences with a total of " << lengthSum(store.readSeqStore) << " nucleotides.\n";

/*	//DEBUG Corrrection Indel Pos
	//check the new CorrectionIndelPos struct and the linked list
	String<CorrectionIndelPos> correctionList;
	String<unsigned int> firstCorrectionForRead;
 	// append reverse complements of reads
	unsigned readCount = length(store.readSeqStore);
	Dna5String tmp;
	for (unsigned i = 0; i < readCount; ++i)
		{
			tmp = store.readSeqStore[i];
			reverseComplement(tmp);
			appendValue(store.readSeqStore, tmp);
		}

	_testCorrectionStruct(correctionList,firstCorrectionForRead,store);
/// for debugging output vorziehen:
 // write in file all input reads with the corrected one
        std::ofstream out(toCString(getArgumentValue(parser, 1)));
        int numCorrected = 0;
        for (unsigned i = 0; i < length(store.readNameStore); ++i)
        {
                // to give the number of reads corrected for several iteration
                if (strContains(toCString(store.readNameStore[i]), "corrected"))
                        ++numCorrected;

                out << '>' << store.readNameStore[i] << std::endl;
                out << store.readSeqStore[i] << std::endl;
        }
///DEBUG END	
	exit(0);	
*/

	// initialise the top and down level by using the log4 from the total number of reads
	// Done in the CorrectRead function at every round now (in case the error rate would change automatically)
//	if (options.fromLevel == 0)
//	{
//		int logRation = static_cast<int>(log10(static_cast<double>(length(store.readSeqStore))) / log10(4.0));
//		options.fromLevel = logRation + 2;
//		options.toLevel   = options.fromLevel + 10;
//		std::cerr << "The estimated top level is " << options.fromLevel << " and the down level is " << options.toLevel << std::endl;
//	}
	
#ifndef FIONA_NOERROROPTIMIZATION
    if (options.verbosity >= 1)
        std::cerr << "Use normal full optimization"<<std::endl;
#endif
    if (options.verbosity >= 1)
        std::cerr << "Building external index with Qgram length " << QGRAM_LENGTH << std::endl;
    options.autolevel = (options.fromLevel <= 1);

    //if (autoCycles) options.cycles = 20;
    String<double> logCorrections;
    String<double> roundsDone;
    double lastAdjRSquare = 0;
    unsigned nfamprev = 0;
    if (options.verbosity >= 1)
    {
        std::cerr << "number iters: " << options.cycles << std::endl;
        std::cerr << "\n"
                     "__RUNNING READ CORRECTION________________________________________\n";
    }

	for (options.cycle = 1; options.cycle <= options.cycles; ++options.cycle)
	{
        if (options.verbosity >= 1)
        {
            std::cerr << std::endl << "Cycle "<< options.cycle;
            if (!autoCycles) std::cerr << " of " << options.cycles;
            std::cerr << std::endl;
        }
		unsigned numCorrected = 0;
		nfamprev = nfamilies;
		nfamilies = 0;
        switch (options.method)
        {
            case CONTROL_FP:
                // use of p-value like a limit
                numCorrected = correctReads(store, options, FionaPoisson());
                break;
            case COUNT:
                numCorrected = correctReads(store, options, FionaCount());
                break;
            case CONTROL_FN:
                numCorrected = correctReads(store, options, FionaPoissonSens());
                break;
            case CLASSIFIER:
                numCorrected = correctReads(store, options, FionaPoissonClassif());
                break;
            case EXPECTED:
                // use an expected value for a certain level
                numCorrected = correctReads(store, options, FionaExpected());
                break;
        }
		//Todo (Hugues) adjust error rate estimate between rounds
		//
        if (options.verbosity >= 1)
            std::cerr << std::endl << "Number of families at nodes:" << nfamilies << std::endl;
		if (options.verbosity >= 1 && options.cycle > 1)
			std::cerr << std::endl << "Relative change: " << ((float) (nfamprev - nfamilies)) / (nfamprev == 0u ? 1 : nfamprev) << std::endl;

//		if (options.acceptedMismatches > 0) --options.acceptedMismatches;

		// TODO maybe to stop if there is not reads corrected in the cycle before
		// if so after each iteration must save the ID for the reads which are corrected
		// thus we can also show the total number of reads that are corrected at the final stage
		//if (autoCycles)
		//{
            resize(logCorrections, options.cycle);
            resize(roundsDone, options.cycle);
            logCorrections[options.cycle-1] = (double)log((double)numCorrected);
            roundsDone[options.cycle-1]     = (double)options.cycle;

            if (options.cycle >= 1)
            {
                //compute adjusted R-Square after fitting model
                LinearModel linearModel;
                linearRegression(linearModel,roundsDone,logCorrections);
                double adjRSquare =  adjustedRSquare(linearModel,roundsDone,logCorrections);
                if (options.verbosity >= 2)
                    std::cerr << "The adjusted R^2 in cycle " << options.cycle << " is " << adjRSquare << "\n";
                if (autoCycles)
                {
                    //do another round if adjusted R square value is better than 0.95
                    if (!bestExpFit && options.cycle > 3)
                    {
                        if (adjRSquare <= 0.95)
                        {
                            if (options.verbosity >= 2)
                                std::cerr <<std::endl<<"Stopped at cycle: "<< options.cycle <<" with adjustedRSquare : "<< adjRSquare <<std::endl;
                            ++options.cycle;
                            break;
                        }
                    }
                    else
                    {

                        if (adjRSquare < lastAdjRSquare)
                        {
                            if (options.verbosity >= 2)
                                std::cerr <<std::endl<<"Stopped at cycle: "<< options.cycle <<" with adjustedRSquare (previous): "<< adjRSquare <<" ("<<lastAdjRSquare <<")"<<std::endl;
                            ++options.cycle;
                            break;
                        }
                        else
                        {
                            lastAdjRSquare=adjRSquare;
                        }
                    }
                }
            }
        //}
	}

    // Write out corrected reads.
    unsigned numCorrected = 0;
    int res = writeOutput(numCorrected, store, options);
    if (res != 0)
        return res;

    if (options.verbosity >= 1 && options.cycles > 1)
        std::cerr << "Total number reads corrected for " << options.cycle-1 << " cycles is " << numCorrected << std::endl;

//	struct rusage usage;
//	getrusage(RUSAGE_SELF, &usage);
    if (options.verbosity >= 1)
    {
        std::cerr << std::endl;
        std::cerr << "Time required for execution: " << SEQAN_PROTIMEDIFF(correction) << " seconds." << std::endl;
//	    std::cerr << "Peak resident memory usage:  " << usage.ru_maxrss / (1024*1024) << " Mb." << std::endl;
    }

	return 0;
}
