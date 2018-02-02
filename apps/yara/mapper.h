// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef APP_YARA_MAPPER_H_
#define APP_YARA_MAPPER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Mapper Options
// ----------------------------------------------------------------------------

struct Options
{
    typedef std::string                     TString;
    typedef std::vector<TString>            TList;
    typedef FileFormat<BamFileOut>::Type    TOutputFormat;

    uint64_t            contigsSize;
    uint64_t            contigsMaxLength;
    uint64_t            contigsSum;

    CharString          contigsIndexFile;
    Pair<CharString>    readsFile;
    CharString          outputFile;
    TOutputFormat       outputFormat;
    SecondaryAlignments secondaryAlignments;
    TList               secondaryAlignmentsList;
    bool                uncompressedBam;
    CharString          readGroup;

    MappingMode         mappingMode;
    float               errorRate;
    float               indelRate;
    float               strataRate;
    Sensitivity         sensitivity;
    TList               sensitivityList;

    bool                singleEnd;
    unsigned            libraryLength;
    unsigned            libraryDev;
//    LibraryOrientation  libraryOrientation;
//    TList               libraryOrientationList;
    bool                verifyMatches;

    unsigned            readsCount;
    unsigned            threadsCount;
    unsigned            hitsThreshold;
    bool                rabema;
    unsigned            verbose;

    CharString          commandLine;
    CharString          version;

    Options() :
        contigsSize(),
        contigsMaxLength(),
        contigsSum(),
        secondaryAlignments(TAG),
        uncompressedBam(false),
        readGroup("none"),
        mappingMode(STRATA),
        errorRate(0.05f),
        indelRate(0.25f),
        strataRate(0.00f),
        sensitivity(HIGH),
        singleEnd(true),
        libraryLength(),
        libraryDev(),
//        libraryOrientation(FWD_REV),
        verifyMatches(true),
        readsCount(100000),
        threadsCount(1),
        hitsThreshold(300),
        rabema(false),
        verbose(0)
    {
#ifdef _OPENMP
        threadsCount = std::thread::hardware_concurrency();
#endif
        appendValue(secondaryAlignmentsList, "tag");
        appendValue(secondaryAlignmentsList, "record");
        appendValue(secondaryAlignmentsList, "omit");

        appendValue(sensitivityList, "low");
        appendValue(sensitivityList, "high");
        appendValue(sensitivityList, "full");

//        appendValue(libraryOrientationList, "fwd-rev");
//        appendValue(libraryOrientationList, "fwd-fwd");
//        appendValue(libraryOrientationList, "rev-rev");
    }
};

// ----------------------------------------------------------------------------
// Mapper Configuration
// ----------------------------------------------------------------------------

template <typename TThreading_       = Parallel,
          typename TSequencing_      = SingleEnd,
          typename TStrategy_        = Strata,
          typename TContigsSize_     = uint8_t,
          typename TContigsLen_      = uint32_t,
          typename TContigsSum_      = uint32_t,
          typename TAlloc_           = MMap<>,
          unsigned BUCKETS_          = 3>
struct ReadMapperConfig
{
    typedef TThreading_         TThreading;
    typedef TSequencing_        TSequencing;
    typedef TStrategy_          TStrategy;
    typedef TContigsSize_       TContigsSize;
    typedef TContigsLen_        TContigsLen;
    typedef TContigsSum_        TContigsSum;
    typedef TAlloc_             TAlloc;

    static const unsigned BUCKETS = BUCKETS_;
};

// ----------------------------------------------------------------------------
// Mapper Traits
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct MapperTraits
{
    typedef typename TConfig::TThreading                            TThreading;
    typedef typename TConfig::TSequencing                           TSequencing;
    typedef typename TConfig::TStrategy                             TStrategy;
    typedef typename TConfig::TContigsSize                          TContigsSize;
    typedef typename TConfig::TContigsLen                           TContigsLen;
    typedef typename TConfig::TContigsSum                           TContigsSum;
    typedef typename TConfig::TAlloc                                TAlloc;

    typedef SeqStore<void, YaraContigsConfig<TAlloc> >              TContigs;
    typedef typename TContigs::TSeqs                                TContigSeqs;
    typedef typename TContigs::TSeqNames                            TContigNames;
    typedef typename Value<TContigSeqs>::Type                       TContig;
    typedef typename StringSetPosition<TContigSeqs>::Type           TContigsPos;

    typedef YaraFMConfig<TContigsSize, TContigsLen, TContigsSum, TAlloc> TIndexConfig;
    typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
    typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;
    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef typename Fibre<TIndex, FibreSA>::Type                   TSA;

    typedef SeqStore<void, YaraReadsConfig>                         TReads;
    typedef typename If<IsSameType<TSequencing, PairedEnd>,
                        Pair<SeqFileIn>, SeqFileIn>::Type           TReadsFileIn;
    typedef PrefetchedFile<TReadsFileIn, TReads, TThreading>        TReadsFile;
    typedef FormattedFile<Bam, Output, TContigNames>                TOutputFile;

    typedef typename TReads::TSeqs                                  TReadSeqs;
    typedef typename Value<TReadSeqs>::Type                         TReadSeq;
    typedef typename Size<TReadSeqs>::Type                          TReadSeqsSize;
    typedef String<TReadSeqsSize>                                   TSeedsCount;

    typedef ReadsContext<TSpec, TConfig>                            TReadsContext;

    typedef StringSet<TReadSeqs, Segment<TReadSeqs> >               TSeeds;
    typedef Tuple<TSeeds, TConfig::BUCKETS>                         TSeedsBuckets;

    typedef Hit<TIndexSize, HammingDistance>                        THit;
    typedef String<THit>                                            THits;
    typedef Tuple<THits, TConfig::BUCKETS>                          THitsBuckets;
    typedef String<TIndexSize>                                      THitsCounts;
    typedef ConcurrentAppender<THits>                               THitsAppender;

    typedef StringSet<TSeedsCount, Owner<ConcatDirect<> > >         TRanks;
    typedef Tuple<TRanks, TConfig::BUCKETS>                         TRanksBuckets;

    typedef Limits<TContigsSize, TContigsLen, TContigsSum>          TMatchSpec;
    typedef Match<TMatchSpec>                                       TMatch;
    typedef String<TMatch>                                          TMatches;
    typedef ConcurrentAppender<TMatches>                            TMatchesAppender;
    typedef StringSet<TMatches, Segment<TMatches> >                 TMatchesSet;
    typedef typename Suffix<TMatches>::Type                         TMates;
    typedef StringSet<TMates, Segment<TMates> >                     TMatesSet;

    typedef typename Position<TMatches>::Type                       TMatchesPos;
    typedef String<TMatchesPos>                                     TMatchesPositions;
    typedef ModifiedString<TMatches, ModPos<TMatchesPositions> >    TMatchesView;
    typedef StringSet<TMatchesView, Segment<TMatchesView> >         TMatchesViewSet;
    typedef ModifiedString<TMatchesView, ModPos<TMatchesPositions> > TMatchesViewView;
    typedef String<double>                                          TMatchesProbs;

    typedef String<CigarElement<> >                                 TCigar;
    typedef StringSet<TCigar, Segment<TCigar> >                     TCigarSet;
    typedef StringSetLimits<TCigarSet>::Type                        TCigarLimits;
};

// ----------------------------------------------------------------------------
// Mapper Stats
// ----------------------------------------------------------------------------

template <typename TValue>
struct Stats
{
    TValue loadContigs;
    TValue loadReads;
    TValue collectSeeds;
    TValue findSeeds;
    TValue classifyReads;
    TValue rankSeeds;
    TValue extendHits;
    TValue sortMatches;
    TValue compactMatches;
    TValue selectPairs;
    TValue verifyMatches;
    TValue alignMatches;
    TValue writeMatches;

    unsigned long loadedReads;
    unsigned long mappedReads;
    unsigned long pairedReads;
    unsigned long rescuedReads;

    Stats() :
        loadContigs(0),
        loadReads(0),
        collectSeeds(0),
        findSeeds(0),
        classifyReads(0),
        rankSeeds(0),
        extendHits(0),
        sortMatches(0),
        compactMatches(0),
        selectPairs(0),
        verifyMatches(0),
        alignMatches(0),
        writeMatches(0),
        loadedReads(0),
        mappedReads(0),
        pairedReads(0),
        rescuedReads(0)
    {}
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig = void>
struct Mapper
{
    typedef MapperTraits<TSpec, TConfig>    Traits;

    Options const &                     options;
    Timer<double>                       timer;
    Stats<double>                       stats;

    unsigned                            libraryLength;
    unsigned                            libraryDev;

    typename Traits::TContigs           contigs;
    typename Traits::TIndex             index;
    typename Traits::TReads             reads;

    typename Traits::TReadsFile         readsFile;
    typename Traits::TOutputFile        outputFile;

    typename Traits::TReadsContext      ctx;
    typename Traits::TSeedsBuckets      seeds;
    typename Traits::THitsBuckets       hits;
    typename Traits::TRanksBuckets      ranks;

    typename Traits::TMatches           matchesByCoord;
    typename Traits::TMatchesSet        matchesSetByCoord;

    typename Traits::TMatchesPositions  matchesPositions;
    typename Traits::TMatchesPositions  primaryMatchesPositions;
    typename Traits::TMatchesView       matchesByErrors;
    typename Traits::TMatchesViewSet    matchesSetByErrors;
    typename Traits::TMatchesViewSet    optimalMatchesSet;
    typename Traits::TMatchesViewSet    suboptimalMatchesSet;
    typename Traits::TMatchesViewView   primaryMatches;
    typename Traits::TMatchesProbs      primaryMatchesProbs;

    typename Traits::TMates             matesByCoord;
    typename Traits::TMatesSet          matesSetByCoord;

    typename Traits::TCigar             cigars;
    typename Traits::TCigarSet          cigarSet;

    Mapper(Options const & options) :
        options(options),
        libraryLength(),
        libraryDev(),
        readsFile(options.readsCount)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function configureThreads()
// ----------------------------------------------------------------------------
// Sets the number of threads that OpenMP can spawn.

template <typename TSpec, typename TConfig>
inline void configureThreads(Mapper<TSpec, TConfig> & me)
{
    omp_set_num_threads(me.options.threadsCount);

    if (me.options.verbose > 0)
        std::cerr << "Threads count:\t\t\t" << omp_get_max_threads() << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadContigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadContigs(Mapper<TSpec, TConfig> & me)
{
    start(me.timer);
    try
    {
        if (!open(me.contigs, toCString(me.options.contigsIndexFile), OPEN_RDONLY))
            throw RuntimeError("Error while opening reference file.");
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }
    stop(me.timer);
    me.stats.loadContigs += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Loading reference:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadContigsIndex()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadContigsIndex(Mapper<TSpec, TConfig> & me)
{
    start(me.timer);
    try
    {
        if (!open(me.index, toCString(me.options.contigsIndexFile), OPEN_RDONLY))
            throw RuntimeError("Error while opening reference index file.");
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference index.");
    }
    stop(me.timer);
    me.stats.loadContigs += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Loading reference index:\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function openReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void openReads(Mapper<TSpec, TConfig> & me)
{
    _openReadsImpl(me, typename TConfig::TSequencing());
}

template <typename TSpec, typename TConfig>
inline void _openReadsImpl(Mapper<TSpec, TConfig> & me, SingleEnd)
{
    if (!open(me.readsFile, toCString(me.options.readsFile.i1)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void _openReadsImpl(Mapper<TSpec, TConfig> & me, PairedEnd)
{
    if (!open(me.readsFile, toCString(me.options.readsFile.i1), toCString(me.options.readsFile.i2)))
        throw RuntimeError("Error while opening reads file.");
}

// ----------------------------------------------------------------------------
// Function closeReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void closeReads(Mapper<TSpec, TConfig> & me)
{
    close(me.readsFile);
}

// ----------------------------------------------------------------------------
// Function loadReads()
// ----------------------------------------------------------------------------
// Loads one block of reads.

template <typename TSpec, typename TConfig>
inline void loadReads(Mapper<TSpec, TConfig> & me)
{
    typedef typename MapperTraits<TSpec, TConfig>::TMatch   TMatch;

    start(me.timer);

    readRecords(me.reads, me.readsFile);

    if (maxLength(me.reads.seqs, typename TConfig::TThreading()) > MemberLimits<TMatch, ReadSize>::VALUE)
        throw RuntimeError("Maximum read length exceeded.");

    // Append reverse complemented reads.
    appendReverseComplement(me.reads);

    stop(me.timer);

    me.stats.loadReads += getValue(me.timer);
    me.stats.loadedReads += getReadsCount(me.reads.seqs);

    if (me.options.verbose > 1)
    {
        std::cerr << "Loading reads:\t\t\t" << me.timer << std::endl;
        std::cerr << "Reads count:\t\t\t" << getReadsCount(me.reads.seqs) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function clearReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clearReads(Mapper<TSpec, TConfig> & me)
{
    clear(me.reads);
}

// ----------------------------------------------------------------------------
// Function openOutputFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void openOutputFile(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef typename TTraits::TContigSeqs           TContigSeqs;
    typedef typename Value<TContigSeqs>::Type       TContigSeq;

    bool opened = false;

    if (empty(me.options.outputFile))
    {
        // Output to cout.
        if (me.options.uncompressedBam)
        {
            // Turn off BAM compression.
            setFormat(me.outputFile, me.options.outputFormat);
            opened = _open(me.outputFile, std::cout, Nothing(), False());
        }
        else
        {
            opened = open(me.outputFile, std::cout, me.options.outputFormat);
        }
    }
    else
    {
        // Output to file.
        opened = open(me.outputFile, toCString(me.options.outputFile), OPEN_WRONLY | OPEN_CREATE);
    }

    if (!opened) throw RuntimeError("Error while opening output file.");

    setContigNames(context(me.outputFile), me.contigs.names);

    // Fill contig lengths.
    resize(contigLengths(context(me.outputFile)), length(me.contigs.seqs));
    transform(contigLengths(context(me.outputFile)), me.contigs.seqs, [](TContigSeq const & seq) { return length(seq); });

    // Write header.
    BamHeader header;
    fillHeader(header, me.options);
    writeHeader(me.outputFile, header);
}

// ----------------------------------------------------------------------------
// Function closeOutputFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void closeOutputFile(Mapper<TSpec, TConfig> & me)
{
    close(me.outputFile);
}

// ----------------------------------------------------------------------------
// Function initSeeds()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void initSeeds(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        setHost(me.seeds[bucketId], readSeqs);
}

// ----------------------------------------------------------------------------
// Function clearSeeds()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clearSeeds(Mapper<TSpec, TConfig> & me)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
    {
        clear(me.seeds[bucketId]);
        clear(me.ranks[bucketId]);
        shrinkToFit(me.seeds[bucketId]);
        shrinkToFit(me.ranks[bucketId]);
    }
}

// ----------------------------------------------------------------------------
// Function initReadsContext()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void initReadsContext(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    clear(me.ctx);
    resize(me.ctx, readSeqs);
}

// ----------------------------------------------------------------------------
// Function collectSeeds()
// ----------------------------------------------------------------------------
// Collects seeds from all reads.

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TReadSeqs>
inline void collectSeeds(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>                TTraits;
    typedef SeedsCollector<Counter, TTraits>            TCounter;
    typedef SeedsCollector<void, TTraits>               TFiller;

    typename TTraits::TSeedsCount seedsCounts;

    start(me.timer);
    TCounter counter(me.ctx, me.seeds[ERRORS], seedsCounts, ERRORS, readSeqs, me.options);
    TFiller filler(me.ctx, me.seeds[ERRORS], seedsCounts, ERRORS, readSeqs, me.options);
    stop(me.timer);
    me.stats.collectSeeds += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cerr << "Seeding time:\t\t\t" << me.timer << std::endl;
        std::cerr << "Seeds count:\t\t\t" << length(me.seeds[ERRORS]) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function findSeeds()
// ----------------------------------------------------------------------------

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TBucketId>
inline void findSeeds(Mapper<TSpec, TConfig> & me, TBucketId bucketId)
{
    start(me.timer);
    if (ERRORS > 0)
    {
        // Estimate the number of hits.
        reserve(me.hits[bucketId], lengthSum(me.seeds[bucketId]) * Power<ERRORS, 2>::VALUE, Exact());
        if (me.options.sensitivity == FULL)
            _findSeedsImpl(me, me.hits[bucketId], me.seeds[bucketId], ERRORS, EditDistance());
        else
            _findSeedsImpl(me, me.hits[bucketId], me.seeds[bucketId], ERRORS, HammingDistance());
    }
    else
    {
        reserve(me.hits[bucketId], length(me.seeds[bucketId]), Exact());
        _findSeedsImpl(me, me.hits[bucketId], me.seeds[bucketId], ERRORS, Exact());
    }
    stop(me.timer);
    me.stats.findSeeds += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cerr << "Filtering time:\t\t\t" << me.timer << std::endl;
        std::cerr << "Hits count:\t\t\t" <<
               countHits<unsigned long>(me.hits[bucketId], typename TConfig::TThreading()) << std::endl;
    }
}

template <typename TSpec, typename TConfig, typename THits, typename TSeeds, typename TErrors, typename TDistance>
inline void _findSeedsImpl(Mapper<TSpec, TConfig> & me, THits & hits, TSeeds & seeds, TErrors errors, TDistance)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef FilterDelegate<TSpec, TTraits>          TDelegate;
    typedef typename TTraits::THitsAppender         TAppender;

    TAppender appender(hits);
    TDelegate delegate(appender);

    // Find hits.
    find(me.index, seeds, errors, delegate, Backtracking<TDistance>(), typename TConfig::TThreading());

    // Sort the hits by seedId.
    if (IsSameType<typename TConfig::TThreading, Parallel>::VALUE)
        sortHits(hits, typename TConfig::TThreading());
}

// ----------------------------------------------------------------------------
// Function classifyReads()
// ----------------------------------------------------------------------------
// Classifies the reads by hardness.

template <typename TSpec, typename TConfig>
inline void classifyReads(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>                TTraits;
    typedef ReadsClassifier<TSpec, TTraits>             TClassifier;

    start(me.timer);
    TClassifier classifier(me.ctx, me.hits[0], me.seeds[0], me.options);
    stop(me.timer);
    me.stats.classifyReads += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cerr << "Classification time:\t\t" << me.timer << std::endl;
        std::cerr << "Hits count:\t\t\t" <<
               countHits<unsigned long>(me.hits[0], typename TConfig::TThreading()) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function rankSeeds()
// ----------------------------------------------------------------------------
// Rank the seeds in all buckets.

template <typename TSpec, typename TConfig>
inline void rankSeeds(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef SeedsRanker<TSpec, TTraits>     TSeedsRanker;

    typename TTraits::THitsCounts hitsCounts;

    start(me.timer);
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        TSeedsRanker ranker(hitsCounts, me.ranks[bucketId], me.seeds[bucketId], me.hits[bucketId], me.options);
    stop(me.timer);
    me.stats.rankSeeds += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Ranking time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function clearHits()
// ----------------------------------------------------------------------------
// Clears the hits in all buckets.

template <typename TSpec, typename TConfig>
inline void clearHits(Mapper<TSpec, TConfig> & me)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
    {
        clear(me.hits[bucketId]);
        shrinkToFit(me.hits[bucketId]);
    }
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------
// Counts the hits in all buckets.

template <typename TSpec, typename TConfig>
inline unsigned long countHits(Mapper<TSpec, TConfig> const & me)
{
    unsigned long hitsCount = 0;

    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        hitsCount += countHits<unsigned long>(me.hits[bucketId], typename TConfig::TThreading());

    return hitsCount;
}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------
// Extends the hits in a bucket.

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TBucketId>
inline void extendHits(Mapper<TSpec, TConfig> & me, TBucketId bucketId)
{
    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef HitsExtender<TSpec, TTraits>        THitsExtender;
    typedef typename TTraits::TMatchesAppender  TMatchesAppender;

    start(me.timer);
    TMatchesAppender appender(me.matchesByCoord);
    THitsExtender extender(me.ctx, appender, me.contigs.seqs,
                           me.seeds[bucketId], me.hits[bucketId], me.ranks[bucketId], ERRORS,
                           indexSA(me.index), me.options);
    stop(me.timer);
    me.stats.extendHits += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cerr << "Extension time:\t\t\t" << me.timer << std::endl;
        std::cerr << "Matches count:\t\t\t" << length(me.matchesByCoord) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function reserveMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reserveMatches(Mapper<TSpec, TConfig> & me)
{
    // Estimate the number of matches.
    reserve(me.matchesByCoord, countHits(me) / 3);
}

// ----------------------------------------------------------------------------
// Function aggregateMatches()
// ----------------------------------------------------------------------------
// Aggregate matches by readId.

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void aggregateMatches(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef typename TTraits::TMatch        TMatch;

    start(me.timer);
    // Sort matches by readId and bucket them.
    sort(me.matchesByCoord, MatchSorter<TMatch, ReadId>(), typename TConfig::TThreading());
    setHost(me.matchesSetByCoord, me.matchesByCoord);
    bucket(me.matchesSetByCoord, Getter<TMatch, ReadId>(), getReadsCount(readSeqs), typename TConfig::TThreading());
    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Sorting time:\t\t\t" << me.timer << std::endl;

    // Remove duplicate matches (sorts the matches by genomic coordinate).
    start(me.timer);
    removeDuplicates(me.matchesSetByCoord, typename TConfig::TThreading());
    stop(me.timer);
    me.stats.compactMatches += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cerr << "Compaction time:\t\t" << me.timer << std::endl;
        std::cerr << "Matches count:\t\t\t" << lengthSum(me.matchesSetByCoord) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function clearMatches()
// ----------------------------------------------------------------------------
// Clears all matches.

template <typename TSpec, typename TConfig>
inline void clearMatches(Mapper<TSpec, TConfig> & me)
{
    clear(me.matchesSetByCoord);
    clear(me.optimalMatchesSet);
    clear(me.suboptimalMatchesSet);

    clear(me.matchesByCoord);
    shrinkToFit(me.matchesByCoord);
    clear(me.primaryMatches);
//    shrinkToFit(me.primaryMatches);

    clear(me.matesSetByCoord);
}

// ----------------------------------------------------------------------------
// Function rankMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void rankMatches(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>                    TTraits;
    typedef typename TTraits::TMatch                        TMatch;
    typedef typename TTraits::TMatchesPos                   TMatchesPos;
    typedef typename TTraits::TMatchesSet                   TMatchesSet;
    typedef typename TTraits::TMatchesViewSet               TMatchesViewSet;
    typedef typename Value<TMatchesSet const>::Type         TMatchesSetValue;
    typedef typename Iterator<TMatchesSetValue const, Standard>::Type TMatchesSetValueIt;
    typedef typename Value<TMatchesViewSet const>::Type     TMatchesViewSetValue;
    typedef typename Iterator<TMatchesViewSet const, Standard>::Type TMatchesViewSetIt;
    typedef typename Size<TReadSeqs>::Type                  TReadId;
    typedef typename Size<TMatchesSetValue>::Type           TMatchesSize;
    typedef std::uniform_int_distribution<TMatchesSize>     TMatchesRnd;
    typedef String<unsigned>                                TLibraryLengths;

    start(me.timer);
    // Create a position modifier of the matches from the identity permutation.
    assign(me.matchesPositions, seqan::Range<TMatchesSize>(0, length(me.matchesByCoord)), Exact());
    setHost(me.matchesByErrors, me.matchesByCoord);
    setCargo(me.matchesByErrors, me.matchesPositions);

    // Bucket matches in the position modifier.
    setHost(me.matchesSetByErrors, me.matchesByErrors);
    assign(stringSetLimits(me.matchesSetByErrors), stringSetLimits(me.matchesSetByCoord), Exact());
    assign(stringSetPositions(me.matchesSetByErrors), stringSetPositions(me.matchesSetByCoord), Exact());

    // Sort matches by errors.
    forEach(me.matchesSetByErrors, sortMatches<TMatchesViewSetValue, Errors>, typename TTraits::TThreading());

    // Select all co-optimal matches.
    assign(me.optimalMatchesSet, me.matchesSetByErrors);
    clipMatches(me.optimalMatchesSet, countMatchesInBestStratum<TMatchesViewSetValue>, typename TTraits::TThreading());

    // Select all sub-optimal matches.
    assign(me.suboptimalMatchesSet, me.matchesSetByErrors);
    clipMatches(me.suboptimalMatchesSet, [&](TMatchesViewSetValue const & matches)
    {
        if (empty(matches)) return TMatchesSize(0);

        TReadId readId = getMember(front(matches), ReadId());

        return countMatchesInStrata(matches, getReadStrata<TMatch>(me.options, length(readSeqs[readId])));
    },
    typename TTraits::TThreading());

    // Append an invalid match to matches by coord.
    resize(me.matchesByCoord, length(me.matchesByCoord) + 1, Exact());
    setInvalid(back(me.matchesByCoord));
    // Update matches by errors.
    resize(me.matchesPositions, length(me.matchesPositions) + 1, Exact());
    setPosition(me.matchesByErrors, length(me.matchesByErrors) - 1, length(me.matchesByCoord) - 1);

    // Initialize primary matches.
    setHost(me.primaryMatches, me.matchesByErrors);
    assign(me.primaryMatchesPositions, stringSetPositions(me.matchesSetByErrors), Exact());
    setCargo(me.primaryMatches, me.primaryMatchesPositions);

    // Choose primary matches among best matches.
    iterate(me.optimalMatchesSet, [&](TMatchesViewSetIt const & matchesIt)
    {
        // Use one generator per thread.
        std::default_random_engine generator;

        TReadId readId = position(matchesIt, me.optimalMatchesSet);
        TMatchesViewSetValue const & matches = value(matchesIt);

        // Set unmapped reads as invalid.
        if (empty(matches))
        {
            setPosition(me.primaryMatches, readId, length(me.matchesByErrors) - 1);
        }
        // Choose match at random.
        else
        {
            TMatchesRnd rnd(0, length(matches) - 1);
            setPosition(me.primaryMatches, readId, position(me.primaryMatches, readId) + rnd(generator));
        }
    },
    Standard(), typename TTraits::TThreading());

    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);
    if (me.options.verbose > 1)
        std::cerr << "Sorting time:\t\t\t" << me.timer << std::endl;

    // Update mapped reads.
    transform(me.ctx.mapped, me.primaryMatches, isValid<typename TTraits::TMatchSpec>, typename TTraits::TThreading());

    if (me.options.verbose > 0)
    {
        unsigned long mappedReads = count(me.ctx.mapped, true, typename TTraits::TThreading());
        me.stats.mappedReads += mappedReads;

        if (me.options.verbose > 1)
            std::cerr << "Mapped reads:\t\t\t" << mappedReads << std::endl;
    }

    if (IsSameType<typename TConfig::TSequencing, SingleEnd>::VALUE) return;

    start(me.timer);
    // Estimate library mean length and deviation if one of them was not provided.
    if (!me.options.libraryLength || !me.options.libraryDev)
    {
        // Collect library lengths from unique optimal pairs.
        TLibraryLengths libraryLengths;
        reserve(libraryLengths, getPairsCount(readSeqs), Exact());
        ConcurrentAppender<TLibraryLengths> libraryLengthsAppender(libraryLengths);
        forAllMatchesPairs(me.optimalMatchesSet, readSeqs, [&](TMatchesViewSetValue const & firstMatches, TMatchesViewSetValue const & secondMatches)
        {
            if (length(firstMatches) == 1 && length(secondMatches) == 1)
            {
                TMatch const & firstMatch = front(firstMatches);
                TMatch const & secondMatch = front(secondMatches);

                if (contigEqual(firstMatch, secondMatch) && orientationProper(firstMatch, secondMatch))
                    appendValue(libraryLengthsAppender, getLibraryLength(firstMatch, secondMatch), Insist(), typename TTraits::TThreading());
            }
        },
        typename TTraits::TThreading());

        // If library mean length and deviation cannot be estimated proceed as single-ended.
        if (empty(libraryLengths)) return;

        // Remove library outliers > 6 * median.
        unsigned libraryMedian = nthElement(libraryLengths, length(libraryLengths) / 2, typename TTraits::TThreading());
        removeIf(libraryLengths, std::bind2nd(std::greater<unsigned>(), 6.0 * libraryMedian), typename TTraits::TThreading());

        // If library mean length and deviation cannot be estimated proceed as single-ended.
        if (empty(libraryLengths)) return;

        // Compute library mean.
        unsigned librarySum = accumulate(libraryLengths, 0u, typename TTraits::TThreading());
        float libraryMean = std::max(librarySum / static_cast<float>(length(libraryLengths)), 1.0f);

        // Compute library standard deviation.
        String<float> libraryDiffs;
        resize(libraryDiffs, length(libraryLengths), Exact());
        transform(libraryDiffs, libraryLengths, std::bind2nd(std::minus<float>(), libraryMean), typename TTraits::TThreading());
        float librarySqSum = innerProduct(libraryDiffs, 0.0f, typename TTraits::TThreading());
        float libraryDev = std::max(std::sqrt(librarySqSum / static_cast<float>(length(libraryLengths))), 1.0f);

        if (me.options.verbose > 1)
        {
            std::cerr << "Library median:\t\t\t" << libraryMedian << std::endl;
            std::cerr << "Library mean:\t\t\t" << libraryMean << std::endl;
            std::cerr << "Library stddev:\t\t\t" << libraryDev << std::endl;
        }

        // Set library mean and error as just computed.
        me.libraryLength = libraryMean;
        me.libraryDev = libraryDev;
    }

    // Overwrite library mean and error if provided in input.
    if (me.options.libraryLength)
        me.libraryLength = me.options.libraryLength;
    if (me.options.libraryDev)
        me.libraryDev = me.options.libraryDev;

    resize(me.primaryMatchesProbs, getReadsCount(readSeqs), 0.0, Exact());

    // Enumerate feasible pairs.
    forAllMatchesPairs(me.matchesSetByCoord, readSeqs, [&](TMatchesSetValue const & firstMatches, TMatchesSetValue const & secondMatches)
    {
        TReadId firstId = getMember(front(firstMatches), ReadId());
        TReadId secondId = getMember(front(secondMatches), ReadId());

        SEQAN_ASSERT(isMapped(me.ctx, firstId));
        SEQAN_ASSERT(isMapped(me.ctx, secondId));

        double firstMatchOptimalRate = toErrorRate(readSeqs, firstId, getMinErrors(me.ctx, firstId));
        double secondMatchOptimalRate = toErrorRate(readSeqs, secondId, getMinErrors(me.ctx, secondId));

        auto firstBestCount = countMatchesInBestStratum(me.optimalMatchesSet[firstId]);
        auto firstSubCount = length(me.suboptimalMatchesSet[firstId]) - firstBestCount;

        auto secondBestCount = countMatchesInBestStratum(me.optimalMatchesSet[secondId]);
        auto secondSubCount = length(me.suboptimalMatchesSet[secondId]) - secondBestCount;

        // First mate match with all second mate matches.
        Pair<TMatchesSetValueIt, double> firstPrimary =
        findPrimaryMatch(firstMatches, secondMatches,
                         firstMatchOptimalRate, secondMatchOptimalRate,
                         secondBestCount, secondSubCount,
                         readSeqs, me.contigs.seqs,
                         me.libraryLength, me.libraryDev);

        // Second mate match with all first mate matches.
        Pair<TMatchesSetValueIt, double> secondPrimary =
        findPrimaryMatch(secondMatches, firstMatches,
                         secondMatchOptimalRate, firstMatchOptimalRate,
                         firstBestCount, firstSubCount,
                         readSeqs, me.contigs.seqs,
                         me.libraryLength, me.libraryDev);

        // No feasible pair found.
        if (atEnd(getValueI1(firstPrimary), firstMatches) || atEnd(getValueI1(secondPrimary), secondMatches))
            return;

        // Get matches by coords positions.
        TMatchesPos firstPosByCoord = stringSetPositions(me.matchesSetByCoord)[firstId] +
                                      position(getValueI1(firstPrimary), firstMatches);
        TMatchesPos secondPosByCoord = stringSetPositions(me.matchesSetByCoord)[secondId] +
                                       position(getValueI1(secondPrimary), secondMatches);

        // Translate matches by coords positions into in matches by errors positions.
        auto firstPosBegin = begin(cargo(me.matchesByErrors), Standard()) + stringSetPositions(me.matchesSetByErrors)[firstId];
        auto firstPosEnd = firstPosBegin + length(me.matchesSetByErrors[firstId]);
        auto firstPos = std::find(firstPosBegin, firstPosEnd, firstPosByCoord);
        auto firstPosByErrors = position(firstPos, cargo(me.matchesByErrors));

        auto secondPosBegin = begin(cargo(me.matchesByErrors), Standard()) + stringSetPositions(me.matchesSetByErrors)[secondId];
        auto secondPosEnd = secondPosBegin + length(me.matchesSetByErrors[secondId]);
        auto secondPos = std::find(secondPosBegin, secondPosEnd, secondPosByCoord);
        auto secondPosByErrors = position(secondPos, cargo(me.matchesByErrors));

        // Set primary matches positions.
        setPosition(me.primaryMatches, firstId, firstPosByErrors);
        setPosition(me.primaryMatches, secondId, secondPosByErrors);
        SEQAN_ASSERT(isEqual(value(getValueI1(firstPrimary)), me.primaryMatches[firstId]));
        SEQAN_ASSERT(isEqual(value(getValueI1(secondPrimary)), me.primaryMatches[secondId]));

        // Set primary matches probabilities.
        me.primaryMatchesProbs[firstId] = getValueI2(firstPrimary);
        me.primaryMatchesProbs[secondId] = getValueI2(secondPrimary);

        // Set reads as properly paired.
//        if (isProper(me.primaryMatches[firstId], me.primaryMatches[secondId], me.libraryLength, me.libraryDev))
//        {
            setPaired(me.ctx, firstId);
            setPaired(me.ctx, secondId);
//        }
    },
    typename TTraits::TThreading());

    stop(me.timer);
    me.stats.selectPairs += getValue(me.timer);

    // Update paired reads.
    if (me.options.verbose > 0)
    {
        unsigned long pairedReads = count(me.ctx.paired, true, typename TTraits::TThreading());
        me.stats.pairedReads += pairedReads;

        if (me.options.verbose > 1)
        {
            std::cerr << "Pairing time:\t\t\t" << me.timer << std::endl;
            std::cerr << "Paired reads:\t\t\t" << pairedReads << std::endl;
        }
    }
}

// ----------------------------------------------------------------------------
// Function verifyMatches()
// ----------------------------------------------------------------------------
// Verifies all mates in within the insert window of their matches.

template <typename TSpec, typename TConfig>
inline void verifyMatches(Mapper<TSpec, TConfig> & me)
{
    _verifyMatchesImpl(me, typename TConfig::TSequencing());
}

template <typename TSpec, typename TConfig, typename TSequencing>
inline void _verifyMatchesImpl(Mapper<TSpec, TConfig> & /* me */, TSequencing) {}

template <typename TSpec, typename TConfig>
inline void _verifyMatchesImpl(Mapper<TSpec, TConfig> & me, PairedEnd)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef AnchorsVerifier<TSpec, TTraits>         TMatchesVerifier;
    typedef typename TTraits::TMatchesAppender      TMatchesAppender;
    typedef typename TTraits::TMatch                TMatch;
    typedef typename TTraits::TMatesSet             TMatesSet;
    typedef typename Iterator<TMatesSet>::Type      TMatesSetIt;

    start(me.timer);
    unsigned long anchorsCount = length(me.matchesByCoord);
    TMatchesAppender appender(me.matchesByCoord);
    TMatchesVerifier verifier(me.ctx, appender,
                              me.contigs.seqs, me.reads.seqs,
                              me.optimalMatchesSet,
                              me.libraryLength, me.libraryDev,
                              me.options);

    // Sort matches by readId and bucket them.
    me.matesByCoord = suffix(me.matchesByCoord, anchorsCount);
    sort(me.matesByCoord, MatchSorter<TMatch, ReadId>(), typename TConfig::TThreading());
    setHost(me.matesSetByCoord, me.matesByCoord);
    bucket(me.matesSetByCoord, Getter<TMatch, ReadId>(), getReadsCount(me.reads.seqs), typename TConfig::TThreading());

    resize(cargo(me.matchesByErrors), length(host(me.matchesByErrors)), Exact());
    iota(suffix(cargo(me.matchesByErrors), anchorsCount), anchorsCount);

    // Update primary matches with mates.
    iterate(me.matesSetByCoord, [&](TMatesSetIt const & matesSetIt)
    {
        auto mateId = position(matesSetIt, me.matesSetByCoord);
        auto anchorId = getMateId(me.reads.seqs, mateId);
        auto const & mates = value(matesSetIt);

        if (!empty(mates))
        {
            SEQAN_ASSERT(isMapped(me.ctx, anchorId));
            SEQAN_ASSERT_NOT(isPaired(me.ctx, anchorId));
            SEQAN_ASSERT(isValid(me.primaryMatches[anchorId]));

            SEQAN_ASSERT_NOT(isMapped(me.ctx, mateId));
            SEQAN_ASSERT_NOT(isPaired(me.ctx, anchorId));
            SEQAN_ASSERT_NOT(isValid(me.primaryMatches[mateId]));

            setPosition(me.primaryMatches, mateId, stringSetPositions(me.matesSetByCoord)[mateId] + beginPosition(me.matesByCoord));
            SEQAN_ASSERT(isEqual(me.primaryMatches[mateId], front(mates)));

            setMapped(me.ctx, mateId);
            setPaired(me.ctx, mateId);
            setPaired(me.ctx, anchorId);

            // Set primary matches probabilities.
            double errorRate = getErrorRate(me.primaryMatches[anchorId], me.reads.seqs);
            auto bestCount = countMatchesInBestStratum(me.optimalMatchesSet[anchorId]);
            auto subCount = length(me.suboptimalMatchesSet[anchorId]) - bestCount;
            me.primaryMatchesProbs[anchorId] = getMatchProb(errorRate, errorRate, bestCount, subCount);
            me.primaryMatchesProbs[mateId] = me.primaryMatchesProbs[anchorId];

            SEQAN_ASSERT(isMapped(me.ctx, mateId));
            SEQAN_ASSERT(isPaired(me.ctx, anchorId));
            SEQAN_ASSERT(isValid(me.primaryMatches[mateId]));
        }
    },
    Standard(), Serial());

    stop(me.timer);
    me.stats.verifyMatches += getValue(me.timer);

    if (me.options.verbose > 0)
    {
        me.stats.rescuedReads += length(me.matchesByCoord) - anchorsCount;
    }
    if (me.options.verbose > 1)
    {
        std::cerr << "Rescued reads:\t\t\t" << length(me.matchesByCoord) - anchorsCount << std::endl;
        std::cerr << "Verification time:\t\t" << me.timer << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function alignMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void alignMatches(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef MatchesAligner<LinearGaps, TTraits>     TLinearAligner;
    typedef MatchesAligner<AffineGaps , TTraits>    TAffineAligner;

    start(me.timer);
    setHost(me.cigarSet, me.cigars);
    typename TTraits::TCigarLimits cigarLimits;

    if (me.options.rabema)
        TLinearAligner aligner(me.cigarSet, cigarLimits, me.primaryMatches, me.contigs.seqs, me.reads.seqs, me.options);
    else
        TAffineAligner aligner(me.cigarSet, cigarLimits, me.primaryMatches, me.contigs.seqs, me.reads.seqs, me.options);

    stop(me.timer);
    me.stats.alignMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Alignment time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function clearAlignments()
// ----------------------------------------------------------------------------
// Clears all cigars.

template <typename TSpec, typename TConfig>
inline void clearAlignments(Mapper<TSpec, TConfig> & me)
{
    clear(me.cigars);
    clear(me.cigarSet);
    shrinkToFit(me.cigarSet);
}

// ----------------------------------------------------------------------------
// Function writeMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void writeMatches(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef MatchesWriter<TSpec, TTraits>       TMatchesWriter;

    start(me.timer);
    TMatchesWriter writer(me.outputFile,
                          me.suboptimalMatchesSet,
                          me.primaryMatches, me.primaryMatchesProbs, me.cigarSet,
                          me.ctx, me.reads,
                          me.options);
    stop(me.timer);
    me.stats.writeMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Output time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void mapReads(Mapper<TSpec, TConfig> & me)
{
    _mapReadsImpl(me, me.reads.seqs, typename TConfig::TStrategy());
}

// ----------------------------------------------------------------------------
// Function _mapReadsImpl(); All
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs, All)
{
    initReadsContext(me, readSeqs);
    initSeeds(me, readSeqs);

    collectSeeds<0>(me, readSeqs);
    findSeeds<0>(me, 0);
    classifyReads(me);
    collectSeeds<1>(me, readSeqs);
    collectSeeds<2>(me, readSeqs);
    findSeeds<1>(me, 1);
    if (me.options.sensitivity == LOW)
        findSeeds<1>(me, 2);
    else
        findSeeds<2>(me, 2);
    reserveMatches(me);
    extendHits<0>(me, 0);
    extendHits<1>(me, 1);
    extendHits<2>(me, 2);
    clearSeeds(me);
    clearHits(me);
    aggregateMatches(me, readSeqs);
    rankMatches(me, readSeqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
    writeMatches(me);
    clearMatches(me);
    clearAlignments(me);
}

// ----------------------------------------------------------------------------
// Function _mapReadsImpl(); Strata
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs, Strata)
{
    initReadsContext(me, readSeqs);
    initSeeds(me, readSeqs);

    collectSeeds<0>(me, readSeqs);
    findSeeds<0>(me, 0);
    classifyReads(me);
    collectSeeds<1>(me, readSeqs);
    collectSeeds<2>(me, readSeqs);
    findSeeds<0>(me, 1);
    findSeeds<0>(me, 2);
    rankSeeds(me);
    reserveMatches(me);
    extendHits<0>(me, 0);
    extendHits<0>(me, 1);
    extendHits<0>(me, 2);
    clearSeeds(me);
    clearHits(me);

    initSeeds(me, readSeqs);
    collectSeeds<1>(me, readSeqs);
    findSeeds<1>(me, 1);
    collectSeeds<2>(me, readSeqs);
    findSeeds<1>(me, 2);
    rankSeeds(me);
    // TODO(esiragusa): filter out hits with distance < 1.
    extendHits<1>(me, 1);
    extendHits<1>(me, 2);
    clearSeeds(me);
    clearHits(me);

    if (me.options.sensitivity > LOW)
    {
        initSeeds(me, readSeqs);
        collectSeeds<2>(me, readSeqs);
        findSeeds<2>(me, 2);
        rankSeeds(me);
        // TODO(esiragusa): filter out hits with distance < 2.
        extendHits<2>(me, 2);
        clearHits(me);
        clearSeeds(me);
    }

    aggregateMatches(me, readSeqs);
    rankMatches(me, readSeqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
    writeMatches(me);
    clearMatches(me);
    clearAlignments(me);
}

// ----------------------------------------------------------------------------
// Function printStats()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TValue>
inline void printStats(Mapper<TSpec, TConfig> const & me, Timer<TValue> const & timer)
{
    printRuler(std::cerr);

    TValue total = getValue(timer) / 100.0;

    std::cerr << "Total time:\t\t\t" << getValue(timer) << " sec" << std::endl;
    std::cerr << "Genome loading time:\t\t" << me.stats.loadContigs << " sec" << "\t\t" << me.stats.loadContigs / total << " %" << std::endl;
    std::cerr << "Reads loading time:\t\t" << me.stats.loadReads << " sec" << "\t\t" << me.stats.loadReads / total << " %" << std::endl;
    std::cerr << "Seeding time:\t\t\t" << me.stats.collectSeeds << " sec" << "\t\t" << me.stats.collectSeeds / total << " %" << std::endl;
    std::cerr << "Filtering time:\t\t\t" << me.stats.findSeeds << " sec" << "\t\t" << me.stats.findSeeds / total << " %" << std::endl;
    std::cerr << "Classification time:\t\t" << me.stats.classifyReads << " sec" << "\t\t" << me.stats.classifyReads / total << " %" << std::endl;
    if (IsSameType<typename TConfig::TStrategy, Strata>::VALUE)
        std::cerr << "Ranking time:\t\t\t" << me.stats.rankSeeds << " sec" << "\t\t" << me.stats.rankSeeds / total << " %" << std::endl;
    std::cerr << "Extension time:\t\t\t" << me.stats.extendHits << " sec" << "\t\t" << me.stats.extendHits / total << " %" << std::endl;
    std::cerr << "Sorting time:\t\t\t" << me.stats.sortMatches << " sec" << "\t\t" << me.stats.sortMatches / total << " %" << std::endl;
    std::cerr << "Compaction time:\t\t" << me.stats.compactMatches << " sec" << "\t\t" << me.stats.compactMatches / total << " %" << std::endl;
    if (IsSameType<typename TConfig::TSequencing, PairedEnd>::VALUE)
    {
        std::cerr << "Pairing time:\t\t\t" << me.stats.selectPairs << " sec" << "\t\t" << me.stats.selectPairs / total << " %" << std::endl;
        std::cerr << "Verification time:\t\t" << me.stats.verifyMatches << " sec" << "\t\t" << me.stats.verifyMatches / total << " %" << std::endl;
    }
    std::cerr << "Alignment time:\t\t\t" << me.stats.alignMatches << " sec" << "\t\t" << me.stats.alignMatches / total << " %" << std::endl;
    std::cerr << "Output time:\t\t\t" << me.stats.writeMatches << " sec" << "\t\t" << me.stats.writeMatches / total << " %" << std::endl;

    printRuler(std::cerr);

    double totalReads = me.stats.loadedReads / 100.0;
    std::cerr << "Total reads:\t\t\t" << me.stats.loadedReads << std::endl;
    std::cerr << "Mapped reads:\t\t\t" << me.stats.mappedReads << "\t\t" << me.stats.mappedReads / totalReads << " %" << std::endl;
    if (IsSameType<typename TConfig::TSequencing, PairedEnd>::VALUE)
    {
        std::cerr << "Paired reads:\t\t\t" << me.stats.pairedReads << "\t\t" << me.stats.pairedReads / totalReads << " %" << std::endl;
        std::cerr << "Rescued reads:\t\t\t" << me.stats.rescuedReads << "\t\t" << me.stats.rescuedReads / totalReads << " %" << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void runMapper(Mapper<TSpec, TConfig> & me)
{
    Timer<double> timer;

    start(timer);

    configureThreads(me);

    if (me.options.verbose > 1) printRuler(std::cerr);

    loadContigs(me);
    loadContigsIndex(me);

    // Open output file and write header.
    openOutputFile(me);
    openReads(me);

    // Process reads in blocks.
    while (true)
    {
        if (me.options.verbose > 1) printRuler(std::cerr);
        loadReads(me);
        if (empty(me.reads.seqs)) break;
        mapReads(me);
        clearReads(me);
    }

    closeReads(me);
    closeOutputFile(me);

    stop(timer);

    if (me.options.verbose > 0)
        printStats(me, timer);
}

// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------

template <typename TContigsSize, typename TContigsLen, typename TContigsSum,
          typename TThreading, typename TSequencing, typename TStrategy>
inline void spawnMapper(Options const & options,
                        TThreading const & /* tag */,
                        TSequencing const & /* tag */,
                        TStrategy const & /* tag */)
{
    typedef ReadMapperConfig<TThreading, TSequencing, TStrategy, TContigsSize, TContigsLen, TContigsSum>  TConfig;

    Mapper<void, TConfig> mapper(options);
    runMapper(mapper);
}

#endif  // #ifndef APP_YARA_MAPPER_H_
