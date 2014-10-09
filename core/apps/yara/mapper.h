// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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

    CharString          contigsIndexFile;
    Pair<CharString>    readsFile;
    CharString          outputFile;
    TOutputFormat       outputFormat;
    bool                outputSecondary;
    bool                uncompressedBam;
    CharString          readGroup;

    MappingMode         mappingMode;
    float               errorRate;
    float               strataRate;
    bool                quick;

    bool                singleEnd;
    unsigned            libraryLength;
    unsigned            libraryError;
    LibraryOrientation  libraryOrientation;
    TList               libraryOrientationList;
//    bool                anchorOne;

    unsigned            readsCount;
    unsigned            threadsCount;
    unsigned            hitsThreshold;
    bool                rabema;
    unsigned            verbose;

    CharString          commandLine;
    CharString          version;

    Options() :
        outputSecondary(false),
        uncompressedBam(false),
        readGroup("none"),
        mappingMode(STRATA),
        errorRate(0.05f),
        strataRate(0.00f),
        quick(false),
        singleEnd(true),
        libraryLength(200),
        libraryError(200),
        libraryOrientation(FWD_REV),
//        anchorOne(false),
        readsCount(100000),
        threadsCount(1),
        hitsThreshold(300),
        rabema(false),
        verbose(0)
    {
        appendValue(libraryOrientationList, "fwd-rev");
        appendValue(libraryOrientationList, "fwd-fwd");
        appendValue(libraryOrientationList, "rev-rev");
    }
};

// ----------------------------------------------------------------------------
// Mapper Configuration
// ----------------------------------------------------------------------------

template <typename TThreading_      = Parallel,
          typename TSequencing_     = SingleEnd,
          typename TStrategy_       = Strata,
          unsigned BUCKETS_         = 3>
struct ReadMapperConfig
{
    typedef TThreading_     TThreading;
    typedef TSequencing_    TSequencing;
    typedef TStrategy_      TStrategy;
//    typedef TAnchoring_     TAnchoring;

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
//    typedef typename TConfig::TAnchoring                            TAnchoring;

    typedef SeqStore<void, YaraContigsConfig>                       TContigs;
    typedef typename TContigs::TSeqs                                TContigSeqs;
    typedef typename Value<TContigSeqs>::Type                       TContig;
    typedef typename StringSetPosition<TContigSeqs>::Type           TContigsPos;

    typedef Index<YaraContigsFM, YaraIndexSpec>                     TIndex;
    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef typename Fibre<TIndex, FibreSA>::Type                   TSA;

    typedef SeqStore<void, YaraReadsConfig>                         TReads;
    typedef typename If<IsSameType<TSequencing, PairedEnd>,
                        Pair<SeqFileIn>, SeqFileIn>::Type           TReadsFileIn;
    typedef PrefetchedFile<TReadsFileIn, TReads, TThreading>        TReadsFile;
    typedef SmartFile<Bam, Output, YaraContigs>                     TOutputFile;

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

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;
    typedef StringSet<TMatches, Segment<TMatches> >                 TMatchesSet;
    typedef ConcurrentAppender<TMatches>                            TMatchesAppender;

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
    TValue alignMatches;
    TValue writeMatches;

    unsigned long loadedReads;
    unsigned long mappedReads;
    unsigned long pairedReads;

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
        alignMatches(0),
        writeMatches(0),
        loadedReads(0),
        mappedReads(0),
        pairedReads(0)
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

    typename Traits::TContigs           contigs;
    typename Traits::TIndex             index;
    typename Traits::TReads             reads;

    typename Traits::TReadsFile         readsFile;
    typename Traits::TOutputFile        outputFile;

    typename Traits::TReadsContext      ctx;
    typename Traits::TSeedsBuckets      seeds;
    typename Traits::THitsBuckets       hits;
    typename Traits::TRanksBuckets      ranks;

    typename Traits::TMatches           matches;
    typename Traits::TMatchesSet        matchesSet;
    typename Traits::TMatchesSet        bestMatchesSet;
    typename Traits::TMatchesSet        suboptimalMatchesSet;
    typename Traits::TMatches           primaryMatches;

    typename Traits::TCigar             cigars;
    typename Traits::TCigarSet          cigarSet;

    Mapper(Options const & options) :
        options(options),
        readsFile(options.readsCount)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getReadErrors()
// ----------------------------------------------------------------------------
// Returns the absolute number of errors for a given read sequence.

template <typename TReadSeqSize>
inline TReadSeqSize getReadErrors(Options const & options, TReadSeqSize readSeqLength)
{
    return std::min((TReadSeqSize)(readSeqLength * options.errorRate),
                    (TReadSeqSize)MemberLimits<Match<void>, Errors>::VALUE);
}

// ----------------------------------------------------------------------------
// Function getReadStrata()
// ----------------------------------------------------------------------------
// Returns the absolute number of strata for a given read sequence.

template <typename TReadSeqSize>
inline TReadSeqSize getReadStrata(Options const & options, TReadSeqSize readSeqLength)
{
    return std::min((TReadSeqSize)(readSeqLength * options.strataRate),
                    (TReadSeqSize)MemberLimits<Match<void>, Errors>::VALUE);
}

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
        if (!open(me.contigs, toCString(me.options.contigsIndexFile)))
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
        if (!open(me.index, toCString(me.options.contigsIndexFile)))
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
    start(me.timer);

    readRecords(me.reads, me.readsFile);

    if (maxLength(me.reads.seqs, typename TConfig::TThreading()) > MemberLimits<Match<void>, ReadSize>::VALUE)
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

    setNameStore(context(me.outputFile), me.contigs.names);

    // Fill contig lengths.
    resize(sequenceLengths(context(me.outputFile)), length(me.contigs.seqs));
    transform(sequenceLengths(context(me.outputFile)), me.contigs.seqs, [&](TContigSeq const & seq) { return length(seq); });

    // Write header.
    BamHeader header;
    fillHeader(header, me.options);
    writeRecord(me.outputFile, header);
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
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef HitsExtender<TSpec, TTraits>    THitsExtender;

    typename TTraits::TMatchesAppender appender(me.matches);

    start(me.timer);
    THitsExtender extender(me.ctx, appender, me.contigs.seqs,
                           me.seeds[bucketId], me.hits[bucketId], me.ranks[bucketId], ERRORS,
                           indexSA(me.index), me.options);
    stop(me.timer);
    me.stats.extendHits += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cerr << "Extension time:\t\t\t" << me.timer << std::endl;
        std::cerr << "Matches count:\t\t\t" << length(me.matches) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function reserveMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reserveMatches(Mapper<TSpec, TConfig> & me)
{
    // Estimate the number of matches.
    reserve(me.matches, countHits(me) / 3);
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

    // Bucket sort matches by readId.
    start(me.timer);
    setHost(me.matchesSet, me.matches);
    sort(me.matches, MatchSorter<TMatch, ReadId>(), typename TConfig::TThreading());
    bucket(me.matchesSet, Getter<TMatch, ReadId>(), getReadsCount(readSeqs), typename TConfig::TThreading());
    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Sorting time:\t\t\t" << me.timer << std::endl;

    start(me.timer);
    removeDuplicates(me.matchesSet, typename TConfig::TThreading());
    stop(me.timer);
    me.stats.compactMatches += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cerr << "Compaction time:\t\t" << me.timer << std::endl;
        std::cerr << "Matches count:\t\t\t" << lengthSum(me.matchesSet) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function verifyMatches()
// ----------------------------------------------------------------------------
// Verifies all mates in within the insert window of their matches.

//template <typename TSpec, typename TConfig, typename TReadSeqs>
//inline void verifyMatches(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs)
//{
//    _verifyMatchesImpl(me, readSeqs, typename TConfig::TAnchoring());
//}
//
//template <typename TSpec, typename TConfig, typename TReadSeqs, typename TAnchoring>
//inline void _verifyMatchesImpl(Mapper<TSpec, TConfig> & /* me */, TReadSeqs & /* readSeqs */, TAnchoring) {}
//
//template <typename TSpec, typename TConfig, typename TReadSeqs>
//inline void _verifyMatchesImpl(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs, AnchorOne)
//{
//    typedef MapperTraits<TSpec, TConfig>    TTraits;
//    typedef AnchorsVerifier<TSpec, TTraits> TMatchesVerifier;
//
//    start(me.timer);
//    TMatchesVerifier verifier(me.ctx, me.pairs,
//                              me.contigs.seqs, readSeqs,
//                              me.matchesSet, me.options);
//    stop(me.timer);
//
//    if (me.options.verbose > 1)
//    {
//        std::cerr << "Verification time:\t\t" << me.timer << std::endl;
//        std::cerr << "Mates count:\t\t\t" << length(me.pairs) << std::endl;
//        std::cerr << "Mapped pairs:\t\t\t" <<
//                countMappedReads(readSeqs, me.pairs, typename TConfig::TThreading()) << std::endl;
//    }
//}

// ----------------------------------------------------------------------------
// Function clearMatches()
// ----------------------------------------------------------------------------
// Clears all matches.

template <typename TSpec, typename TConfig>
inline void clearMatches(Mapper<TSpec, TConfig> & me)
{
    clear(me.matchesSet);
    clear(me.bestMatchesSet);
    clear(me.suboptimalMatchesSet);

    clear(me.matches);
    shrinkToFit(me.matches);
    clear(me.primaryMatches);
    shrinkToFit(me.primaryMatches);
}

// ----------------------------------------------------------------------------
// Function rankMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void rankMatches(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>                    TTraits;
    typedef typename TTraits::TMatchesSet                   TMatchesSet;
    typedef typename Iterator<TMatchesSet, Standard>::Type  TMatchesIt;
    typedef typename Value<TMatchesSet>::Type               TMatches;
    typedef PairsSelector<TSpec, TTraits>                   TPairsSelector;
    typedef typename Size<TReadSeqs>::Type                  TReadId;

    // Sort matches by errors.
    start(me.timer);
    iterate(me.matchesSet, sortMatches<TMatchesIt, Errors>, Standard(), typename TTraits::TThreading());
//    forEach(me.matchesSet, sortMatches<TMatches, Errors>, typename TTraits::TThreading());
    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);
    if (me.options.verbose > 1)
        std::cerr << "Sorting time:\t\t\t" << me.timer << std::endl;

    // Select all co-optimal matches.
    assign(me.bestMatchesSet, me.matchesSet);
    clipMatches(me.bestMatchesSet, countMatchesInBestStratum<TMatches>, typename TTraits::TThreading());

    // Select all sub-optimal matches.
    assign(me.suboptimalMatchesSet, me.matchesSet);
    clipMatches(me.suboptimalMatchesSet, [&](TMatches const & matches)
    {
        if (empty(matches)) return typename Size<TMatches>::Type(0);

        TReadId readId = getMember(front(matches), ReadId());

        return countMatchesInStrata(matches, getReadStrata(me.options, length(readSeqs[readId])));
    },
    typename TTraits::TThreading());

    // Initialize primary matches.
    resize(me.primaryMatches, getReadsCount(readSeqs), Exact());
    forEach(me.primaryMatches, setInvalid<void>, typename TTraits::TThreading());

    // Try to pair mates.
    if (IsSameType<typename TConfig::TSequencing, PairedEnd>::VALUE)
    {
        start(me.timer);

        // Concordant pairs of first co-optimal match with second sub-optimal match.
        TPairsSelector selectorOptSubConcordant(me.primaryMatches, me.ctx, readSeqs, me.bestMatchesSet, me.suboptimalMatchesSet, me.options);
        // Concordant pairs of first sub-optimal match with second co-optimal match.
        TPairsSelector selectorSubOptConcordant(me.primaryMatches, me.ctx, readSeqs, me.suboptimalMatchesSet, me.bestMatchesSet, me.options);

        // Mark paired mates as properly paired.
        iterate(me.primaryMatches, [&](typename Iterator<TMatches, Standard>::Type & matchesIt)
        {
            if (isValid(*matchesIt)) setPaired(me.ctx, getMember(*matchesIt, ReadId()));
        },
        Standard(), typename TTraits::TThreading());

        // Concordant co-optimal matches on the same chromosome outside of the expected insert size.
        Options pairing = me.options;
        pairing.libraryError = MaxValue<unsigned>::VALUE;
        TPairsSelector selectorOptOptConcordant(me.primaryMatches, me.ctx, readSeqs, me.bestMatchesSet, me.bestMatchesSet, pairing);

        // Any pair of co-optimal matches on the same chromosome.
        pairing.libraryOrientation = ANY;
        pairing.libraryError = MaxValue<unsigned>::VALUE;
        TPairsSelector selectorOptOptAny(me.primaryMatches, me.ctx, readSeqs, me.bestMatchesSet, me.bestMatchesSet, pairing);

        stop(me.timer);
        me.stats.selectPairs += getValue(me.timer);
    }

    // Randomly choose primary matches among co-optimal ones.
    MatchesPicker<TMatches> picker;
    iterate(me.primaryMatches, [&](typename Iterator<TMatches, Standard>::Type & matchesIt)
    {
        if (!isValid(*matchesIt)) *matchesIt = picker(me.bestMatchesSet[position(matchesIt, me.primaryMatches)]);
    },
    Standard(), Serial());

    unsigned long mappedReads = 0;
    if (me.options.verbose > 0)
    {
        transform(me.ctx.mapped, me.primaryMatches, isValid<void>, typename TTraits::TThreading());
        mappedReads = count(me.ctx.mapped, true, typename TTraits::TThreading());
        me.stats.mappedReads += mappedReads;
    }
    if (me.options.verbose > 1)
        std::cerr << "Mapped reads:\t\t\t" << mappedReads << std::endl;

    unsigned long pairedReads = 0;
    if (me.options.verbose > 0)
    {
        pairedReads = count(me.ctx.paired, true, typename TTraits::TThreading());
        me.stats.pairedReads += pairedReads;
    }
    if (me.options.verbose > 1)
    {
        std::cerr << "Pairing time:\t\t\t" << me.timer << std::endl;
        std::cerr << "Paired reads:\t\t\t" << pairedReads << std::endl;
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
                          me.suboptimalMatchesSet, me.primaryMatches, me.cigarSet,
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
    if (me.options.quick)
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
//    verifyMatches(me, readSeqs);
    rankMatches(me, readSeqs);
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

    if (!me.options.quick)
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
//    verifyMatches(me, readSeqs);
    rankMatches(me, readSeqs);
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
        std::cerr << "Pairing time:\t\t\t" << me.stats.selectPairs << " sec" << "\t\t" << me.stats.selectPairs / total << " %" << std::endl;
    std::cerr << "Alignment time:\t\t\t" << me.stats.alignMatches << " sec" << "\t\t" << me.stats.alignMatches / total << " %" << std::endl;
    std::cerr << "Output time:\t\t\t" << me.stats.writeMatches << " sec" << "\t\t" << me.stats.writeMatches / total << " %" << std::endl;

    printRuler(std::cerr);

    double totalReads = me.stats.loadedReads / 100.0;
    std::cerr << "Total reads:\t\t\t" << me.stats.loadedReads << std::endl;
    std::cerr << "Mapped reads:\t\t\t" << me.stats.mappedReads << "\t\t" << me.stats.mappedReads / totalReads << " %" << std::endl;
    if (IsSameType<typename TConfig::TSequencing, PairedEnd>::VALUE)
        std::cerr << "Paired reads:\t\t\t" << me.stats.pairedReads << "\t\t" << me.stats.pairedReads / totalReads << " %" << std::endl;
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

template <typename TThreading, typename TSequencing, typename TStrategy>
inline void spawnMapper(Options const & options,
                        TThreading const & /* tag */,
                        TSequencing const & /* tag */,
                        TStrategy const & /* tag */)
{
    typedef ReadMapperConfig<TThreading, TSequencing, TStrategy>    TConfig;

    Mapper<void, TConfig> mapper(options);
    runMapper(mapper);
}

#endif  // #ifndef APP_YARA_MAPPER_H_
