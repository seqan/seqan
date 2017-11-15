#ifndef RAZERS_RAZERS_PARALLEL_H_
#define RAZERS_RAZERS_PARALLEL_H_

// TODO(holtgrew): Ideally, we do not need any locks.

#include "parallel_misc.h"
#include "parallel_job_queue.h"
#include "razers_match_filter.h"

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

template <typename TSpec>
class Lock;

struct Omp_;
typedef Tag<Omp_> Omp;

template <>
class Lock<Omp>
{
public:
    omp_lock_t lock_;

    Lock() { omp_init_lock(&lock_); }

    ~Lock() { omp_destroy_lock(&lock_); }
};

template <typename TMatches>
struct SingleVerificationResult
{
    std::shared_ptr<TMatches> matches;
    unsigned hitGroupId;
    unsigned windowNo;

    SingleVerificationResult() :
        hitGroupId(0), windowNo(0)
    {}

    SingleVerificationResult(std::shared_ptr<TMatches> & matches_, unsigned hitGroupId_, unsigned windowNo_) :
        matches(matches_), hitGroupId(hitGroupId_), windowNo(windowNo_)
    {}
};

// Stores the results of the verification.
//
// Put into its own class so it can be locked independently of other class
// members.
template <typename TMatches>
class SingleVerificationResults
{
public:
    String<SingleVerificationResult<TMatches> > localMatches;
    Lock<Omp> * lock;

    SingleVerificationResults() :
        lock(new Lock<Omp>()) {}

    SingleVerificationResults(SingleVerificationResults const & other) :
        localMatches(other.localMatches), lock(new Lock<Omp>())
    {
        // Not thread-safe copying since this is only used at the beginning when resizing block local storages string.
    }

    SingleVerificationResults & operator=(SingleVerificationResults const & other)
    {
        if (this == &other)
            return *this;

        localMatches = other.localMatches;
        return *this;
    }

    ~SingleVerificationResults()
    {
        delete lock;
    }

};

template <
    typename TMatches,
    typename TFragmentStore,
    typename TFilterFinder,
    typename TFilterPattern,
    typename TShape /*TODO(holtgrew): Superflous.*/,
    typename TOptions,
    typename TCounts,
    typename TRazerSMode>
struct MapSingleReads {};

template <typename TSpec>
class ThreadLocalStorage;

// ThreadLocalStorage specialization for single-end read mapping in RazerS.
template <
    typename TMatches_,
    typename TFragmentStore,
    typename TFilterFinder_,
    typename TFilterPattern_,
    typename TShape /*TODO(holtgrew): Superflous.*/,
    typename TOptions,
    typename TCounts,
    typename TRazerSMode
    >
class ThreadLocalStorage<
    MapSingleReads<
        TMatches_,
        TFragmentStore,
        TFilterFinder_,
        TFilterPattern_,
        TShape,
        TOptions,
        TCounts,
        TRazerSMode> >
{
public:
    typedef TFilterPattern_ TFilterPattern;
    typedef TFilterFinder_ TFilterFinder;

    typedef TMatches_ TMatches;
#ifdef RAZERS_EXTERNAL_MATCHES
    typedef TMatches TLargeMatches;
#else // #ifdef RAZERS_EXTERNAL_MATCHES
    typedef typename Value<TMatches>::Type TMatchRecord;
    typedef String<TMatchRecord, MMap<ExternalConfigLarge<> > > TLargeMatches;
#endif // #ifdef RAZERS_EXTERNAL_MATCHES

    // The id of this thread.
    unsigned threadId;

    // Each thread needs its local options since the compactionThreshold is changed.
    // TODO(holtgrew): Change overall program structure so this is factorized out of the options struct.
    TOptions options;
    TOptions /*const*/ * globalOptions;

    // Each thread has its own SWIFT finder and pattern object.
    TFilterFinder filterFinder;
    TFilterPattern filterPattern;

    TCounts counts;  // TODO(holtgrew): Artifact?

    TLargeMatches matches;  // TODO(holtgrew): However, not used in global store since global reads-to-reference alignment requires everything to be in memory
    TFragmentStore /*const*/ * globalStore;

    TShape shape;

    String<String<SingleVerificationResult<TMatches> > > verificationResultBuckets;
    String<unsigned> missingInBucket;
    unsigned nextBucketToWriteBack;

    typedef MatchVerifier<TFragmentStore, TMatches, TOptions, TRazerSMode, TFilterPattern, TCounts> TMatchVerifier;
    TMatchVerifier verifier;

    // Mailbox for the verification results.
    SingleVerificationResults<TMatches> verificationResults;

    typedef MatchFilter<typename Spec<TOptions>::Type, typename TFragmentStore::TReadSeqStore, ThreadLocalStorage> TMatchFilter;
    std::shared_ptr<TMatchFilter> matchFilter;

    String<unsigned> splitters;

    ThreadLocalStorage() {}
};

template <typename TSpec>
inline void limitRead(ThreadLocalStorage<TSpec> & tls, unsigned readId, int newLimit)
{
    setMaxErrors(tls, readId, newLimit);
}

template <typename TSpec>
inline void disableRead(ThreadLocalStorage<TSpec> & tls, unsigned readId)
{
    setMaxErrors(tls, readId, -1);
}

template <typename TMatches, typename TFragmentStore, typename THitString, typename TOptions, typename TFilterPattern>
struct SingleVerification;

template <typename TMatches, typename TFragmentStore, typename THitString_, typename TOptions, typename TFilterPattern>
class Job<SingleVerification<TMatches, TFragmentStore, THitString_, TOptions, TFilterPattern> >
{
public:
    typedef SingleVerificationResults<TMatches> TVerificationResults;
    typedef THitString_ THitString;

    int threadId;
    TVerificationResults * verificationResults;
    TFragmentStore * globalStore;
    unsigned contigId;
    unsigned windowNo;
    std::shared_ptr<THitString> hitsPtr;
    unsigned hitGroupId;
    unsigned hitBegin;
    unsigned hitEnd;
    TOptions * options;
    TFilterPattern * filterPattern;

    Job() {}

    Job(int threadId_, TVerificationResults & verificationResults_, TFragmentStore & globalStore_, unsigned contigId_, unsigned windowNo_, std::shared_ptr<THitString> & hitsPtr_, unsigned hitGroupId_, unsigned hitBegin_, unsigned hitEnd_, TOptions & options_, TFilterPattern & filterPattern_) :
        threadId(threadId_), verificationResults(&verificationResults_), globalStore(&globalStore_), contigId(contigId_), windowNo(windowNo_), hitsPtr(hitsPtr_), hitGroupId(hitGroupId_), hitBegin(hitBegin_), hitEnd(hitEnd_), options(&options_), filterPattern(&filterPattern_)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TMatches>
inline void
appendToVerificationResults(SingleVerificationResults<TMatches> & verificationResults, SingleVerificationResult<TMatches> const & results)
{
    omp_set_lock(&verificationResults.lock->lock_);
//     if (length(verificationResults.localMatches) > 0u)
// SEQAN_OMP_PRAGMA(critical)
//         std::cerr << "BEFORE: " << &verificationResults << " length(*front(verificationResults.localMatches).matches) == " << length(*front(verificationResults.localMatches).matches) << ", " << (void*)front(verificationResults.localMatches).matches.get() << std::endl;
    appendValue(verificationResults.localMatches, results);
// SEQAN_OMP_PRAGMA(critical)
//     std::cerr << "AFTER:  " << &verificationResults << " length(*front(verificationResults.localMatches).matches) == " << length(*front(verificationResults.localMatches).matches) << ", " << (void*)front(verificationResults.localMatches).matches.get() << std::endl;
    omp_unset_lock(&verificationResults.lock->lock_);
}

// Uses the read ID to find the correct SWIFT pattern of the TLS, and
// the correct local ID within this SWIFT pattern to update the max
// errors.
//
// We do not disable the read right
template <typename TSpec, typename TReadNo, typename TMaxErrors>
inline void
setMaxErrors(ThreadLocalStorage<TSpec> & tls, TReadNo readNo, TMaxErrors maxErrors)
{
    int localReadNo = readNo - tls.splitters[tls.threadId];
    setMaxErrors(tls.filterPattern, localReadNo, maxErrors);
}

template <
    typename TMatches,
    typename TFragmentStore,
    typename TFilterFinder,
    typename TFilterPattern,
    typename TShape /*TODO(holtgrew): Superflous.*/,
    typename TOptions,
    typename TCounts,
    typename TRazerSMode,
    typename THitString>
void workVerification(ThreadLocalStorage<MapSingleReads<TMatches, TFragmentStore, TFilterFinder, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls,
                      Job<SingleVerification<TMatches, TFragmentStore, THitString, TOptions, TFilterPattern> > & job,
                      String<unsigned> const & splitters)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE
    double start = sysTime();

    typedef typename Iterator<THitString, Standard>::Type THitStringIterator;
    typedef typename TFragmentStore::TContigSeq TContigSeq;

    TContigSeq & contigSeq = job.globalStore->contigStore[job.contigId].seq;
#ifdef RAZERS_BANDED_MYERS
    int64_t contigLength = length(contigSeq);
#endif

    std::shared_ptr<TMatches> localMatches(new TMatches());
    resize(*localMatches, 1);
    clear(*localMatches);

    // Initialize verifier.
    tls.verifier.matches = localMatches.get();
    tls.verifier.options = job.options;
    tls.verifier.filterPattern = job.filterPattern;
    tls.verifier.cnts = 0;

    unsigned offset = splitters[job.threadId];
    for (THitStringIterator it = iter(*job.hitsPtr, job.hitBegin), itEnd = iter(*job.hitsPtr, job.hitEnd); it != itEnd; ++it)
    {
//        if (length(swiftInfix(value(it), job.globalStore->contigStore[job.contigId].seq)) < length(tls.globalStore->readSeqStore[value(it).ndlSeqNo]))
//            continue;  // Skip if hit length < read length.  TODO(holtgrew): David has to fix something in banded myers to make this work.

        unsigned absReadId = offset + value(it).ndlSeqNo;
        tls.verifier.m.readId = absReadId;

#ifdef RAZERS_BANDED_MYERS
        tls.verifier.patternState.leftClip = (it->hstkPos >= 0) ? 0 : -it->hstkPos;   // left clip if match begins left of the genome
        tls.verifier.rightClip = (it->hstkPos + it->bucketWidth <= contigLength) ? 0 : it->hstkPos + it->bucketWidth - contigLength;  // right clip if match end right of the genome
#endif
        matchVerify(tls.verifier, swiftInfix(value(it), contigSeq), absReadId, tls.globalStore->readSeqStore[absReadId], TRazerSMode());
    }

// SEQAN_OMP_PRAGMA(critical)
//     {
//         std::cerr << "thread " << omp_get_thread_num() << "; window " << job.windowNo << " hit group id " << job.hitGroupId << " thread id " << job.threadId << std::endl;
//         std::cerr << "localMatches.data_begin == " << (void*)(localMatches->data_begin) << ", localMatches.data_end == " << (void*)(localMatches->data_end) << std::endl;
//         std::cerr << "length(*localMatches) == " << length(*localMatches) << std::endl;
//     }

    appendToVerificationResults(*job.verificationResults, SingleVerificationResult<TMatches>(localMatches, job.hitGroupId, job.windowNo));
    tls.options.timeVerification += sysTime() - start;

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE
}

template <
    typename TMatches,
    typename TFragmentStore,
    typename TFilterFinder,
    typename TFilterPattern,
    typename TShape /*TODO(holtgrew): Superflous.*/,
    typename TOptions,
    typename TCounts,
    typename TRazerSMode>
void
writeBackToLocal(ThreadLocalStorage<MapSingleReads<TMatches, TFragmentStore, TFilterFinder, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, String<SingleVerificationResult<TMatches> > & verificationHits, bool dontCompact)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
    if (tls.threadId == 0u && tls.options._debugLevel >= 3)
        fprintf(stderr, "[writeback]");
    typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;

    TAlignedReadStoreSize oldSize = length(tls.matches);
    TAlignedReadStoreSize newSize = oldSize;

#ifdef RAZERS_DEFER_COMPACTION
    (void) dontCompact;  // unused

    // (1) Write back verification results into bucket.
    for (unsigned i = 0; i < length(verificationHits); ++i)
    {
        // verificationHits[i].windowNo;
// SEQAN_OMP_PRAGMA(critical)
//         std::cerr << "thread " << omp_get_thread_num() << " got (" << verificationHits[i].windowNo << ", " << verificationHits[i].hitGroupId << ")" << std::endl;
        tls.verificationResultBuckets[verificationHits[i].windowNo][verificationHits[i].hitGroupId] = verificationHits[i];
        tls.missingInBucket[verificationHits[i].windowNo] -= 1;
// SEQAN_OMP_PRAGMA(critical)
//         std::cerr << "thread " << omp_get_thread_num() << " {windowNo == " << verificationHits[i].windowNo << "--(" << tls.missingInBucket[verificationHits[i].windowNo] << ")}" << std::endl;
    }
// SEQAN_OMP_PRAGMA(critical)
//     std::cerr << "thread " << omp_get_thread_num() << " [wrote " << length(verificationHits) << " matches to buckets]" << std::flush;

    unsigned const DELTA = getMaxDeviationOfOrder(tls.filterPattern);
    // std::cerr << "(DELTA=" << DELTA << ")";
    //std::cerr << "[DELTA=" << DELTA << std::endl;
    size_t firstBeginPos = std::numeric_limits<size_t>::max();  // Leftmost sort position, required later for masking.
    size_t firstWindowBegin = std::numeric_limits<size_t>::max();  // Leftmost sort position, required later for masking.
    unsigned bucketsWrittenBack = 0;
    // (2) Write back the longest contiguous sequence of full buckets.
    for (; tls.nextBucketToWriteBack < length(tls.missingInBucket) && tls.missingInBucket[tls.nextBucketToWriteBack] == 0u; ++tls.nextBucketToWriteBack, ++bucketsWrittenBack)
    {
        // std::cerr << "(((WRITING BACK BUCKET " << tls.nextBucketToWriteBack << ")))" << std::endl;
        // (2 a) Compute new size, reserve memory, copy data.
        size_t originalSize = length(tls.matches);
        unsigned idx = tls.nextBucketToWriteBack;
// SEQAN_OMP_PRAGMA(critical)
//         std::cerr << "\n";
        for (unsigned i = 0; i < length(tls.verificationResultBuckets[idx]); ++i)
        {
// SEQAN_OMP_PRAGMA(critical)
//             {
//             std::cerr << "thread " << omp_get_thread_num() << " accessing (" << idx << ", " << i << ")" << std::endl;
//             // unsigned len = length(*tls.verificationResultBuckets[idx][i].matches);
//             std::cerr << "thread " << omp_get_thread_num() << " len=" << length(*tls.verificationResultBuckets[idx][i].matches) << "|" << std::endl;
//             std::cerr << "thread " << omp_get_thread_num() << " newSize=" << newSize << "|" << std::endl;
//             std::cerr << "thread " << omp_get_thread_num() << " matches.data_begin == " << (void*)(tls.verificationResultBuckets[idx][i].matches->data_begin) << ", matches.data_end == " << (void*)(tls.verificationResultBuckets[idx][i].matches->data_end) << std::endl;
//             }
            SEQAN_ASSERT_NEQ(tls.verificationResultBuckets[idx][i].matches.get(), static_cast<TMatches *>(0));
            newSize += length(*tls.verificationResultBuckets[idx][i].matches);
        }
        reserve(tls.matches, newSize);
        for (unsigned i = 0; i < length(tls.verificationResultBuckets[idx]); ++i)
        {
            if (!empty(*tls.verificationResultBuckets[idx][i].matches))
                // std::cerr << "BUCKET FROM WINDOW " << tls.nextBucketToWriteBack << "\t" << front(*tls.verificationResultBuckets[idx][i].matches).beginPos << "\t" << back(*tls.verificationResultBuckets[idx][i].matches).endPos << "\t" << length(*tls.verificationResultBuckets[idx][i].matches) << std::endl;
                append(tls.matches, *tls.verificationResultBuckets[idx][i].matches);
        }

        // (2 b) Get begin position.
        size_t beginPos = originalSize;
        // std::cerr << "originalSize = " << originalSize << std::endl;
        if (beginPos > 0u)
            beginPos -= 1;
        size_t dPos = 1;
        // Exponential search backwards.  After masking, reads are sorted by begin position.
        size_t windowBegin = tls.options.windowSize * idx;
        if (firstWindowBegin == std::numeric_limits<size_t>::max())
            firstWindowBegin = windowBegin;
        while (beginPos > 0u &&
               static_cast<size_t>(tls.matches[beginPos].beginPos) < windowBegin &&
               static_cast<size_t>(tls.matches[beginPos].beginPos + 10 * DELTA) > windowBegin)
        {
            if (beginPos > dPos)
                beginPos -= dPos;
            else
                beginPos = 0;
            dPos *= 2;
        }
        // // Binary search forwards.
        // typedef typename Iterator<TMatches>::Type TIterator;
        // typedef typename Value<TMatches>::Type TMatch;
        // TMatch m;
        // m.beginPos = windowBegin;
        // LessBeginPos<TMatch> cmp;
        // TIterator it = std::lower_bound(begin(tls.matches, Standard()) + beginPos, end(tls.matches, Standard()), m, cmp);
        // beginPos = it - begin(tls.matches, Standard());
        if (firstBeginPos == std::numeric_limits<size_t>::max())
            firstBeginPos = beginPos;

// SEQAN_OMP_PRAGMA(critical)
//         if (length(tls.matches) > 0u)
//             std::cerr << "((MASKING FROM " << tls.matches[beginPos].beginPos << " TO " << windowBegin + tls.options.windowSize << " ~ " << back(tls.matches).beginPos << "))" << std::endl
//                       << "  ->> " << length(tls.matches) - beginPos << " of " << length(tls.matches) << " elements , beginPos == " << beginPos << std::endl;
// SEQAN_OMP_PRAGMA(critical)
//         if (length(tls.matches) > 0u)
//             std::cerr << "((MASKING FROM " << tls.matches[beginPos].beginPos << " TO " << windowBegin + tls.options.windowSize << "))" << std::endl;
// SEQAN_OMP_PRAGMA(critical)
//         std::cerr << "thread " << omp_get_thread_num() << " (masking from " << beginPos << " to end=" << length(tls.matches) << ")" << std::endl;
        // Do not mask duplicates if not in edit distance mode and not using pigeonhole filter
        if (!IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE && tls.options.threshold != 0)
            continue;

        // (2 c) Mask duplicates from beginPos to the end position.
        maskDuplicates(tls.matches, begin(tls.matches, Standard()) + beginPos, end(tls.matches, Standard()), tls.options, TRazerSMode());
    }

    // std::cerr << "[wrote back " << bucketsWrittenBack << " buckets (" << tls.nextBucketToWriteBack << "/" << length(tls.missingInBucket) << ")]" << std::endl;
    if (bucketsWrittenBack > 0u)
    {
        // (3) Update match filter data structure for disabling reads.
        size_t nextWindowBegin = tls.options.windowSize * (tls.nextBucketToWriteBack);
        if (tls.nextBucketToWriteBack == length(tls.missingInBucket)) // is last?
            nextWindowBegin += DELTA;
        // std::cerr << "((REGISTERING/PROCESSING FROM " << firstWindowBegin << " TO " << nextWindowBegin << "))" << std::endl;
        typedef typename Iterator<TMatches>::Type TIterator;
        TIterator itBegin = begin(tls.matches, Standard()) + firstBeginPos;
        TIterator itEnd = end(tls.matches, Standard());
        TIterator it = itBegin;
        // SEQAN_ASSERT_LT(itBegin->beginPos, nextWindowBegin);
        for (; it != itEnd && static_cast<size_t>(it->beginPos + DELTA) <= nextWindowBegin; ++it)
        {
            if (it->orientation == '-')
                continue;                          // Skip masked reads.
            if (it->isRegistered)
                continue;
            it->isRegistered = true;
            registerRead(*tls.matchFilter, it->readId, it->score);
        }
        itEnd = it;
        it = itBegin;
        unsigned disabled = 0;
        for (; it != itEnd; ++it)
        {
            if (it->orientation == '-')
                continue;                          // Skip masked reads.
            disabled += processRead(*tls.matchFilter, it->readId);
        }
        if (tls.options._debugLevel >= 2 && disabled > 0)
            fprintf(stderr, " [%u reads disabled]", disabled);
    }
#else  // #ifdef RAZERS_DEFER_COMPACTION
    for (unsigned i = 0; i < length(verificationHits); ++i)
        newSize += length(*verificationHits[i].matches);

    reserve(tls.matches, newSize);

    // Write back all matches from verification to the block local store.
    for (unsigned i = 0; i < length(verificationHits); ++i)
        append(tls.matches, *verificationHits[i].matches);

    // Possibly compact matches.
    if (!dontCompact && length(tls.matches) > tls.options.compactThresh)
    {
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_COMPACT);
#endif  // #ifdef RAZERS_PROFILE
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typename Size<TAlignedReadStore>::Type oldSize = length(tls.matches);

        // if (tls.threadId == 0u && tls.options._debugLevel >= 3)
        //     fprintf(stderr, "[compact]");
        if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE || tls.options.threshold == 0)
            maskDuplicates(tls.matches, tls.options, TRazerSMode());  // overlapping parallelograms cause duplicates

        compactMatches(tls.matches, tls.counts, tls.options, TRazerSMode(), tls, COMPACT);

        if (length(tls.matches) * 4 > oldSize)       // the threshold should not be raised if too many matches were removed
        {
            while (tls.options.compactThresh < oldSize)
                tls.options.compactThresh *= tls.options.compactMult;
            // tls.options.compactThresh += (tls.options.compactThresh >> 1);  // if too many matches were removed
            if (tls.threadId == 0u && tls.options._debugLevel >= 3)
                fprintf(stderr, "[raising threshold to %u]", unsigned(tls.options.compactThresh));
        }
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_COMPACT);
#endif  // #ifdef RAZERS_PROFILE
    }
#endif  // #ifdef RAZERS_DEFER_COMPACTION

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
}

template <typename TMatches>
void clearLocalMatches(String<TMatches *> & localMatches)
{
    for (unsigned i = 0; i < length(localMatches); ++i)
        delete localMatches[i];
    clear(localMatches);
}

template <typename TMatches>
void clearLocalMatches(String<SingleVerificationResult<TMatches> > & localMatches)
{
    clear(localMatches);
}

// Find read matches in one genome sequence.
//
// The parallelization is simple.  We perform the filtering window-wise.
// Then, we generate verification jobs from the swift hits.  The SWIFT hits
// are distributed by a simple static load balancing of adjacent hits.  These
// are then processed by all "leading" threads, where leading threads are
// those with the largest number of processed windows.
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TThreadLocalStorages,
    typename TContigId,
    typename TCounts,
    typename TSpec,
    typename TShape,
    typename TAlignMode,
    typename TGapMode,
    typename TScoreMode,
    typename TMatchNPolicy,
    typename TFilterSpec>
void _mapSingleReadsParallelToContig(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TThreadLocalStorages & threadLocalStorages,
    String<unsigned> const & splitters,
    TContigId const & contigId,
    TCounts & /*cnts*/,
    char                                                      orientation,
    RazerSCoreOptions<TSpec> & options,
    TShape const & /*shape*/,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & /*mode*/,
    TFilterSpec)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE
    typedef FragmentStore<TFSSpec, TFSConfig>                       TFragmentStore;
    typedef typename TFragmentStore::TContigSeq                     TContigSeq;
    typedef typename TFragmentStore::TReadSeqStore                  TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type const               TRead;
    typedef StringSet<TRead>                                        TReadSet;
    typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >    TIndex;         // q-gram index
    //typedef typename Size<TReadSeqStore>::Type                      TSize;

    //typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
    typedef RazerSCoreOptions<TSpec> TOptions;

    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TMatches TMatches;

    // TODO(holtgrew): What about cnts, mode?

    typedef Finder<TContigSeq, TFilterSpec>                         TFilterFinder;
    typedef Pattern<TIndex, TFilterSpec>                            TFilterPattern;

    typedef RazerSCoreOptions<TSpec> TOptions;

    typedef typename WindowFindResult<TFilterFinder>::Type THitString;
    typedef Job<SingleVerification<TMatches, TFragmentStore, THitString, TOptions, TFilterPattern> > TVerificationJob;

    // Debug output...
    if (options._debugLevel >= 1)
    {
        std::cerr << std::endl << "Process genome seq #" << contigId;
        if (orientation == 'F')
            std::cerr << "[fwd]";
        else
            std::cerr << "[rev]";
    }

    // -----------------------------------------------------------------------
    // Reverse-complement contig if necessary.
    // -----------------------------------------------------------------------
    TContigSeq & contigSeq = store.contigStore[contigId].seq;
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE
    if (orientation == 'R')
        reverseComplement(contigSeq);
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE

    // -----------------------------------------------------------------------
    // Per-contig initialization of thread local storage objects.
    // -----------------------------------------------------------------------
    // TODO(holtgrew): Maybe put into its own function?
    for (unsigned i = 0; i < options.threadCount; ++i)
    {
        // Initialize verifier object.
        threadLocalStorages[i].verifier.onReverseComplement = (orientation == 'R');
        threadLocalStorages[i].verifier.genomeLength = length(contigSeq);
        threadLocalStorages[i].verifier.oneMatchPerBucket = false;
        threadLocalStorages[i].verifier.m.contigId = contigId;
    }

    // -----------------------------------------------------------------------
    // Perform filtration.
    // -----------------------------------------------------------------------
    TaskQueue<TVerificationJob, OmpLock> taskQueue;
    volatile unsigned leaderWindowsDone = 0;  // Number of windows done in leaders.
    volatile unsigned threadsFiltering = options.threadCount;

    // We will create the swift finder for thread 0 first.  This will trigger parallel repeat finding in the SWIFT
    // finder construction.  Then, we copy over the finder to all threads.
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_COPY_FINDER);
#endif  // #ifdef RAZERS_PROFILE

    threadLocalStorages[0].filterFinder = TFilterFinder(store.contigStore[contigId].seq, threadLocalStorages[0].options.repeatLength, 1);

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_COPY_FINDER);
#endif  // #ifdef RAZERS_PROFILE

    SEQAN_OMP_PRAGMA(parallel)
    {
        unsigned windowsDone = 0;

        // Initialization.
        TThreadLocalStorage & tls = threadLocalStorages[omp_get_thread_num()];
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_COPY_FINDER);
#endif  // #ifdef RAZERS_PROFILE
        // tls.filterFinder = TFilterFinder(store.contigStore[contigId].seq, tls.options.repeatLength, 1);
        if (omp_get_thread_num() != 0)
            tls.filterFinder = threadLocalStorages[0].filterFinder;

        // wait until everyone copied the filter of thread 0 before it is changed
        SEQAN_OMP_PRAGMA(barrier)

#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_COPY_FINDER);
        timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

        if (!windowFindBegin(tls.filterFinder, tls.filterPattern, tls.options.errorRate))
            std::cerr << "ERROR: windowFindBegin() failed in thread " << tls.threadId << std::endl;

#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

        // Pre-allocate buckets.
        tls.nextBucketToWriteBack = 0;
        resize(tls.verificationResultBuckets, unsigned(ceil(1.0 * (length(contigSeq) + 1000) / options.windowSize)));  // +1000 for the case where contig seq length is multiple of window size
// SEQAN_OMP_PRAGMA(critical)
//         std::cerr << "window count: " << length(tls.verificationResultBuckets) << std::endl;
        clear(tls.missingInBucket);
        resize(tls.missingInBucket, length(tls.verificationResultBuckets), std::numeric_limits<unsigned>::max());

        // For each filtration window...
        bool hasMore = !empty(host(tls.filterFinder));
        while (hasMore)
        {
#ifdef RAZERS_PROFILE
            timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE
            double filterStart = sysTime();
            // fprintf(stderr, "[filter]");
            hasMore = windowFindNext(tls.filterFinder, tls.filterPattern, tls.options.windowSize);
            // std::cerr << "FILTERING WINDOW " << windowsDone << std::endl << "\t" << tls.options.windowSize;

            windowsDone += 1;  // Local windows done count.
            atomicMax(leaderWindowsDone, windowsDone);

            std::shared_ptr<THitString> hitsPtr(new THitString()); //TODO (weese:) Could we reuse memory here?
            resize(*hitsPtr, 1);
            clear(*hitsPtr);
            using std::swap;
            swap(*hitsPtr, getWindowFindHits(tls.filterFinder));
            THitString & hits = *hitsPtr;
            tls.options.countFiltration += length(hits);
            // std::cerr << "  HITS: " << length(hits) << std::endl;
            // if (length(hits) > 0u)
            //     std::cerr << "  RANGE " << front(hits).hstkPos << "\t" << back(hits).bucketWidth << std::endl;

            // Enqueue verification jobs.
            if (length(hits) > 0u)
            {
                // Compute splitters, given a verification package size and a
                // bound on the package count.
                String<unsigned> splitters;
                unsigned packageCount = tls.options.maxVerificationPackageCount * omp_get_max_threads();
                computeSplittersBySlotSize(splitters, length(hits), tls.options.verificationPackageSize, packageCount);

                // Push verification jobs to the job queue.
                String<TVerificationJob> jobs;
                reserve(jobs, length(splitters) - 1);
                for (unsigned i = 1; i < length(splitters); ++i)
                {
// SEQAN_OMP_PRAGMA(critical)
//                     std::cerr << "\n";
                    appendValue(jobs, TVerificationJob(tls.threadId, tls.verificationResults, store, contigId, windowsDone - 1, hitsPtr, i - 1, splitters[i - 1], splitters[i], *tls.globalOptions, tls.filterPattern));
// SEQAN_OMP_PRAGMA(critical)
//                     std::cerr << "new job(" << tls.threadId << ", tls.verificationResults, store, " << contigId << ", " << windowsDone - 1 << ", hitsPtr, " << i - 1 << ", " << splitters[i - 1] << ", " << splitters[i] << ", *tls.globalOptions, tls.filterPattern)" << std::endl;
                }
                pushFront(taskQueue, jobs);

                // Preallocate space in bucket and initialize "to do" counter.
                clear(tls.verificationResultBuckets[windowsDone - 1]);
                resize(tls.verificationResultBuckets[windowsDone - 1], length(splitters) - 1);
                tls.missingInBucket[windowsDone - 1] = length(splitters) - 1;
                // for (unsigned i = 0; i < length(splitters) - 1; ++i)
                //     tls.missingInBucket[windowsDone - 1] -= splitters[i] == splitters[i + 1];
            }
            else
            {
                tls.missingInBucket[windowsDone - 1] = 0;
            }
            tls.options.timeFiltration += sysTime() - filterStart;
#ifdef RAZERS_PROFILE
            timelineEndTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

            // Perform verification as long as we are a leader and there are filtration jobs to perform.
            while (leaderWindowsDone == windowsDone)
            {
                TVerificationJob job;
                if (!popFront(job, taskQueue))
                    break;
                // fprintf(stderr, "[verify]");
                workVerification(tls, job, splitters);
            }

            // Write back verification results for this thread so far.
            //
            // First, swap out the current set of local stores from the verification results.
            omp_set_lock(&tls.verificationResults.lock->lock_);
            String<SingleVerificationResult<TMatches> > localMatches;
// SEQAN_OMP_PRAGMA(critical)
//             std::cerr << "thread " << omp_get_thread_num() << " SWAPPING " << &tls.verificationResults.localMatches << "\n";
            swap(localMatches, tls.verificationResults.localMatches);
            omp_unset_lock(&tls.verificationResults.lock->lock_);
            // Don't compact matches if in configured 'block fraction' of genome.
            size_t hstckLen = tls.filterFinder.endPos - tls.filterFinder.startPos;
            size_t hstckLeft = tls.filterFinder.endPos - tls.filterFinder.curPos;
            double fracTodo = 1.0  * hstckLeft / hstckLen;
            bool dontCompact = tls.options.noCompactFrac >= fracTodo;
            // Write back the contents of these stores to the thread-local store.
            writeBackToLocal(tls, localMatches, dontCompact);
// #ifndef RAZERS_DEFER_COMPACTION
            clearLocalMatches(localMatches);
// #endif  // #ifndef RAZERS_DEFER_COMPACTION
        }

        // Finalization
        windowFindEnd(tls.filterFinder, tls.filterPattern);

        SEQAN_OMP_PRAGMA(atomic)
        threadsFiltering -= 1;

        // Continue to try to help verify.
        while (threadsFiltering > 0u)
        {
            TVerificationJob job;
            if (popFront(job, taskQueue))
                workVerification(tls, job, splitters);
        }

        // After every thread is done with everything, write back once more.
        SEQAN_OMP_PRAGMA(barrier)
        writeBackToLocal(tls, tls.verificationResults.localMatches, true);
// #ifndef RAZERS_DEFER_COMPACTION
        clearLocalMatches(tls.verificationResults.localMatches);
// #endif  // #ifndef RAZERS_DEFER_COMPACTION
        // std::cerr << "AT END OF CONTIG #alignments " << length(tls.matches) << std::endl;
    }


    // NOTE:
    // We never undo the reverse-complementing!
    // It is not necessary as the contigs are freed by unlockAndFreeContig
    // They are loaded again in dumpMatches if necessary
    //
    // if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
    //     if (orientation == 'R')	reverseComplement(contigSeq);	// we have to restore original orientation
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE
}

// Global initialization of block local storages.
template <typename TThreadLocalStorages, typename TFragmentStore, typename TSplitters, typename TShape, typename TOptions>
void initializeThreadLocalStoragesSingle(TThreadLocalStorages & threadLocalStorages,
                                         TFragmentStore /*const*/ & store,
                                         TSplitters const & splitters,
                                         TShape /*const*/ & shape,
                                         TOptions /*const*/ & options)
{
    SEQAN_ASSERT_GT(length(splitters), 1u);
    int threadCount = length(splitters) - 1;

    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TFilterPattern TFilterPattern;
    typedef typename Host<TFilterPattern>::Type TIndex;
    typedef typename Position<typename TFragmentStore::TContigStore>::Type TPosition;

    resize(threadLocalStorages, threadCount);
    SEQAN_OMP_PRAGMA(parallel for schedule(static, 1))
    for (int i = 0; i < threadCount; ++i)
    {
        TThreadLocalStorage & tls = threadLocalStorages[i];

        tls.threadId = i;
        tls.globalStore = &store;
        tls.shape = shape;
        tls.options = options;  // TODO(holtgrew): Copy for stats and threshold, really good?
        tls.globalOptions = &options;
        tls.splitters = splitters;

#ifdef RAZERS_DEFER_COMPACTION
        typedef typename TThreadLocalStorage::TMatchFilter TMatchFilter;
        double READ_FRAC_WITH_HISTO = 0.01;
        tls.matchFilter.reset(new TMatchFilter(tls.splitters[tls.threadId + 1] - tls.splitters[tls.threadId], options.matchHistoStartThreshold, READ_FRAC_WITH_HISTO, tls, tls.splitters[tls.threadId], tls.globalStore->readSeqStore, tls.options));
        tls.options.compactThresh = std::numeric_limits<unsigned>::max();
#endif // #ifdef RAZERS_DEFER_COMPACTION

        // Clear pattern and set parameters.
        TFilterPattern & filterPattern = tls.filterPattern;
        clear(filterPattern);

        // Initialize the index.
        TIndex & index = host(tls.filterPattern);
        clear(index);
//        indexText(index).limitsValid = false;
//        assign(indexText(index).strings, infix(store.readSeqStore, splitters[i], splitters[i + 1]), Exact());

        clear(indexText(index));
        for (TPosition j = splitters[i]; j < splitters[i + 1]; ++j)
            appendValue(indexText(index), store.readSeqStore[j]);
        //unsigned x = length(indexText(index));
        //fprintf(stderr, "Index #%d has %u entries.\n", i, x);
        index.shape = shape;

#ifdef RAZERS_OPENADDRESSING
        index.alpha = options.loadFactor;
#endif
        cargo(index).abundanceCut = options.abundanceCut;
        cargo(index)._debugLevel = options._debugLevel;

        // Configure filter pattern
        // (if this is a pigeonhole filter, all sequences must be appended first)
        _applyFilterOptions(filterPattern, options);

        tls.filterPattern.params.printDots = (tls.threadId == 0) && (tls.options._debugLevel > 0);
    }
}

// Write back from thread local storages to global store.
template <typename TFragmentStore,
          typename TThreadLocalStorages>
void
writeBackToGlobalStore(
    TFragmentStore & target,
    TThreadLocalStorages /*const*/ & threadLocalStorages,
    bool isSingleEnd)      // begin/end already swapped for paired-end reads
{
    typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
    typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TMatches TMatches;
    typedef typename Iterator<TMatches, Standard>::Type TMatchesIterator;

    // Update the IDs and calculate new size so the prefix increment can be
    // used in the loops.
    TAlignedReadStoreSize oldSize = length(target.alignedReadStore);
    TAlignedReadStoreSize newSize = oldSize;

    for (unsigned i = 0; i < length(threadLocalStorages); ++i)
        newSize += length(threadLocalStorages[i].matches);

    // Resize first so copying happens at most once and not every for each
    // block in the worst case
    resize(target.alignedReadStore, newSize);
    resize(target.alignQualityStore, newSize);

    // Append single block stores.
    // TODO(holtgrew): Do in parallel!
    for (unsigned i = 0; i < length(threadLocalStorages); ++i)
    {
        TMatchesIterator it = begin(threadLocalStorages[i].matches, Standard());
        TMatchesIterator itEnd = end(threadLocalStorages[i].matches, Standard());
        for (; it != itEnd; ++it, ++oldSize)
        {
            using std::swap;
            if (isSingleEnd && it->orientation == 'R')
                swap(it->beginPos, it->endPos);
            target.alignedReadStore[oldSize] = TAlignedReadStoreElem(oldSize, it->readId, it->contigId, it->beginPos, it->endPos);
            if (!isSingleEnd)
                SEQAN_ASSERT_NEQ(it->pairMatchId, Value<TMatches>::Type::INVALID_ID);
            target.alignedReadStore[oldSize].pairMatchId = it->pairMatchId;
            target.alignQualityStore[oldSize] = TAlignedQualStoreElem(it->pairScore, it->score, -it->score);
        }
    }
}

template <typename TFragmentStore,
          typename TThreadLocalStorages>
void
writeBackToGlobalStore(
    TFragmentStore & target,
    TThreadLocalStorages /*const*/ & threadLocalStorages)
{
    writeBackToGlobalStore(target, threadLocalStorages, true);
}

// Performs splitting of reads, initialization of OpenMP and the calls
// mapSingleReadsParallelToContig for each contig.
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
int _mapSingleReadsParallel(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & mode,
    TFilterSpec)
{
    typedef FragmentStore<TFSSpec, TFSConfig>                       TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore                  TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type const               TRead;
    typedef StringSet<TRead>                                        TReadSet;
    typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >    TIndex;         // q-gram index
    //typedef typename Size<TReadSeqStore>::Type                      TSize;

    typedef Pattern<TIndex, TFilterSpec>                            TFilterPattern;
    //typedef Pattern<TRead, MyersUkkonen>                            TMyersPattern;  // verifier
    typedef RazerSCoreOptions<TSpec> TOptions;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
    typedef typename TFragmentStore::TContigSeq                     TContigSeq;
    typedef Finder<TContigSeq, TFilterSpec>                         TFilterFinder;

    typedef typename TFragmentStore::TContigSeq TContigSeq;
    typedef typename Position<TContigSeq>::Type TContigPos;
    typedef MatchRecord<TContigPos> TMatchRecord;
    typedef String<TMatchRecord> TMatches;

    // -----------------------------------------------------------------------
    // Initialize global information.
    // -----------------------------------------------------------------------

    // Save OpenMP maximal thread count so we can restore it below, then set
    // from options.
    if (options._debugLevel >= 1)
        std::cerr << std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_INIT);
    double beginInit = sysTime();
#endif  // #ifdef RAZERS_PROFILE

    // Clear/initialize global stats.
    options.countFiltration = 0;
    options.countVerification = 0;
    options.timeMapReads = 0;
    options.timeDumpResults = 0;

    // -----------------------------------------------------------------------
    // Initialize thread local storages.
    // -----------------------------------------------------------------------
    SEQAN_PROTIMESTART(initTime);
    String<unsigned> splitters;
    computeSplittersBySlotCount(splitters, length(store.readNameStore), options.threadCount);
    typedef ThreadLocalStorage<MapSingleReads<TMatches, TFragmentStore, TFilterFinder, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > TThreadLocalStorage;
    String<TThreadLocalStorage> threadLocalStorages;
    initializeThreadLocalStoragesSingle(threadLocalStorages, store, splitters, shape, options);

#ifdef RAZERS_PROFILE
    double endInit = sysTime();
    std::cerr << "TIME initialization: " << (endInit - beginInit) << " s";
    timelineEndTask(TASK_INIT);
#endif  // #ifdef RAZERS_PROFILE
    double timeInitialization = SEQAN_PROTIMEDIFF(initTime);
    if (options._debugLevel >= 1)
        std::cerr << std::endl << "Initialization took              \t" << timeInitialization << " seconds" << std::endl;

    // -----------------------------------------------------------------------
    // Perform parallel mapping.
    // -----------------------------------------------------------------------

    // Save compaction threshold and set global threshold to infinity, so matchVerify does not compact!
    int oldThreshold = options.compactThresh;
    options.compactThresh = std::numeric_limits<unsigned>::max();

    // For each contig: Map reads in parallel.
    SEQAN_PROTIMESTART(findTime);
    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId)
    {
        lockContig(store, contigId);
        if (options.forward)
            _mapSingleReadsParallelToContig(store, threadLocalStorages, splitters, contigId, cnts, 'F', options, shape, mode, TFilterSpec());
        if (options.reverse)
            _mapSingleReadsParallelToContig(store, threadLocalStorages, splitters, contigId, cnts, 'R', options, shape, mode, TFilterSpec());
        unlockAndFreeContig(store, contigId);
    }
    double endMapping = sysTime();
#ifdef RAZERS_PROFILE
    std::cerr << std::endl << "TIME mapping: " << (endMapping - endInit) << " s" << std::endl;
#endif  // #ifdef RAZERS_PROFILE

#ifdef RAZERS_EXTERNAL_MATCHES
    // Compute whether to use slow, sequential sorting or parallel in-memory sorting.
    uint64_t totalMatchCount = 0;
    uint64_t maxMatchCount = 0;
    for (unsigned i = 0; i < length(threadLocalStorages); ++i)
    {
        totalMatchCount += length(threadLocalStorages[i].matches);
        maxMatchCount = _max(maxMatchCount, length(threadLocalStorages[i].matches));
    }
    bool useExternalSort = false;
    bool useSequentialCompaction = false;
    if (options.availableMatchesMemorySize == -1)
    {
        useExternalSort = true;
    }
    else if (options.availableMatchesMemorySize != 0)
    {
        typedef typename Value<TMatches>::Type TMatch;
        int64_t totalMemoryRequired = sizeof(TMatch) * totalMatchCount;
        int64_t maxMemoryRequired = sizeof(TMatch) * maxMatchCount;
        if (options.availableMatchesMemorySize < totalMemoryRequired)
        {
            if (options.availableMatchesMemorySize < maxMemoryRequired)
            {
                useExternalSort = true;
            }
            else
            {
                useSequentialCompaction = true;
            }
        }
    }

    // std::cerr << "useExternalSort == " << useExternalSort << "\n"
    //           << "useSequentialCompaction == " << useSequentialCompaction << "\n";

    // Switch between using parallel compaction, sequential compaction, and
    // sequential compaction with external sorting.  The actual switch for the
    // sorting is in function compactMatches().
    if (useSequentialCompaction || useExternalSort)
    {
        for (unsigned i = 0; i < length(threadLocalStorages); ++i)
        {
            // remove duplicates when mapping with gaps or using pigeonhole filter
            if (IsSameType<TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
                maskDuplicates(threadLocalStorages[omp_get_thread_num()].matches, options, mode);
            Nothing nothing;
            CompactMatchesMode compactMode = useSequentialCompaction ? COMPACT_FINAL : COMPACT_FINAL_EXTERNAL;
            // std::cerr << "BEFORE FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
            compactMatches(threadLocalStorages[omp_get_thread_num()].matches, cnts, options, mode, nothing, compactMode);
            // std::cerr << "AFTER FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
        }
    }
    else
    {
#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
    SEQAN_OMP_PRAGMA(parallel)
    {
// TODO(holtgrew): We would really like to stop the additional masking step, the incremental masking SHOULD have taken care of this. Thus, the following should be ifndef and not ifdef.
#ifndef RAZERS_DEFER_COMPACTION
        // remove duplicates when mapping with gaps or using pigeonhole filter
        if (IsSameType<TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
            maskDuplicates(threadLocalStorages[omp_get_thread_num()].matches, options, mode);
#endif  // #ifndef RAZERS_DEFER_COMPACTION
        Nothing nothing;
        // std::cerr << "BEFORE FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
// SEQAN_OMP_PRAGMA(critical)
//             std::cerr << "BEFORE FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
        compactMatches(threadLocalStorages[omp_get_thread_num()].matches, cnts, options, mode, nothing, COMPACT_FINAL);
        // std::cerr << "AFTER FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
// SEQAN_OMP_PRAGMA(critical)
//             std::cerr << "AFTER FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
    }
    SEQAN_OMP_PRAGMA(barrier)
#ifdef RAZERS_EXTERNAL_MATCHES
}

#endif // #ifdef RAZERS_EXTERNAL_MATCHES

    // Write back local stores to global stores.
    writeBackToGlobalStore(store, threadLocalStorages);
    // std::cerr << "length(threadLocalStorages[0].matches) == " << length(threadLocalStorages[0].matches) << std::endl;
    // std::cerr << "length(store.alignedReadStore) == " << length(store.alignedReadStore) << std::endl;
    double endWriteback = sysTime();
    options.timeFsCopy = endWriteback - endMapping;

    // TODO(holtgrew): Sum up cnts?!

    // -----------------------------------------------------------------------
    // Collect global statistics, cleanup.
    // -----------------------------------------------------------------------

    // Restore old compaction threshold.
    options.compactThresh = oldThreshold;

    // Add up thread-local filtration and verification counters and print totals.
    for (unsigned i = 0; i < length(threadLocalStorages); ++i)
    {
        options.countFiltration += threadLocalStorages[i].options.countFiltration;
        options.countVerification += threadLocalStorages[i].options.countVerification;
    }

    if (options._debugLevel >= 1)
    {
        for (unsigned i = 0; i < length(threadLocalStorages); ++i)
        {
            std::cerr << "Thread #" << i << std::endl;
            std::cerr << "  Masking duplicates took        \t" << threadLocalStorages[i].options.timeMaskDuplicates << " seconds" << std::endl;
            std::cerr << "  Compacting matches took        \t" << threadLocalStorages[i].options.timeCompactMatches << " seconds" << std::endl;
            std::cerr << "  Time for filtration            \t" << threadLocalStorages[i].options.timeFiltration << " seconds" << std::endl;
            std::cerr << "  Time for verifications         \t" << threadLocalStorages[i].options.timeVerification << " seconds" << std::endl;
        }
        std::cerr << "Time for copying back            \t" << options.timeFsCopy << " seconds" << std::endl;
    }

    options.timeMapReads = SEQAN_PROTIMEDIFF(findTime);
    if (options._debugLevel >= 1)
        std::cerr << std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << std::endl;
    if (options._debugLevel >= 1)
    {
        std::cerr << std::endl;
        std::cerr << "___FILTRATION_STATS____" << std::endl;
        std::cerr << "Filtration counter:      " << options.countFiltration << std::endl;
        std::cerr << "Successful verifications: " << options.countVerification << std::endl;
    }

    // Restore global state.
    omp_set_num_threads(oldMaxThreads);

    return 0;
}

// Performs splitting of reads, initialization of OpenMP and the calls
// mapSingleReadsParallelToContig for each contig.
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
int _mapSingleReadsParallel(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & mode)
{
    if (options.threshold > 0)
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
        return _mapSingleReadsParallel(store, cnts, options, shape, mode, Swift<TSwiftSpec>());
    }
    else
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, void, Hamming_>::Type TPigeonholeSpec;
        return _mapSingleReadsParallel(store, cnts, options, Shape<Dna, OneGappedShape>(), mode, Pigeonhole<TPigeonholeSpec>());
    }
}

}  // namespace seqan

#endif  // #ifndef RAZERS_RAZERS_PARALLEL_H_
