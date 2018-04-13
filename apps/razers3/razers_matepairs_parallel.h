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
// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.

// TODO(holtgrew): aligned read id can probably go away...

#ifndef SEQAN_HEADER_RAZERS_MATEPAIRS_PARALLEL_H
#define SEQAN_HEADER_RAZERS_MATEPAIRS_PARALLEL_H

#include <numeric>

#include <seqan/misc/dequeue.h>

#include "razers_parallel.h"
#include "razers_matepairs.h"
#include "razers_paired_match_filter.h"

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// Stores the results of the verification.
//
// Put into its own class so it can be locked independently of other class
// members.
template <typename TMatches>
class PairedVerificationResults
{
public:
    String<TMatches *> localMatches;
    Lock<Omp> * lock;

    PairedVerificationResults() :
        lock(new Lock<Omp>()) {}

    PairedVerificationResults(PairedVerificationResults const & other) :
        localMatches(other.localMatches), lock(new Lock<Omp>())
    {
        // Not thread-safe copying since this is only used at the beginning when resizing block local storages string.
    }

    PairedVerificationResults & operator=(PairedVerificationResults const & other)
    {
        if (this == &other)
            return *this;

        localMatches = other.localMatches;
        return *this;
    }

    ~PairedVerificationResults()
    {
        delete lock;
    }

};

template <typename TMatches, typename TFragmentStore, typename TFilterFinderL, typename TFilterFinderR, typename TFilterPattern, typename TShape /*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode>
struct MapPairedReads {};

// ThreadLocalStorage specialization for single-end read mapping in RazerS.
template <typename TMatches_, typename TFragmentStore, typename TFilterFinderL_, typename TFilterFinderR_, typename TFilterPattern_, typename TShape_ /*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode>
class ThreadLocalStorage<MapPairedReads<TMatches_, TFragmentStore, TFilterFinderL_, TFilterFinderR_, TFilterPattern_, TShape_, TOptions, TCounts, TRazerSMode> >
{
public:
    typedef typename TFragmentStore::TReadSeqStore                  TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type                     TRead;
    typedef StringSet<TRead>                                        TReadSet;
    typedef TShape_ TShape;
    typedef TMatches_ TMatches;

    typedef TFilterPattern_ TFilterPattern;
    typedef TFilterFinderL_ TFilterFinderL;
    typedef TFilterFinderR_ TFilterFinderR;

    typedef typename TFragmentStore::TContigSeq             TGenome;
    typedef typename Infix<TGenome>::Type TGenomeInfix;

    typedef typename Value<TMatches>::Type TMatch;

    // The id of this thread.
    unsigned threadId;

    // Each thread needs its local options since the compactionThreshold is changed.
    // TODO(holtgrew): Change overall program structure so this is factorized out of the options struct.
    TOptions options;
    TOptions /*const*/ * globalOptions;

    // Split the read seq store from fragment store into two parts.
    TReadSet readSetL, readSetR;

    // Each thread has its own SWIFT finder and pattern object.
    TFilterFinderL filterFinderL;
    TFilterFinderR filterFinderR;
    TFilterPattern filterPatternL, filterPatternR;

    TCounts counts;  // TODO(holtgrew): Artifact?

#ifdef RAZERS_EXTERNAL_MATCHES
    typedef TMatches TLargeMatches;
#else // #ifdef RAZERS_EXTERNAL_MATCHES
    typedef typename Value<TMatches>::Type TMatchRecord;
    typedef String<TMatchRecord, MMap<ExternalConfigLarge<> > > TLargeMatches;
#endif // #ifdef RAZERS_EXTERNAL_MATCHES

    TLargeMatches matches;
    TFragmentStore /*const*/ * globalStore;

    TShape shape;

    typedef MatchVerifier<TFragmentStore, TMatches, TOptions, TRazerSMode, TFilterPattern, TCounts> TMatchVerifier;
    TMatchVerifier verifierL, verifierR;

    // Mailbox for the verification results.
    PairedVerificationResults<TMatches> verificationResults;

    typedef PairedMatchFilter<typename Spec<TOptions>::Type, typename TFragmentStore::TReadSeqStore, ThreadLocalStorage> TMatchFilter;
    std::shared_ptr<TMatchFilter> matchFilter;

    String<unsigned> splitters;

    unsigned completeWindows;

    TGenomeInfix genomeInf;

    // The last potential match fifo from the verification can be thread
    // local.
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedRead;
    typedef typename Value<TAlignQualityStore>::Type        TAlignQuality;
    typedef Pair<int64_t, TMatch>   TDequeueValue;
    typedef Dequeue<TDequeueValue>                          TDequeue;
    typedef typename TDequeue::TIter                        TDequeueIterator;
    TDequeue fifo;                      // stores left-mate potential matches
    String<int64_t> fifoLastPotMatchNo; // last number of a left-mate potential
    int64_t fifoLastNo;                 // last number over all left-mate pot. matches in the queue
    int64_t fifoFirstNo;                // first number over all left-mate pot. match in the queue

    ThreadLocalStorage() :
        fifoLastNo(0), fifoFirstNo(0) {}
};

template <typename TThreadLocalStorage>
class FilterPatternLSetMaxErrorsWrapper
{
public:
    TThreadLocalStorage & tls;
    FilterPatternLSetMaxErrorsWrapper(TThreadLocalStorage & tls_) :
        tls(tls_) {}
};

template <typename TThreadLocalStorage>
class FilterPatternRSetMaxErrorsWrapper
{
public:
    TThreadLocalStorage & tls;
    FilterPatternRSetMaxErrorsWrapper(TThreadLocalStorage & tls_) :
        tls(tls_) {}
};

template <typename TMatches, typename TFragmentStore, typename THitString, typename TOptions, typename TFilterPattern>
struct PairedVerification;

template <typename THitString>
void
buildHitSplittersAndPartitionHits(String<size_t> & splitters, THitString & hitString, unsigned packageCount, unsigned matePairCount)
{
    typedef typename Iterator<THitString>::Type THitStringIterator;

    (void)matePairCount;

    // TODO(holtgrew): Optimize with packageCount power of 2 and/or use libdiv?
    // TODO(holtgrew): Or maybe logarithmic search faster than modulo?

    unsigned const a = 1664525;
    unsigned const c = 1013904223;

    // Partition hitString into buckets by ndlSeqNo % packageCount.
    //
    // First, build counters.
    clear(splitters);
    resize(splitters, packageCount + 1, 0);
    splitters[0] = 0;
    for (THitStringIterator it = begin(hitString, Standard()), itEnd = end(hitString, Standard()); it != itEnd; ++it)
    {
        SEQAN_ASSERT_LEQ(it->ndlSeqNo, matePairCount);
        unsigned idx = (a * it->ndlSeqNo + c) % packageCount + 1;
        splitters[idx] += 1;
        SEQAN_ASSERT_LEQ(splitters[it->ndlSeqNo % packageCount + 1], length(hitString));
    }
    std::partial_sum(begin(splitters, Standard()), end(splitters, Standard()), begin(splitters, Standard()));
    // Second, copy into temporary buffer.
    String<size_t> offsets(splitters);
    THitString buffer;
    resize(buffer, length(hitString));
    for (THitStringIterator it = begin(hitString, Standard()), itEnd = end(hitString, Standard()); it != itEnd; ++it)
    {
        unsigned idx = (a * it->ndlSeqNo + c) % packageCount;
        buffer[offsets[idx]] = *it;
        offsets[idx] += 1;
        if (idx < length(offsets) - 1)
            SEQAN_ASSERT_LEQ(offsets[idx], offsets[idx + 1]);
        if (idx > 0u)
            SEQAN_ASSERT_LEQ(offsets[idx - 1], offsets[idx]);
    }

    // Finally, write out results.
    using std::swap;
    swap(hitString, buffer);

    // std::cout << "SPLITTERS: ";
    // for (unsigned i = 0; i < length(splitters); ++i) {
    //     std::cout << splitters[i] << " ";
    // }
    // std::cout << " last should be " << length(hitString) << std::endl;
}

template <typename TMatches, typename TFragmentStore, typename THitString_, typename TOptions, typename TFilterPattern>
class Job<PairedVerification<TMatches, TFragmentStore, THitString_, TOptions, TFilterPattern> >
{
public:
    typedef PairedVerificationResults<TMatches> TVerificationResults;
    typedef THitString_ THitString;
    typedef std::shared_ptr<THitString> THitStringPtr;

    int threadId;
    TVerificationResults * verificationResults;
    TFragmentStore * globalStore;
    unsigned contigId;
    char orientation;
    // Hit strings of previous window for the left finder.
    THitStringPtr prevHitsPtrL;
    size_t prevHitsLBegin, prevHitsLEnd;
    // Current hit string of current finder.
    THitStringPtr hitsPtrL;
    size_t hitsLBegin, hitsLEnd;
    THitStringPtr hitsPtrR;
    size_t hitsRBegin, hitsREnd;
    int64_t rightWindowBegin;
    TOptions * options;
    TFilterPattern * filterPatternL;
    TFilterPattern * filterPatternR;

    Job() {}

    Job(int threadId_, TVerificationResults & verificationResults_, TFragmentStore & globalStore_, unsigned contigId_, char orientation_, THitStringPtr & prevHitsPtrL_, size_t prevHitsLBegin_, size_t prevHitsLEnd_, THitStringPtr & hitsPtrL_, size_t hitsLBegin_, size_t hitsLEnd_, THitStringPtr & hitsPtrR_, size_t hitsRBegin_, size_t hitsREnd_, int64_t rightWindowBegin_, TOptions & options_, TFilterPattern & filterPatternL_, TFilterPattern & filterPatternR_) :
        threadId(threadId_), verificationResults(&verificationResults_), globalStore(&globalStore_), contigId(contigId_), orientation(orientation_), prevHitsPtrL(prevHitsPtrL_), prevHitsLBegin(prevHitsLBegin_), prevHitsLEnd(prevHitsLEnd_), hitsPtrL(hitsPtrL_), hitsLBegin(hitsLBegin_), hitsLEnd(hitsLEnd_), hitsPtrR(hitsPtrR_), hitsRBegin(hitsRBegin_), hitsREnd(hitsREnd_), rightWindowBegin(rightWindowBegin_), options(&options_), filterPatternL(&filterPatternL_), filterPatternR(&filterPatternR_)
    {
        SEQAN_ASSERT_LEQ(hitsLBegin, length(*hitsPtrL));
        SEQAN_ASSERT_LEQ(hitsLEnd, length(*hitsPtrL));
        SEQAN_ASSERT_LEQ(hitsRBegin, length(*hitsPtrR));
        SEQAN_ASSERT_LEQ(hitsREnd, length(*hitsPtrR));
        if (prevHitsPtrL.get() != 0)
        {
            SEQAN_ASSERT_LEQ(prevHitsLBegin, length(*prevHitsPtrL));
            SEQAN_ASSERT_LEQ(prevHitsLEnd, length(*prevHitsPtrL));
        }
    }

};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================


// Allow disabling reads in compactPairMatches() for left-mate read set.
//
// We do not disable the read right
template <typename TThreadLocalStorage, typename TReadNo, typename TMaxErrors>
void
setMaxErrors(FilterPatternLSetMaxErrorsWrapper<TThreadLocalStorage> & wrapper,
             TReadNo pairNo,
             TMaxErrors maxErrors)
{
    // std::cerr << std::endl << "SET MAX ERRORS LEFT" << std::endl;
    SEQAN_ASSERT_LT(pairNo, wrapper.tls.splitters[wrapper.tls.threadId + 1]);
    SEQAN_ASSERT_GEQ(pairNo, wrapper.tls.splitters[wrapper.tls.threadId]);
    setMaxErrors(wrapper.tls.filterPatternL, pairNo - wrapper.tls.splitters[wrapper.tls.threadId], maxErrors);
}

// Allow disabling reads in compactPairMatches() for left-mate read set.
//
// We do not disable the read right
template <typename TThreadLocalStorage, typename TReadNo, typename TMaxErrors>
void
setMaxErrors(FilterPatternRSetMaxErrorsWrapper<TThreadLocalStorage> & wrapper,
             TReadNo pairNo,
             TMaxErrors maxErrors)
{
    // std::cerr << std::endl << "SET MAX ERRORS RIGHT" << std::endl;
    SEQAN_ASSERT_LT(pairNo, wrapper.tls.splitters[wrapper.tls.threadId + 1]);
    SEQAN_ASSERT_GEQ(pairNo, wrapper.tls.splitters[wrapper.tls.threadId]);
    setMaxErrors(wrapper.tls.filterPatternR, pairNo - wrapper.tls.splitters[wrapper.tls.threadId], maxErrors);
}

template <typename TFragmentStore>
inline
void
appendToVerificationResults(PairedVerificationResults<TFragmentStore> & verificationResults, TFragmentStore * storePtr)
{
    omp_set_lock(&verificationResults.lock->lock_);
    appendValue(verificationResults.localMatches, storePtr);
    omp_unset_lock(&verificationResults.lock->lock_);
}

template <typename TThreadLocalStorages, typename TFragmentStore, typename TSplitters, typename TShape, typename TOptions>
void initializeThreadLocalStoragesPaired(TThreadLocalStorages & threadLocalStorages,
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
    typedef typename TThreadLocalStorage::TReadSet TReadSet;

    resize(threadLocalStorages, threadCount);
    SEQAN_OMP_PRAGMA(parallel for schedule(static, 1))
    for (int i = 0; i < threadCount; ++i)
    {
        TThreadLocalStorage & tls = threadLocalStorages[i];

        // Initialize properties that are simple to set.
        tls.threadId = i;
        tls.globalStore = &store;
        tls.shape = shape;
        tls.options = options;  // TODO(holtgrew): Copy for stats and threshold, really good?
        tls.globalOptions = &options;
        tls.splitters = splitters;

        // Split mate-pairs over two indices.
        unsigned pairCount = splitters[i + 1] - splitters[i];
        TReadSet & readSetL = tls.readSetL;
        TReadSet & readSetR = tls.readSetR;
        resize(readSetL, pairCount, Exact());
        resize(readSetR, pairCount, Exact());

#ifdef RAZERS_DEFER_COMPACTION
        typedef typename TThreadLocalStorage::TMatchFilter TMatchFilter;
        double READ_FRAC_WITH_HISTO = 0.01;
        tls.matchFilter.reset(new TMatchFilter(tls.splitters[i + 1] - tls.splitters[i], options.matchHistoStartThreshold, READ_FRAC_WITH_HISTO, tls, tls.splitters[i], tls.globalStore->readSeqStore, tls.options));
        tls.options.compactThresh = std::numeric_limits<unsigned>::max();
#endif // #ifdef RAZERS_DEFER_COMPACTION

        unsigned offset = splitters[i];
        for (unsigned j = 0; j < pairCount; ++j)
        {
            assign(readSetL[j], store.readSeqStore[store.matePairStore[offset + j].readId[0]]);
            assign(readSetR[j], store.readSeqStore[store.matePairStore[offset + j].readId[1]]);
        }
        reverseComplement(readSetR);

        // Clear patterns and set parameters.
        TFilterPattern & filterPatternL = tls.filterPatternL;
        TFilterPattern & filterPatternR = tls.filterPatternR;
        clear(filterPatternL);
        clear(filterPatternR);

        // Initialize the indices.
        // TODO(holtgrew): Necessary to split into readSetL and readSetR if we assign separately anyway?
        TIndex & indexL = host(tls.filterPatternL);
        clear(indexL);
        clear(indexText(indexL));
        reserve(indexText(indexL), length(readSetL), Exact());
        for (TPosition j = 0, jEnd = length(readSetL); j < jEnd; ++j)
            appendValue(indexText(indexL), readSetL[j]);
        indexL.shape = shape;
#ifdef RAZERS_OPENADDRESSING
        indexL.alpha = options.loadFactor;
#endif
        cargo(indexL).abundanceCut = options.abundanceCut;
        cargo(indexL)._debugLevel = options._debugLevel;
        indexRequire(indexL, QGramSADir());

        TIndex & indexR = host(tls.filterPatternR);
        clear(indexR);
        clear(indexText(indexR));
        reserve(indexText(indexR), length(readSetR), Exact());
        for (TPosition j = 0, jEnd = length(readSetR); j < jEnd; ++j)
            appendValue(indexText(indexR), readSetR[j]);
        indexR.shape = shape;
#ifdef RAZERS_OPENADDRESSING
        indexR.alpha = options.loadFactor;
#endif
        cargo(indexR).abundanceCut = options.abundanceCut;
        cargo(indexR)._debugLevel = options._debugLevel;
        indexRequire(indexR, QGramSADir());

        // Configure filter pattern
        // (if this is a pigeonhole filter, all sequences must be appended first)
        _applyFilterOptions(filterPatternL, options);
        _applyFilterOptions(filterPatternR, options);
        filterPatternL.params.printDots = false;
        filterPatternR.params.printDots = (tls.threadId == 0) && (tls.options._debugLevel > 0);
    }
}

template <typename TMatches, typename TFragmentStore, typename TFilterFinderL, typename TFilterFinderR, typename TFilterPattern, typename TShape /*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename THitString>
void workVerification(ThreadLocalStorage<MapPairedReads<TMatches, TFragmentStore, TFilterFinderL, TFilterFinderR, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls,
                      Job<PairedVerification<TMatches, TFragmentStore, THitString, TOptions, TFilterPattern> > & job,
                      String<unsigned> const & splitters)
{
    typedef ThreadLocalStorage<MapPairedReads<TMatches, TFragmentStore, TFilterFinderL, TFilterFinderR, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > TThreadLocalStorage;

    typedef typename TThreadLocalStorage::TReadSet TReadSet;

    typedef typename TFragmentStore::TContigSeq             TGenome;

    typedef typename TFragmentStore::TMatePairStore         TMatePairStore;
    // typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
    // typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
    typedef typename Value<TMatePairStore>::Type            TMatePair;
    // typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
    // typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
    typedef typename TThreadLocalStorage::TReadSet TReadSet;
    //typedef Index<TReadSet, IndexQGram<TShape>  >   TReadIndex;
    // typedef typename Id<TAlignedRead>::Type					TId;

    typedef typename TFragmentStore::TContigSeq             TGenome;
    typedef typename Size<TGenome>::Type                    TSize;
    typedef typename Position<TGenome>::Type                TGPos;
    typedef typename MakeSigned_<TGPos>::Type               TSignedGPos;
    //typedef typename Infix<TGenome>::Type                   TGenomeInf;

    typedef typename Iterator<THitString, Standard>::Type THitStringIter;

    typedef typename Value<TMatches>::Type TMatch;

    typedef MatchVerifier<
        TFragmentStore,
        TMatches,
        TOptions,
        TRazerSMode,
        TFilterPattern,
        TCounts>                                            TVerifier;

    // MATE-PAIR FILTRATION
    typedef Pair<int64_t, TMatch>   TDequeueValue;
    typedef Dequeue<TDequeueValue>                          TDequeue;
    typedef typename TDequeue::TIter                        TDequeueIterator;

    double startVerify = sysTime();

    // buffer variable to extend window in previous left hits string by.
    // TODO(holtgrew): DELTA has to be set to a better value, probably.
    unsigned const DELTA = getMaxDeviationOfOrder(tls.filterPatternL);

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE

    // Allocate fragmentstore for job's results.
    TMatches * localMatches = new TMatches();
    resize(*localMatches, 1);
    clear(*localMatches);

    // Thread-wide offset for reads.
    unsigned threadIdOffset = splitters[job.threadId];

    // Initialize verifiers.
    tls.verifierL.matches = localMatches;
    tls.verifierL.options = job.options;
    tls.verifierL.filterPattern = job.filterPatternL;
    tls.verifierL.cnts = 0;

    tls.verifierR.matches = localMatches;
    tls.verifierR.options = job.options;
    tls.verifierR.filterPattern = job.filterPatternR;
    tls.verifierR.cnts = 0;

    const unsigned NOT_VERIFIED = 1u << (8 * sizeof(unsigned) - 1);

    TOptions & options = tls.options;

    TFilterPattern & filterPatternL = *job.filterPatternL;
    TFilterPattern & filterPatternR = *job.filterPatternR;
    (void)filterPatternR;
    // TFilterFinderL & filterFinderL = *job.filterFinderL;
    // TFilterFinderR & filterFinderR = *job.filterFinderR;
    TReadSet & readSetL = indexText(host(*job.filterPatternL));
    TReadSet & readSetR = indexText(host(*job.filterPatternR));
    TVerifier & verifierL = tls.verifierL;
    TVerifier & verifierR = tls.verifierR;

    if (empty(readSetL))
        return;

    // distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
    // distance >= libLen - libErr - 2*parWidth + shapeLen
    TSize readLength = length(readSetL[0]);
    TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(filterPatternL)));
    TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(filterPatternL)));
    TGPos scanShift = (minDistance < 0) ? 0 : minDistance;

    Pair<TGPos> gPair;

    // Make sure the queue is empty and all entries in the last potential
    // match no map point before the first no in this queue.
    clear(tls.fifo);
    tls.fifoFirstNo = 0;
    tls.fifoLastNo = 0;
    // TODO(holtgrew): Could do length(...)/job.stride but would have to do so everywhere else below, too.
    // TODO(holtgrew): Get around the clear() and resize-with-fill somehow.
    clear(tls.fifoLastPotMatchNo);
    resize(tls.fifoLastPotMatchNo, length(indexText(host(filterPatternL))), (int64_t) - 2, Exact());

    TGenome & genome = tls.globalStore->contigStore[job.contigId].seq;

    TSize gLength = length(genome);

    TMatch mR;
    TDequeueValue fL(-2, mR);   // to supress uninitialized warnings

    //	unsigned const preFetchMatches = 2048;

    // -----------------------------------------------------------------------
    // Enqueue hits from the previous left window.
    // -----------------------------------------------------------------------
    // if (job.prevHitsPtrL.get() == 0)
    //     std::cerr << "\nLENGTH NO PREV HITS" << std::endl;
    // else
    //     std::cerr << "\nLENGTH length(*job.prevHitsPtrL) == " << length(*job.prevHitsPtrL) << std::endl;
    if (job.prevHitsPtrL.get() != 0 && job.prevHitsLEnd > job.prevHitsLBegin)
    {
        THitStringIter it = iter(*job.prevHitsPtrL, job.prevHitsLEnd, Standard());
        THitStringIter itBegin = iter(*job.prevHitsPtrL, job.prevHitsLBegin, Standard());

        do
        {
            --it;
            if (it->hstkPos + maxDistance + DELTA < job.rightWindowBegin)
                break;
        }
        while (it != itBegin);

        // it now points either to beginning or left of the actual entry we
        // want to point to.  Move it right if it is not actually left of the entry.
        if (it->hstkPos + maxDistance + DELTA == job.rightWindowBegin)  // TODO(holtgrew): Keep in sync with above expression, except for </==
            ++it;

        // Now enqueue all hits in the overlap that belong to this job.
        for (THitStringIter itEnd = iter(*job.prevHitsPtrL, job.prevHitsLEnd, Standard()); it != itEnd; ++it)
        {
            // std::cerr << " [from left window]" << std::flush;
            gPair = Pair<TGPos, TGPos>(_max(static_cast<TSignedGPos>(0), static_cast<TSignedGPos>(it->hstkPos)), _min(it->hstkPos + it->bucketWidth, static_cast<TSignedGPos>(length(genome))));

            fL.i1 = tls.fifoLastPotMatchNo[it->ndlSeqNo];
            tls.fifoLastPotMatchNo[it->ndlSeqNo] = tls.fifoLastNo++;

            fL.i2.readId = tls.globalStore->matePairStore[threadIdOffset + it->ndlSeqNo].readId[0] | NOT_VERIFIED;
            fL.i2.beginPos = static_cast<TSignedGPos>(it->hstkPos);
            fL.i2.endPos = gPair.i2;

            // std::cerr << "\nPREV SWIFT\tL\t" << tls.globalStore->matePairStore[it->ndlSeqNo].readId[0] << "\t" << tls.globalStore->readNameStore[tls.globalStore->matePairStore[it->ndlSeqNo].readId[0]] << "\t" << fL.i2.beginPos << "\t" << fL.i2.endPos << std::endl;

            pushBack(tls.fifo, fL);
        }
    }

    // -----------------------------------------------------------------------
    // Perform paired-end verification.
    // -----------------------------------------------------------------------
    THitStringIter itL = iter(*job.hitsPtrL, job.hitsLBegin, Standard());
    THitStringIter itEndL = iter(*job.hitsPtrL, job.hitsLEnd, Standard());
    // Iterate over all filtration results are returned by SWIFT.
    for (THitStringIter itR = iter(*job.hitsPtrR, job.hitsRBegin, Standard()), itEndR = iter(*job.hitsPtrR, job.hitsREnd, Standard()); itR != itEndR; ++itR)
    {
        ++options.countFiltration;
#if 0
        CharString pref = prefix(tls.globalStore->readNameStore[2 * (threadIdOffset + itR->ndlSeqNo) + 1], length("EAS20_8_6_1_248_1397"));
        CharString s = "EAS20_8_6_1_248_1397";
        if (pref == s)
            std::cerr << "GOTCHA" << std::endl;
        // SEQAN_OMP_PRAGMA(critical)
        // {
        //     std::cerr << tls.globalStore->readNameStore[2 * (threadIdOffset + itR->ndlSeqNo) + 1] << std::endl;
        // }
#endif
#ifdef RAZERS_DEBUG_MATEPAIRS
        std::cerr << "\nSWIFT\tR\t" << itR->ndlSeqNo << "\t" << tls.globalStore->readNameStore[2 * (threadIdOffset + itR->ndlSeqNo) + 1] << "\t" << scanShift + itR->hstkPos << "\t" << scanShift + itR->hstkPos + itR->bucketWidth << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS

        unsigned matePairId = itR->ndlSeqNo;
        TGPos rEndPos = itR->hstkPos + itR->bucketWidth + scanShift;
        TGPos doubleParWidth = 2 * itR->bucketWidth;

        // (1) Remove out-of-window left mates from fifo.
        while (!empty(tls.fifo) && (TSignedGPos)front(tls.fifo).i2.endPos + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
        {
#ifdef RAZERS_DEBUG_MATEPAIRS
            if ((front(tls.fifo).i2.readId & ~NOT_VERIFIED) > length(tls.globalStore->readNameStore))
                std::cerr << "\nPOP\tL\t" << "[bad read #" << (front(tls.fifo).i2.readId & ~NOT_VERIFIED) << "]\t" << front(tls.fifo).i2.beginPos << "\t" << front(tls.fifo).i2.endPos << std::endl;
            else
                std::cerr << "\nPOP\tL\t" << tls.globalStore->readNameStore[front(tls.fifo).i2.readId & ~NOT_VERIFIED] << "\t" << front(tls.fifo).i2.beginPos << "\t" << front(tls.fifo).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
            popFront(tls.fifo);
            ++tls.fifoFirstNo;
        }

        // (2) Add within-window left mates to fifo.
        while (empty(tls.fifo) || (TSignedGPos)back(tls.fifo).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
        {
            // Get next left-mate hits, if any and go to next.  This
            // corresponds to a find() in the sequential non-window case.
            if (itL != itEndL)
            {
                ++options.countFiltration;
#ifdef RAZERS_DEBUG_MATEPAIRS
                std::cerr << "\nSWIFT\tL\t" << itL->ndlSeqNo << "\t" << tls.globalStore->readNameStore[2 * (threadIdOffset + itL->ndlSeqNo)] << "\t" << itL->hstkPos << "\t" << itL->hstkPos + itL->bucketWidth << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                gPair = Pair<TGPos, TGPos>(_max(static_cast<TSignedGPos>(0), static_cast<TSignedGPos>(itL->hstkPos)), _min(itL->hstkPos + itL->bucketWidth, static_cast<TSignedGPos>(length(genome))));
                if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
                {
                    // link in
                    fL.i1 = tls.fifoLastPotMatchNo[itL->ndlSeqNo];
                    tls.fifoLastPotMatchNo[itL->ndlSeqNo] = tls.fifoLastNo++;

                    // Translate from thread-local to global id.
                    fL.i2.readId = (tls.globalStore->matePairStore[(threadIdOffset + itL->ndlSeqNo)].readId[0]) | NOT_VERIFIED;
                    fL.i2.beginPos = gPair.i1;
                    fL.i2.endPos = gPair.i2;

                    pushBack(tls.fifo, fL);
                }

                ++itL;
            }
            else
            {
                break;
            }
        }

        int bestLeftScore = std::numeric_limits<int>::min();
        int bestLibSizeError = std::numeric_limits<int>::max();
        TDequeueIterator bestLeft = TDequeueIterator();

        bool rightVerified = false;
        TDequeueIterator it;
        unsigned leftReadId = tls.globalStore->matePairStore[(threadIdOffset + matePairId)].readId[0];
        int64_t last = (int64_t) - 1;
        int64_t lastValid = (int64_t) - 1;
        int64_t i;
        for (i = tls.fifoLastPotMatchNo[matePairId]; tls.fifoFirstNo <= i; last = i, i = (*it).i1)
        {
            it = &value(tls.fifo, i - tls.fifoFirstNo);

            // search left mate
            //			if (((*it).i2.readId & ~NOT_VERIFIED) == leftReadId)
            //			        ^== we need not to test anymore, as only corr. left mates are traversed
            //						via the linked list beginning from lastPotMatchNo[matePairId]
            {
                // verify left mate (equal seqNo), if not done already
                if ((*it).i2.readId & NOT_VERIFIED)
                {
                    if ((TSignedGPos)(*it).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
                    {
#ifdef RAZERS_BANDED_MYERS
                        verifierL.patternState.leftClip = ((*it).i2.beginPos >= 0) ? 0 : -(*it).i2.beginPos;  // left clip if match begins left of the genome
#endif
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "\nVERIFY\tL\t" << matePairId << "\t" << tls.globalStore->readNameStore[2 * (threadIdOffset + matePairId)] << "\t" << (TSignedGPos)(*it).i2.beginPos << "\t" << (*it).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        ++options.countVerification;
//                        if (i==0)
//                        std::cout<<"here"<<std::endl;

                        // adjust sink position according to insert size
                        if (!rightVerified)
                            verifierL.sinkPos = (TSignedGPos)(itR->hstkPos + itR->bucketWidth) - options.libraryLength;

                        if (matchVerify(verifierL, infix(genome, ((*it).i2.beginPos >= 0) ? (TSignedGPos)(*it).i2.beginPos : (TSignedGPos)0, (TSignedGPos)(*it).i2.endPos),
                                        leftReadId, readSetL[matePairId], TRazerSMode()))
                        {
#ifdef RAZERS_DEBUG_MATEPAIRS
                            std::cerr << "  YES: " << verifierL.m.beginPos << "\t" << verifierL.m.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
//                            std::cerr << "  Left+ " << verifierL.m.endPos << std::endl;
                            verifierL.m.readId = (*it).i2.readId & ~NOT_VERIFIED;       // has been verified positively
                            (*it).i2 = verifierL.m;
                        }
                        else
                        {
                            (*it).i2.readId = ~NOT_VERIFIED;                // has been verified negatively
#ifdef RAZERS_DEBUG_MATEPAIRS
                            std::cerr << "  NO" << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                            continue;                                       // we intentionally do not set lastPositive to i
                        }                                                   // to remove i from linked list
                    }
                    else
                    {
                        lastValid = i;
                        continue;                                           // left pot. hit is out of tolerance window
                    }
                } //else {}													// left match is verified already

                // short-cut negative matches
                if (last != lastValid)
                {
                    SEQAN_ASSERT_NEQ(lastValid, i);
                    if (lastValid == (int64_t) - 1)
                        tls.fifoLastPotMatchNo[matePairId] = i;
                    else
                        value(tls.fifo, lastValid - tls.fifoFirstNo).i1 = i;
                }
                lastValid = i;

                if (!rightVerified)                                         // here a verfied left match is available
                {
#ifdef RAZERS_DEBUG_MATEPAIRS
                    std::cerr << "\nVERIFY\tR\t" << matePairId << "\t" << tls.globalStore->readNameStore[2 * (threadIdOffset + matePairId) + 1] << "\t" << itR->hstkPos << "\t" << itR->hstkPos + itR->bucketWidth << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                    ++options.countVerification;
//                    if (pref == s)
//                        std::cerr << "\nVERIFY\tR\t" << matePairId << "\t" << tls.globalStore->readNameStore[2 * (threadIdOffset + matePairId) + 1] << "\t" << itR->hstkPos << "\t" << itR->hstkPos + itR->bucketWidth << std::endl;
                    if (matchVerify(verifierR, swiftInfix(*itR, tls.genomeInf), 2 * (threadIdOffset + matePairId) + 1, readSetR[matePairId], TRazerSMode()))
                    {
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "  YES: " << verifierR.m.beginPos << "\t" << verifierR.m.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        rightVerified = true;
                        mR = verifierR.m;
//                        if (pref == s)
//                            std::cerr << "BREAK HERE" << std::endl;
                    }
                    else
                    {
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "  NO" << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        // Break out of lastPotMatch loop, rest of find(right SWIFT results loop will not
                        // be executed since bestLeftScore remains untouched.
                        i = (*it).i1;
                        break;
                    }
                }

                /*
                  if ((*it).i2.readId == leftReadId)
                  {
                  bestLeft = it;
                  bestLeftScore = (*it).i3.score;
                  break;
                  }
                */
                if ((*it).i2.readId == leftReadId)
                {
                    int score = (*it).i2.score;
                    if (bestLeftScore <= score)
                    {
                        // distance between left mate beginning and right mate end
                        int64_t dist = (int64_t)verifierR.m.endPos - (int64_t)(*it).i2.beginPos;
                        int libSizeError = options.libraryLength - dist;
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "    libSizeError = " << libSizeError << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        if (libSizeError < 0)
                            libSizeError = -libSizeError;
                        if (libSizeError > options.libraryError)
                            continue;
                        if (bestLeftScore == score)
                        {
                            if (bestLibSizeError > libSizeError)
                            {
                                bestLibSizeError = libSizeError;
                                bestLeft = it;
                            }
                        }
                        else
                        {
                            bestLeftScore = score;
                            bestLibSizeError = libSizeError;
                            bestLeft = it;
                            // if (bestLeftScore == 0) break;	// TODO: replace if we have real qualities
                        }
                    }
                }
            }
        }

        // (3) Short-cut negative matches.
        if (last != lastValid)
        {
            SEQAN_ASSERT_NEQ(lastValid, i);
            if (lastValid == (int64_t) - 1)
                tls.fifoLastPotMatchNo[matePairId] = i;
            else
                value(tls.fifo, lastValid - tls.fifoFirstNo).i1 = i;
        }

        // (4) Verify right mate, if left mate matches.
        if (bestLeftScore != std::numeric_limits<int>::min())
        {
            fL.i2 = (*bestLeft).i2;

            // transform mate readNo to global readNo
            TMatePair & mp     = tls.globalStore->matePairStore[(threadIdOffset + matePairId)];
            fL.i2.readId      = mp.readId[0];
            mR.readId         = mp.readId[1];
            mR.orientation    = (job.orientation == 'F') ? 'R' : 'F';
            fL.i2.orientation = job.orientation;

            // transform coordinates to the forward strand
            if (job.orientation == 'F')
            {
                TSize temp = mR.beginPos;
                mR.beginPos = mR.endPos;
                mR.endPos = temp;
            }
            else
            {
                fL.i2.beginPos = gLength - fL.i2.beginPos;
                fL.i2.endPos = gLength - fL.i2.endPos;
                TSize temp = mR.beginPos;
                mR.beginPos = gLength - mR.endPos;
                mR.endPos = gLength - temp;
                // dist = -dist;
            }

            // set a unique pair id
            fL.i2.pairMatchId = mR.pairMatchId = options.nextPairMatchId * options.threadCount + tls.threadId;
            if (++options.nextPairMatchId == TMatch::INVALID_ID)
                options.nextPairMatchId = 0;

            // score the whole match pair
            fL.i2.pairScore = mR.pairScore = fL.i2.score + mR.score;
            fL.i2.libDiff = mR.libDiff = bestLibSizeError;

            // both mates match with correct library size
            /*								std::cout << "found " << matePairId << " on " << orientation << contigId;
                                            std::cout << " dist:" << dist;
                                            if (orientation=='F')
                                            std::cout << " \t_" << fL.i2.beginPos+1 << "_" << mR.endPos;
                                            else
                                            std::cout << " \t_" << mR.beginPos+1 << "_" << mL.endPos;
                                            //							std::cout << " L_" << (*bestLeft).beginPos << "_" << (*bestLeft).endPos << "_" << (*bestLeft).editDist;
                                            //							std::cout << " R_" << mR.beginPos << "_" << mR.endPos << "_" << mR.editDist;
                                            std::cout << std::endl;
            */
            if (!options.spec.DONT_DUMP_RESULTS)
            {
                appendValue(*localMatches, fL.i2, Generous());
                appendValue(*localMatches, mR, Generous());
                SEQAN_ASSERT_EQ(length(*localMatches) % 2, 0u);
//                        if (pref == s) {
//                            std::cerr << "ADDED" << std::endl;
//                            std::cerr << "  (" << mR.beginPos << ", " << mR.endPos << ")" << std::endl;
//                            std::cerr << "  (" << fL.i2.beginPos << ", " << fL.i2.endPos << ")" << std::endl;
//                        }

#ifdef RAZERS_DEBUG_MATEPAIRS
                std::cerr << "\nHIT\tL\t" << fL.i2.readId << "\t" << tls.globalStore->readNameStore[fL.i2.readId] << "\t" << fL.i2.beginPos << "\t" << fL.i2.endPos << std::endl;
                std::cerr << "\nHIT\tR\t" << mR.readId << "\t" << tls.globalStore->readNameStore[mR.readId] << "\t" << mR.beginPos << "\t" << mR.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
            }
            // }
        }
        // }
    }

    appendToVerificationResults(*job.verificationResults, localMatches);

    tls.options.timeVerification += sysTime() - startVerify;
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE
}

template <typename TMatches, typename TFragmentStore, typename TFilterFinderL, typename TFilterFinderR, typename TFilterPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
void
writeBackToLocal(ThreadLocalStorage<MapPairedReads<TMatches, TFragmentStore, TFilterFinderL, TFilterFinderR, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, String<TMatches *> & verificationHits, bool dontCompact)
{
    // TODO(holtgrew): The same as single read writeback! Really? Combine?
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
    if (tls.threadId == 0u && tls.options._debugLevel >= 3)
        fprintf(stderr, "[writeback]");
    typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
    //typedef ThreadLocalStorage<MapPairedReads<TMatches, TFragmentStore, TFilterFinderL, TFilterFinderR, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > TThreadLocalStorage;

    // Compute new length for the matches.
    TAlignedReadStoreSize oldSize = length(tls.matches);
    TAlignedReadStoreSize newSize = oldSize;

    for (unsigned i = 0; i < length(verificationHits); ++i)
    {
        SEQAN_ASSERT_EQ(length(*verificationHits[i]) % 2, 0u);
        TMatches * bucket = verificationHits[i];
        newSize += length(*bucket);
    }

    // Write back all matches from verification to the block local store.
    reserve(tls.matches, newSize, Generous());
    SEQAN_ASSERT_EQ(length(tls.matches) % 2, 0u);
    for (unsigned i = 0; i < length(verificationHits); ++i)
        append(tls.matches, *verificationHits[i]);
    SEQAN_ASSERT_EQ(length(tls.matches) % 2, 0u);

#ifdef RAZERS_DEFER_COMPACTION
    (void) dontCompact;

    // Register found mate pairs in match filter.
    typedef typename Iterator<TMatches>::Type TIterator;
    TIterator itBegin = begin(tls.matches, Standard()) + oldSize;
    TIterator itEnd = end(tls.matches, Standard());
    TIterator it = itBegin;
    for (; it != itEnd; ++it)
    {
        if (it->orientation == '-')
            continue;  // Skip masked reads.
        if (it->isRegistered)
            continue;
        it->isRegistered = true;
        registerRead(*tls.matchFilter, it->readId / 2, it->pairScore);

#if SEQAN_ENABLE_DEBUG
        TIterator it2 = it;
#endif  // #if SEQAN_ENABLE_DEBUG
        ++it;  // Skip second mate.
        SEQAN_ASSERT_EQ(it->readId / 2, it2->readId / 2);
    }
    // Limit/distable reads.
    unsigned disabled = 0;
    for (it = itBegin; it != itEnd; ++it, ++it)
    {
        if (it->orientation == '-')
            continue;  // Skip masked reads.
        disabled += processRead(*tls.matchFilter, it->readId / 2);
    }
    if (tls.options._debugLevel >= 2)
        fprintf(stderr, " [%u reads disabled]", disabled);

#else  // #ifdef RAZERS_DEFER_COMPACTION
       // Possibly compact matches.
    if (!dontCompact && length(tls.matches) > tls.options.compactThresh)
    {
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_COMPACT);
#endif  // #ifdef RAZERS_PROFILE
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typename Size<TAlignedReadStore>::Type oldSize = length(tls.matches);

        if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE || tls.options.threshold == 0)
            maskDuplicates(tls.matches, tls.options, TRazerSMode());    // overlapping parallelograms cause duplicates

        FilterPatternLSetMaxErrorsWrapper<TThreadLocalStorage> wrapperL(tls);
        FilterPatternRSetMaxErrorsWrapper<TThreadLocalStorage> wrapperR(tls);
        compactPairMatches(*tls.globalStore, tls.matches, tls.counts, tls.options, wrapperL, wrapperR, COMPACT);

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

// Find read matches in one genome sequence.
//
// The parallelization is similar to the single-end mode.  We perform the
// filtering window-wise.  Then, we generate verification jobs from the swift
// hits.  Here, we distribute them to the jobs by a hash function (num % len),
// constructable and reconstructable by the stride (result of modulo) and the
// module number.  These are then processed by all "leading" threads, where
// leading threads are those with the largest number of processed windows.
template <
    typename TFSSpec,
    typename TFSConfig,
    typename TThreadLocalStorages,
    typename TCounts,
    typename TRazerSOptions,
    typename TRazerSMode>
void _mapMatePairReadsParallel(
    FragmentStore<TFSSpec, TFSConfig> & store,
    unsigned                                  contigId,             // ... and its sequence number
    TThreadLocalStorages & threadLocalStorages,
    String<unsigned> const & splitters,
    TCounts & /*cnts*/,                                  // TODO(holtgrew): What about this?
    char                                      orientation,          // q-gram index of reads
    TRazerSOptions & options,
    TRazerSMode                       const & /*mode*/)
{
    typedef FragmentStore<TFSSpec, TFSConfig>               TFragmentStore;
    //typedef typename TFragmentStore::TMatePairStore         TMatePairStore;
    //typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    //typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;
    //typedef typename Value<TMatePairStore>::Type            TMatePair;
    //typedef typename Value<TAlignedReadStore>::Type         TAlignedRead;
    //typedef typename Value<TAlignQualityStore>::Type        TAlignQuality;
    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    //typedef typename TThreadLocalStorage::TReadSet TReadSet;
    //typedef typename TThreadLocalStorage::TShape TShape;
    //typedef Index<TReadSet, IndexQGram<TShape>  >   TReadIndex;
    //typedef typename Id<TAlignedRead>::Type                 TId;

    typedef typename TFragmentStore::TContigSeq             TGenome;
    //typedef typename Size<TGenome>::Type                    TSize;
    typedef typename Position<TGenome>::Type                TGPos;
    typedef typename MakeSigned_<TGPos>::Type               TSignedGPos;
    //typedef typename Infix<TGenome>::Type                   TGenomeInf;

    typedef TRazerSOptions TOptions;

    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TMatches TMatches;
    //typedef typename Value<TMatches>::Type TMatch;

    // FILTRATION
    typedef typename TThreadLocalStorage::TFilterFinderL    TFilterFinderL;
    typedef typename TThreadLocalStorage::TFilterFinderR    TFilterFinderR;
    typedef typename TThreadLocalStorage::TFilterPattern    TFilterPattern;

    // MATE-PAIR FILTRATION
    //typedef Pair<int64_t, TMatch>   TDequeueValue;
    //typedef Dequeue<TDequeueValue>                          TDequeue;
    //typedef typename TDequeue::TIter                        TDequeueIterator;

    // VERIFICATION
    //typedef MatchVerifier<
    //    TFragmentStore,
    //    TMatches,
    //    TRazerSOptions,
    //    TRazerSMode,
    //    TFilterPattern,
    //    TCounts
    //    > TVerifier;

    typedef typename TFilterFinderL::THitString THitString;
    typedef Job<PairedVerification<TMatches, TFragmentStore, THitString, TOptions, TFilterPattern> > TVerificationJob;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE

    // iterate all genomic sequences
    if (options._debugLevel >= 1)
    {
        std::cerr << std::endl << "Process genome seq #" << contigId;
        if (orientation == 'F')
            std::cerr << "[fwd]";
        else
            std::cerr << "[rev]";
    }

    // -----------------------------------------------------------------------
    // Guard against too small contigs.
    // -----------------------------------------------------------------------

    // distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
    // distance >= libLen - libErr - 2*parWidth + shapeLen
    // TSize readLength = length(threadLocalStorages[0].readSetL[0]); // XXX
    // TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(threadLocalStorages[0].filterPatternL))); // XXX
    TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(threadLocalStorages[0].filterPatternL)));
    TGPos scanShift = (minDistance < 0) ? 0 : minDistance;

    // exit if contig is shorter than library size
    TGenome & genome = store.contigStore[contigId].seq;
    if (length(genome) <= scanShift)
        return;

    // -----------------------------------------------------------------------
    // Reverse-complement contig if necessary.
    // -----------------------------------------------------------------------
    // lockContig(store, contigId);
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE
    if (orientation == 'R')
        reverseComplement(genome);
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE

    // -----------------------------------------------------------------------
    // Per-contig initialization of thread local storage objects.
    // -----------------------------------------------------------------------
    // TODO(holtgrew): Maybe put into its own function?
    SEQAN_OMP_PRAGMA(parallel for schedule(static, 1))
    for (int i = 0; i < static_cast<int>(options.threadCount); ++i)
    {
        // Initialize verifier objects.
        threadLocalStorages[i].verifierL.onReverseComplement = (orientation == 'R');
        threadLocalStorages[i].verifierL.genomeLength = length(genome);
        threadLocalStorages[i].verifierL.oneMatchPerBucket = true;
        threadLocalStorages[i].verifierL.m.contigId = contigId;

        threadLocalStorages[i].verifierR.onReverseComplement = (orientation == 'R');
        threadLocalStorages[i].verifierR.genomeLength = length(genome);
        threadLocalStorages[i].verifierR.oneMatchPerBucket = true;
        threadLocalStorages[i].verifierR.m.contigId = contigId;

        // threadLocalStorages[i].filterFinderL  = TFilterFinderL(genome, options.repeatLength, 1);
        typedef typename TFragmentStore::TContigSeq TGenome;
        typedef typename Infix<TGenome>::Type       TGenomeInfix;
        TGenomeInfix inf(genome, scanShift, length(genome));
        set(threadLocalStorages[i].genomeInf, inf);
        // threadLocalStorages[i].filterFinderR  = TFilterFinderR(threadLocalStorages[i].genomeInf, options.repeatLength, 1);
    }

    // -----------------------------------------------------------------------
    // Perform filtration.
    // -----------------------------------------------------------------------
    TaskQueue<TVerificationJob, OmpLock> taskQueue;
    volatile unsigned leaderWindowsDone = 0;  // Number of windows done in leaders.
    volatile unsigned threadsFiltering = options.threadCount;

    // We will create the swift finders for thread 0 first.  This will trigger parallel repeat finding in the SWIFT
    // finder construction.  Then, we copy over the finder to all threads.
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_COPY_FINDER);
#endif  // #ifdef RAZERS_PROFILE

    threadLocalStorages[0].filterFinderL  = TFilterFinderL(store.contigStore[contigId].seq, options.repeatLength, 1);
    threadLocalStorages[0].filterFinderR  = TFilterFinderR(threadLocalStorages[0].genomeInf, options.repeatLength, 1);

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
        {
            tls.filterFinderL = threadLocalStorages[0].filterFinderL;
            tls.filterFinderR = threadLocalStorages[0].filterFinderR;
        }
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_COPY_FINDER);
#endif  // #ifdef RAZERS_PROFILE
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE


        TFilterPattern & filterPatternL = tls.filterPatternL;
        TFilterPattern & filterPatternR = tls.filterPatternR;
        TFilterFinderL & filterFinderL = tls.filterFinderL;
        TFilterFinderR & filterFinderR = tls.filterFinderR;

#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE
        if (!windowFindBegin(filterFinderL, filterPatternL, tls.options.errorRate))
            std::cerr << "ERROR: windowFindBegin() failed for left reads in thread " << tls.threadId << std::endl;
        if (!windowFindBegin(filterFinderR, filterPatternR, tls.options.errorRate))
            std::cerr << "ERROR: windowFindBegin() failed for right reads in thread " << tls.threadId << std::endl;
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

        // Previous left hits will be stored here.
        std::shared_ptr<THitString> previousLeftHits;
        std::shared_ptr<THitString> leftHits;
        std::shared_ptr<THitString> rightHits;
        // Declare hits splitters and initialize for left since we need it as "previous left" below.
        String<size_t> leftHitsSplitters;
        resize(leftHitsSplitters, options.maxVerificationPackageCount + 1, 0);
        String<size_t> rightHitsSplitters;
        String<size_t> previousLeftHitsSplitters;

        // For each filtration window...
        bool hasMore = !empty(tls.readSetL);
        while (hasMore)
        {
#ifdef RAZERS_PROFILE
            timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE
            double filterStart = sysTime();

            // TDequeue fifo;						// stores left-mate potential matches  // XXX
            // String<int64_t> lastPotMatchNo;		// last number of a left-mate potential // XXX
            // int64_t lastNo = 0;					// last number over all left-mate pot. matches in the queue // XXX
            // int64_t firstNo = 0;				// first number over all left-mate pot. match in the queue // XXX
            // Pair<TGPos> gPair;

            // resize(lastPotMatchNo, length(host(filterPatternL)), (int64_t)-1, Exact());  // XXX

            // TSize gLength = length(genome);  // XXX

            // TAlignedRead mR;
            // TAlignQuality qR;
            // TDequeueValue fL(-1, mR, qR);	// to supress uninitialized warnings

            // Search for hits from next window.
            int delta = windowsDone == 0 ? scanShift : 0;  // First window of right finder is smaller.
            hasMore = windowFindNext(tls.filterFinderR, tls.filterPatternR, tls.options.windowSize - delta);
            bool ret = windowFindNext(tls.filterFinderL, tls.filterPatternL, tls.options.windowSize);
            (void) ret;
            SEQAN_ASSERT(ret == hasMore);
            windowsDone += 1;  // Local windows done count.
            atomicMax(leaderWindowsDone, windowsDone);

            size_t rightWindowBegin = beginPosition(tls.genomeInf) + tls.options.windowSize * (windowsDone - 1);

            // Create verification jobs.
            // std::cerr << "\nSWIFT HIT COUNT\t" << length(getSwiftHits(tls.filterFinderL)) << "\t" << length(getSwiftHits(tls.filterFinderR)) << std::endl;
            // if (windowsDone > 1)
            //     std::cerr << "\nPOSITION       \t" << tls.filterFinderL.curPos << "\t" << tls.filterFinderR.curPos << std::endl;
            if (length(getWindowFindHits(tls.filterFinderL)) > 0u || length(getWindowFindHits(tls.filterFinderR)) > 0u)
            {
                using std::swap;
                String<TVerificationJob> jobs;

                // Update previous left hits and splitters.
                previousLeftHits = leftHits;
                swap(previousLeftHitsSplitters, leftHitsSplitters);
                // Update new left hits and splitters.
                leftHits.reset(new THitString());
                resize(*leftHits, 1);
                clear(*leftHits);
                swap(*leftHits, getWindowFindHits(tls.filterFinderL));
                buildHitSplittersAndPartitionHits(leftHitsSplitters, *leftHits, options.maxVerificationPackageCount, length(indexText(host(filterPatternL))));
                rightHits.reset(new THitString());
                swap(*rightHits, getWindowFindHits(tls.filterFinderR));
                buildHitSplittersAndPartitionHits(rightHitsSplitters, *rightHits, options.maxVerificationPackageCount, length(host(filterPatternL)));

                for (unsigned i = 0; i < options.maxVerificationPackageCount; ++i)
                {
                    // std::cerr << "i == " << i << " length(leftHitsSplitters) == " << length(leftHitsSplitters) << ", length(rightHitsSplitters) == " << length(rightHitsSplitters) << ", length(previousLeftHitsSplitters) == " << length(previousLeftHitsSplitters) << ", verificationPackageCount == " << options.maxVerificationPackageCount << std::endl;
                    // std::cerr << "length(leftHits) == " << length(leftHits) << ", length(rightHits) == " << length(rightHits) << std::endl;
                    SEQAN_ASSERT_LEQ(previousLeftHitsSplitters[i], previousLeftHitsSplitters[i + 1]);
                    SEQAN_ASSERT_LEQ(rightHitsSplitters[i], rightHitsSplitters[i + 1]);
                    SEQAN_ASSERT_LEQ(leftHitsSplitters[i], leftHitsSplitters[i + 1]);
                    if (rightHitsSplitters[i] == rightHitsSplitters[i + 1] || (previousLeftHitsSplitters[i] == previousLeftHitsSplitters[i + 1] && leftHitsSplitters[i] == leftHitsSplitters[i + 1]))
                        continue;
                    appendValue(jobs, TVerificationJob(tls.threadId, tls.verificationResults, store, contigId, orientation, previousLeftHits, previousLeftHitsSplitters[i], previousLeftHitsSplitters[i + 1], leftHits, leftHitsSplitters[i], leftHitsSplitters[i + 1], rightHits, rightHitsSplitters[i], rightHitsSplitters[i + 1], rightWindowBegin, *tls.globalOptions, tls.filterPatternL, tls.filterPatternR));
                }

                pushFront(taskQueue, jobs);
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
                workVerification(tls, job, splitters);
            }

            // Write back verification results for this thread so far.
            //
            // First, swap out the current set of local stores from the verification results.
            omp_set_lock(&tls.verificationResults.lock->lock_);
            String<TMatches *> localMatches;
            using std::swap;
            swap(localMatches, tls.verificationResults.localMatches);
            omp_unset_lock(&tls.verificationResults.lock->lock_);
            // Don't compact matches if in configured 'block fraction' of genome.
            size_t hstckLen = filterFinderR.endPos - filterFinderR.startPos;
            size_t hstckLeft = filterFinderR.endPos - filterFinderR.curPos;
            double fracTodo = 1.0  * hstckLeft / hstckLen;
            bool dontCompact = tls.options.noCompactFrac >= fracTodo;
            // Write back the contents of these stores to the thread-local store.
            writeBackToLocal(tls, localMatches, dontCompact);
            clearLocalMatches(localMatches);
        }

        // Finalization
        windowFindEnd(filterFinderL, filterPatternL);
        windowFindEnd(filterFinderR, filterPatternR);

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
        clearLocalMatches(tls.verificationResults.localMatches);
    }

    // NOTE: We never re-reverse complement since this function is only called
    // twice, in the right order regarding the orientation parameters.
    // if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
    //  if (orientation == 'R')	reverseComplement(genome);	// we have to restore original orientation
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE
}

//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)

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
int _mapMatePairReadsParallel(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy>  const & mode,
    TFilterSpec)
{
    typedef FragmentStore<TFSSpec, TFSConfig>           TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
    typedef RazerSCoreOptions<TSpec> TOptions;

    typedef typename Value<TReadSeqStore>::Type         TRead;
    typedef StringSet<TRead>                            TReadSet;
#ifndef RAZERS_OPENADDRESSING
    typedef Index<TReadSet, IndexQGram<TShape> >                    TIndex;         // q-gram index
#else
    typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >    TIndex;
#endif

    typedef Pattern<TIndex, TFilterSpec>                TFilterPattern; // filter
    //typedef Pattern<TRead const, MyersUkkonen>          TMyersPattern;  // verifier

    typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef Finder<TContigSeq, TFilterSpec>             TFilterFinderL;
    typedef typename Infix<TContigSeq>::Type            TContigInf;
    typedef Finder<TContigInf, TFilterSpec>             TFilterFinderR;

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

//TODO:REVCOMP right mates
/*
    options.compMask[4] = (options.matchN)? 15: 0;
    if (options.gapMode == RAZERS_GAPPED)
    {
        int pairCount = length(store.matePairStore);
        resize(forwardPatternsL, pairCount, Exact());
        resize(forwardPatternsR, pairCount, Exact());
        String<Dna5String> tmps;
        resize(tmps, omp_get_max_threads());
        SEQAN_OMP_PRAGMA(parallel for schedule(static))
        for (int i = 0; i < pairCount; ++i)
        {
#ifdef RAZERS_NOOUTERREADGAPS
            if (!empty(store.readSeqStore[store.matePairStore[i].readId[0]]) && !empty(store.readSeqStore[store.matePairStore[i].readId[1]]))
            {
                setHost(forwardPatternsL[i], prefix(store.readSeqStore[store.matePairStore[i].readId[0]], length(store.readSeqStore[store.matePairStore[i].readId[0]]) - 1));
                tmps[omp_get_thread_num()] = prefix(store.readSeqStore[store.matePairStore[i].readId[1]], length(store.readSeqStore[store.matePairStore[i].readId[1]]) - 1);
                reverseComplement(tmps[omp_get_thread_num()]);
                setHost(forwardPatternsR[i], tmps[omp_get_thread_num()]);
            }
#else
            setHost(forwardPatternsL[i], store.readSeqStore[store.matePairStore[i].readId[0]]);
            tmps[omp_get_thread_num()] = store.readSeqStore[store.matePairStore[i].readId[1]];
            reverseComplement(tmps[omp_get_thread_num()]);
            setHost(forwardPatternsR[i], tmps[omp_get_thread_num()]);
#endif
            _patternMatchNOfPattern(forwardPatternsL[i], options.matchN);
            _patternMatchNOfPattern(forwardPatternsR[i], options.matchN);
            _patternMatchNOfFinder(forwardPatternsL[i], options.matchN);
            _patternMatchNOfFinder(forwardPatternsR[i], options.matchN);
        }
    }
#endif  // #ifdef RAZERS_BANDED_MYERS
*/

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
    computeSplittersBySlotCount(splitters, length(store.matePairStore), options.threadCount);
    typedef ThreadLocalStorage<MapPairedReads<TMatches, TFragmentStore, TFilterFinderL, TFilterFinderR, TFilterPattern, TShape, TOptions, TCounts, TRazerSMode> > TThreadLocalStorage;
    String<TThreadLocalStorage> threadLocalStorages;
    initializeThreadLocalStoragesPaired(threadLocalStorages, store, splitters, shape, options);

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

    SEQAN_PROTIMESTART(findTime);
    for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId)
    {
        // lock to prevent releasing and loading the same contig twice
        // (once per _mapSingleReadsToContig call)
        lockContig(store, contigId);
        if (options.forward)
            _mapMatePairReadsParallel(store, contigId, threadLocalStorages, splitters, cnts, 'F', options, mode);
        if (options.reverse)
            _mapMatePairReadsParallel(store, contigId, threadLocalStorages, splitters, cnts, 'R', options, mode);
        unlockAndFreeContig(store, contigId);
    }
#ifdef RAZERS_PROFILE
    double endMapping = sysTime();
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

    // Switch between using parallel compaction, sequential compaction, and
    // sequential compaction with external sorting.  The actual switch for the
    // sorting is in function compactMatches().
    if (useSequentialCompaction || useExternalSort)
    {
        for (unsigned i = 0; i < length(threadLocalStorages); ++i)
        {
            if (IsSameType<TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
                maskDuplicates(threadLocalStorages[omp_get_thread_num()].matches, options, mode);
            Nothing nothing;
            CompactMatchesMode compactMode = useSequentialCompaction ? COMPACT_FINAL : COMPACT_FINAL_EXTERNAL;
            // std::cerr << "BEFORE FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
            compactPairMatches(store, threadLocalStorages[omp_get_thread_num()].matches, cnts, options, mode, nothing, compactMode);
            // std::cerr << "AFTER FINAL COMPACTION " << length(threadLocalStorages[omp_get_thread_num()].matches) << std::endl;
        }
    }
    else
    {
#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
        // Write back local stores to global stores.
    if (options._debugLevel >= 2)
        for (unsigned i = 0; i < length(threadLocalStorages); ++i)
            std::cerr << "thread " << i << " has " << length(threadLocalStorages[i].matches) << " aligned reads." << std::endl;
    // Mask duplicates and compact matches in parallel in-memory if not using external
    // matches (queried through macros) and or when using external matches and using
    // parallel compaction (!useSequentialCompaction) and using internal sorting
    // (!useExternalSort).
    SEQAN_OMP_PRAGMA(parallel)
    {
        if (IsSameType<TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
            maskDuplicates(threadLocalStorages[omp_get_thread_num()].matches, options, mode);
        Nothing nothing;
        compactPairMatches(store, threadLocalStorages[omp_get_thread_num()].matches, cnts, options, nothing, nothing, COMPACT_FINAL);
    }
#ifdef RAZERS_EXTERNAL_MATCHES
}

#endif // #ifdef RAZERS_EXTERNAL_MATCHES

    writeBackToGlobalStore(store, threadLocalStorages, false);
#ifdef RAZERS_PROFILE
    double endWriteback = sysTime();
    std::cerr << "TIME back to global: " << (endWriteback - endMapping) << " s" << std::endl;
#endif  // #ifdef RAZERS_PROFILE

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

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE

    // restore original orientation (R-reads are infixes of ConcatDirect StringSet)
    SEQAN_OMP_PRAGMA(parallel for schedule(static, 1))
    for (int i = 0; i < (int)length(threadLocalStorages); ++i)
        reverseComplement(threadLocalStorages[i].readSetR);

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE

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
        std::cerr << "Filtration counter:  " << options.countFiltration << std::endl;
        std::cerr << "Verification counter: " << options.countVerification << std::endl;
    }

    // Restore global state.
    omp_set_num_threads(oldMaxThreads);

    return 0;
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
int _mapMatePairReadsParallel(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & mode)
{
    if (options.threshold > 0)
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
        return _mapMatePairReadsParallel(store, cnts, options, shape, mode, Swift<TSwiftSpec>());
    }
    else
    {
        typedef typename If<IsSameType<TGapMode, RazerSGapped>, void, Hamming_>::Type TPigeonholeSpec;
        return _mapMatePairReadsParallel(store, cnts, options, Shape<Dna, OneGappedShape>(), mode, Pigeonhole<TPigeonholeSpec>());
    }
}

} // End namespace

#endif
