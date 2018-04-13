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

#ifndef SEQAN_HEADER_RAZERS_MATEPAIRS_H
#define SEQAN_HEADER_RAZERS_MATEPAIRS_H

#include <seqan/misc/dequeue.h>

namespace seqan {

// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.


//////////////////////////////////////////////////////////////////////////////
// Definitions

#ifdef RAZERS_MEMOPT

template <typename TMPReadSet, typename TShape>
struct SAValue<Index<TMPReadSet, IndexQGram<TShape> > >
{
    typedef Pair<
        unsigned,
        unsigned,
        BitCompressed<24, 8>    // max. 16M reads of length < 256
        > Type;
};

#else

template <typename TMPReadSet, typename TShape>
struct SAValue<Index<TMPReadSet, IndexQGram<TShape> > >
{
    typedef Pair<
        unsigned,               // many reads
        unsigned                // of arbitrary length
        > Type;
};

#endif


// template <typename TMPReadSet, typename TShape, typename TSpec>
// struct Cargo< Index<TMPReadSet, IndexQGram<TShape, TSpec> > > {
//  typedef struct {
//      double		abundanceCut;
//      int			_debugLevel;
//  } Type;
// };

#ifdef RAZERS_PRUNE_QGRAM_INDEX

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TMPReadSet, typename TShape>
inline bool _qgramDisableBuckets(Index<TMPReadSet, IndexQGram<TShape> > & index)
{
    typedef Index<TMPReadSet, IndexQGram<TShape>    >   TReadIndex;
    typedef typename Fibre<TReadIndex, QGramDir>::Type  TDir;
    typedef typename Iterator<TDir, Standard>::Type     TDirIterator;
    typedef typename Value<TDir>::Type                  TSize;

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

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TFSSpec, typename TFSConfig, typename TRazerSOptions>
bool loadReads(
    FragmentStore<TFSSpec, TFSConfig>   & store,
	SeqFileIn                           & leftMates,        // left mates file
    const char                          * fileNameR,        // right mates file
    TRazerSOptions & options)
{
    bool countN = !(options.matchN || options.outputFormat == 1);

    SeqFileIn rightMates;

    bool success;
    if (!isEqual(fileNameR, "-"))
        success = open(rightMates, fileNameR);
    else
        success = open(rightMates, std::cin);

    if (!success)
        return false;


    String<uint64_t>    qualSum;
    String<Dna5Q>       seq[2];
    CharString          qual[2];
    CharString          seqId[2];

    unsigned seqCount = 0;
    unsigned kickoutcount = 0;
    unsigned maxReadLength = 0;

	while (!atEnd(leftMates) && !atEnd(rightMates))
    {
        ++seqCount;
        
        readRecord(seqId[0], seq[0], qual[0], leftMates);
        readRecord(seqId[1], seq[1], qual[1], rightMates);

        if (options.readNaming == 0 || options.readNaming == 3)
        {
            if (!options.fullFastaId)
            {
                cropAfterFirst(seqId[0], IsWhitespace());   // read Fasta id up to the first whitespace
                cropAfterFirst(seqId[1], IsWhitespace());   // read Fasta id up to the first whitespace
            }
            if (options.readNaming == 0)
            {
                append(seqId[0], "/L", Exact());
                append(seqId[1], "/R", Exact());
            }
        }
        else
        {
            clear(seqId[0]);
            clear(seqId[1]);
        }

        if (countN)
        {
            for (int j = 0; j < 2; ++j)
            {
                int maxBase = (int)(0.8 * length(seq[j]));
                int allowed[5] =
                { maxBase, maxBase, maxBase, maxBase, (int)(options.errorRate * length(seq[j]))};
                for (unsigned k = 0; k < length(seq[j]); ++k)
                    if (--allowed[ordValue(getValue(seq[j], k))] == 0)
                    {
//						std::cout << "Ignoring mate-pair: " << seq[0] << " " << seq[1] << std::endl;
                        clear(seq[0]);
                        clear(seq[1]);
                        clear(seqId[0]);
                        clear(seqId[1]);
                        clear(qual[0]);
                        clear(qual[1]);
                        ++kickoutcount;
                        break;
                    }
            }
        }

        for (int j = 0; j < 2; ++j)
        {
            // store dna and quality together
            assignQualities(seq[j], qual[j]);

            if (options.trimLength > 0 && length(seq[j]) > (unsigned)options.trimLength)
                resize(seq[j], options.trimLength);

            unsigned len = length(seq[j]);
            if (length(qualSum) <= len)
            {
                resize(qualSum, len, 0u);
                resize(options.readLengths, len + 1, 0u);
            }
            ++options.readLengths[len];

            for (unsigned i = 0; i < len; ++i)
                qualSum[i] += getQualityValue(seq[j][i]);
        }
        appendMatePair(store, seq[0], seq[1], seqId[0], seqId[1]);
        if (maxReadLength < length(seq[0]))
            maxReadLength = length(seq[0]);
        if (maxReadLength < length(seq[1]))
            maxReadLength = length(seq[1]);
    }


    if (atEnd(leftMates) != atEnd(rightMates))
    {
        if (options._debugLevel > 1)
        {
            std::cerr << "Warning: Unexpected end in one of both paired-end files.\n";
            return false;

        }
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
    if ((unsigned)sa.i2 < maxReadLength - 1)
    {
        std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << std::endl;
        seqCount = 0;
    }

    if (options._debugLevel > 1 && kickoutcount > 0)
        std::cerr << "Ignoring " << kickoutcount << " low quality mate-pairs.\n";
    return seqCount > 0;
}

template <typename TFragmentStore>
struct LessPairScore :
    public std::binary_function<
        typename Value<typename TFragmentStore::TAlignedReadStore>::Type,
        typename Value<typename TFragmentStore::TAlignedReadStore>::Type,
        bool>
{
    TFragmentStore & mainStore;
    TFragmentStore & threadStore;

    LessPairScore(TFragmentStore & _mainStore, TFragmentStore & _threadStore) :
        mainStore(_mainStore), threadStore(_threadStore) {}

    inline bool operator()(
        typename Value<typename TFragmentStore::TAlignedReadStore>::Type const & a,
        typename Value<typename TFragmentStore::TAlignedReadStore>::Type const & b) const
    {
        typedef typename TFragmentStore::TReadStore         TReadStore;
        typedef typename TFragmentStore::TAlignedReadStore  TAlignedReadStore;
        typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;
        typedef typename Value<TReadStore>::Type            TRead;
        typedef typename Value<TAlignedReadStore>::Type     TAlignedRead;
        typedef typename Value<TAlignQualityStore>::Type    TQual;
        //typedef typename Id<TRead>::Type                    TId;

        // pair number
        if (b.readId == TAlignedRead::INVALID_ID) return false;

        if (a.readId == TAlignedRead::INVALID_ID) return true;

        TRead const & ra = mainStore.readStore[a.readId];
        TRead const & rb = mainStore.readStore[b.readId];
        if (ra.matePairId < rb.matePairId) return true;

        if (ra.matePairId > rb.matePairId) return false;

        // quality
        if (a.id == TAlignedRead::INVALID_ID) return false;

        if (b.id == TAlignedRead::INVALID_ID) return true;

        TQual const & qa = threadStore.alignQualityStore[a.id];
        TQual const & qb = threadStore.alignQualityStore[b.id];
        if (qa.pairScore > qb.pairScore) return true;

        if (qa.pairScore < qb.pairScore) return false;

        if (a.libDiff < b.libDiff) return true;

        if (a.libDiff > b.libDiff) return false;

        return a.pairMatchId < b.pairMatchId;
    }

};

template <typename TFragmentStore, typename TReadMatch>
struct LessPairErrors :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    typedef typename TFragmentStore::TReadStore         TReadStore;
    typedef typename Value<TReadStore>::Type            TRead;

    TFragmentStore const & mainStore;

    LessPairErrors(TFragmentStore const & mainStore_) :
        mainStore(mainStore_)
    {}

    inline bool operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (b.readId == TReadMatch::INVALID_ID) return false;

        if (a.readId == TReadMatch::INVALID_ID) return true;

        unsigned matePairIdA = a.readId >> 1;
        unsigned matePairIdB = b.readId >> 1;
        if (matePairIdA < matePairIdB) return true;

        if (matePairIdA > matePairIdB) return false;

        // quality
        if (a.pairScore > b.pairScore) return true;

        if (a.pairScore < b.pairScore) return false;

        if (a.libDiff < b.libDiff) return true;

        if (a.libDiff > b.libDiff) return false;

        if (a.pairMatchId < b.pairMatchId) return true;

        if (a.pairMatchId > b.pairMatchId) return false;

        return a.readId < b.readId;
    }

};


template <typename TFragmentStore, typename TReadMatch>
struct LessPairErrors3Way :
    public std::binary_function<TReadMatch, TReadMatch, bool>
{
    typedef typename TFragmentStore::TReadStore         TReadStore;
    typedef typename Value<TReadStore>::Type            TRead;

    TFragmentStore const & mainStore;

    LessPairErrors3Way(TFragmentStore const & mainStore_) :
        mainStore(mainStore_)
    {}

    inline int operator()(TReadMatch const & a, TReadMatch const & b) const
    {
        // read number
        if (b.readId == TReadMatch::INVALID_ID) return -1;
        if (a.readId == TReadMatch::INVALID_ID) return 1;

        unsigned matePairIdA = a.readId >> 1;
        unsigned matePairIdB = b.readId >> 1;
        if (matePairIdA < matePairIdB) return -1;
        if (matePairIdA > matePairIdB) return 1;

        // quality
        if (a.pairScore > b.pairScore) return -1;
        if (a.pairScore < b.pairScore) return 1;

        if (a.libDiff < b.libDiff) return -1;
        if (a.libDiff > b.libDiff) return 1;

        if (a.pairMatchId < b.pairMatchId) return -1;
        if (a.pairMatchId > b.pairMatchId) return 1;

        if (a.readId < b.readId) return -1;
        if (a.readId > b.readId) return 1;

        return 0;
    }

};

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template <typename TFragmentStore, typename TMatches, typename TCounts, typename TSpec, typename TFilterL, typename TFilterR>
void compactPairMatches(
    TFragmentStore & store,            // all but aligned reads up to the global writeback
    TMatches & matches,                // aligned read only
    TCounts &,
    RazerSCoreOptions<TSpec> & options,
    TFilterL & filterL,
    TFilterR & filterR,
    CompactMatchesMode        compactMode)
{
    typedef typename Value<TMatches>::Type                          TMatch;
    typedef typename Iterator<TMatches, Standard>::Type             TIterator;

    SEQAN_ASSERT_EQ(length(store.alignedReadStore) % 2, 0u);

    // fprintf(stderr, "[pair-compact]");
    double beginTime = sysTime();
    unsigned matePairId = -2;
    unsigned hitCount = 0;
    unsigned hitCountCutOff = options.maxHits;
    int scoreDistCutOff = std::numeric_limits<int>::min();
    int scoreRangeBest = (options.scoreDistanceRange == 0u) ? std::numeric_limits<int>::min() : -(int)options.scoreDistanceRange;

    TIterator it = begin(matches, Standard());
    TIterator itEnd = end(matches, Standard());
    TIterator dit = it;
    TIterator ditBeg = it;
    unsigned disabled = 0;

    // sort
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
#ifdef RAZERS_EXTERNAL_MATCHES
    if (compactMode == COMPACT_FINAL_EXTERNAL)
    {
        typedef Pipe<TMatches, Source<> > TSource;
        typedef LessPairErrors3Way<TFragmentStore, TMatch> TLess;
        typedef Pool<TMatch, SorterSpec<SorterConfigSize<TLess, typename Size<TSource>::Type> > > TSorterPool;

        TLess cmp(store);

        TSource source(matches);
        TSorterPool sorter(cmp);
        sorter << source;
        matches << sorter;

        for (unsigned i = 1; i < length(matches); ++i)
            SEQAN_ASSERT_LEQ(cmp(matches[i - 1], matches[i]), 0);
    }
    else
    {
#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
    std::sort(it, itEnd, LessPairErrors<TFragmentStore, TMatch>(store));
//	sortAlignedReads(threadStore, LessPairScore<TFragmentStore>(mainStore, threadStore));
#ifdef RAZERS_EXTERNAL_MATCHES
}

#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

    for (; it != itEnd; ++it)
    {
        SEQAN_ASSERT_EQ(it->pairMatchId, (it + 1)->pairMatchId);

        // ignore pair alignments if one of the mates is marked as deleted (<=> orientation is '-')
        if (it->orientation == '-' || (it + 1)->orientation == '-')
        {
            ++it;
            continue;
        }

        // std::cerr << *it << std::endl;
        // std::cerr << *(it + 1) << std::endl;
        // SEQAN_ASSERT(it->pairMatchId == (it + 1)->pairMatchId);
        // if (*it).readId == TMatch::INVALID_ID || (*it).pairMatchId == TMatch::INVALID_ID) continue;
        if (matePairId == store.readStore[it->readId].matePairId)
        {
            if (it->pairScore <= scoreDistCutOff)
            {
                typename Iterator<String<unsigned char>, Standard>::Type p = begin(options.errorCutOff, Standard()) + (2 * matePairId);

                int maxErrors = -scoreDistCutOff;
                if (*p > (unsigned)maxErrors)
                    *p = maxErrors;
                ++p;
                if (*p > (unsigned)maxErrors)
                    *p = maxErrors;

                setMaxErrors(filterL, matePairId, maxErrors - 1);
                setMaxErrors(filterR, matePairId, maxErrors - 1);

                while (it != itEnd && matePairId == store.readStore[it->readId].matePairId)
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
                    int maxErrors = -(*it).pairScore;
                    if (options.purgeAmbiguous && (*it).pairScore > scoreRangeBest)
                        maxErrors = 0;

                    setMaxErrors(filterL, matePairId, maxErrors - 1);
                    setMaxErrors(filterR, matePairId, maxErrors - 1);

                    typename Iterator<String<unsigned char>, Standard>::Type p = begin(options.errorCutOff, Standard()) + (2 * matePairId);
                    if (*p > (unsigned)maxErrors)
                        *p = maxErrors;
                    ++p;
                    if (*p > (unsigned)maxErrors)
                        *p = maxErrors;

                    if (maxErrors == 0 && options._debugLevel >= 2)
                        disabled += 1;
                    // std::cerr << "(pair #" << matePairId << " disabled)";

                    if (options.purgeAmbiguous)
                    {
                        if (it->pairScore > scoreRangeBest || compactMode == COMPACT_FINAL || compactMode == COMPACT_FINAL_EXTERNAL)
                        {
                            dit = ditBeg;
                        }
                        else
                        {
                            *dit = *it; ++dit; ++it;
                            *dit = *it; ++dit;
                            continue;
                        }
                    }
                }
#endif
                ++it;
                continue;
            }
        }
        else
        {
            matePairId = store.readStore[(*it).readId].matePairId;
            hitCount = 0;
            if (options.scoreDistanceRange > 0)
                scoreDistCutOff = (*it).pairScore - options.scoreDistanceRange;
            ditBeg = dit;
        }
        *dit = *it; ++dit; ++it;
        *dit = *it; ++dit;
    }
    unsigned origSize = length(matches);
    resize(matches, dit - begin(matches, Standard()));
//	compactAlignedReads(matches);

    options.timeCompactMatches += sysTime() - beginTime;

    if (options._debugLevel >= 2)
    {
        fprintf(stderr, "[%u reads disabled]", disabled);
        unsigned newSize = length(matches);
        fprintf(stderr, "[%u of %u alignments removed]", unsigned(origSize - newSize), unsigned(origSize));
    }
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
    typename TMatches,
    typename TFSSpec,
    typename TFSConfig,
    typename TReadIndex,
    typename TFilterSpec,
    typename TCounts,
    typename TRazerSOptions,
    typename TRazerSMode>
void _mapMatePairReads(
    TMatches & matches,
    FragmentStore<TFSSpec, TFSConfig> & store,
    unsigned                                  contigId,             // ... and its sequence number
    Pattern<TReadIndex, TFilterSpec> & filterPatternL,
    Pattern<TReadIndex, TFilterSpec> & filterPatternR,
    TCounts & cnts,
    char                                      orientation,          // q-gram index of reads
    TRazerSOptions & options,
    TRazerSMode                       const & mode)
{
    typedef FragmentStore<TFSSpec, TFSConfig>               TFragmentStore;
    typedef typename TFragmentStore::TMatePairStore         TMatePairStore;
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    //typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;
    typedef typename Value<TMatePairStore>::Type            TMatePair;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedRead;
    //typedef typename Value<TAlignQualityStore>::Type        TAlignQuality;
    typedef typename Fibre<TReadIndex, FibreText>::Type TReadSet;
    //typedef typename Id<TAlignedRead>::Type                 TId;

    typedef typename TFragmentStore::TContigSeq             TGenome;
    typedef typename Size<TGenome>::Type                    TSize;
    typedef typename Position<TGenome>::Type                TGPos;
    typedef typename MakeSigned_<TGPos>::Type               TSignedGPos;
    typedef typename Infix<TGenome>::Type                   TGenomeInf;

    typedef typename Value<TMatches>::Type TMatch;

    // FILTRATION
    typedef Finder<TGenome, TFilterSpec>                    TFilterFinderL;
    typedef Finder<TGenomeInf, TFilterSpec>                 TFilterFinderR;
    typedef Pattern<TReadIndex, TFilterSpec>                TFilterPattern;

    // MATE-PAIR FILTRATION
    typedef Pair<int64_t, TMatch>                           TDequeueValue;
    typedef Dequeue<TDequeueValue>                          TDequeue;
    typedef typename TDequeue::TIter                        TDequeueIterator;

    // VERIFICATION
    typedef MatchVerifier<
        TFragmentStore,
        TMatches,
        TRazerSOptions,
        TRazerSMode,
        TFilterPattern,
        TCounts
        > TVerifier;

    const unsigned NOT_VERIFIED = 1u << (8 * sizeof(unsigned) - 1);

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
    TGenome & genome = store.contigStore[contigId].seq;
    if (orientation == 'R')
        reverseComplement(genome);

    TReadSet & readSetL = host(host(filterPatternL));
    TReadSet & readSetR = host(host(filterPatternR));
    TVerifier   verifierL(matches, options, filterPatternL, cnts);
    TVerifier   verifierR(matches, options, filterPatternR, cnts);

    verifierL.oneMatchPerBucket = true;
    verifierR.oneMatchPerBucket = true;
    verifierL.m.contigId = contigId;
    verifierR.m.contigId = contigId;

    if (empty(readSetL))
        return;

    // distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
    // distance >= libLen - libErr - 2*parWidth + shapeLen
    TSize readLength = length(readSetL[0]);
    TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(filterPatternL)));
    TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(filterPatternL)));
    TGPos scanShift = (minDistance < 0) ? 0 : minDistance;

    // exit if contig is shorter than library size
    if (length(genome) <= scanShift)
        return;

    TGenomeInf genomeInf = infix(genome, scanShift, length(genome));
    TFilterFinderL filterFinderL(genome, options.repeatLength, 1);
    TFilterFinderR filterFinderR(genomeInf, options.repeatLength, 1);

    TDequeue fifo;                      // stores left-mate potential matches
    String<int64_t> lastPotMatchNo;     // last number of a left-mate potential
    int64_t lastNo = 0;                 // last number over all left-mate pot. matches in the queue
    int64_t firstNo = 0;                // first number over all left-mate pot. match in the queue
    Pair<TGPos> gPair;

    resize(lastPotMatchNo, length(host(filterPatternL)), (int64_t) - 2, Exact());

    TSize gLength = length(genome);

    TMatch mR;
    TDequeueValue fL(-2, mR);   // to supress uninitialized warnings

    // Iterate over all filtration results are returned by SWIFT.
    while (find(filterFinderR, filterPatternR, options.errorRate))
    {
        ++options.countFiltration;

#ifdef RAZERS_DEBUG_MATEPAIRS
        std::cerr << "\nSWIFT\tR\t" << filterPatternR.curSeqNo << "\t" << store.readNameStore[2 * filterPatternR.curSeqNo + 1];
        std::cerr << "\t" << scanShift + beginPosition(filterFinderR) << "\t" << scanShift + endPosition(filterFinderR) << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS

        unsigned matePairId = filterPatternR.curSeqNo;
        TGPos rEndPos = endPosition(filterFinderR) + scanShift;
        TGPos doubleParWidth = 2 * (*filterFinderR.curHit).bucketWidth;

        // (1) Remove out-of-window left mates from fifo.
        while (!empty(fifo) && (TSignedGPos)front(fifo).i2.endPos + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
        {
#ifdef RAZERS_DEBUG_MATEPAIRS
            if (front(fifo).i2.readId > length(store.readNameStore))
                std::cerr << "\nPOP\tL\t" << "[bad read]" << "\t" << front(fifo).i2.beginPos << "\t" << front(fifo).i2.endPos << std::endl;
            else
                std::cerr << "\nPOP\tL\t" << store.readNameStore[front(fifo).i2.readId & ~NOT_VERIFIED] << "\t" << front(fifo).i2.beginPos << "\t" << front(fifo).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
//            std::cerr << "  -Left [" << front(fifo).i2.endPos << "\t" << front(fifo).i2.beginPos << ')' << std::endl;
            popFront(fifo);
            ++firstNo;
        }

        // (2) Add within-window left mates to fifo.
        while (empty(fifo) || (TSignedGPos)back(fifo).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
        {
            if (find(filterFinderL, filterPatternL, options.errorRate))
            {
                ++options.countFiltration;
#ifdef RAZERS_DEBUG_MATEPAIRS
                std::cerr << "\nSWIFT\tL\t" << filterPatternL.curSeqNo << "\t" << store.readNameStore[2 * filterPatternL.curSeqNo];
                std::cerr << "\t" << beginPosition(filterFinderL) << "\t" << endPosition(filterFinderL) << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                gPair = positionRange(filterFinderL);
                if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
                {
                    // link in
                    fL.i1 = lastPotMatchNo[filterPatternL.curSeqNo];
                    lastPotMatchNo[filterPatternL.curSeqNo] = lastNo++;

                    fL.i2.readId = store.matePairStore[filterPatternL.curSeqNo].readId[0] | NOT_VERIFIED;
                    fL.i2.beginPos = beginPosition(filterFinderL);
                    fL.i2.endPos = gPair.i2;

//            std::cerr << "  +Left \t" << firstNo + length(fifo) << ":\t[" << fL.i2.endPos << "\t" << fL.i2.beginPos << ')' << std::endl;
                    pushBack(fifo, fL);
                }
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
        unsigned leftReadId = store.matePairStore[matePairId].readId[0];
        int64_t last = (int64_t) - 1;
        int64_t lastValid = (int64_t) - 1;
        int64_t i;
        for (i = lastPotMatchNo[matePairId]; firstNo <= i; last = i, i = (*it).i1)
        {
//            std::cout<< "\t[" << i << "]" << "\t" << fifo[3].i1 << std::endl;
            it = &value(fifo, i - firstNo);

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
                        std::cerr << "\nVERIFY\tL\t" << matePairId << "\t" << store.readNameStore[2 * matePairId] << "\t" << (TSignedGPos)(*it).i2.beginPos << "\t" << (*it).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        ++options.countVerification;
//                        if (i==0)
//                        std::cout<<"here"<<std::endl;

                        // adjust sink position according to insert size
                        if (!rightVerified)
                            verifierL.sinkPos = (TSignedGPos)endPosition(filterFinderR) - options.libraryLength;

                        if (matchVerify(verifierL, infix(genome, ((*it).i2.beginPos >= 0) ? (TSignedGPos)(*it).i2.beginPos : (TSignedGPos)0, (TSignedGPos)(*it).i2.endPos),
                                        leftReadId, readSetL[matePairId], mode))
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
                        lastPotMatchNo[matePairId] = i;
                    else
                        value(fifo, lastValid - firstNo).i1 = i;
                }
                lastValid = i;

                if (!rightVerified)                                         // here a verfied left match is available
                {
#ifdef RAZERS_DEBUG_MATEPAIRS
                    std::cerr << "\nVERIFY\tR\t" << matePairId << "\t" << store.readNameStore[2 * matePairId + 1] << "\t" << beginPosition(filterFinderR) << "\t" << endPosition(filterFinderR) << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                    ++options.countVerification;
                    if (matchVerify(verifierR, infix(filterFinderR), 2 * matePairId + 1, readSetR[matePairId], mode))
                    {
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "  YES: " << verifierR.m.beginPos << "\t" << verifierR.m.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        rightVerified = true;
                        mR = verifierR.m;
                        // adjust sink position according to insert size
                        verifierL.sinkPos = (TSignedGPos)verifierR.m.endPos - options.libraryLength;
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
/*
                        if (orientation == 'F')
                            std::cout << (int64_t)(*it).i2.beginPos << "\t" << (int64_t)verifierR.m.beginPos;
                        else
                            std::cout << (int64_t)(*it).i2.endPos << "\t" << (int64_t)verifierR.m.endPos;
                        std::cout << '\t' << dist << '\t' << libSizeError << std::endl;
*/
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
                lastPotMatchNo[matePairId] = i;
            else
                value(fifo, lastValid - firstNo).i1 = i;
        }

        // verify right mate, if left mate matches
        if (bestLeftScore != std::numeric_limits<int>::min())
        {
            fL.i2 = (*bestLeft).i2;

            // transform mate readNo to global readNo
            TMatePair & mp     = store.matePairStore[matePairId];
            fL.i2.readId      = mp.readId[0];
            mR.readId         = mp.readId[1];
            mR.orientation    = (orientation == 'F') ? 'R' : 'F';
            fL.i2.orientation = orientation;

            // transform coordinates to the forward strand
            if (orientation == 'F')
            {
// TODO (weese:) Manuel, doesn't this violate the invariant begin<end for MatchRecords?
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
            fL.i2.pairMatchId = mR.pairMatchId = options.nextPairMatchId;
            if (++options.nextPairMatchId == TAlignedRead::INVALID_ID)
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
                appendValue(matches, fL.i2, Generous());
                appendValue(matches, mR, Generous());

#ifdef RAZERS_DEBUG_MATEPAIRS
                std::cerr << "\nHIT\tL\t" << fL.i2.readId << "\t" << store.readNameStore[fL.i2.readId] << "\t" << fL.i2.beginPos << "\t" << fL.i2.endPos << std::endl;
                std::cerr << "\nHIT\tR\t" << mR.readId << "\t" << store.readNameStore[mR.readId] << "\t" << mR.beginPos << "\t" << mR.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS

                if ((int64_t)length(store.alignedReadStore) > options.compactThresh)
                {
                    typename Size<TAlignedReadStore>::Type oldSize = length(store.alignedReadStore);
                    if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
                        maskDuplicates(matches, options, mode);         // overlapping parallelograms cause duplicates
                    compactPairMatches(store, matches, cnts, options, filterPatternL, filterPatternR, COMPACT);

                    if (length(store.alignedReadStore) * 4 > oldSize)                   // the threshold should not be raised
                        options.compactThresh = (int64_t)(options.compactThresh * options.compactMult);
                    //options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed

                    if (options._debugLevel >= 2)
                        std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
                }
            }
            // XXX
            // }
            // XXX
        }
        // XXX
        // }
        // XXX
    }

    if (!unlockAndFreeContig(store, contigId))                      // if the contig is still used
        if (orientation == 'R')
            reverseComplement(genome);
    // we have to restore original orientation
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
int _mapMatePairReads(
    FragmentStore<TFSSpec, TFSConfig> & store,
    TCounts & cnts,
    RazerSCoreOptions<TSpec> & options,
    TShape const & shape,
    RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy>  const & mode,
    TFilterSpec)
{
    typedef FragmentStore<TFSSpec, TFSConfig>           TFragmentStore;
    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;

    typedef typename Value<TReadSeqStore>::Type         TRead;
    typedef StringSet<TRead>                            TReadSet;
#ifndef RAZERS_OPENADDRESSING
    typedef Index<TReadSet, IndexQGram<TShape> >    TIndex;         // q-gram index
#else
    typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >    TIndex;
#endif

    //typedef typename If<
    //    IsSameType<TGapMode, RazerSGapped>,
    //    SwiftSemiGlobal,
    //    SwiftSemiGlobalHamming>::Type           TSwiftSpec;
    typedef Pattern<TIndex, TFilterSpec>                TFilterPattern; // filter

    typedef typename TFragmentStore::TContigSeq TContigSeq;
    typedef typename Position<TContigSeq>::Type TContigPos;
    typedef MatchRecord<TContigPos> TMatchRecord;

//	std::cout << "SA-TYPE:" <<sizeof(typename SAValue<TIndex>::Type)<<std::endl;

    // split mate-pairs over two indices
    TReadSet readSetL, readSetR;
    unsigned pairCount = length(store.matePairStore);
    resize(readSetL, pairCount, Exact());
    resize(readSetR, pairCount, Exact());

    for (unsigned i = 0; i < pairCount; ++i)
    {
        assign(readSetL[i], store.readSeqStore[store.matePairStore[i].readId[0]]);
        assign(readSetR[i], store.readSeqStore[store.matePairStore[i].readId[1]]);
    }
    reverseComplement(readSetR);

    // configure q-gram index
    TIndex filterIndexL(readSetL, shape);
    TIndex filterIndexR(readSetR, shape);
#ifdef RAZERS_OPENADDRESSING
    filterIndexL.alpha = options.loadFactor;
    filterIndexR.alpha = options.loadFactor;
#endif

    cargo(filterIndexL).abundanceCut = options.abundanceCut;
    cargo(filterIndexR).abundanceCut = options.abundanceCut;
    cargo(filterIndexL)._debugLevel = options._debugLevel;
    cargo(filterIndexR)._debugLevel = options._debugLevel;

    // configure Filter
    TFilterPattern filterPatternL(filterIndexL);
    TFilterPattern filterPatternR(filterIndexR);

    // right mate qualities are reversed -> reverse right shape
    reverse(indexShape(filterIndexR));

    // Configure filter pattern
    // (if this is a pigeonhole filter, all sequences must be appended first)
    _applyFilterOptions(filterPatternL, options);
    _applyFilterOptions(filterPatternR, options);
    filterPatternL.params.printDots = false; // only one should print the dots
    filterPatternR.params.printDots = options._debugLevel > 0;

    // clear stats
    options.countFiltration = 0;
    options.countVerification = 0;
    options.timeMapReads = 0;
    options.timeDumpResults = 0;

    options.timeDumpResults = 0;
    SEQAN_PROTIMESTART(find_time);

    // We collect the matches in a more compact data structure than the
    // AlignedReadStoreElement from FragmentStore.
    String<TMatchRecord> matches;

    for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId)
    {
        // lock to prevent releasing and loading the same contig twice
        // (once per _mapSingleReadsToContig call)
        lockContig(store, contigId);

//		std::cout<<"contigLen: "<<length(store.contigStore[contigId].seq)<<std::endl;

        if (options.forward)
            _mapMatePairReads(matches, store, contigId, filterPatternL, filterPatternR, cnts, 'F', options, mode);

        if (options.reverse)
            _mapMatePairReads(matches, store, contigId, filterPatternL, filterPatternR, cnts, 'R', options, mode);

        unlockAndFreeContig(store, contigId);
    }

    double beginCopyTime = sysTime();
    // Final compact matches
    if (IsSameType<TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
        maskDuplicates(matches, options, mode);  // overlapping parallelograms cause duplicates
    compactPairMatches(store, matches, cnts, options, filterPatternL, filterPatternR, COMPACT_FINAL);
    // Write back to store.
    reserve(store.alignedReadStore, length(matches), Exact());
    reserve(store.alignQualityStore, length(matches), Exact());
    typedef typename Iterator<String<TMatchRecord>, Standard>::Type TIterator;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
    typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
    for (TIterator it = begin(matches), itEnd = end(matches); it != itEnd; ++it)
    {
        SEQAN_ASSERT_NEQ(it->orientation, '-');
        SEQAN_ASSERT(!(it->orientation == 'F') || (it->beginPos <= it->endPos));  // implication
        SEQAN_ASSERT(!(it->orientation == 'R') || (it->beginPos >= it->endPos));  // implication
        // if (it->orientation == 'R')
        //     std::swap(it->beginPos, it->endPos);
        appendValue(store.alignedReadStore, TAlignedReadStoreElem(length(store.alignQualityStore), it->readId, it->contigId, it->beginPos, it->endPos));
        back(store.alignedReadStore).pairMatchId = it->pairMatchId;
        appendValue(store.alignQualityStore, TAlignedQualStoreElem(it->pairScore, it->score, -it->score));
    }
    options.timeFsCopy = sysTime() - beginCopyTime;

    // restore original orientation (R-reads are infixes of ConcatDirect StringSet)
    reverseComplement(readSetR);

    options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
    if (options._debugLevel >= 1)
        std::cerr << std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << std::endl;
    if (options._debugLevel >= 1)
    {
        std::cerr << "Time for copying back            \t" << options.timeFsCopy << " seconds" << std::endl;
        std::cerr << std::endl;
        std::cerr << "___FILTRATION_STATS____" << std::endl;
        std::cerr << "Filtration counter:  " << options.countFiltration << std::endl;
        std::cerr << "Verification counter: " << options.countVerification << std::endl;
    }

    return 0;
}

} // End namespace

#endif
