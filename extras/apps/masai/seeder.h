// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
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

#ifndef SEQAN_EXTRAS_MASAI_SEEDER_H_
#define SEQAN_EXTRAS_MASAI_SEEDER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TReadsDelegate, typename THitsDelegate, typename TSpec = MultipleBacktracking>
struct Seeder
{
    TFragmentStore & store;
    TReadsDelegate & readsDelegate;
    THitsDelegate & hitsDelegate;

    TReadsWotd          readsWotd;
    TReadsQGram         readsQGram;

    unsigned            readsCount;
    unsigned            seedsCount;
    unsigned long       hitsCount;

    Seeder(TFragmentStore & store, TReadsDelegate & readsDelegate, THitsDelegate & hitsDelegate) :
        store(store),
        readsDelegate(readsDelegate),
        hitsDelegate(hitsDelegate),
        readsCount(0),
        seedsCount(0),
        hitsCount(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TReadsDelegate, typename THitsDelegate, typename TSpec>
void indexSeedsExact(Seeder<TReadsDelegate, THitsDelegate, TSpec> & seeder,
                     TReadSeqSize seedsLength, TReadSeqSize firstSeed, TReadSeqSize lastSeed)
{
    typedef typename Fibre<TReadsQGram, QGramSA>::Type          TReadsIndexSAFibre;
    typedef typename Fibre<TReadsQGram, QGramDir>::Type         TReadsIndexDirFibre;
    typedef typename Fibre<TReadsQGram, QGramShape>::Type       TReadsIndexShape;
    typedef typename Fibre<TReadsQGram, QGramBucketMap>::Type   TReadsIndexBucketMap;

    typedef typename Value<TReadsIndexDirFibre>::Type           TSize;
    typedef Iterator<TReadSeq, Standard>::Type                  TReadSeqIterator;

    TReadSeqStoreSize readsCount = seeder.readsCount * 2;

    seeder.readsQGram = TReadsQGram(seeder.store.readSeqStore);

    setStepSize(seeder.readsQGram, seedsLength);

    TReadsIndexSAFibre & sa     = indexSA(seeder.readsQGram);
    TReadsIndexDirFibre & dir    = indexDir(seeder.readsQGram);
    TReadsIndexShape & shape  = indexShape(seeder.readsQGram);
    TReadsIndexBucketMap & bucketMap = indexBucketMap(seeder.readsQGram);

    // Resize suffix array and directory.
    resize(sa, (lastSeed - firstSeed) * readsCount, Exact());
    resize(dir, _fullDirLength(seeder.readsQGram), Exact());

    // Clear directory.
    _qgramClearDir(dir, bucketMap);

    // Count qgrams.
    for (TReadSeqStoreSize readId = 0; readId < readsCount; ++readId)
    {
        // Skip disabled reads.
        if (isDisabled(seeder.readsDelegate, readId))
            continue;

        TReadSeq & read          = seeder.store.readSeqStore[readId];
        TReadSeqIterator itText = begin(read, Standard());

        itText += seedsLength * firstSeed;
        for (TSize i = firstSeed; i < lastSeed; ++i)
        {
            ++dir[requestBucket(bucketMap, hash(shape, itText))];
            itText += seedsLength;
        }
    }

    // Compute cumulative sum.
    _qgramCummulativeSum(dir, False());

    // Fill suffix array.
    for (TReadSeqStoreSize readId = 0; readId < readsCount; ++readId)
    {
        // Skip disabled reads.
        if (isDisabled(seeder.readsDelegate, readId))
            continue;

        TReadSeq & read          = seeder.store.readSeqStore[readId];
        TReadSeqIterator itText = begin(read, Standard());

        typename Value<TReadsIndexSAFibre>::Type localPos;
        assignValueI1(localPos, readId);
        assignValueI2(localPos, 0);

        itText += seedsLength * firstSeed;
        for (TSize i = firstSeed; i < lastSeed; ++i)
        {
            assignValueI2(localPos, seedsLength * i);

            sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;

            itText += seedsLength;
        }
    }

    // Refine buckets.
//    _refineQGramIndex(sa, dir, indexText(seeder.readsQGram), weight(shape), seedsLength);
//    _setHost(seeder.readsQGram);
}

// ============================================================================

template <typename TReadsDelegate, typename THitsDelegate, typename TSpec>
void indexSeedsApproximate(Seeder<TReadsDelegate, THitsDelegate, TSpec> & seeder,
                           TReadSeqSize seedsLength, TReadSeqSize firstSeed, TReadSeqSize lastSeed)
{
    typedef typename Fibre<TReadsWotd, FibreSA>::Type           TReadsIndexSAFibre;
    typedef typename Value<TReadsIndexSAFibre>::Type            TReadsIndexSAPos;

//    clear(seeder.readsWotd);
    seeder.readsWotd = TReadsWotd(seeder.store.readSeqStore);
    TReadSeqStoreSize readsCount = seeder.readsCount * 2;

    TReadsIndexSAFibre & sa = indexSA(seeder.readsWotd);

    reserve(sa, (lastSeed - firstSeed) * readsCount, Exact());

    for (TReadSeqStoreSize readId = 0; readId < readsCount; ++readId)
    {
        // Skip disabled reads.
        if (isDisabled(seeder.readsDelegate, readId))
            continue;

        for (TReadSeqSize seed = firstSeed; seed < lastSeed; ++seed)
        {
            TReadsIndexSAPos localPos;
            assignValueI1(localPos, readId);
            assignValueI2(localPos, seed * seedsLength);
            appendValue(sa, localPos, Exact());
        }
    }
}

template <typename TReadsDelegate, typename THitsDelegate, typename TSpec, typename TDepth>
void visitSeedsApproximate(Seeder<TReadsDelegate, THitsDelegate, TSpec> & seeder, TDepth depth)
{
    typedef typename Iterator<TReadsWotd, TopDown<ParentLinks<> > >::Type  TReadsIndexIterator;
    TReadsIndexIterator readsIt(seeder.readsWotd);

    do
    {
        if (repLength(readsIt) >= depth || !goDown(readsIt))
            if (!goRight(readsIt))
                while (goUp(readsIt) && !goRight(readsIt))
                    ;
    }
    while (!isRoot(readsIt));
}

// ============================================================================

template <typename TReadsDelegate, typename THitsDelegate, typename TSpec, typename TGenomeIndex, typename TDistance>
void findSeedsExact(Seeder<TReadsDelegate, THitsDelegate, TSpec> & seeder,
                    TGenomeIndex & genomeIndex, TReadSeqSize seedsLength, TReadSeqSize errorsPerSeed, TDistance)
{
    typedef Backtracking<TDistance, Stretched<> >           TBacktracking;
//	typedef Backtracking<TDistance>                         TBacktracking;
    typedef Finder<TGenomeIndex, TBacktracking>             TFinder;
    typedef Pattern<TReadsQGram, TBacktracking>             TPattern;

    TFinder finder(genomeIndex);
    TPattern pattern(seeder.readsQGram, seedsLength);

//    seeder.hitsCount = count_exact_iterative(finder, pattern, 0);

    while (find(finder, pattern, errorsPerSeed))
    {
        // Skip disabled reads.
        if (isDisabled(seeder.readsDelegate, position(pattern).i1))
            continue;

        ++seeder.hitsCount;

        onSeedHit(seeder.hitsDelegate,
                  position(finder).i1, beginPosition(finder).i2,
                  position(pattern).i1, beginPosition(pattern).i2,
                  0);
    }
}

template <typename TReadsDelegate, typename THitsDelegate, typename TSpec, typename TGenomeIndex, typename TDistance>
void findSeedsApproximate(Seeder<TReadsDelegate, THitsDelegate, TSpec> & seeder,
                          TGenomeIndex & genomeIndex, TReadSeqSize seedsLength, TReadSeqSize errorsPerSeed, TDistance)
{
    typedef Backtracking<TDistance>                         TBacktracking;
    typedef Finder<TGenomeIndex, TBacktracking>             TFinder;
    typedef Pattern<TReadsWotd, TBacktracking>              TPattern;

    TFinder finder(genomeIndex);
    TPattern pattern(seeder.readsWotd, seedsLength);

//    seeder.hitsCount = count_iterative(finder, pattern, errorsPerSeed);

    while (find(finder, pattern, errorsPerSeed))
    {
        // Skip disabled reads.
        if (isDisabled(seeder.readsDelegate, position(pattern).i1))
            continue;

        ++seeder.hitsCount;

        onSeedHit(seeder.hitsDelegate,
                  position(finder).i1, beginPosition(finder).i2,
                  position(pattern).i1, beginPosition(pattern).i2,
                  pattern.prefix_aligner.errors);
    }
}

// ============================================================================

template <typename TReadsDelegate, typename THitsDelegate, typename TSpec, typename TGenomeIndex, typename TDistance>
void find(Seeder<TReadsDelegate, THitsDelegate, TSpec> & seeder, TGenomeIndex & genomeIndex,
          TReadSeqSize seedsLength, TReadSeqSize errorsPerSeed, TReadSeqSize firstSeed, TReadSeqSize lastSeed, TDistance)
{
    if (errorsPerSeed > 0)
    {
        indexSeedsApproximate(seeder, seedsLength, firstSeed, lastSeed);
        findSeedsApproximate(seeder, genomeIndex, seedsLength, errorsPerSeed, TDistance());
    }
    else
    {
        indexSeedsExact(seeder, seedsLength, firstSeed, lastSeed);
        findSeedsExact(seeder, genomeIndex, seedsLength, errorsPerSeed, TDistance());
    }

//    indexSeedsExact(seeder, seedsLength, firstSeed, lastSeed);
//    findSeedsExact(seeder, genomeIndex, seedsLength, errorsPerSeed, TDistance());
}

// ============================================================================

template <typename TReadsDelegate, typename THitsDelegate, typename TGenomeIndex, typename TDistance>
void findSeedsExact(Seeder<TReadsDelegate, THitsDelegate, SingleBacktracking> & seeder, TGenomeIndex & genomeIndex,
                    TReadSeqSize seedsLength, TReadSeqSize firstSeed, TReadSeqSize lastSeed, TDistance)
{
    typedef Finder<TGenomeIndex, FinderSTree>       TFinder;
    typedef Pattern<TReadSeq>                       TPattern;

    TFinder finder(genomeIndex);

    TReadSeqStoreSize readsCount = seeder.readsCount * 2;

    for (TReadSeqStoreSize readId = 0; readId < readsCount; ++readId)
    {
        // Skip disabled reads.
        if (isDisabled(seeder.readsDelegate, readId))
            continue;

        TReadSeq & read = seeder.store.readSeqStore[readId];

        for (TReadSeqSize seed = firstSeed; seed < lastSeed; ++seed)
        {
            TPattern pattern(infix(read, seedsLength * seed, seedsLength * (seed + 1)));

            clear(finder);

            while (find(finder, pattern))
            {
                ++seeder.hitsCount;

                onSeedHit(seeder.hitsDelegate,
                          position(finder).i1, beginPosition(finder).i2,
                          readId, seedsLength * seed,
                          0);
            }
        }
    }
}

template <typename TReadsDelegate, typename THitsDelegate, typename TGenomeIndex, typename TDistance>
void findSeedsApproximate(Seeder<TReadsDelegate, THitsDelegate, SingleBacktracking> & seeder, TGenomeIndex & genomeIndex,
                          TReadSeqSize seedsLength, TReadSeqSize errorsPerSeed, TReadSeqSize firstSeed, TReadSeqSize lastSeed,
                          TDistance)
{
    typedef Backtracking<TDistance>                 TBacktracking;
    typedef Finder<TGenomeIndex, TBacktracking>     TFinder;
    typedef Pattern<TReadSeq, TBacktracking>        TPattern;

    TFinder finder(genomeIndex);

    TReadSeqStoreSize readsCount = seeder.readsCount * 2;

    for (TReadSeqStoreSize readId = 0; readId < readsCount; ++readId)
    {
        // Skip disabled reads.
        if (isDisabled(seeder.readsDelegate, readId))
            continue;

        TReadSeq & read = seeder.store.readSeqStore[readId];

        for (TReadSeqSize seed = firstSeed; seed < lastSeed; ++seed)
        {
            TPattern pattern(infix(read, seedsLength * seed, seedsLength * (seed + 1)));

            clear(finder);

            while (find(finder, pattern, errorsPerSeed))
            {
                ++seeder.hitsCount;

                onSeedHit(seeder.hitsDelegate,
                          position(finder).i1, beginPosition(finder).i2,
                          readId, seedsLength * seed,
                          pattern.prefix_aligner.errors);
            }
        }
    }
}

// ============================================================================

template <typename TReadsDelegate, typename THitsDelegate, typename TGenomeIndex, typename TDistance>
void find(Seeder<TReadsDelegate, THitsDelegate, SingleBacktracking> & seeder, TGenomeIndex & genomeIndex,
          TReadSeqSize seedsLength, TReadSeqSize errorsPerSeed, TReadSeqSize firstSeed, TReadSeqSize lastSeed, TDistance)
{
    if (errorsPerSeed > 0)
        findSeedsApproximate(seeder, genomeIndex, seedsLength, errorsPerSeed, firstSeed, lastSeed, TDistance());
    else
        findSeedsExact(seeder, genomeIndex, seedsLength, firstSeed, lastSeed, TDistance());
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_SEEDER_H_
