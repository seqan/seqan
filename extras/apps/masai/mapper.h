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
// This file contains the Mapper class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_MAPPER_H_
#define SEQAN_EXTRAS_MASAI_MAPPER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find.h>

#include "tags.h"
#include "store.h"
#include "index.h"
#include "indexer.h"
#include "manager.h"
#include "writer.h"
#include "extender.h"
#include "seeder.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

// TODO(esiragusa):Use typename TGenomeIndex instead of TSpec.
template <typename TSpec = void>
struct Mapper
{
    typedef Indexer<TSpec>  TIndexer;

    TFragmentStore      store;
    TIndexer            indexer;

    unsigned            readsCount;
    TReadSeqSize        seedLength;

    bool                writeCigar;
    bool                verifyHits;
    bool                dumpResults;

    // TODO(esiragusa):Remove writeCigar from Mapper members.
    Mapper(TReadSeqSize seedLength, bool writeCigar = true, bool verifyHits = true, bool dumpResults = true) :
        indexer(store),
        readsCount(0),
        seedLength(seedLength),
        writeCigar(writeCigar),
        verifyHits(verifyHits),
        dumpResults(dumpResults)
    {}
};

// ----------------------------------------------------------------------------
// Class Seeding
// ----------------------------------------------------------------------------

// TODO(esiragusa): refactor class Seeding.
template <typename TErrors = unsigned char, typename TSpec = void>
struct Seeding
{
    typedef Pair<TReadSeqSize, TErrors> TSeed;
    typedef String<TSeed>               TSeeds;

    TSeeds          seeds;
    TReadSeqSize    seedLength;

    Seeding(TReadSeqSize readLength, TErrors errors, TReadSeqSize seedLength) :
        seedLength(seedLength)
    {
        _computeSeeds(*this, readLength, errors);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TErrors, typename TSpec>
void _computeSeeds(Seeding<TErrors, TSpec> & seeding, TReadSeqSize readLength, TErrors errors)
{
    typedef Seeding<TErrors, TSpec>   TSeeding;
    typedef typename TSeeding::TSeed  TSeed;

    clear(seeding.seeds);

    TReadSeqSize seedsCount = readLength / seeding.seedLength;
    TErrors errorsPerSeed = errors / seedsCount;

    TReadSeqSize firstSeeds = (errors % seedsCount) + 1;

    std::cout << "Seeds:\t\t\t\t";

    if (errorsPerSeed > 0)
    {
        // Remaining seeds get errorsPerSeed - 1 errors.
        for (unsigned seed = 0; seed < seedsCount - firstSeeds; ++seed)
        {
            std::cout << "(" << seeding.seedLength << "," << (unsigned)(errorsPerSeed - 1) << ") ";
            appendValue(seeding.seeds, TSeed(seeding.seedLength, errorsPerSeed - 1));
        }
    }

    // First seeds get errorsPerSeed errors.
    for (unsigned seed = seedsCount - firstSeeds; seed < seedsCount; ++seed)
    {
        std::cout << "(" << seeding.seedLength << "," << (unsigned)errorsPerSeed << ") ";
        appendValue(seeding.seeds, TSeed(seeding.seedLength, errorsPerSeed));
    }

    std::cout << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadReads()                                                [Mapper]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString>
bool loadReads(Mapper<TSpec> & mapper, TString const & readsFile)
{
    if (!loadReads(mapper.store, readsFile))
        return false;

    mapper.readsCount = length(mapper.store.readSeqStore);

    _loadReadsRC(mapper);

    return true;
}

template <typename TSpec>
bool _loadReadsRC(Mapper<TSpec> & mapper)
{
    reserve(mapper.store.readSeqStore, mapper.readsCount * 2, Exact());

    for (TReadSeqStoreSize readId = 0; readId < mapper.readsCount; ++readId)
    {
        TReadSeq const & read = mapper.store.readSeqStore[readId];
        appendValue(mapper.store.readSeqStore, read);
        reverseComplement(back(mapper.store.readSeqStore));
    }

    return true;
}

template <typename TSpec>
TReadSeqSize _readsLength(Mapper<TSpec> & mapper)
{
    return length(mapper.store.readSeqStore[0]);
}

// ----------------------------------------------------------------------------
// Function mapReads()                                                 [Mapper]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString, typename TErrors, typename TDistance, typename TStrategy, typename TMultiple>
bool mapReads(Mapper<TSpec> & mapper,
              TString const & mappedReadsFile,
              TErrors errors,
              TDistance const & /*tag*/,
              TStrategy const & /*tag*/,
              TMultiple const & /*tag*/,
              Sam const & /*tag*/)
{
    typedef Match<>                                         TMatch;
    typedef String<char, External<> >                       TStream;
    typedef MatchWriter<TStream, TDistance, Sam>            TMatchWriter;
    typedef MatchManager<TMatch, TMatchWriter, TStrategy>   TMatchManager;
    typedef Extender<TMatchManager, TDistance>              TExtender;
    typedef Seeder<TMatchManager, TExtender, TMultiple>     TSeeder;

    TStream file;
    if (mapper.dumpResults)
        open(file, toCString(mappedReadsFile), OPEN_RDWR | OPEN_CREATE);

    TMatchWriter writer(file, mapper.store, mapper.readsCount, mapper.dumpResults);
    TMatchManager manager(writer, mapper.readsCount);
    TExtender extender(mapper.store, manager, mapper.readsCount, mapper.seedLength, mapper.verifyHits);
    TSeeder seeder(mapper.store, manager, extender);
    seeder.readsCount = mapper.readsCount;

    // TODO(esiragusa):Remove writeCigar from mapper members.
    writer.writeCigar = mapper.writeCigar;

    _mapReads(mapper, seeder, manager, errors, TStrategy());

    return true;
}

template <typename TSpec, typename TString, typename TErrors, typename TDistance, typename TStrategy, typename TMultiple>
bool mapReads(Mapper<TSpec> & mapper,
              TString const & mappedReadsFile,
              TErrors errors,
              TDistance const & /*tag*/,
              TStrategy const & /*tag*/,
              TMultiple const & /*tag*/,
              Raw const & /*tag*/)
{
    typedef Match<>                                         TMatch;
    typedef String<TMatch, External<> >                     TStream;
    typedef MatchWriter<TStream, TDistance, Raw>            TMatchWriter;
    typedef MatchManager<TMatch, TMatchWriter, TStrategy>   TMatchManager;
    typedef Extender<TMatchManager, TDistance>              TExtender;
    typedef Seeder<TMatchManager, TExtender, TMultiple>     TSeeder;

    TStream file;
    if (mapper.dumpResults)
        open(file, toCString(mappedReadsFile), OPEN_RDWR | OPEN_CREATE);

    TMatchWriter writer(file, mapper.store, mapper.readsCount, mapper.dumpResults);
    TMatchManager manager(writer, mapper.readsCount);
    TExtender extender(mapper.store, manager, mapper.readsCount, mapper.seedLength, mapper.verifyHits);
    TSeeder seeder(mapper.store, manager, extender);
    seeder.readsCount = mapper.readsCount;

    _mapReads(mapper, seeder, manager, errors, TStrategy());

    return true;
}

// ============================================================================

template <typename TSpec, typename TSeeder, typename TMatches, typename TErrors>
bool _mapReads(Mapper<TSpec> & mapper, TSeeder & seeder, TMatches & matches, TErrors errors, All const & /*tag*/)
{
    typedef Seeding<>                           TSeeding;
    typedef TSeeding::TSeeds                    TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TSeeding seeding(_readsLength(mapper), errors, mapper.seedLength);

    seeder.hitsDelegate.minErrorsPerRead = 0;
    seeder.hitsDelegate.maxErrorsPerRead = errors;

    unsigned position = 0;
    TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

    std::cout << "Errors:\t\t\t\t" << errors << std::endl;

    for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
    {
        find(seeder, mapper.indexer.genomeIndex,
             getValueI1(*seedsIt), getValueI2(*seedsIt),
             position, position + 1,
             HammingDistance());

        // TODO(esiragusa):Compute minErrorsPerRead from seeds.
//        seeder.hitsDelegate.minErrorsPerRead += getValueI2(*seedsIt);
        seeder.hitsDelegate.minErrorsPerRead++;

        ++position;
    }

    std::cout << "Hits:\t\t\t\t" << seeder.hitsCount << std::endl;
    std::cout << "Matches:\t\t\t" << matches.matchesCount << std::endl;

    return true;
}

template <typename TSpec, typename TSeeder, typename TMatches, typename TErrors>
bool _mapReads(Mapper<TSpec> & mapper, TSeeder & seeder, TMatches & matches, TErrors errors, AllBest const & /*tag*/)
{
    _mapReadsByStratum(mapper, seeder, matches, std::min(errors, (TErrors)1), AllBest());

    if (errors > 1)
    {
        matches.errors = 2;
        _mapReadsBySeed(mapper, seeder, matches, errors, AllBest());
    }

    return true;
}

template <typename TSpec, typename TSeeder, typename TMatches, typename TErrors>
bool _mapReadsBySeed(Mapper<TSpec> & mapper, TSeeder & seeder, TMatches & matches, TErrors errors, AllBest const & /*tag*/)
{
    typedef Seeding<>                           TSeeding;
    typedef TSeeding::TSeeds                    TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TSeeding seeding(_readsLength(mapper), errors, mapper.seedLength);

    seeder.hitsDelegate.minErrorsPerRead = 0;
    seeder.hitsDelegate.maxErrorsPerRead = errors;

    unsigned position = 0;
    TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

    for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
    {
        std::cout << "Errors:\t\t\t\t" << seeder.hitsDelegate.minErrorsPerRead << std::endl;

        find(seeder, mapper.indexer.genomeIndex,
             getValueI1(*seedsIt), getValueI2(*seedsIt),
             position, position + 1,
             HammingDistance());

        ++position;

        // TODO(esiragusa):Compute minErrorsPerRead from seeds.
//        seeder.hitsDelegate.minErrorsPerRead += getValueI2(*seedsIt);
        seeder.hitsDelegate.minErrorsPerRead++;

        std::cout << "Hits:\t\t\t\t" << seeder.hitsCount << std::endl;
        std::cout << "Matches:\t\t\t" << matches.matchesCount << std::endl;

//        matches.errors += getValueI2(*seedsIt);
//        raiseErrorThreshold(matches);
        matches.errors = std::max((TErrors)matches.errors, (TErrors)position);
    }

    return true;
}

template <typename TSpec, typename TSeeder, typename TMatches, typename TErrors>
bool _mapReadsByStratum(Mapper<TSpec> & mapper, TSeeder & seeder, TMatches & matches, TErrors errors, AllBest const & /*tag*/)
{
    typedef Seeding<>                           TSeeding;
    typedef TSeeding::TSeeds                    TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    seeder.hitsDelegate.minErrorsPerRead = 0;
    seeder.hitsDelegate.maxErrorsPerRead = 0;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TReadSeqSize readsLength_ = _readsLength(mapper);
        TReadSeqSize seedLength_ = std::max(mapper.seedLength, readsLength_ / (errors_ + 1));

        TSeeding seeding(readsLength_, errors_, seedLength_);

        unsigned position = 0;
        TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

        for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
        {
            find(seeder, mapper.indexer.genomeIndex,
                 getValueI1(*seedsIt), getValueI2(*seedsIt),
                 position, position + 1,
                 HammingDistance());

            ++position;
        }

        seeder.hitsDelegate.minErrorsPerRead++;
        seeder.hitsDelegate.maxErrorsPerRead++;

        raiseErrorThreshold(matches);

        std::cout << "Hits:\t\t\t\t" << seeder.hitsCount << std::endl;
        std::cout << "Matches:\t\t\t" << matches.matchesCount << std::endl;
    }

    return true;
}

template <typename TSpec, typename TSeeder, typename TMatches, typename TErrors>
bool _mapReads(Mapper<TSpec> & mapper, TSeeder & seeder, TMatches & matches, TErrors errors, AnyBest const & /*tag*/)
{
    typedef Seeding<>                           TSeeding;
    typedef TSeeding::TSeeds                    TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TReadSeqSize readsLength = _readsLength(mapper);
    TReadSeqSize seedLength  = mapper.seedLength;
    TReadSeqSize seedCount   = readsLength / seedLength;
    TErrors stratumDelta     = seedCount - 1;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TErrors seedErrors_ = errors_ / seedCount;
        TReadSeqSize seed = errors_ % seedCount;

        std::cout << "Seed:\t\t\t\t(" << seedLength << "," << seedErrors_ << ")" << std::endl;

        seeder.hitsDelegate.minErrorsPerRead = seed;
        seeder.hitsDelegate.maxErrorsPerRead = std::min(errors, errors_ + stratumDelta);

        find(seeder, mapper.indexer.genomeIndex, seedLength, seedErrors_, seed, seed + 1, HammingDistance());

        std::cout << "Hits:\t\t\t\t" << seeder.hitsCount << std::endl;
        std::cout << "Matches:\t\t\t" << matches.matchesCount << std::endl;

        raiseErrorThreshold(matches);
    }

    return true;
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_MAPPER_H_
