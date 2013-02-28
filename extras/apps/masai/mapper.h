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

#include "tags.h"
#include "store.h"
#include "seeder.h"
#include "extender.h"
#include "manager.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadMapperConfig
// ----------------------------------------------------------------------------

template <typename TDistance_       = EditDistance,
          typename TStrategy_       = AnyBest,
          typename TBacktracking_   = MultipleBacktracking>
struct ReadMapperConfig
{
	typedef TDistance_          TDistance;
	typedef TStrategy_          TStrategy;
    typedef TBacktracking_      TBacktracking;
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec = void, typename TConfig = ReadMapperConfig<> >
struct Mapper
{
    typedef MatchManager<TDelegate, typename TConfig::TStrategy>            TManager;
    typedef Extender<TManager, typename TConfig::TDistance>                 TExtender;
    typedef Seeder<TManager, TExtender, typename TConfig::TBacktracking>    TSeeder;

    Holder<TReads>      reads;
    TDelegate           & delegate;
    TManager            _manager;
    TExtender           _extender;
    TSeeder             _seeder;
    TReadSeqSize        _seedLength;

    Mapper(TReads & reads, TDelegate & delegate, bool disableExtender = false) :
        reads(reads),
        delegate(delegate),
        _manager(delegate, reads.readsCount),
        _extender(value(reads._store), _manager, reads.readsCount, disableExtender),
        _seeder(value(reads._store), _manager, _extender),
        _seedLength(0)
    {
        _seeder.readsCount = reads.readsCount;
    }
};

// ----------------------------------------------------------------------------
// Class Seeding_
// ----------------------------------------------------------------------------

// TODO(esiragusa): Remove class Seeding_
template <typename TErrors = unsigned char, typename TSpec = void>
struct Seeding_
{
    typedef Pair<TReadSeqSize, TErrors> TSeed;
    typedef String<TSeed>               TSeeds;

    TSeeds          seeds;
    TReadSeqSize    seedLength;

    Seeding_(TReadSeqSize readLength, TErrors errors, TReadSeqSize seedLength) :
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

// ----------------------------------------------------------------------------
// Function _computeSeeds()                                          [Seeding_]
// ----------------------------------------------------------------------------

template <typename TErrors, typename TSpec>
void _computeSeeds(Seeding_<TErrors, TSpec> & seeding, TReadSeqSize readLength, TErrors errors)
{
    typedef Seeding_<TErrors, TSpec>   TSeeding_;
    typedef typename TSeeding_::TSeed  TSeed;

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
// Function setSeedLength()                                            [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TSize>
void setSeedLength(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TSize seedLength)
{
    mapper._seedLength = seedLength;
    mapper._extender.seedLength = seedLength;
}

// ----------------------------------------------------------------------------
// Function setReads()                                                 [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
void setReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TReads & reads)
{
    setValue(mapper.reads, reads);

    mapper._seeder.readsCount = reads.readsCount;
    mapper._extender.readsCount = reads.readsCount;
    mapper._manager.readsCount = reads.readsCount;

//    setReads(mapper._seeder, reads);
//    setReads(mapper._extender, reads);
//    setReads(mapper._manager, reads);
}

// ----------------------------------------------------------------------------
// Function getReads()                                                 [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
inline typename Reference<TReads>::Type
getReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper)
{
    return value(mapper.reads);
}

// ----------------------------------------------------------------------------
// Function unmappedReads()                                            [Mapper]
// ----------------------------------------------------------------------------

//template <typename TReads, typename TDelegate, typename TSpec, typename TConfig>
//void unmappedReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper)
//{
//}

// ----------------------------------------------------------------------------
// Function mapReads()                                                 [Mapper]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors)
{
//    mapper._extender.genome = genomeIndex.genome;
    _mapReads(mapper, genomeIndex, errors, typename TConfig::TStrategy());
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                           [Mapper<All>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               All const & /*tag*/)
{
    typedef Seeding_<>                          TSeeding_;
    typedef TSeeding_::TSeeds                   TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TSeeding_ seeding(avgSeqLength(getReads(mapper)), errors, mapper._seedLength);

    mapper._extender.minErrorsPerRead = 0;
    mapper._extender.maxErrorsPerRead = errors;

    unsigned position = 0;
    TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

    std::cout << "Errors:\t\t\t\t" << errors << std::endl;

    for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
    {
        find(mapper._seeder, genomeIndex.index,
             getValueI1(*seedsIt), getValueI2(*seedsIt),
             position, position + 1,
             HammingDistance());

        // TODO(esiragusa):Compute minErrorsPerRead from seeds.
//        mapper._extender.minErrorsPerRead += getValueI2(*seedsIt);
        mapper._extender.minErrorsPerRead++;

        ++position;
    }

    std::cout << "Hits:\t\t\t\t" << mapper._seeder.hitsCount << std::endl;
    std::cout << "Matches:\t\t\t" << mapper._manager.matchesCount << std::endl;
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                       [Mapper<AllBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               AllBest const & /*tag*/)
{
    _mapReadsByStratum(mapper, genomeIndex, std::min(errors, (TErrors)1), AllBest());

    if (errors > 1)
    {
        mapper._manager.errors = 2;
        _mapReadsBySeed(mapper, genomeIndex, errors, AllBest());
    }
}

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReadsBySeed(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
                     AllBest const & /*tag*/)
{
    typedef Seeding_<>                          TSeeding_;
    typedef TSeeding_::TSeeds                   TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TSeeding_ seeding(avgSeqLength(getReads(mapper)), errors, mapper._seedLength);

    mapper._extender.minErrorsPerRead = 0;
    mapper._extender.maxErrorsPerRead = errors;

    unsigned position = 0;
    TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

    for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
    {
        std::cout << "Errors:\t\t\t\t" << mapper._extender.minErrorsPerRead << std::endl;

        find(mapper._seeder, genomeIndex.index,
             getValueI1(*seedsIt), getValueI2(*seedsIt),
             position, position + 1,
             HammingDistance());

        ++position;

        // TODO(esiragusa):Compute minErrorsPerRead from seeds.
//        mapper._extender.minErrorsPerRead += getValueI2(*seedsIt);
        mapper._extender.minErrorsPerRead++;

        std::cout << "Hits:\t\t\t\t" << mapper._seeder.hitsCount << std::endl;
        std::cout << "Matches:\t\t\t" << mapper._manager.matchesCount << std::endl;

//        mapper._manager.errors += getValueI2(*seedsIt);
//        raiseErrorThreshold(mapper._manager);
        mapper._manager.errors = std::max((TErrors)(mapper._manager.errors), (TErrors)position);
    }
}

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReadsByStratum(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
                        AllBest const & /*tag*/)
{
    typedef Seeding_<>                          TSeeding_;
    typedef TSeeding_::TSeeds                   TSeeds;
    typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    mapper._extender.minErrorsPerRead = 0;
    mapper._extender.maxErrorsPerRead = 0;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TReadSeqSize readsLength_ = avgSeqLength(getReads(mapper));
        TReadSeqSize seedLength_ = std::max(mapper._seedLength, readsLength_ / (errors_ + 1));

        TSeeding_ seeding(readsLength_, errors_, seedLength_);

        unsigned position = 0;
        TSeedsIterator seedsEnd = end(seeding.seeds, Standard());

        for (TSeedsIterator seedsIt = begin(seeding.seeds, Standard()); seedsIt != seedsEnd; ++seedsIt)
        {
            find(mapper._seeder, genomeIndex.index,
                 getValueI1(*seedsIt), getValueI2(*seedsIt),
                 position, position + 1,
                 HammingDistance());

            ++position;
        }

        mapper._extender.minErrorsPerRead++;
        mapper._extender.maxErrorsPerRead++;

        raiseErrorThreshold(mapper._manager);

        std::cout << "Hits:\t\t\t\t" << mapper._seeder.hitsCount << std::endl;
        std::cout << "Matches:\t\t\t" << mapper._manager.matchesCount << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                         [Mapper<KBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               KBest const & /*tag*/)
{
    //typedef Seeding_<>                          TSeeding_;
    //typedef TSeeding_::TSeeds                   TSeeds;
    //typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TReadSeqSize readsLength = avgSeqLength(getReads(mapper));
    TReadSeqSize seedLength  = mapper._seedLength;
    TReadSeqSize seedCount   = readsLength / seedLength;
    TErrors stratumDelta     = seedCount - 1;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TErrors seedErrors_ = errors_ / seedCount;
        TReadSeqSize seed = errors_ % seedCount;

        std::cout << "Seed:\t\t\t\t(" << seedLength << "," << seedErrors_ << ")" << std::endl;

        mapper._extender.minErrorsPerRead = seed;
        mapper._extender.maxErrorsPerRead = std::min(errors, errors_ + stratumDelta);

        find(mapper._seeder, genomeIndex.index, seedLength, seedErrors_, seed, seed + 1, HammingDistance());

        std::cout << "Hits:\t\t\t\t" << mapper._seeder.hitsCount << std::endl;
        std::cout << "Matches:\t\t\t" << mapper._manager.matchesCount << std::endl;

        raiseErrorThreshold(mapper._manager);
    }
}

// ----------------------------------------------------------------------------
// Function _mapReads()                                       [Mapper<AnyBest>]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TConfig, typename TGenomeIndex, typename TErrors>
void _mapReads(Mapper<TReads, TDelegate, TSpec, TConfig> & mapper, TGenomeIndex & genomeIndex, TErrors errors,
               AnyBest const & /*tag*/)
{
    //typedef Seeding_<>                          TSeeding_;
    //typedef TSeeding_::TSeeds                   TSeeds;
    //typedef Iterator<TSeeds, Standard>::Type    TSeedsIterator;

    TReadSeqSize readsLength = avgSeqLength(getReads(mapper));
    TReadSeqSize seedLength  = mapper._seedLength;
    TReadSeqSize seedCount   = readsLength / seedLength;
    TErrors stratumDelta     = seedCount - 1;

    for (TErrors errors_ = 0; errors_ <= errors; ++errors_)
    {
        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;

        TErrors seedErrors_ = errors_ / seedCount;
        TReadSeqSize seed = errors_ % seedCount;

        std::cout << "Seed:\t\t\t\t(" << seedLength << "," << seedErrors_ << ")" << std::endl;

        mapper._extender.minErrorsPerRead = seed;
        mapper._extender.maxErrorsPerRead = std::min(errors, errors_ + stratumDelta);

//        if (mapper._extender.maxErrorsPerRead == errors)
//            mapper._extender.maxErrorsPerRead = errorsLossy;

        find(mapper._seeder, genomeIndex.index, seedLength, seedErrors_, seed, seed + 1, HammingDistance());

        std::cout << "Hits:\t\t\t\t" << mapper._seeder.hitsCount << std::endl;
        std::cout << "Matches:\t\t\t" << mapper._manager.matchesCount << std::endl;

        raiseErrorThreshold(mapper._manager);
    }

//    for (TErrors errors_ = errors + 1; errors_ <= errorsLossy; ++errors_)
//    {
//        std::cout << "Errors:\t\t\t\t" << errors_ << std::endl;
//        std::cout << "Matches:\t\t\t" << mapper._manager.matchesCount << std::endl;
//
//        raiseErrorThreshold(mapper._manager);
//    }
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_MAPPER_H_
