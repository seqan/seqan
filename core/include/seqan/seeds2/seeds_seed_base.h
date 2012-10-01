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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// The class Seed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_BASE_H_
#define SEQAN_SEEDS_SEEDS_SEED_BASE_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// Mixin member for mixing in scores into seeds.
template <typename TScore>
struct ScoreMixin_
{
    TScore _score;
};

/**
.Tag.Seed Specs
..summary:Specialization tags for @Class.Seed@.
..cat:Seed Handling
..tag.Simple:Simple seed that only stores coordinates of the begin end end point and the diagonals.
..tag.Chained:Seed that stores the dot-plot diagonals (i.e. matches) that the seed consists of.
..see:Spec.SimpleSeed
..see:Spec.ChainedSeed
..include:seqan/seeds.h
*/

// Default configuration for seeds without score.
struct DefaultSeedConfig
{
    typedef size_t TPosition;
    typedef size_t TSize;
    typedef MakeSigned_<size_t>::Type TDiagonal;
    typedef False THasScore;
    typedef Nothing TScoreValue;
    typedef Nothing TScoreMixin;
};

// Default configuration for seeds with score.
struct DefaultSeedConfigScore
{
    typedef size_t TPosition;
    typedef size_t TSize;
    typedef MakeSigned_<size_t>::Type TDiagonal;
    typedef True THasScore;
    typedef int TScoreValue;
    typedef ScoreMixin_<int> TScoreMixin;
};

/**
.Class.Seed:
..summary:Describe a seed.
..cat:Seed Handling
..signature:Seed<TSpec, TConfig>
..param.TSpec:The seed specialization type.
..param.TConfig:The configuration object to use for this seed.
..include:seqan/seeds.h
 */
template <typename TSpec, typename TConfig = DefaultSeedConfig>
class Seed;

// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.Position.param.T:type:Class.Seed
///.Metafunction.Position.class:Class.Seed
template <typename TSpec, typename TConfig>
struct Position<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TPosition Type;
};

template <typename TSpec, typename TConfig>
struct Position<Seed<TSpec, TConfig> const>
        : Position<Seed<TSpec, TConfig> > {};

///Metafunction.Size.param.T:type:Class.Seed
///Metafunction.Size.class:Class.Seed
template <typename TSpec, typename TConfig>
struct Size<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TSize Type;
};

template <typename TSpec, typename TConfig>
struct Size<Seed<TSpec, TConfig> const>
        : Size<Seed<TSpec, TConfig> > {};

/**
.Metafunction.Diagonal:
..cat:Seed Handling
..summary:Returns type of the value for the diagonal of a seed.
..signature:Diagonal<T>::Type
..class:Class.Seed
..param.T:Type of the seed to retrieve the diagonal for.
...type:Class.Seed
..include:seqan/seeds2.h
 */
template <typename T>
struct Diagonal;

template <typename TSpec, typename TConfig>
struct Diagonal<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TDiagonal Type;
};

template <typename TSpec, typename TConfig>
struct Diagonal<Seed<TSpec, TConfig> const>
        : Diagonal<Seed<TSpec, TConfig> > {};

/**
.Metafunction.HasScore:
..cat:Seed Handling
..summary:Returns True if the seed stores a score, False otherwise.
..signature:HasScore<T>::Type
..class:Class.Seed
..param.T:Type of the seed to retrieve whether it has a score for.
...type:Class.Seed
..include:seqan/seeds2.h
 */
template <typename T>
struct HasScore;

template <typename TSpec, typename TConfig>
struct HasScore<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::THasScore Type;
};

template <typename TSpec, typename TConfig>
struct HasScore<Seed<TSpec, TConfig> const>
        : HasScore<Seed<TSpec, TConfig> > {};

/**
.Metafunction.SeedScore:
..cat:Seed Handling
..summary:Returns type of the value for the score of a seed.
..signature:SeedScore<T>::Type
..class:Class.Seed
..param.T:Type of the seed to retrieve the score for.
...type:Class.Seed
..include:seqan/seeds2.h
 */
template <typename T>
struct SeedScore;

template <typename TSpec, typename TConfig>
struct SeedScore<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TScoreValue Type;
};

template <typename TSpec, typename TConfig>
struct SeedScore<Seed<TSpec, TConfig> const>
        : SeedScore<Seed<TSpec, TConfig> > {};

// ===========================================================================
// Functions
// ===========================================================================

// TODO(holtgrew): COULD introduce {get,set}{Begin,End}(dim, value), but probably only necessary to make consistent with multi dimensional chaining interface.

/**
.Function.assign.param.source.type:Class.Seed
.Function.assign.class:Class.Seed
.Function.assign.param.target.type:Class.Seed
.Function.assign.class:Class.Seed
.Function.move.param.source.type:Class.Seed
.Function.move.class:Class.Seed
.Function.move.param.target.type:Class.Seed
.Function.move.class:Class.Seed
..include:seqan/seeds2.h
*/

/**
.Function.getBeginDim0:
..summary: Returns the first position of the seed in the query.
..cat:Seed Handling
..signature:beginDim0(seed)
..class:Class.Seed
..param.seed:The seed whose query position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getEndDim0:
..summary: Returns the last position of the seed in the query.
..cat:Seed Handling
..signature:endDim0(seed)
..class:Class.Seed
..param.seed:The seed whose last in the query position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getBeginDim1:
..summary: Returns the first position of the seed in the database.
..cat:Seed Handling
..signature:beginDim1(seed)
..class:Class.Seed
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getEndDim1:
..summary: Returns the last position of the seed in the database.
..cat:Seed Handling
..signature:endDim1(seed)
..class:Class.Seed
..param.seed:The seed whose last in the database position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getLowerDiagonal:
..summary: Returns the most left diagonal of the seed (maximum diagonal value).
..cat:Seed Handling
..signature:getLowerDiagonal(seed)
..class:Class.Seed
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns:The most left diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getLowerDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._lowerDiagonal;
}

/**
.Function.setLowerDiagonal:
..summary: Sets a new value for the most left diagonal.
..cat:Seed Handling
..signature:setLowerDiagonal(seed, diag)
..class:Class.Seed
..param.seed:The seed whose left diagonal value should be updated.
...type:Class.Seed
..param.diag:The new value for the most left diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig, typename TPosition>
inline void 
setLowerDiagonal(Seed<TSpec, TConfig> & seed, 
				 TPosition newDiag)
{
	SEQAN_CHECKPOINT;
	seed._lowerDiagonal = newDiag;
}

/**
.Function.getUpperDiagonal:
..summary: Returns the most right diagonal of the seed (minimum diagonal value).
..cat:Seed Handling
..signature:getUpperDiagonal(seed)
..class:Class.Seed
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns:The most right diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getUpperDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._upperDiagonal;
}

/**
.Function.setUpperDiagonal:
..summary: Sets a new value for the most right diagonal.
..cat:Seed Handling
..signature:setUpperDiagonal(seed, diag)
..class:Class.Seed
..param.seed:The seed whose right diagonal value should be updated.
...type:Class.Seed
..param.diag:The new value for the most right diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig, typename TPosition>
inline void 
setUpperDiagonal(Seed<TSpec, TConfig> & seed, 
				 TPosition newDiag)
{
	SEQAN_CHECKPOINT;
	seed._upperDiagonal = newDiag;
}


/**
.Function.getSeedSize
..summary:Returns the number of matches and mismatches of the seed.  This is the longest true diagonal fitting into its dimensions.
..signature:getSeedSize(seed)
..class:Class.Seed
..remark:"Seed size" is mostly called "seed length" in the literature.  However, in SeqAn reverse "length" is reserved to be the size of a container.
..include:seqan/seeds2.h
*/
template <typename TSpec, typename TConfig>
inline typename Size<Seed<TSpec, TConfig> >::Type
getSeedSize(Seed<TSpec, TConfig> & seed)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): What if reverse seed?
    return _max(getEndDim0(seed) - getBeginDim0(seed), getEndDim1(seed) - getBeginDim1(seed));
}


template <typename TSpec, typename TConfig>
inline typename Size<Seed<TSpec, TConfig> >::Type
getSeedSize(Seed<TSpec, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): What if reverse seed?
    return _max(getEndDim0(seed) - getBeginDim0(seed), getEndDim1(seed) - getBeginDim1(seed));
}


// Computed values, based on properties returned by getters.

// TODO(holtgrew): Rename to getBeginDiagonal.
/**
.Function.getStartDiagonal:
..summary: Returns the diagonal of the start point.
..cat:Seed Handling
..signature:startDiagonal(seed)
..class:Class.Seed
..param.seed:The seed whose start diagonal should be returned.
...type:Class.Seed
..returns:The diagonal of the start point.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getStartDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return getBeginDim1(seed) - getBeginDim0(seed);
}

/**
.Function.getEndDiagonal:
..summary: Returns the diagonal of the end point.
..cat:Seed Handling
..signature:endDiagonal(seed)
..class:Class.Seed
..param.seed:The seed whose end diagonal should be returned.
...type:Class.Seed
..returns:The diagonal of the end point.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getEndDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return getEndDim1(seed) - getEndDim0(seed);
}

// Functions for seeds with the ScoreMixin_.

// Case: No score, do not update anything.
template <typename TSpec, typename TConfig>
inline void
_updateSeedsScoreMergeHelper(Seed<TSpec, TConfig> & /*seed*/, Seed<TSpec, TConfig> const & /*other*/, False const &)
{
    SEQAN_CHECKPOINT;
}

// Case: Seed has score.  The assumption is that the score is
// proportional to the size of the seed.  Each seed contributes a
// fraction of its score that is proportional to the fraction of the
// score it contributes.
template <typename TSpec, typename TConfig>
inline void
_updateSeedsScoreMergeHelper(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other, True const &)
{
    SEQAN_CHECKPOINT;
    typedef Seed<TSpec, TConfig> TSeed;
    typedef typename Size<TSeed>::Type TSize;

    // Compute new size.
    TSize newBegin0 = _min(getBeginDim0(seed), getBeginDim0(other));
    TSize newEnd0 = _max(getEndDim0(seed), getEndDim0(other));
    TSize newBegin1 = _min(getBeginDim1(seed), getBeginDim1(other));
    TSize newEnd1 = _max(getEndDim1(seed), getEndDim1(other));
    TSize newSize = _max(newEnd0 - newBegin0, newEnd1 - newBegin1);
    // New seed should be larger than either old one and overlap should be > 0.
    SEQAN_ASSERT_GEQ(newSize, getSeedSize(seed));
    SEQAN_ASSERT_GEQ(newSize, getSeedSize(other));
    SEQAN_ASSERT_LEQ(newSize, getSeedSize(seed) + getSeedSize(other));
    TSize overlap = getSeedSize(seed) + getSeedSize(other) - newSize;
    // Overlap should be smaller than or equal to either seed size.
    SEQAN_ASSERT_GEQ(getSeedSize(seed), overlap);
    SEQAN_ASSERT_GEQ(getSeedSize(other), overlap);

    // Compute fraction each seed contributes.
    TSize total = getSeedSize(seed) + getSeedSize(other) - overlap;
    double fracSeed = static_cast<double>(getSeedSize(seed) - 0.5 * overlap) / static_cast<double>(total);
    double fracOther = static_cast<double>(getSeedSize(other) - 0.5 * overlap) / static_cast<double>(total);
    typedef typename SeedScore<TSeed>::Type TScoreValue;
	TScoreValue newScore = static_cast<TScoreValue>(round(fracSeed * getScore(seed) + fracOther * getScore(other)));
    setScore(seed, newScore);
}


// Update the score for merging.  If the seeds do not have scores, nothing is done.
template <typename TSpec, typename TConfig>
inline void
_updateSeedsScoreMerge(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other)
{
    SEQAN_CHECKPOINT;
    typedef Seed<TSpec, TConfig> TSeed;
    _updateSeedsScoreMergeHelper(seed, other, typename HasScore<TSeed>::Type());
}

// Case: No score, do not update anything.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreSimpleChainHelper(Seed<TSpec, TConfig> & /*seed*/, Seed<TSpec, TConfig> const & /*other*/, Score<TScoreValue, Simple> const & /*scoringScheme*/, False const &)
{
    SEQAN_CHECKPOINT;
}

// Case: Seeds have scores.  Update the score of seed according to the
// gap size.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreSimpleChainHelper(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other, Score<TScoreValue, Simple> const & scoringScheme, True const &)
{
    SEQAN_CHECKPOINT;
    typedef Seed<TSpec, TConfig> TSeed;
    typedef typename Size<TSeed>::Type TSize;

    // TODO(holtgrew): seed must be the one to the upper left.

    // Only linear gap costs are supported.
    SEQAN_ASSERT_EQ(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme));
    // Scores for gaps and mismatches must be penalties.
    SEQAN_ASSERT_LT(scoreGap(scoringScheme), static_cast<TScoreValue>(0));
    SEQAN_ASSERT_LT(scoreMismatch(scoringScheme), static_cast<TScoreValue>(0));

    // We use a simple heuristic for updating the seed's score the
    // systematically overestimates the gap costs: We close the gap
    // with a maximal diagonal and then remaining indels or just
    // indels, whatever yields a lower score.
    TSize maxDist = _max(getBeginDim0(other) - getEndDim0(seed), getBeginDim1(other) - getEndDim1(seed));
    TSize minDist = _max(getBeginDim0(other) - getEndDim0(seed), getBeginDim1(other) - getEndDim1(seed));
    TSize diagLen = minDist;
    TSize indelLen = maxDist - minDist;
    TScoreValue gapScore1 = diagLen * scoreMismatch(scoringScheme) + indelLen * scoreGap(scoringScheme);
    TScoreValue gapScore2 = (maxDist + minDist) * scoreGap(scoringScheme);
    TScoreValue gapScore = _max(gapScore1, gapScore2);

    // The new score is the sum of the seed scores and the better of
    // the gap scores computed above.
    setScore(seed, getScore(seed) + getScore(other) + gapScore);
}


// Update the score for simple chaining.  If the seeds do not have scores, nothing is done.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreSimpleChain(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other, Score<TScoreValue, Simple> const & scoringScheme)
{
    SEQAN_CHECKPOINT;
    typedef Seed<TSpec, TConfig> TSeed;
    _updateSeedsScoreSimpleChainHelper(seed, other, scoringScheme, typename HasScore<TSeed>::Type());
}


// Case: No score, do not update anything.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreChaosHelper(Seed<TSpec, TConfig> & /*seed*/, Seed<TSpec, TConfig> const & /*other*/, TScoreValue const & /*scoreDelta*/, False const &)
{
    SEQAN_CHECKPOINT;
}

// Case: Seed has score.  The new score is the sum of the other seed's
// score and the delta as computed by the chaining routine.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreChaosHelper(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other, TScoreValue const & scoreDelta, True const &)
{
    SEQAN_CHECKPOINT;
    setScore(seed, getScore(seed) + getScore(other) + scoreDelta);
}


// Update the score for CHAOS chaining.  If the seeds do not have scores, nothing is done.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreChaos(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other, TScoreValue const & scoreDelta)
{
    SEQAN_CHECKPOINT;
    typedef Seed<TSpec, TConfig> TSeed;
    _updateSeedsScoreChaosHelper(seed, other, scoreDelta, typename HasScore<TSeed>::Type());
}


template <typename TSeed>
inline typename SeedScore<TSeed>::Type
getScore(TSeed const & seed)
{
    SEQAN_CHECKPOINT;
    return seed._score;
}


template <typename TSeed, typename TScore>
inline void
setScore(TSeed & seed, TScore const & score)
{
    SEQAN_CHECKPOINT;
    seed._score = score;
}


// The part of assign() for seeds that have no score mixin.
template <typename TSpec, typename TConfig>
inline void
_assignScoreMixin(Seed<TSpec, TConfig> & /*seed*/, Seed<TSpec, TConfig> const & /*other*/, False const &)
{
    SEQAN_CHECKPOINT;
    // Do nothing, no mixin.
}


// The part of assign() for seeds that have a score mixin.
template <typename TSpec, typename TConfig>
inline void
_assignScoreMixin(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other, True const &)
{
    SEQAN_CHECKPOINT;
    // Update score.
    seed._score = other._score;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_BASE_H_
