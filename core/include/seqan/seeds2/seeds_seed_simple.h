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
// Specialization "Simple" for class Seed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
#define SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Spec.SimpleSeed:
..summary:Describes a seed with start and end position and diagonal upper and lower bounds.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, SimpleSeed>
..param.TPosition:The type number that should be used. Must have negative numbers (e.g. int/long).
.Memfunc.SimpleSeed#Seed:
..class:Spec.SimpleSeed
..summary:Constructor
..signature: Seed<TPosition, SimpleSeed> ()
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, length)
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, qEndPos, dEndPos)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.qEndPos: End in query sequence.
..param.dEndPos: End in database sequence.
..param.length: Length of the seed.
..include:seqan/seeds.h
*/
template <typename TConfiguration>
class Seed<Simple, TConfiguration>
        : public TConfiguration::TScoreMixin
{
  public:
    typedef typename TConfiguration::TPosition TPosition;
    typedef typename TConfiguration::TSize TSize;
    typedef typename TConfiguration::TDiagonal TDiagonal;

    typedef typename TConfiguration::TScoreMixin TScoreMixin;

    TPosition _beginDim0;
    TPosition _beginDim1;
    TPosition _endDim0;
    TPosition _endDim1;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;

    Seed()
            : TScoreMixin(),
              _beginDim0(0), _beginDim1(0), _endDim0(0), _endDim1(0),
              _lowerDiagonal(0), _upperDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition beginDim0, TPosition beginDim1, TPosition seedLength)
            : TScoreMixin(),
              _beginDim0(beginDim0),
              _beginDim1(beginDim1),
              _endDim0(beginDim0 + seedLength),
              _endDim1(beginDim1 + seedLength),
              _lowerDiagonal(static_cast<TDiagonal>(beginDim1 - beginDim0)),
              _upperDiagonal(static_cast<TDiagonal>(beginDim1 - beginDim0))
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_GEQ(_upperDiagonal, _lowerDiagonal);
    }

    Seed(TPosition beginDim0, TPosition beginDim1, TPosition endDim0,
         TPosition endDim1)
            : TScoreMixin(),
              _beginDim0(beginDim0),
              _beginDim1(beginDim1),
              _endDim0(endDim0),
              _endDim1(endDim1),
              _lowerDiagonal(_min(static_cast<TDiagonal>(beginDim1 - beginDim0), static_cast<TDiagonal>(endDim1 - endDim0))),
              _upperDiagonal(_max(static_cast<TDiagonal>(beginDim1 - beginDim0), static_cast<TDiagonal>(endDim1 - endDim0)))
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_GEQ(_upperDiagonal, _lowerDiagonal);
    }

    template <typename TSeed2>
    Seed(TSeed2 const & other)
            : TScoreMixin(),
              _beginDim0(getBeginDim0(other)),
              _beginDim1(getBeginDim1(other)),
              _endDim0(getEndDim0(other)),
              _endDim1(getEndDim1(other)),
              _lowerDiagonal(getLowerDiagonal(other)),
              _upperDiagonal(getUpperDiagonal(other))
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_GEQ(_upperDiagonal, _lowerDiagonal);
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, Seed<Simple, TConfig> const & seed)
{
    return stream << "Seed<Simple, TConfig>(" << getBeginDim0(seed)
                  << ", " << getBeginDim1(seed) << ", "
                  << getEndDim0(seed) << ", "
                  << getEndDim1(seed) << ", lower diag = " << getLowerDiagonal(seed) << ", upper diag = " << getUpperDiagonal(seed) << ")";
}

template <typename TConfig>
inline bool
operator==(Seed<Simple, TConfig> const & a, Seed<Simple, TConfig> const & b)
{
    SEQAN_CHECKPOINT;
    return a._beginDim0 == b._beginDim0 &&
            a._beginDim1 == b._beginDim1 &&
            a._endDim0 == b._endDim0 &&
            a._endDim1 == b._endDim1 &&
            a._upperDiagonal == b._upperDiagonal &&
            a._lowerDiagonal == b._lowerDiagonal;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getBeginDim0(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._beginDim0;
}

/**
.Function.setBeginDim0:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setBeginDim0(seed, start)
..class:Spec.SimpleSeed
..param.seed:The seed whose start position should be updated.
...type:Spec.SimpleSeed
..param.start:The query position where the seed should start.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setBeginDim0(Seed<Simple, TConfig> & seed, 
             TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._beginDim0 = newLeftPosition;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getEndDim0(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._endDim0;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getBeginDim1(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._beginDim1;
}

/**
.Function.setBeginDim1:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setBeginDim1(seed, start)
..class:Spec.SimpleSeed
..param.seed:The seed whose start position should be updated.
...type:Spec.SimpleSeed
..param.start:The database position where the seed should start.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setBeginDim1(Seed<Simple, TConfig> & seed, 
             TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._beginDim1 = newLeftPosition;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getEndDim1(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._endDim1;
}

/**
.Function.setEndDim0:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setEndDim0(seed, end)
..class:Spec.SimpleSeed
..param.seed:The seed whose end position should be updated.
...type:Spec.SimpleSeed
..param.end:The query position where the seed should end.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setEndDim0(Seed<Simple, TConfig> & seed, 
           TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._endDim0 = newRightPosition;
}

/**
.Function.setEndDim1:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setEndDim1(seed, end)
..class:Spec.SimpleSeed
..param.seed:The seed whose end position should be updated.
...type:Spec.SimpleSeed
..param.end:The database position where the seed should end.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setEndDim1(Seed<Simple, TConfig> & seed, 
           TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._endDim1 = newRightPosition;
}

// Basic Functions

template <typename TConfig>
void
move(Seed<Simple, TConfig> & target, Seed<Simple, TConfig> & source)
{
    SEQAN_CHECKPOINT;
    target._beginDim0 = source._beginDim0;
    target._beginDim1 = source._beginDim1;
    target._endDim0 = source._endDim0;
    target._endDim1 = source._endDim1;
    target._lowerDiagonal = source._lowerDiagonal;
    target._upperDiagonal = source._upperDiagonal;
    _assignScoreMixin(target, source, typename HasScore<Seed<Simple, TConfig> >::Type());
}

template <typename TConfig>
void
assign(Seed<Simple, TConfig> & target, Seed<Simple, TConfig> const & source)
{
    SEQAN_CHECKPOINT;
    target._beginDim0 = source._beginDim0;
    target._beginDim1 = source._beginDim1;
    target._endDim0 = source._endDim0;
    target._endDim1 = source._endDim1;
    target._lowerDiagonal = source._lowerDiagonal;
    target._upperDiagonal = source._upperDiagonal;
    _assignScoreMixin(target, source, typename HasScore<Seed<Simple, TConfig> >::Type());
}

template <typename TConfig>
void
assign(Seed<Simple, TConfig> & target, Seed<Simple, TConfig> & source)
{
    SEQAN_CHECKPOINT;
    typedef Seed<Simple, TConfig> TSeed;
    assign(target, const_cast<TSeed const &>(source));
}

// Debug Output

struct Tikz_ {};

template <typename TStream, typename TConfig>
inline void
__write(TStream & stream, Seed<Simple, TConfig> const & seed, Tikz_ const &)
{
//IOREV _nodoc_ specialization not documented
    stream << "\\draw[seed] (" << getBeginDim1(seed) << ", -" << getBeginDim0(seed) << ") -- (" << (getEndDim1(seed) - 1) << ", -" << (getEndDim0(seed) - 1) << ");" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
