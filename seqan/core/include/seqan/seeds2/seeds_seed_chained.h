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
// Specialization "Chained" for class Seed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
#define SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

struct Chained_;
typedef Tag<Chained_> ChainedSeed;  // TODO(holtgrew): Chained already taken as template in file. Maybe prefer non-parameterized types for simpler names.

/**
.Spec.ChainedSeed
..summary:Describes a seed with start and end position2 and diagonal upper and lower bounds. Additionaly diagonal segments
between start and end position2 are stored.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, ChainedSeed>
..param.TPosition:The type of number that schuld be used. Must have negative numbers (e.g. int/long).
.Memfunc.ChainedSeed#Seed:
..class:Spec.ChainedSeed
..summary:Constructor
..signature: Seed<TPosition, ChainedSeed> ()
..signature: Seed<TPosition, ChainedSeed> (qStartPos, dStartPos, length)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.length: Length of the seed.
..include:seqan/seeds.h
*/
template <typename TConfig>
class Seed<ChainedSeed, TConfig>
        : public TConfig::TScoreMixin
{
public:
    typedef typename TConfig::TPosition TPosition;
    typedef typename TConfig::TSize TSize;
    typedef typename TConfig::TDiagonal TDiagonal;

    typedef typename TConfig::TScoreMixin TScoreMixin;

    typedef SeedDiagonal<TPosition, TSize> TSeedDiagonal;

    ::std::list<TSeedDiagonal> _seedDiagonals;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;

    Seed() : TScoreMixin(), _lowerDiagonal(0), _upperDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition beginDim0, TPosition beginDim1, TPosition seedLength)
            : TScoreMixin(),
              _lowerDiagonal(beginDim1 - beginDim0),
              _upperDiagonal(beginDim1 - beginDim0)
    {
        SEQAN_CHECKPOINT;
        appendValue(_seedDiagonals, TSeedDiagonal(beginDim0, beginDim1, seedLength));
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.Value.param.T:Spec.ChainedSeed
///.Metafunction.Value.class:Spec.ChainedSeed
template <typename TConfig>
struct Value<Seed<ChainedSeed, TConfig> >
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Position<TSeed_>::Type TPosition_;
    typedef typename Size<TSeed_>::Type TSize_;

    typedef SeedDiagonal<TPosition_, TSize_> Type;
};

template <typename TConfig>
struct Value<Seed<ChainedSeed, TConfig> const>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Position<TSeed_>::Type TPosition_;
    typedef typename Size<TSeed_>::Type TSize_;

    typedef SeedDiagonal<TPosition_, TSize_> const Type;
};

///.Metafunction.Size.param.T:Spec.ChainedSeed
///.Metafunction.Size.class:Spec.ChainedSeed
template <typename TConfig>
struct Size<Seed<ChainedSeed, TConfig> >
{
    typedef typename TConfig::TSize Type;
};

template <typename TConfig>
struct Size<Seed<ChainedSeed, TConfig> const>
        : Size<Seed<ChainedSeed, TConfig> > {};

///.Metafunction.Iterator.param.T:Spec.ChainedSeed
///.Metafunction.Iterator.class:Spec.ChainedSeed
template <typename TConfig>
struct Iterator<Seed<ChainedSeed, TConfig>, Standard>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef typename ::std::list<TSeedDiagonal_>::iterator Type;
};

template <typename TConfig>
struct Iterator<Seed<ChainedSeed, TConfig> const, Standard>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef typename ::std::list<TSeedDiagonal_>::const_iterator Type;
};

///.Metafunction.Reference.param.T:Spec.ChainedSeed
///.Metafunction.Reference.class:Spec.ChainedSeed
template <typename TConfig>
struct Reference<Seed<ChainedSeed, TConfig> >
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef TSeedDiagonal_ & Type;
};

template <typename TConfig>
struct Reference<Seed<ChainedSeed, TConfig> const>
{
    typedef Seed<ChainedSeed, TConfig> TSeed_;
    typedef typename Value<TSeed_>::Type TSeedDiagonal_;
    typedef TSeedDiagonal_ const & Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, Seed<ChainedSeed, TConfig> const & seed)
{
    typedef Seed<ChainedSeed, TConfig> const TSeed;
    typedef typename Iterator<TSeed>::Type TIterator;

    stream << "Seed<ChainedSeed, TConfig>([";
    for (TIterator it = begin(seed); it != end(seed); ++it) {
        if (it != begin(seed))
            stream << ", ";
        stream << *it;
    }
    stream << "])";
    return stream;
}

template <typename TConfig>
inline bool
operator==(Seed<ChainedSeed, TConfig> const & a, Seed<ChainedSeed, TConfig> const & b)
{
    SEQAN_CHECKPOINT;
    return a._seedDiagonals == b._seedDiagonals &&
            a._upperDiagonal == b._upperDiagonal &&
            a._lowerDiagonal == b._lowerDiagonal;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getBeginDim0(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return front(seed._seedDiagonals).beginDim0;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getEndDim0(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return back(seed._seedDiagonals).beginDim0 + back(seed._seedDiagonals).length;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getBeginDim1(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return front(seed._seedDiagonals).beginDim1;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getEndDim1(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return back(seed._seedDiagonals).beginDim1 + back(seed._seedDiagonals).length;
}

/**
.Function.length.param.object.type:Spec.ChainedSeed
.Function.length.class:Spec.ChainedSeed
..include:seqan/seeds2.h
*/
template <typename TConfig>
inline typename Size<Seed<ChainedSeed, TConfig> >::Type
length(Seed<ChainedSeed, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    return length(seed._seedDiagonals);
}

/**
.Function.appendDiagonal
..summary: Adds diagonal to the Chained Seed.
..cat:Seed Handling
..signature:appendDiag(seed, diagonal)
..class:Spec.ChainedSeed
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeed
..param.diag: The diagonal to add.
...type:Class.SeedDiagonal
...remarks: A diagonal consists of three values: 1: start in 1. sequence, 2: start in 2. sequence, 3: length of match
..include:seqan/seeds2.h
*/
template <typename TConfig>
inline void
appendDiagonal(Seed<ChainedSeed, TConfig> & seed,
               typename Value<Seed<ChainedSeed, TConfig> >::Type const & diagonal)
{
    SEQAN_CHECKPOINT;

    if (length(seed) > 0) {
        SEQAN_ASSERT_LEQ(back(seed._seedDiagonals).beginDim0 + back(seed._seedDiagonals).length, diagonal.beginDim0);
        SEQAN_ASSERT_LEQ(back(seed._seedDiagonals).beginDim1 + back(seed._seedDiagonals).length, diagonal.beginDim1);
    }

    appendValue(seed._seedDiagonals, diagonal);
}

/**
.Function.truncateDiagonals
..summary:Removes diagonals from the given first one to the end of the seed's diagonals.
..cat:Seed Handling
..signature:truncateDiagonals(seed, first)
..class:Spec.ChainedSeed
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeed
..param.first: Iterator the first diagonal to remove.
..include:seqan/seeds2.h
*/
template <typename TConfig>
inline void
truncateDiagonals(Seed<ChainedSeed, TConfig> & seed,
                  typename Iterator<Seed<ChainedSeed, TConfig> >::Type const & first)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Add erase() to std::list adaptors?
    seed._seedDiagonals.erase(first, seed._seedDiagonals.end());
}

/**
.Function.begin.param.object.type:Spec.ChainedSeed
..class:Spec.ChainedSeed
..include:seqan/seeds2.h
*/
template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> >::Type
begin(Seed<ChainedSeed, TConfig> & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.begin();
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
begin(Seed<ChainedSeed, TConfig> const & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.begin();
}

/**
.Function.front.param.object.type:Spec.ChainedSeed
..class:Spec.ChainedSeed
..include:seqan/seeds2.h
*/
template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> >::Type
front(Seed<ChainedSeed, TConfig> & seed)
{
    SEQAN_CHECKPOINT;
    return front(seed._seedDiagonals);
}

template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> const>::Type
front(Seed<ChainedSeed, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    return front(seed._seedDiagonals);
}

/**
.Function.back.param.container.type:Spec.ChainedSeed
..class:Spec.ChainedSeed
..include:seqan/seeds2.h
*/
template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> >::Type
back(Seed<ChainedSeed, TConfig> & seed)
{
    SEQAN_CHECKPOINT;
    return back(seed._seedDiagonals);
}

template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> const>::Type
back(Seed<ChainedSeed, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    return back(seed._seedDiagonals);
}

/**
.Function.end.param.object.type:Spec.ChainedSeed
..class:Spec.ChainedSeed
..include:seqan/seeds2.h
*/
template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> >::Type
end(Seed<ChainedSeed, TConfig> & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.end();
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
end(Seed<ChainedSeed, TConfig> const & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.end();
}

// Basic Functions

template <typename TConfig>
void
move(Seed<ChainedSeed, TConfig> & target, Seed<ChainedSeed, TConfig> & source)
{
    SEQAN_CHECKPOINT;
    std::swap(target._seedDiagonals, source._seedDiagonals);
    target._lowerDiagonal = source._lowerDiagonal;
    target._upperDiagonal = source._upperDiagonal;
    _assignScoreMixin(target, source, typename HasScore<Seed<Simple, TConfig> >::Type());
}

template <typename TConfig>
void
assign(Seed<ChainedSeed, TConfig> & target, Seed<ChainedSeed, TConfig> const & source)
{
    SEQAN_CHECKPOINT;
    assign(target._seedDiagonals, source._seedDiagonals);
    target._lowerDiagonal = source._lowerDiagonal;
    target._upperDiagonal = source._upperDiagonal;
    _assignScoreMixin(target, source, typename HasScore<Seed<Simple, TConfig> >::Type());
}

template <typename TConfig>
void
assign(Seed<ChainedSeed, TConfig> & target, Seed<ChainedSeed, TConfig> & source)
{
    SEQAN_CHECKPOINT;
    typedef Seed<ChainedSeed, TConfig> TSeed;
    assign(target, const_cast<TSeed const &>(source));
}

// Debug Output

template <typename TStream, typename TConfig>
inline void
__write(TStream & stream, Seed<ChainedSeed, TConfig> const & seed, Tikz_ const &)
{
//IOREV _nodoc_ specialization not documented
    // Overall seed.
    stream << "\\draw[seed] (" << getBeginDim1(seed) << ", -" << getBeginDim0(seed) << ") -- (" << (getEndDim1(seed) - 1) << ", -" << (getEndDim0(seed) - 1) << ");" << std::endl;
    // Diagonals.
    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Iterator<TSeed const, Standard>::Type TIterator;
    for (TIterator it = begin(seed); it != end(seed); ++it)
        stream << "\\draw[seed diagonal] (" << it->beginDim1 << ", -" << it->beginDim0 << ") -- (" << (it->beginDim1 + it->length - 1) << ", -" << (it->beginDim0 + it->length - 1) << ");" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
