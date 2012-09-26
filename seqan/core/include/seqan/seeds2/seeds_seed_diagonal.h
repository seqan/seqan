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
// SeedDiagonal objects store information about seed parts/segments of a
// ChainedSeed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_
#define SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Class.SeedDiagonal
..summary:Store the information about a seed segment.
..signature:SeedDiagonal<TPosition, TSize>
..cat:Seed Handling
..param.TPosition:The type to use for positions.
..param.TSize:The type to use for the seed length.
..include:seqan/seeds2.h
*/
template <typename TPosition, typename TSize>
class SeedDiagonal
{
public:
/**
.Memvar.SeedDiagonal#beginDim0:The position in the query (dimension 0).
*/
    TPosition beginDim0;
/**
 .Memvar.SeedDiagonal#beginDim0:The position in the database sequence (dimension 1).
*/
    TPosition beginDim1;
/**
.Memvar.SeedDiagonal#length:The length of the diagonal.
*/
    TSize length;
    
    SeedDiagonal()
            : beginDim0(0), beginDim1(0), length(0)
    { SEQAN_CHECKPOINT; }

    SeedDiagonal(TPosition _beginDim0, TPosition _beginDim1, TSize _length)
            : beginDim0(_beginDim0), beginDim1(_beginDim1), length(_length)
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================
  
///.Metafunction.Position.param.T:SeedDiagonal
///.Metafunction.Position.class:SeedDiagonal
template <typename TPosition, typename TSize>
struct Position<SeedDiagonal<TPosition, TSize> >
{
    typedef TPosition Type;
};

template <typename TPosition, typename TSize>
struct Position<SeedDiagonal<TPosition, TSize> const>
  : Position<SeedDiagonal<TPosition, TSize> > {};

///.Metafunction.Size.param.T:SeedDiagonal
///.Metafunction.Size.class:SeedDiagonal
template <typename TPosition, typename TSize>
struct Size<SeedDiagonal<TPosition, TSize> >
{
  typedef TSize Type;
};

template <typename TPosition, typename TSize>
struct Size<SeedDiagonal<TPosition, TSize> const>
  : Size<SeedDiagonal<TPosition, TSize> > {};
  
// ===========================================================================
// Functions
// ===========================================================================

template <typename TStream, typename TPosition, typename TSize>
inline TStream &
operator<<(TStream & stream, SeedDiagonal<TPosition, TSize> const & seedDiagonal)
{
    return stream << "SeedDiagonal(" << seedDiagonal.beginDim0 << ", " << seedDiagonal.beginDim1 << ", "  << seedDiagonal.length << ")";
}


template <typename TPosition, typename TSize>
inline bool
operator==(SeedDiagonal<TPosition, TSize> const & a, SeedDiagonal<TPosition, TSize> const & b)
{
    SEQAN_CHECKPOINT;
    return a.beginDim0 == b.beginDim0 && a.beginDim1 == b.beginDim1 && a.length == b.length;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_
