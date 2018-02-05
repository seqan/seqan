// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: René Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BLOCK_BAND_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BLOCK_BAND_H_

#include <cstdlib>

namespace seqan
{

template <typename TValue, typename TTag>
class NamedType
{
public:

    NamedType() = default;

    explicit NamedType(TValue const & _v) : val(_v)
    {}

    template <typename TValue_ = TValue, typename = std::enable_if_t<std::is_reference<TValue_>::value>>
    explicit NamedType(TValue && _v) : val(std::move(_v))
    {}

    TValue & get()
    {
        return val;
    }

    TValue const & get() const
    {
        return val;
    }

    operator TValue() const
    {
        return val;
    }

private:
    TValue val{};
};

namespace impl
{

using HorizontalBandPos = NamedType<size_t, struct DPStaticBandGeneratorHorizontalBandPos_>;
using VerticalBandPos   = NamedType<size_t, struct DPStaticBandGeneratorVerticalBandPos_>;

using GridColumnIndex   = NamedType<size_t, struct DPStaticBandGeneratorGridColumnIndex_>;
using GridRowIndex      = NamedType<size_t, struct DPStaticBandGeneratorGridRowIndex_>;
using GridColumnSize    = NamedType<size_t, struct DPStaticBandGeneratorGridColumnSize_>;
using GridRowSize       = NamedType<size_t, struct DPStaticBandGeneratorGridRowSize_>;

class DPWavefrontBand
{
public:
    using TGrid = String<bool, Packed<>>;

    TGrid    _grid{};         // blocks to create the task graph on.
    uint32_t _hBandPos{0};   // offset for the column size
    uint32_t _vBandPos{0};   // offset for the column size
    uint32_t _colSize{0};   // offset for the column size
    uint32_t _rowSize{0};   // offset for the column size

    DPWavefrontBand() = delete;

    DPWavefrontBand(HorizontalBandPos const & hBandPos,
                    VerticalBandPos   const & vBandPos,
                    GridColumnSize    const & colSize,
                    GridRowSize       const & rowSize) :
        _hBandPos(hBandPos.get()),
        _vBandPos(vBandPos.get()),
        _colSize(colSize.get()),
        _rowSize(rowSize.get())
    {
        resize(_grid, _colSize * _rowSize, false, Exact());
    }
};

struct DPWavefrontBandIterSpec_;
using DPWavefrontBandIterSpec = Tag<DPWavefrontBandIterSpec_>;

}  // namespace impl

template <typename TBand>
class Iter<TBand, impl::DPWavefrontBandIterSpec>
{
public:

    using TUnderlyingIterator = typename Iterator<std::conditional_t<std::is_same<TBand, TBand const>::value,
                                                                     typename TBand::TGrid const,
                                                                     typename TBand::TGrid>,
                                                  Standard>::Type;

    TBand *             _contPtr{nullptr};
    TUnderlyingIterator _underlyingIter{};

    Iter() = default;

    template <typename TOtherBand, typename = std::enable_if_t<IsConstructibleFrom<TBand, TOtherBand>::VALUE>>
    Iter(Iter<TOtherBand, impl::DPWavefrontBandIterSpec> const & other) :
        _contPtr{other._contPtr},
        _underlyingIter{other._underlyingIter}
    {}

    template <typename TOtherBand, typename = std::enable_if_t<IsConstructibleFrom<TBand, TOtherBand>::VALUE>>
    Iter(Iter<TOtherBand, impl::DPWavefrontBandIterSpec> && other) :
        _contPtr{other._contPtr},
        _underlyingIter{std::move(other._underlyingIter)}
    {}

    // Create from underlying iterator.
    Iter(TBand & _cont, TUnderlyingIterator const & _iter) :
        _contPtr{&_cont},
        _underlyingIter{_iter}
    {}

    template <typename TOtherBand, typename = std::enable_if_t<IsConstructibleFrom<TBand, TOtherBand>::VALUE>>
    Iter & operator=(Iter<TOtherBand, impl::DPWavefrontBandIterSpec> other)
    {
        using std::swap;
        swap(_contPtr, other._contPtr);
        swap(_underlyingIter, other._underlyingIter);
        return *this;
    }

    ~Iter() = default;
};

// --------------------------------------------------------------------------
// Function operator*()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto operator*(Iter<TBand, impl::DPWavefrontBandIterSpec> & me)
{
    return *me._underlyingIter;
}

// --------------------------------------------------------------------------
// Function operator*()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto operator*(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return *me._underlyingIter;
}

// --------------------------------------------------------------------------
// Function operator++()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto & operator++(Iter<TBand, impl::DPWavefrontBandIterSpec> & me)
{
    ++me._underlyingIter;
    return me;
}

template <typename TBand>
inline auto operator++(Iter<TBand, impl::DPWavefrontBandIterSpec> & me, int)
{
    Iter<TBand, impl::DPWavefrontBandIterSpec> tmp{me};
    ++me._underlyingIter;
    return tmp;
}

// --------------------------------------------------------------------------
// Function operator+=()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto & operator+=(Iter<TBand, impl::DPWavefrontBandIterSpec> & me, int const off)
{
    me._underlyingIter += off;
    return me;
}

// --------------------------------------------------------------------------
// Function operator+()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto operator+(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me, int const off)
{
    Iter<TBand, impl::DPWavefrontBandIterSpec> tmp{me};
    tmp += off;
    return tmp;
}

// --------------------------------------------------------------------------
// Function operator--()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto & operator--(Iter<TBand, impl::DPWavefrontBandIterSpec> & me)
{
    --me._underlyingIter;
    return me;
}

template <typename TBand>
inline auto operator--(Iter<TBand, impl::DPWavefrontBandIterSpec> & me, int)
{
    Iter<TBand, impl::DPWavefrontBandIterSpec> tmp{me};
    --me._underlyingIter;
    return tmp;
}

// --------------------------------------------------------------------------
// Function operator-=()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto & operator-=(Iter<TBand, impl::DPWavefrontBandIterSpec> & me, int const off)
{
    me._underlyingIter -= off;
    return me;
}

// --------------------------------------------------------------------------
// Function operator-()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto operator-(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me, int const off)
{
    Iter<TBand, impl::DPWavefrontBandIterSpec> tmp{me};
    tmp -= off;
    return tmp;
}

// --------------------------------------------------------------------------
// Function operator-()
// --------------------------------------------------------------------------

template <typename TBand>
inline size_t
operator-(Iter<TBand, impl::DPWavefrontBandIterSpec> const & lhs,
          Iter<TBand, impl::DPWavefrontBandIterSpec> const & rhs)
{
    return lhs._underlyingIter - rhs._underlyingIter;
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename TBandLhs, typename TBandRhs,
          typename = std::enable_if_t<IsConstructibleFrom<TBandLhs, TBandRhs>::VALUE>>
inline bool
operator==(Iter<TBandLhs, impl::DPWavefrontBandIterSpec> const & lhs,
           Iter<TBandRhs, impl::DPWavefrontBandIterSpec> const & rhs)
{
    return lhs._underlyingIter == rhs._underlyingIter;
}

// --------------------------------------------------------------------------
// Function operator!=()
// --------------------------------------------------------------------------

template <typename TBandLhs, typename TBandRhs,
          typename = std::enable_if_t<IsConstructibleFrom<TBandLhs, TBandRhs>::VALUE>>
inline bool
operator!=(Iter<TBandLhs, impl::DPWavefrontBandIterSpec> const & lhs,
           Iter<TBandRhs, impl::DPWavefrontBandIterSpec> const & rhs)
{
    return !(lhs._underlyingIter == rhs._underlyingIter);
}

// --------------------------------------------------------------------------
// Function container()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto &
container(Iter<TBand, impl::DPWavefrontBandIterSpec> & me)
{
    return *me._contPtr;
}

template <typename TBand>
inline auto const &
container(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return *me._contPtr;
}

// --------------------------------------------------------------------------
// Function rowIndex()
// --------------------------------------------------------------------------

template <typename TBand>
inline impl::GridRowIndex
rowIndex(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return impl::GridRowIndex{position(me) % container(me)._colSize};
}

// --------------------------------------------------------------------------
// Function columnIndex()
// --------------------------------------------------------------------------

template <typename TBand>
inline impl::GridColumnIndex
columnIndex(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return impl::GridColumnIndex{position(me) / container(me)._colSize};
}

// --------------------------------------------------------------------------
// Function coordinates()
// --------------------------------------------------------------------------

template <typename TBand>
inline std::pair<impl::GridColumnIndex, impl::GridRowIndex>
coordinates(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return {columnIndex(me), rowIndex(me)};
}

// --------------------------------------------------------------------------
// Function position()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto
position(Iter<TBand const, impl::DPWavefrontBandIterSpec> const & me)
{
    return me - begin(container(me), Standard{});
}

template <typename TBand>
inline auto
position(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return me - begin(*const_cast<TBand*>(&container(me)), Standard{});
}

// --------------------------------------------------------------------------
// Function hasPredecessorLeft()
// --------------------------------------------------------------------------

template <typename TBand>
inline bool
hasPredecessorLeft(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return (columnIndex(me).get() > 0) ? container(me)._grid[position(me) - container(me)._colSize] : false;
}

// --------------------------------------------------------------------------
// Function hasPredecessorAbove()
// --------------------------------------------------------------------------

template <typename TBand>
inline bool
hasPredecessorAbove(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return (rowIndex(me).get() > 0) ? container(me)._grid[position(me) - 1] : false;
}

// --------------------------------------------------------------------------
// Function hasSuccessorRight()
// --------------------------------------------------------------------------

template <typename TBand>
inline bool
hasSuccessorRight(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return (columnIndex(me).get() < container(me)._rowSize - 1) ?
        container(me)._grid[position(me) + container(me)._colSize] :
        false;
}

// --------------------------------------------------------------------------
// Function hasSuccessorBelow()
// --------------------------------------------------------------------------

template <typename TBand>
inline bool
hasSuccessorBelow(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return (rowIndex(me).get() < container(me)._colSize - 1) ? container(me)._grid[position(me) + 1] : false;
}

// --------------------------------------------------------------------------
// Function successorRight()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto
successorRight(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return me + static_cast<int>(container(me)._colSize);
}

// --------------------------------------------------------------------------
// Function successorBelow()
// --------------------------------------------------------------------------

template <typename TBand>
inline auto
successorBelow(Iter<TBand, impl::DPWavefrontBandIterSpec> const & me)
{
    return me + 1;
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

template <typename TIterSpec>
inline auto begin(impl::DPWavefrontBand & band, Tag<TIterSpec> const /*tag*/)
{
    return Iter<impl::DPWavefrontBand, impl::DPWavefrontBandIterSpec>{band, begin(band._grid, Standard{})};
}

template <typename TIterSpec>
inline auto begin(impl::DPWavefrontBand const & band, Tag<TIterSpec> const /*tag*/)
{
    return Iter<impl::DPWavefrontBand const, impl::DPWavefrontBandIterSpec>{band, begin(band._grid, Standard{})};
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

template <typename TIterSpec>
inline auto end(impl::DPWavefrontBand & band, Tag<TIterSpec> const /*tag*/)
{
    return Iter<impl::DPWavefrontBand, impl::DPWavefrontBandIterSpec>{band, end(band._grid, Standard{})};
}

template <typename TIterSpec>
inline auto end(impl::DPWavefrontBand const & band, Tag<TIterSpec> const /*tag*/)
{
    return Iter<impl::DPWavefrontBand const, impl::DPWavefrontBandIterSpec>{band, end(band._grid, Standard{})};
}

// --------------------------------------------------------------------------
// Function transformToGrid()
// --------------------------------------------------------------------------

template <typename TSeqH, typename TSeqV, typename TBandSpec>
inline auto transformToGrid(TSeqH const & seqH,
                            TSeqV const & seqV,
                            size_t const blockSize,
                            DPBandConfig<TBandSpec> const & bandConfig)
{
    // temp variables for lower and upper diagonal.
    uint32_t _lowerDiag{};
    uint32_t _upperDiag{};

    SEQAN_IF_CONSTEXPR (std::is_same<TBandSpec, BandOff>::value)
    {  // If no band is available, we will use a band that spans over the full matrix.
        _lowerDiag = length(seqV);
        _upperDiag = length(seqH);
    }
    else
    {
        using TSizeType = decltype(length(seqH));
        _lowerDiag = std::min(static_cast<TSizeType>(std::abs(lowerDiagonal(bandConfig))), length(seqV));
        _upperDiag = std::min(static_cast<TSizeType>(std::abs(upperDiagonal(bandConfig))), length(seqH));
    }

    SEQAN_ASSERT_GT(_lowerDiag, 0u);
    SEQAN_ASSERT_GT(_upperDiag, 0u);

    // Number of blocks in either dimension.
    impl::GridColumnSize cSize{(length(seqV) + blockSize) / blockSize};
    impl::GridRowSize    rSize{(length(seqH) + blockSize) / blockSize};

    impl::DPWavefrontBand band{impl::HorizontalBandPos{_upperDiag}, impl::VerticalBandPos{_lowerDiag}, cSize, rSize};

    impl::GridColumnIndex currCol;
    impl::GridRowIndex    currRow;

    /*     0     1     2    3     4     5
       |-----|-----|---*-|-----|-----|-----|
       |     |     |    *|     |     |     |
    0  |     |     |     *     |     |     |
       |-----|-----|-----|*----|-----|-----|
       |*    |     |     | *   |     |     |
    1  | *   |     |     |  *  |     |     |
       |--*--|-----|-----|---*-|-----|-----|
       |   * |     |     |    *|     |     |
    2  |    *|     |     |     *     |     |
       |-----*-----|-----|-----|*----|-----|
       |     |*    |     |     | *   |     |
    3  |     | *   |     |     |  *  |     |
       |-----|--*--|-----|-----|---*-|-----|
       |     |   * |     |     |    *|     |
    4  |     |    *|     |     |     *     |
       |-----|-----*-----|-----|-----|*----|
    * Blocks part of the band:
    * {0,0}, {0,1}, {0,2}
    * {1,0}, {1,1}, {1,2}, {1,3}, {1,4}
    * {2,0}, {2,1}, {2,2}, {2,3}, {2,4}
    * {3,0}, {3,1}, {3,2}, {3,3}, {3,4}
    * {4,2}, {4,3}, {4,4}
    * {5,4}
    */

    using seqan::begin;
    using seqan::end;
    // Move along the matrix and set all blocks active that are part of the band.
    for (auto it = begin(band, Standard()); it != end(band, Standard()); ++it)
    {
       std::tie(currCol, currRow) = coordinates(it);
       if (currCol.get() * blockSize <= (_upperDiag + (blockSize * currRow.get() + (blockSize - 1))) &&
           currRow.get() * blockSize <= (_lowerDiag + (blockSize * currCol.get() + (blockSize - 1))))
            *it = true;
    }
    return band;
}

} // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BLOCK_BAND_H_
