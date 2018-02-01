
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// The navigator for the sparse score dp-matrix. This class also provides an
// iterator for the active and the previous column. It stores the neighbouring
// cells needed for the recursion formula.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_SPARSE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_SPARSE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPMatrixNavigator                        [SparseDPMatrix, ScoreMatrix]
// ----------------------------------------------------------------------------

// Specialization of the score matrix navigator for a sparse dp matrix.
template <typename TValue, typename THost, typename TNavigationSpec>
class DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>, DPScoreMatrix, TNavigationSpec>
{
public:
    typedef  DPMatrix_<TValue, SparseDPMatrix, THost> TDPMatrix_;
    typedef typename Pointer_<TDPMatrix_>::Type TDPMatrixPointer_;
    typedef typename Iterator<TDPMatrix_, Standard>::Type TDPMatrixIterator;

    template <typename TBandSpec,
              std::enable_if_t<std::is_same<TBandSpec, BandOff>::value, int> = 0>
    DPMatrixNavigator_(TDPMatrix_ & matrix,
                       DPBandConfig<TBandSpec> const & /*band*/)
    {
        _ptrDataContainer = &matrix;
        _activeColIterator = begin(matrix, Standard());
        _prevColIterator = _activeColIterator;
        _laneLeap = 1 - _dataLengths(matrix)[DPMatrixDimension_::VERTICAL];
        *_activeColIterator = TValue();
    }

    template <typename TBandSpec,
              std::enable_if_t<!std::is_same<TBandSpec, BandOff>::value, int> = 0>
    DPMatrixNavigator_(TDPMatrix_ & matrix,
                       DPBandConfig<TBandSpec> const & band)
    {
        typedef typename Size<TDPMatrix_>::Type TSize;
        typedef std::make_signed_t<TSize> TSignedSize;
        _ptrDataContainer = &matrix;

        // Band begins within the first row.
        if (lowerDiagonal(band) >= 0)
        {
            _laneLeap = 0;
            _activeColIterator = begin(matrix, Standard()) + length(matrix, DPMatrixDimension_::VERTICAL) - 1;
        }
        else if (upperDiagonal(band) <= 0) // Band begins within the first column
        {
            _laneLeap = 1 - _dataLengths(matrix)[DPMatrixDimension_::VERTICAL];
            _activeColIterator = begin(matrix, Standard());
        }
        else  // Band intersects with the point of origin.
        {
            _laneLeap = _max(lowerDiagonal(band), 1 - static_cast<TSignedSize>(length(matrix, DPMatrixDimension_::VERTICAL)));
            _activeColIterator = begin(matrix, Standard()) + length(matrix, DPMatrixDimension_::VERTICAL) + _laneLeap - 1;
        }
        _prevColIterator = _activeColIterator;
        *_activeColIterator = TValue();
    }

    TDPMatrixPointer_ _ptrDataContainer{nullptr};   // Pointer to the underlying matrix to navigate on.
    int               _laneLeap{0};                 // The distance to leap when going to the next column.
    size_t            _prevColIteratorOffset{0};    // Offset to reset the previous column iterator when going to the next cell.
    TDPMatrixIterator _activeColIterator{};         // The iterator over the active column.
    TDPMatrixIterator _prevColIterator{};           // The iterator over the previous column. Only needed in the banded case.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _goNextCell                                   [unbanded, FirstCell]
// ----------------------------------------------------------------------------

// specialized for initalization column.
template <typename TValue, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            FirstCell const &)
{
    // no-op
}

// all other column types.
template <typename TValue, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            FirstCell const &)
{
    // Set to begin of column.
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                   [unbanded, InnerCell]
// ----------------------------------------------------------------------------

// specialized for the initialization column.
template <typename TValue, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            InnerCell const &)
{
    ++dpNavigator._activeColIterator;
}

// version for the all other column types.
template <typename TValue, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            InnerCell const &)
{
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                    [unbanded, LastCell]
// ----------------------------------------------------------------------------

// specilaized for initialization column.
template <typename TValue, typename THost>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            LastCell const &)
{
    ++dpNavigator._activeColIterator;
}

// version for all other column types.
template <typename TValue, typename THost,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            LastCell const &)
{
    ++dpNavigator._activeColIterator; // go to next cell.
}

// ----------------------------------------------------------------------------
// Function previousCellHorizontal()
// ----------------------------------------------------------------------------

// unbanded.
template <typename TValue, typename THost>
inline typename Reference<DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                             DPScoreMatrix,
                                             NavigateColumnWise>
                         >::Type
previousCellHorizontal(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                          DPScoreMatrix,
                                          NavigateColumnWise> & dpNavigator)
{
    return *dpNavigator._activeColIterator;
}

template <typename TValue, typename THost, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                             DPScoreMatrix,
                                             NavigateColumnWise> const
                         >::Type
previousCellHorizontal(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                          DPScoreMatrix,
                                          NavigateColumnWise> const & dpNavigator)
{
    return *dpNavigator._activeColIterator;
}

// banded.
template <typename TValue, typename THost>
inline typename Reference<DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                             DPScoreMatrix,
                                             NavigateColumnWiseBanded>
                         >::Type
previousCellHorizontal(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                          DPScoreMatrix,
                                          NavigateColumnWiseBanded> & dpNavigator)
{
    return *dpNavigator._prevColIterator;
}

template <typename TValue, typename THost>
inline typename Reference<DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                             DPScoreMatrix,
                                             NavigateColumnWiseBanded> const
                         >::Type
previousCellHorizontal(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix, THost>,
                                          DPScoreMatrix,
                                          NavigateColumnWiseBanded> const & dpNavigator)
{
    return *(dpNavigator._prevColIterator);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_SPARSE_H_
