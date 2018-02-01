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
// Author: Renï¿½ Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// This file implements the class DPMatrix and its specialization
// FullDPMatrix. The DPMatrix is a wrapper class around the Matrix<TValue,2>
// class. Thus we can implement different specializations for the dp-matrix
// that are used through the different dp-algorithms. For example, we need
// a full dp matrix to store the traceback or the score for the Waterman-
// Eggert algorithm, while for the other dp-algorithms we only need one
// column vector to compute the scores. The default dp-matrix specialization
// can be selected using a special meta-function.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TAlgorithm>
struct DefaultScoreMatrixSpec_;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag MatrixMember
// ----------------------------------------------------------------------------

struct DPMatrixMember_;
typedef Tag<DPMatrixMember_> DPMatrixMember;

// ----------------------------------------------------------------------------
// Tag SparseDPMatrix
// ----------------------------------------------------------------------------

struct SparseDPMatrix_;
typedef Tag<SparseDPMatrix_> SparseDPMatrix;

// ----------------------------------------------------------------------------
// Tag FullDPMatrix
// ----------------------------------------------------------------------------

struct FullDPMatrix_;
typedef Tag<FullDPMatrix_> FullDPMatrix;

// ----------------------------------------------------------------------------
// Enum DPMatrixDimension
// ----------------------------------------------------------------------------

// Used to globally specify the correct dimension and the correct size of
// dimension for the dp matrix.
struct DPMatrixDimension_
{
    typedef unsigned int TValue;

    static const TValue VERTICAL = 0u;
    static const TValue HORIZONTAL = 1u;
    static const TValue DIMENSION = 2u;
};

// ----------------------------------------------------------------------------
// Class DPMatrix_
// ----------------------------------------------------------------------------

// The dp matrix used as a score matrix and as a trace-back matrix.
template <typename TValue, typename TMatrixSpec, typename THost = String<TValue> >
class DPMatrix_
{};


// Default dp matrix implementation stores all cells of the dp matrix in the
// underlying two-dimensional matrix.
template <typename TValue, typename THost>
class DPMatrix_<TValue, FullDPMatrix, THost>
{
public:

    typedef typename Member<DPMatrix_, DPMatrixMember>::Type TMatrix;

    Holder<TMatrix>   data_host;  // The host containing the actual matrix.

    DPMatrix_() :
        data_host()
    {
        create(data_host);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultScoreMatrixSpec_
// ----------------------------------------------------------------------------

// This meta-function determines the default specialization of dp matrix
// based on the given algorithm tag.
template <typename TAlgorithm>
struct DefaultScoreMatrixSpec_
{
    typedef SparseDPMatrix Type;
};

// TODO(rmaerker): Move to WatermanEggert implementation?
template <>
struct DefaultScoreMatrixSpec_<LocalAlignment_<WatermanEggert> >
{
    typedef FullDPMatrix Type;
};

// ----------------------------------------------------------------------------
// Metafunction DataHost_
// ----------------------------------------------------------------------------

// Returns the type of the underlying matrix.
template <typename TValue, typename TMatrixSpec, typename THost>
struct Member<DPMatrix_<TValue, TMatrixSpec, THost>, DPMatrixMember>
{
    typedef Matrix<TValue, 2, THost> Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct Member<DPMatrix_<TValue, TMatrixSpec, THost> const, DPMatrixMember>
{
    typedef Matrix<TValue, 2, THost> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction SizeArr_
// ----------------------------------------------------------------------------

// Returns the type of the containers to store the dimensions and the factors
// in order to move properly in the matrix.
template <typename TDPMatrix>
struct SizeArr_ {};

template <typename TValue, typename TMatrixSpec, typename THost>
struct SizeArr_<DPMatrix_<TValue, TMatrixSpec, THost> >
{
    typedef DPMatrix_<TValue, TMatrixSpec, THost> TDPMatrix_;
    typedef typename Member<TDPMatrix_, DPMatrixMember>::Type TDataHost_;
    typedef typename SizeArr_<TDataHost_>::Type Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct SizeArr_<DPMatrix_<TValue, TMatrixSpec, THost> const>
{
    typedef DPMatrix_<TValue, TMatrixSpec, THost> TDPMatrix_;
    typedef typename Member<TDPMatrix_, DPMatrixMember>::Type TDataHost_;
    typedef typename SizeArr_<TDataHost_>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
struct Spec<DPMatrix_<TValue, TMatrixSpec, THost> >
{
    typedef TMatrixSpec Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct Spec<DPMatrix_<TValue, TMatrixSpec, THost> const>:
    Spec<DPMatrix_<TValue, TMatrixSpec, THost> >{};


// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
struct Value<DPMatrix_<TValue, TMatrixSpec, THost> >
{
    typedef TValue Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct Value<DPMatrix_<TValue, TMatrixSpec, THost> const>
{
    typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
struct Size<DPMatrix_<TValue, TMatrixSpec, THost> >
{
    typedef typename DPMatrix_<TValue, TMatrixSpec, THost>::TMatrix TDataMatrix_;
    typedef typename Size<TDataMatrix_>::Type Type;
};


// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
struct Host<DPMatrix_<TValue, TMatrixSpec, THost> >
{
    typedef THost Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct Host<DPMatrix_<TValue, TMatrixSpec, THost> const>
{
    typedef THost const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

// There are two iterator types. The standard iterator returns a
// non-rooted iterator to the underlying vector of the hosted two-dimensional
// matrix. The rooted iterator returns the iterator defined by the
// hosted matrix object which is a position iterator.
template <typename TValue, typename TMatrixSpec, typename THost>
struct Iterator<DPMatrix_<TValue, TMatrixSpec, THost>, Standard>
{
    typedef DPMatrix_<TValue, TMatrixSpec, THost> TDPMatrix_;
    typedef typename  Host<TDPMatrix_>::Type THost_;
    typedef typename Iterator<THost_, Standard>::Type Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct Iterator<DPMatrix_<TValue, TMatrixSpec, THost> const, Standard>
{
    typedef DPMatrix_<TValue, TMatrixSpec, THost> const TDPMatrix_;
    typedef typename  Host<TDPMatrix_>::Type THost_;
    typedef typename Iterator<THost_ const, Standard>::Type Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct Iterator<DPMatrix_<TValue, TMatrixSpec, THost>, Rooted>
{
    typedef DPMatrix_<TValue, TMatrixSpec, THost> TDPMatrix_;
    typedef typename Member<TDPMatrix_, DPMatrixMember>::Type TDataMatrix_;
    typedef typename Iterator<TDataMatrix_, Rooted>::Type Type;
};

template <typename TValue, typename TMatrixSpec, typename THost>
struct Iterator<DPMatrix_<TValue, TMatrixSpec, THost> const, Rooted>
{
    typedef DPMatrix_<TValue, TMatrixSpec, THost> TDPMatrix_;
    typedef typename Member<TDPMatrix_, DPMatrixMember>::Type TDataMatrix_;
    typedef typename Iterator<TDataMatrix_ const, Rooted>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _checkCorrectDimension()
// ----------------------------------------------------------------------------

// Checks whether a given value applies to the correct dimension.
inline bool _checkCorrectDimension(DPMatrixDimension_::TValue dim)
{
    return dim < DPMatrixDimension_::DIMENSION;
}

// ----------------------------------------------------------------------------
// Function _dataHost()
// ----------------------------------------------------------------------------

// Returns a reference to the hosted matrix.
template <typename TValue, typename TMatrixSpec, typename THost>
inline Holder<typename Host<DPMatrix_<TValue, TMatrixSpec, THost> >::Type> &
_dataHost(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix)
{
    return _dataHost(value(dpMatrix.data_host));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline Holder<typename Host<DPMatrix_<TValue, TMatrixSpec, THost> >::Type> const &
_dataHost(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix)
{
    return _dataHost(value(dpMatrix.data_host));
}

// ----------------------------------------------------------------------------
// Function _dataLengths()
// ----------------------------------------------------------------------------

// Returns a reference to the _dataLengths container of the hosted matrix.
template <typename TValue, typename TMatrixSpec, typename THost>
inline typename SizeArr_<DPMatrix_<TValue, TMatrixSpec, THost> >::Type &
_dataLengths(DPMatrix_<TValue, TMatrixSpec, THost>&dpMatrix)
{
    return _dataLengths(value(dpMatrix.data_host));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename SizeArr_<DPMatrix_<TValue, TMatrixSpec, THost> const>::Type &
_dataLengths(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix)
{
    return _dataLengths(value(dpMatrix.data_host));
}

// ----------------------------------------------------------------------------
// Function _dataFactors()
// ----------------------------------------------------------------------------

// Returns a reference to the _dataFactors container of the hosted matrix.
template <typename TValue, typename TMatrixSpec, typename THost>
inline typename SizeArr_<DPMatrix_<TValue, TMatrixSpec, THost> >::Type &
_dataFactors(DPMatrix_<TValue, TMatrixSpec, THost>&dpMatrix)
{
    return _dataFactors(value(dpMatrix.data_host));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename SizeArr_<DPMatrix_<TValue, TMatrixSpec, THost> const>::Type &
_dataFactors(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix)
{
    return _dataFactors(value(dpMatrix.data_host));
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

// Returns the value of the matrix at the given host position.
template <typename TValue, typename TMatrixSpec, typename THost, typename TPosition>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec, THost> >::Type
value(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix,
      TPosition const & pos)
{
    return value(value(dpMatrix.data_host), pos);
}

template <typename TValue, typename TMatrixSpec, typename THost, typename TPosition>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec, THost> const>::Type
value(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix,
      TPosition const & pos)
{
    return value(value(dpMatrix.data_host), pos);
}

// Returns the value of the matrix at the two given coordinates.
template <typename TValue, typename TMatrixSpec, typename THost, typename TPositionV, typename TPositionH>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec, THost> >::Type
value(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix,
      TPositionV const & posDimV,
      TPositionH const & posDimH)
{
    return value(value(dpMatrix.data_host), posDimV, posDimH);
}

template <typename TValue, typename TMatrixSpec, typename THost, typename TPositionV, typename TPositionH>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec, THost> const>::Type
value(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix,
      TPositionV const & posDimV,
      TPositionH const & posDimH)
{
    return value(value(dpMatrix.data_host), posDimV, posDimH);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// Returns the length of a given dimension.
template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Size<DPMatrix_<TValue, TMatrixSpec, THost> const>::Type
length(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix,
       unsigned int dimension)
{
    SEQAN_ASSERT(_checkCorrectDimension(dimension));

    return length(value(dpMatrix.data_host), dimension);
}

// Returns the overall length of the underlying vector of the hosted matrix.
template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Size<DPMatrix_<TValue, TMatrixSpec, THost> const>::Type
length(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix)
{
    return length(value(dpMatrix.data_host));  // Note that even if the dimensional lengths are set but the matrix was not resized
    // this function returns 0 or the previous length of the host before the resize.
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
inline void
clear(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix)
{
    clear(_dataLengths(dpMatrix));
    resize(_dataLengths(dpMatrix), 2, 0);
    clear(_dataFactors(dpMatrix));
    resize(_dataFactors(dpMatrix), 2, 0);
    _dataFactors(dpMatrix)[DPMatrixDimension_::VERTICAL] = 1u;
    clear(host(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
inline bool
empty(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix)
{
    return empty(host(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function setLength()
// ----------------------------------------------------------------------------

// Sets the new length of a given dimension.
template <typename TValue, typename TMatrixSpec, typename THost, typename TSize>
inline void
setLength(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix,
          unsigned int dimension,
          TSize const & newLength)
{
    SEQAN_ASSERT(_checkCorrectDimension(dimension));
    setLength(value(dpMatrix.data_host), dimension, newLength);
}

// ----------------------------------------------------------------------------
// Function updateFactors()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Size<DPMatrix_<TValue, TMatrixSpec, THost> >::Type
updateFactors(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix)
{
    typedef typename Size<DPMatrix_<TValue, TMatrixSpec, THost> >::Type TSize;

    TSize factor_ = _dataFactors(dpMatrix)[0] * length(dpMatrix, 0);
    for (unsigned int i = 1; (factor_ > 0) && (i < dimension(value(dpMatrix.data_host))); ++i)
    {
        _dataFactors(dpMatrix)[i] = factor_;
        factor_ *= length(dpMatrix, i);
    }
    return factor_;
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

// Resizes the matrix. Note, the lengths of the dimensions have to be set before.
template <typename TValue, typename TMatrixSpec, typename THost>
inline void
resize(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix)
{
    typedef typename Size<DPMatrix_<TValue, TMatrixSpec, THost> >::Type TSize;

    TSize reqSize = updateFactors(dpMatrix);
    if (reqSize > length(dpMatrix))
        resize(host(dpMatrix), reqSize, Exact());
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline void
resize(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix,
       TValue const & fillValue)
{
    typedef typename Size<DPMatrix_<TValue, TMatrixSpec, THost> >::Type TSize;

    TSize reqSize = updateFactors(dpMatrix);
    if (reqSize > length(dpMatrix))
        resize(host(dpMatrix), reqSize, fillValue, Exact());
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost>, Standard>::Type
begin(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix, Standard)
{
    return begin(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost> const, Standard>::Type
begin(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix, Standard)
{
    return begin(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost>, Rooted>::Type
begin(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix, Rooted)
{
    return begin(value(dpMatrix.data_host));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost> const, Rooted>::Type
begin(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix, Rooted)
{
    return begin(value(dpMatrix.data_host));
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost>, Standard>::Type
end(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix, Standard)
{
    return end(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost> const, Standard>::Type
end(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix, Standard)
{
    return end(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost>, Rooted>::Type
end(DPMatrix_<TValue, TMatrixSpec, THost> & dpMatrix, Rooted)
{
    return end(value(dpMatrix.data_host));
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec, THost>, Rooted>::Type
end(DPMatrix_<TValue, TMatrixSpec, THost> const & dpMatrix, Rooted)
{
    return end(value(dpMatrix.data_host));
}

// ----------------------------------------------------------------------------
// Function coordinate()
// ----------------------------------------------------------------------------

// Returns the coordinate of a host positio in a given dimension.
template <typename TValue, typename THost, typename TPosition>
inline typename Position<DPMatrix_<TValue, FullDPMatrix, THost> >::Type
coordinate(DPMatrix_<TValue, FullDPMatrix, THost> const & dpMatrix,
           TPosition hostPos,
           typename DPMatrixDimension_::TValue dimension)
{
    return coordinate(value(dpMatrix.data_host), hostPos, dimension);
}

// ----------------------------------------------------------------------------
// Function toGlobalPosition()
// ----------------------------------------------------------------------------

// Returns the current position of the navigator within the matrix.
template <typename TValue, typename THost,
          typename TPosH,
          typename TPosV>
inline typename Position<DPMatrix_<TValue, FullDPMatrix, THost> >::Type
toGlobalPosition(DPMatrix_<TValue, FullDPMatrix, THost> const & dpMatrix,
                 TPosH const horizontalCoordinate,
                 TPosV const verticalCoordinate)
{
    return horizontalCoordinate * length(dpMatrix, DPMatrixDimension_::VERTICAL) + verticalCoordinate;
}

} // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_H_
