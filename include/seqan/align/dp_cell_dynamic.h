// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Implements the dynamic gap model published in "Dynamic Gaps Selector:
// A Smith Waterman Sequence Alignment Algorithm with Affine Gap Model
// Optimization" by Gianvito Urgese et al.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_DP_CELL_DYNAMIC_H_
#define INCLUDE_SEQAN_ALIGN_DP_CELL_DYNAMIC_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetFlagMaskType
// ----------------------------------------------------------------------------

namespace impl
{
namespace dp_cell
{
template <typename TScoreValue>
struct FlagMaskType
{
    using Type = typename If<Is<SimdVectorConcept<TScoreValue> >,
                             TScoreValue,
                             uint8_t>::Type;
};
}  // namespace dp_cell
}  // namespace impl

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct DynamicGapExtensionHorizontal_;
typedef Tag<DynamicGapExtensionHorizontal_> DynamicGapExtensionHorizontal;

struct DynamicGapExtensionVertical_;
typedef Tag<DynamicGapExtensionVertical_> DynamicGapExtensionVertical;

enum DynamicGapsMask
{
    MASK_VERTICAL_GAP = 1,
    MASK_HORIZONTAL_GAP = 2
};

// ----------------------------------------------------------------------------
// Class DPCell                                                   [DynamicGaps]
// ----------------------------------------------------------------------------

// The specialization for linear gap cost function.
// It solely stores the maximal score.
template <typename TScoreValue>
class DPCell_<TScoreValue, DynamicGaps>
{
public:
    using TFlagMaskType = typename impl::dp_cell::FlagMaskType<TScoreValue>::Type;

    TScoreValue     _score      = DPCellDefaultInfinity<DPCell_>::VALUE;
    TFlagMaskType   _flagMask   = TFlagMaskType();

    DPCell_() = default;
    
    // Copy c'tor.
    DPCell_(DPCell_ const & other) : _score(other._score), _flagMask(other._flagMask)
    {}

    // Move c'tor.
    DPCell_(DPCell_ && other) : DPCell_()
    {
        swap(*this, other);
    }

    // Construct with score.
    DPCell_(TScoreValue const & pScore) : _score(pScore)
    {}
    
    // Assignment and move operator.
    DPCell_& operator=(DPCell_ other)
    {
        swap(*this, other);
        return *this;
    }

    // Assignment of score.
    DPCell_ &
    operator=(TScoreValue const & score)
    {
        _score = score;
        return *this;
    }
    
    ~DPCell_() = default;
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TCell, typename TBoolV, typename TBoolH>
struct SetGapExtension;

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, False, False>
{
    static const char VALUE = 0;
};

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, False, True>
{
    static const char VALUE = MASK_HORIZONTAL_GAP;
};

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, True, False>
{
    static const char VALUE = MASK_VERTICAL_GAP;
};

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, True, True>
{
    static const char VALUE = MASK_HORIZONTAL_GAP | MASK_VERTICAL_GAP;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TScoreValue, typename TFlag>
inline void _setBit(DPCell_<TScoreValue, DynamicGaps> & cell,
                    TFlag const & /*flag*/,
                    DynamicGapExtensionVertical const & /*tag*/)
{
    if (IsSameType<TFlag, True>::VALUE)
        cell._flagMask |= MASK_VERTICAL_GAP;
    else
        cell._flagMask &= ~MASK_VERTICAL_GAP;
}

template <typename TScoreValue, typename TFlag>
inline void _setBit(DPCell_<TScoreValue, DynamicGaps> & cell,
                    TFlag const & /*flag*/,
                    DynamicGapExtensionHorizontal const & /*tag*/)
{
    if (IsSameType<TFlag, True>::VALUE)
        cell._flagMask |= MASK_HORIZONTAL_GAP;
    else
        cell._flagMask &= ~MASK_HORIZONTAL_GAP;
}

template <typename TScoreValue, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >,bool)
isGapExtension(DPCell_<TScoreValue, DynamicGaps> const & cell,
               TSpec const & /*spec*/)
{
    if (IsSameType<TSpec, DynamicGapExtensionHorizontal>::VALUE)
        return cell._flagMask & MASK_HORIZONTAL_GAP;
    else
        return cell._flagMask & MASK_VERTICAL_GAP;
}

template <typename TScoreValue, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >,TScoreValue)
isGapExtension(DPCell_<TScoreValue, DynamicGaps> const & cell,
               TSpec const & /*spec*/)
{
    return blend(cell._flagMask & createVector<TScoreValue>(MASK_VERTICAL_GAP),
                 cell._flagMask & createVector<TScoreValue>(MASK_HORIZONTAL_GAP),
                 createVector<TScoreValue>(IsSameType<TSpec, DynamicGapExtensionHorizontal>::VALUE));
}

template <typename TScoreValue, typename TFlagV, typename TFlagH>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, void)
setGapExtension(DPCell_<TScoreValue, DynamicGaps> & cell,
                TFlagV const & /*vert*/, TFlagH const & /*hori*/)
{
    cell._flagMask = SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, TFlagV, TFlagH>::VALUE;
}

template <typename TScoreValue, typename TFlagV, typename TFlagH>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, void)
setGapExtension(DPCell_<TScoreValue, DynamicGaps> & cell,
                TFlagV const & /*vert*/, TFlagH const & /*hori*/,
                TScoreValue const & cmp)
{
    cell._flagMask = blend(cell._flagMask, createVector<TScoreValue>(SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, TFlagV, TFlagH>::VALUE), cmp);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

// Needed for banded chain alignment for the std::set.
template <typename TScoreValueLeft, typename TScoreValueRight>
inline bool operator<(DPCell_<TScoreValueLeft, DynamicGaps> const & left,
                      DPCell_<TScoreValueRight, DynamicGaps> const & right)
{
    return left._score < right._score;
}

template <typename TScoreValue>
inline void
swap(DPCell_<TScoreValue, DynamicGaps> & lhs,
     DPCell_<TScoreValue, DynamicGaps> & rhs)
{
    std::swap(lhs._score, rhs._score);
    std::swap(lhs._flagMask, rhs._flagMask);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_DP_CELL_DYNAMIC_H_
