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
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================
// Operations on alignments such as integration
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_ALIGNMENT_OPERATIONS_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_ALIGNMENT_OPERATIONS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TSource0, typename TGapSpec0,
          typename TSource1, typename TGapSpec1,
          typename TPos>
inline void
integrateGaps(Gaps<TSource0, TGapSpec0> & targetRow,
              Gaps<TSource1, TGapSpec1> const & sourceRow,
              TPos const viewPos)
{
    typedef typename Iterator<Gaps<TSource0, TGapSpec0>, Standard>::Type TTargetIt;
    typedef typename Iterator<Gaps<TSource1, TGapSpec1> const, Standard>::Type TSourceIt;

    // This assertion ensures that the number of sequence characters after viewPos is greater than or equal to
    // the number of source characters in the clipped infix row.
    SEQAN_ASSERT_GEQ(endPosition(targetRow) - toSourcePosition(targetRow, viewPos),
                     endPosition(sourceRow) - beginPosition(sourceRow));

    // init iterators
    TTargetIt it = iter(targetRow, viewPos);

    // walk through Gaps containers and copy gaps
    for (TSourceIt sIt = begin(sourceRow, Standard()), sItEnd = end(sourceRow, Standard()); sIt != sItEnd;)
    {
        TPos gapSize = countGaps(sIt);
        insertGaps(it, gapSize);
        goFurther(it, gapSize+1);
        goFurther(sIt, gapSize+1);
    }
}

template <typename TSource0, typename TGapSpec0,
          typename TSource1, typename TGapSpec1>
inline void
integrateGaps(Gaps<TSource0, TGapSpec0> & targetRow,
              Gaps<TSource1, TGapSpec1> const & sourceRow)
{
    typename Position<TSource0>::Type viewPos = beginPosition(source(sourceRow)) // correct for infixes
                                              - beginPosition(source(targetRow)) // ...
                                              + beginPosition(sourceRow);        // respect source clipping

    integrateGaps(targetRow, sourceRow, toViewPosition(targetRow, viewPos));
}

// ----------------------------------------------------------------------------
// Function integrateAlign()
// ----------------------------------------------------------------------------

/*!
 * @fn integrateAlign
 * @headerfile <seqan/align.h>
 * @brief Integrates an alignment into another by copying the gaps.
 *
 * @signature void integrateAlign(align1, align2[, positions]);
 *
 * @param[in,out] align1    Target Alignment object into which align2 is to be integrated.
 * @param[in]     align2    Alignment object that is to be integrated into align1.
 * @param[in]     positions The integration positions in align1 for all rows (view positions), String of positions.
 *
 * @section Examples
 *
 * @include demos/dox/align/integrate_align.cpp
 *
 * The output is as follows:
 *
 * @include demos/dox/align/integrate_align.cpp.stdout
 */

template <typename TSource1, typename TSpec1, typename TSource2, typename TSpec2, typename TPos>
void integrateAlign(Align<TSource1, TSpec1> & align,
                    Align<TSource2, TSpec2> const & infixAlign,
                    String<TPos> const & viewPos)
{
    SEQAN_ASSERT_EQ_MSG(length(rows(infixAlign)), length(rows(align)), "Both align objects need same number of rows.");
    typedef typename Size<Align<TSource1, TSpec1> >::Type TSize;
    //NOTE(h-2): could be parallelized
    for (TSize i = 0; i < length(rows(align)); ++i)
        integrateGaps(row(align, i), row(infixAlign, i), viewPos[i]);
}

template <typename TSource1, typename TSpec1, typename TSource2, typename TSpec2>
void integrateAlign(Align<TSource1, TSpec1> & align,
                    Align<TSource2, TSpec2> const & infixAlign)
{
    SEQAN_ASSERT_EQ_MSG(length(rows(infixAlign)), length(rows(align)), "Both align objects need same number of rows.");
    typedef typename Size<Align<TSource1, TSpec1> >::Type TSize;
    //NOTE(h-2): could be parallelized
    for (TSize i = 0; i < length(rows(align)); ++i)
        integrateGaps(row(align, i), row(infixAlign, i));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_ALIGNMENT_OPERATIONS_H_
