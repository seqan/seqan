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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Split alignment implementation.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_SPLIT_ALIGN_SPLIT_INTERFACE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_SPLIT_ALIGN_SPLIT_INTERFACE_H_

#include "dp_scout_split.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Tag for the split alignment algorithm in DPProfile_.

template <typename TSpec = Default>
struct SplitAlignment_ {};

// Tag for the split alignment algorithm.

struct SplitAlignmentAlgo_;
typedef Tag<SplitAlignmentAlgo_> SplitAlignmentAlgo;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForAlignmentAlgorithm_
// ----------------------------------------------------------------------------

// For the split alignment, we will use our SplitAlignmentScout specialization of DPScout.

template <typename TSpec>
struct ScoutSpecForAlignmentAlgorithm_<SplitAlignment_<TSpec>, DPScoutState_<SplitAlignmentScout> >
{
    typedef SplitAlignmentScout Type;
};

template <typename TSpec>
struct ScoutSpecForAlignmentAlgorithm_<SplitAlignment_<TSpec> const, DPScoutState_<SplitAlignmentScout> >
{
    typedef SplitAlignmentScout Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsSplitAlignment_
// ----------------------------------------------------------------------------

// Convenience function used in the DP configuration below.

template <typename TParam>
struct IsSplitAlignment_ : False {};

template <typename TSpec>
struct IsSplitAlignment_<SplitAlignment_<TSpec> >:
    True {};

template <typename TSpec>
struct IsSplitAlignment_<SplitAlignment_<TSpec> const>:
    True {};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsSplitAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> >:
    IsSplitAlignment_<TAlgoSpec> {};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsSplitAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> const>:
    IsSplitAlignment_<TAlgoSpec> {};

// ----------------------------------------------------------------------------
// Metafunction IsFreeEndGap_
// ----------------------------------------------------------------------------

// We want the same free endgaps configuration as for global alignments.

template <typename TSpec, typename TRow>
struct IsFreeEndGap_<SplitAlignment_<TSpec>, TRow> :
        IsFreeEndGap_<GlobalAlignment_<TSpec>, TRow>
{};

template <typename TSpec, typename TRow>
struct IsFreeEndGap_<SplitAlignment_<TSpec> const, TRow> :
        IsFreeEndGap_<GlobalAlignment_<TSpec> const, TRow>
{};

// ----------------------------------------------------------------------------
// Metafunction IsGlobalAlignment_
// ----------------------------------------------------------------------------

// We use similar functionality as the global alignment.

template <typename TSpec>
struct IsGlobalAlignment_<SplitAlignment_<TSpec> > :
    True
{};

template <typename TSpec>
struct IsGlobalAlignment_<SplitAlignment_<TSpec> const> :
    True
{};

// ----------------------------------------------------------------------------
// Metafunction TraceTail_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct TraceTail_<SplitAlignment_<TSpec> > :
    False
{};

// ----------------------------------------------------------------------------
// Metafunction LastColumnEnabled_
// ----------------------------------------------------------------------------

template <typename TSpec, typename TColumnDescriptor>
struct LastColumnEnabled_<SplitAlignment_<TSpec>, TColumnDescriptor>
{
    typedef typename If<IsSameType<typename TColumnDescriptor::TColumnProperty, DPFinalColumn>,
                        typename IsFreeEndGap_<SplitAlignment_<TSpec>, DPLastColumn>::Type,
                        False>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction LastRowEnabled_
// ----------------------------------------------------------------------------

template <typename TSpec, typename TColumnDescriptor>
struct LastRowEnabled_<SplitAlignment_<TSpec>, LastCell, TColumnDescriptor>
{
    typedef typename  If<Or<IsSameType<typename TColumnDescriptor::TLocation, PartialColumnBottom>,
                            IsSameType<typename TColumnDescriptor::TLocation, FullColumn> >,
                         typename IsFreeEndGap_<SplitAlignment_<TSpec>, DPLastRow>::Type,
                         False>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction DPMetaColumn_
// ----------------------------------------------------------------------------

template <typename TSpec, typename TGapCosts, typename TTraceFlag, typename TColumnType>
struct DPMetaColumn_<DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag>, MetaColumnDescriptor<TColumnType, FullColumn> >
{
    typedef DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag> TDPProfile;
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> Horizontal | Zero, All, All
    // If FinalColumn -> Horizontal | Zero, All, All

    typedef typename If<Or<IsSameType<TColumnType, DPInitialColumn>,
                           IsFreeEndGap_<TDPProfile, DPFirstRow> >, RecursionDirectionZero, RecursionDirectionHorizontal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};


template <typename TSpec, typename TGapCosts, typename TTraceFlag, typename TColumnType>
struct DPMetaColumn_<DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag>, MetaColumnDescriptor<TColumnType, PartialColumnTop> >
{
    typedef DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag> TDPProfile;
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // How does the recursion directions look like?

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> Horizontal | Zero, All, LowerBand
    // If FinalColumn -> Horizontal | Zero, All, LowerBand

    typedef typename If<Or<IsSameType<TColumnType, DPInitialColumn>,
                           IsFreeEndGap_<TDPProfile, DPFirstRow> >, RecursionDirectionZero, RecursionDirectionHorizontal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionLowerDiagonal>::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

template <typename TSpec, typename TGapCosts, typename TTraceFlag, typename TColumnType>
struct DPMetaColumn_<DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag>, MetaColumnDescriptor<TColumnType, PartialColumnMiddle> >
{
    typedef DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag> TDPProfile;
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> UpperDiagonal, All, LowerDiagonal
    // If FinalColumn -> UpperDiagonal, All, LowerDiagonal

    typedef typename If<IsSameType<TColumnType, DPInitialColumn>, RecursionDirectionZero, RecursionDirectionUpperDiagonal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionLowerDiagonal>::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

template <typename TSpec, typename TGapCosts, typename TTraceFlag, typename TColumnType>
struct DPMetaColumn_<DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag>, MetaColumnDescriptor<TColumnType, PartialColumnBottom> >
{
    typedef DPProfile_<SplitAlignment_<TSpec>, TGapCosts, TTraceFlag> TDPProfile;
    typedef typename IsLocalAlignment_<TDPProfile>::Type TIsLocal;

    // If InitialColumn -> Zero, Vertical | Zero, Vertical | Zero  // Within the algorithm we need to define the first row as only one cell if it is no initial column
    // If InnerColumn -> UpperDiagonal, All, All
    // If FinalColumn -> UpperDiagonal, All, All

    typedef typename If<IsSameType<TColumnType, DPInitialColumn>, RecursionDirectionZero, RecursionDirectionUpperDiagonal>::Type TRecursionTypeFirstCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeInnerCell_;
    typedef typename If<IsSameType<TColumnType, DPInitialColumn>,
                        typename If<IsFreeEndGap_<TDPProfile, DPFirstColumn>, RecursionDirectionZero, RecursionDirectionVertical>::Type,
                        RecursionDirectionAll>::Type TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

// ----------------------------------------------------------------------------
// Metafunction SetupAlignmentProfile_
// ----------------------------------------------------------------------------

template <typename TFreeEndGaps, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<SplitAlignmentAlgo, TFreeEndGaps, TGapCosts, TTraceSwitch>
{
    typedef DPProfile_<SplitAlignment_<TFreeEndGaps>, TGapCosts, TTraceSwitch> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _reverseTrace()
// ----------------------------------------------------------------------------

// Reverse a trace string and adapt internal position.
template <typename TPosition, typename TSize, typename TSpec>
void _reverseTrace(String<TraceSegment_<TPosition, TSize>, TSpec> & trace,
                   size_t const lengthH,
                   size_t const lengthV)
{
    typedef String<TraceSegment_<TPosition, TSize>, TSpec> TTrace;
    typedef typename Iterator<TTrace, Rooted>::Type TTraceIter;

    if (empty(trace))
        return;

    for (TTraceIter it = begin(trace, Rooted()); !atEnd(it); goNext(it))
    {
        it->_horizontalBeginPos = lengthH - _getEndHorizontal(*it);
        it->_verticalBeginPos = lengthV - _getEndVertical(*it);
    }
    reverse(trace);
}

// ----------------------------------------------------------------------------
// Function _computeSplitTrace()
// ----------------------------------------------------------------------------

template <typename TTarget,
          typename TSeqH,
          typename TSeqV,
          typename TDPContext,
          typename TMatPos,
          typename TDPType, typename TBandSwitch, typename TFreeEndGaps, typename TTraceConfig>
void _computeSplitTrace(TTarget & target,
                        TSeqH const & seqH,
                        TSeqV const & seqV,
                        TDPContext const & dpContext,
                        TMatPos const matPos,
                        AlignConfig2<TDPType, DPBandConfig<TBandSwitch>, TFreeEndGaps, TTraceConfig> const & config)
{
    typedef typename SetupAlignmentProfile_<TDPType, TFreeEndGaps, LinearGaps, TTraceConfig>::Type TDPProfile;

    typedef typename GetDPTraceMatrix<TDPContext const>::Type TDPTraceMatrixHost;
    typedef typename Value<TDPTraceMatrixHost>::Type TTraceValue;

    typedef DPMatrix_<TTraceValue, FullDPMatrix> TDPTraceMatrix;
    typedef DPMatrixNavigator_<TDPTraceMatrix, DPTraceMatrix<TTraceConfig>, NavigateColumnWise> TDPTraceMatrixNavigator;

    TDPTraceMatrix matrix;
    setLength(matrix, +DPMatrixDimension_::HORIZONTAL, length(seqH) + 1 - std::max(0, lowerDiagonal(config._band)));

    if (IsSameType<TBandSwitch, BandOff>::VALUE)
    {
        setLength(matrix, +DPMatrixDimension_::VERTICAL, length(seqV) + 1);
    }
    else
    {
        int bandSize = _min(static_cast<int>(length(seqH)), upperDiagonal(config._band)) -
                       _max(lowerDiagonal(config._band), -static_cast<int>(length(seqV))) + 1;
        setLength(matrix, +DPMatrixDimension_::VERTICAL, _min(static_cast<int>(length(seqV)) + 1, bandSize));
    }

    setHost(matrix, getDpTraceMatrix(dpContext));
    resize(matrix);
    SEQAN_ASSERT_EQ(length(getDpTraceMatrix(dpContext)), length(matrix));
    TDPTraceMatrixNavigator navi;
    _init(navi, matrix, config._band);
    _computeTraceback(target, navi, matPos, seqH, seqV, config._band, TDPProfile());
}

// ----------------------------------------------------------------------------
// Function _splitAlignmentImpl()
// ----------------------------------------------------------------------------

// We call the long sequence contig and the shorter one read but could be changed roles.
template <typename TContigSeqL,
          typename TReadSeqL,
          typename TContigSeqR,
          typename TReadSeqR,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfigL,
          typename TAlignConfigR,
          typename TGapModel>
auto _splitAlignmentImpl(Gaps<TContigSeqL> & gapsContigL,
                         Gaps<TReadSeqL> & gapsReadL,
                         Gaps<TContigSeqR> & gapsContigR,
                         Gaps<TReadSeqR> & gapsReadR,
                         Score<TScoreValue, TScoreSpec> const & scoringScheme,
                         TAlignConfigL const & alignConfigL,
                         TAlignConfigR const & alignConfigR,
                         TGapModel const & /*gapModel*/)
{
    typedef Gaps<TContigSeqL> TGaps;
    typedef typename Size<TGaps>::Type TSize;
    typedef typename Position<TGaps>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    // Compute trace and split score sequence for the left alignment.
    // We actually need to first compute the scores, than trace from the choosen split position.
    DPContext<TScoreValue, TGapModel> dpContextL;
    DPScoutState_<SplitAlignmentScout> scoutStateL;
    resize(scoutStateL.splitScore, length(source(gapsContigL)) + 1, minValue<TScoreValue>() / 2);
    resize(scoutStateL.splitPos, length(scoutStateL.splitScore));

    String<TTraceSegment> traceL;
    _setUpAndRunAlignment(dpContextL, traceL, scoutStateL, source(gapsContigL), source(gapsReadL), scoringScheme, alignConfigL);

    // Get reversed versions of the right contig and read sequence.
    ModifiedString<TContigSeqR, ModReverse> revContigR(source(gapsContigR));
    ModifiedString<TReadSeqR, ModReverse> revReadR(source(gapsReadR));

    // Compute trace and split score sequence for the right alignment.
    DPContext<TScoreValue, TGapModel> dpContextR;
    DPScoutState_<SplitAlignmentScout> scoutStateR;
    resize(scoutStateR.splitScore, length(source(gapsContigR)) + 1, minValue<TScoreValue>() / 2);
    resize(scoutStateR.splitPos, length(scoutStateR.splitScore));

    String<TTraceSegment> traceR;
    _setUpAndRunAlignment(dpContextR, traceR, scoutStateR, revContigR, revReadR, scoringScheme, alignConfigR);

    // Reverse trace so it fits to the forward right sequences.  Also reverse the trace such that we can directly apply
    // it for the right alignment.
    reverse(scoutStateR.splitScore);
    reverse(scoutStateR.splitPos);

    SEQAN_ASSERT_EQ(length(scoutStateL.splitScore), length(scoutStateR.splitScore));

    // We will split the left and right alignments into two parts such that the alignment score is optimal.  We compute
    // the leftmost best position for a split (equivalent to the best prefix of the first left alignment).  Note that
    // placing the breakpoint at the leftmost position is coherent with the SNPdb semantics but there are other data
    // bases that use rightmost placement.

    // TODO(holtgrew): Make selecting the left/right split position from interface possible? Maybe not necessary.

    auto itBegin = makeZipIterator(begin(scoutStateL.splitScore), begin(scoutStateR.splitScore));
    auto res = std::max_element(itBegin, makeZipIterator(end(scoutStateL.splitScore), end(scoutStateR.splitScore)),
                                [] (auto const & lhs, auto const & rhs)
                                {
                                    return std::get<0>(lhs) + std::get<1>(lhs) < std::get<0>(rhs) + std::get<1>(rhs);
                                });
    auto bestPrefixLength = res - itBegin;

    // std::cerr << "bestPrefixLength = " << bestPrefixLength << "\n";

    // std::cerr << "split store left ";
    // for (unsigned i = 0; i < length(scoutStateL.splitScore); ++i)
    //     fprintf(stderr, " %3d", scoutStateL.splitScore[i]);
    // std::cerr << "\n";
    // std::cerr << "split store right";
    // for (unsigned i = 0; i < length(scoutStateR.splitScore); ++i)
    //     fprintf(stderr, " %3d", scoutStateR.splitScore[i]);
    // std::cerr << "\n";

    // Recompute the best trace starting from the recorded split position for left ...
    clear(traceL);
    _computeSplitTrace(traceL, source(gapsContigL), source(gapsReadL), dpContextL,
                       scoutStateL.splitPos[bestPrefixLength], alignConfigL);
    _adaptTraceSegmentsTo(gapsContigL, gapsReadL, traceL);
    // ... and right anchor.
    clear(traceR);
    _computeSplitTrace(traceR, source(gapsContigR), source(gapsReadR), dpContextR,
                       scoutStateR.splitPos[bestPrefixLength], alignConfigR);
    
    _reverseTrace(traceR, length(source(gapsContigR)), length(source(gapsReadR)));
    _adaptTraceSegmentsTo(gapsContigR, gapsReadR, traceR);

    // We have to correct the clipping position for the left alignment because of the to-right projection. The
    // insertion itself is not part of the alignment.
    TPosition cePosL = toViewPosition(gapsContigL, bestPrefixLength);
    if (bestPrefixLength > 0)
        cePosL = toViewPosition(gapsContigL, bestPrefixLength - 1) + 1;
    setClippedEndPosition(gapsContigL, cePosL);
    setClippedEndPosition(gapsReadL, cePosL);

    return std::make_pair(std::get<0>(*res), std::get<1>(*res));
}

template <typename TContigSeqL, typename TReadSeqL,
          typename TContigSeqR, typename TReadSeqR,
          typename TScoreValue, typename TScoreSpec,
          bool TTop, bool TRight, bool TLeft, bool TBottom, typename TConfigSpec,
          typename TGapModel>
auto _splitAlignmentImpl(Gaps<TContigSeqL> & gapsContigL,
                         Gaps<TReadSeqL> & gapsReadL,
                         Gaps<TContigSeqR> & gapsContigR,
                         Gaps<TReadSeqR> & gapsReadR,
                         Score<TScoreValue, TScoreSpec> const & scoringScheme,
                         AlignConfig<TTop, TRight, TLeft, TBottom, TConfigSpec> const &,
                         int lowerDiagonal,
                         int upperDiagonal,
                         TGapModel const & /*gapModel*/)
{
    typedef typename SubstituteAlignConfig_<AlignConfig<TTop, TRight, TLeft, TBottom> >::Type TFreeEndGaps;
    // Check whether we need to run the banded versions.
    if (lowerDiagonal != minValue<int>() && upperDiagonal != maxValue<int>())
    {
        typedef AlignConfig2<SplitAlignmentAlgo, DPBandConfig<BandOn>, TFreeEndGaps,
        TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TAlignConfigL;
        typedef AlignConfig2<SplitAlignmentAlgo, DPBandConfig<BandOn>, TFreeEndGaps,
        TracebackOn<TracebackConfig_<CompleteTrace, GapsRight> > > TAlignConfigR;
        return _splitAlignmentImpl(gapsContigL, gapsReadL, gapsContigR, gapsReadR, scoringScheme,
                                   TAlignConfigL(lowerDiagonal, upperDiagonal),
                                   TAlignConfigR(lowerDiagonal, upperDiagonal),
                                   TGapModel());
    }
    else
    {
        typedef AlignConfig2<SplitAlignmentAlgo, DPBandConfig<BandOff>, TFreeEndGaps,
        TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TAlignConfigL;
        typedef AlignConfig2<SplitAlignmentAlgo, DPBandConfig<BandOff>, TFreeEndGaps,
        TracebackOn<TracebackConfig_<CompleteTrace, GapsRight> > > TAlignConfigR;
        return _splitAlignmentImpl(gapsContigL, gapsReadL, gapsContigR, gapsReadR, scoringScheme, TAlignConfigL(),
                                   TAlignConfigR(), TGapModel());
    }
}

template <typename TContigSeqL, typename TReadSeqL, typename TContigSeqR, typename TReadSeqR,
          typename TScoreValue, typename TScoreSpec,
          bool TTop, bool TRight, bool TLeft, bool TBottom, typename TConfigSpec>
auto _splitAlignmentImpl(Gaps<TContigSeqL> & gapsContigL,
                         Gaps<TReadSeqL> & gapsReadL,
                         Gaps<TContigSeqR> & gapsContigR,
                         Gaps<TReadSeqR> & gapsReadR,
                         Score<TScoreValue, TScoreSpec> const & scoringScheme,
                         AlignConfig<TTop, TRight, TLeft, TBottom, TConfigSpec> const & config,
                         int lowerDiagonal = minValue<int>(),
                         int upperDiagonal = maxValue<int>())
{
    if (_usesAffineGaps(scoringScheme, source(gapsContigL), source(gapsReadL)))
        return _splitAlignmentImpl(gapsContigL, gapsReadL, gapsContigR, gapsReadR,
                                   scoringScheme, config, lowerDiagonal, upperDiagonal, AffineGaps());
    else
        return _splitAlignmentImpl(gapsContigL, gapsReadL, gapsContigR, gapsReadR,
                                   scoringScheme, config, lowerDiagonal, upperDiagonal, LinearGaps());
}

// ----------------------------------------------------------------------------
// Function splitAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn splitAlignment
 * @headerfile <seqan/align_split.h>
 * @brief Compute split alignments.
 *
 * @signature TScoreValue splitAlignment(alignL,         alignR,         scoringScheme[, config][, lowerDiag, upperDiag]);
 * @signature TScoreValue splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, scoringScheme[, config][, lowerDiag, upperDiag]);
 *
 * @param[in,out] alignL @link Align @endlink object with two rows for the left alignment.
 * @param[in,out] alignR @link Align @endlink object with two rows for the right alignment.
 * @param[in,out] gapsHL @link Gaps @endlink object with the horizontal/contig row for the left alignment.
 * @param[in,out] gapsVL @link Gaps @endlink object with the vertical/read row for the left alignment.
 * @param[in,out] gapsHR @link Gaps @endlink object with the horizontal/contig row for the right alignment.
 * @param[in,out] gapsVR @link Gaps @endlink object with the vertical/read row for the right alignment.
 * @param[in]     scoringScheme The scoring scheme to use for the alignment.
 * @param[in]     config A configuration object of type @link AlignConfig @endlink, to specify free-end-gaps.
 * @param[in]     lowerDiag The lower diagonal.You have to specify the upper and lower diagonals for the left
 *                          alignment.  For the right alignment, the corresponding diagonals are chosen for the
 *                          lower right part of the DP matrix, <tt>int</tt>.
 * @param[in]     upperDiag The lower diagonal.  Also see remark for <tt>lowerDiag</tt>, <tt>int</tt>.
 *
 * @return TScoreValue The sum of the alignment scores of both alignments (Metafunction: @link Score#Value @endlink
 *                     of the type of <tt>scoringScheme</tt>).
 *
 * There are two variants of the split alignment problem.  In the first variant, we want to align two sequences where the
 * first (say the reference) one is shorter than the second (say a read) and the read contains an insertion with respect
 * to the reference.  We now want to align the read agains the reference such that the left part of the read aligns well
 * against the left part of the reference and the right part of the read aligns well against the right part of the
 * reference.  The center gap in the reference is free.
 *
 * For example:
 *
 * @code{.console}
 * reference  AGCATGTTAGATAAGATAGC-----------TGTGCTAGTAGGCAGTCAGCGCCAT
 *            ||||||||||||||||||||           |||||||||||||||||||||||||
 * read       AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT
 * @endcode
 *
 * The second variant is to align two sequences A and B against a reference such that the left part of A aligns well to
 * the left part of the reference and the right part of B aligns well to the right part of the reference.  Together,
 * both reads span the whole reference and overlap with an insertion in the reference.
 *
 * @code{.console}
 * reference  AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
 *            |||||||||||||||||| | ||
 *            AGCATGTTAGATAAGATATCCGTCC
 *            read 1
 *                              ||| |||||||||||||||||||||||
 *                            CCGCTATGCTAGTAGGCAGTCAGCGCCAT
 *                                                   read 2
 * @endcode
 *
 * The resulting alignment of the left/right parts is depicted below. The square brackets indicate clipping positions.
 *
 * @code{.console}
 * reference  AGCATGTTAGATAAGATA    [GCTGTGCTAGTAGGCAGTCAGCGCCAT
 *            ||||||||||||||||||    [ | ||
 *            AGCATGTTAGATAAGATA    [TCCGTCC
 *            read 1
 * reference  AGCATGTTAGATAAGATA]    GTGCTAGTAGGCAGTCAGCGCCAT
 *                              ]     |||||||||||||||||||||||
 *                         CCGCT]    ATGCTAGTAGGCAGTCAGCGCCAT
 *                                                     read 2
 * @endcode
 *
 * In the first case, we want to find the one breakpoint in the reference and the two breakpoints in the reads and the
 * alignment of the left and right well-aligning read parts.  In the second case, we want to find the one breakpoint in
 * the reference and the breakpoint/clipping position in each read.
 *
 * The <tt>splitAlignment()</tt> function takes as the input two alignments.  The sequence in each alignment's first row
 * is the reference and the sequence of the second row is the read.  The sequence has to be the same sequence whereas
 * the reads might differ.  If the reads are the same then this is the same as the first case and if the reads differ
 * then this is the second case.
 *
 * The result is two alignments of the left and right contig path clipped appropriately.  The resulting score is the sum
 * of the scores of both alignments.
 *
 * @section Remarks
 *
 * The DP algorithm is chosen automatically depending on whether the gap open and extension costs are equal.
 *
 * @section Example
 *
 * The following example demonstrates the usage of <tt>splitAlignment</tt> in the first case.  The second case
 * works accordingly.
 *
 * @include demos/dox/align_split/split_alignment.cpp
 *
 * The output is as follows.
 *
 * @include demos/dox/align_split/split_alignment.cpp.stdout
 */

// Variant: unbanded, with Align objects.

template <typename TSequenceL, typename TAlignSpecL, typename TSequenceR, typename TAlignSpecR,
          typename TScoreVal, typename TScoreSpec,
          bool TTop, bool TRight, bool TLeft, bool TBottom, typename TConfigSpec>
int splitAlignment(Align<TSequenceL, TAlignSpecL> & alignL,
                   Align<TSequenceR, TAlignSpecR> & alignR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme,
                   AlignConfig<TTop, TRight, TLeft, TBottom, TConfigSpec> const & config)
{
    SEQAN_ASSERT_EQ_MSG(source(row(alignL, 0)), source(row(alignR, 0)),
                        "Contig must be the same for left and right split alignment.");

    auto tmp = _splitAlignmentImpl(row(alignL, 0), row(alignL, 1), row(alignR, 0), row(alignR, 1), scoringScheme, config);
    return std::get<0>(tmp) + std::get<1>(tmp);
}

template <typename TSequenceL, typename TAlignSpecL, typename TSequenceR, typename TAlignSpecR,
          typename TScoreVal, typename TScoreSpec>
int splitAlignment(Align<TSequenceL, TAlignSpecL> & alignL,
                   Align<TSequenceR, TAlignSpecR> & alignR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    return splitAlignment(alignL, alignR, scoringScheme, AlignConfig<false, false, true, true>());
}

// Variant: unbanded, with Gaps objects.

template <typename TSeqHL, typename TGapSpecHL, typename TSeqVL, typename TGapSpecVL,
          typename TSeqHR, typename TGapSpecHR, typename TSeqVR, typename TGapSpecVR,
          typename TScoreVal, typename TScoreSpec,
          bool TTop, bool TRight, bool TLeft, bool TBottom, typename TConfigSpec>
int splitAlignment(Gaps<TSeqHL, TGapSpecHL> & gapsHL,
                   Gaps<TSeqVL, TGapSpecVL> & gapsVL,
                   Gaps<TSeqHR, TGapSpecHR> & gapsHR,
                   Gaps<TSeqVR, TGapSpecVR> & gapsVR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme,
                   AlignConfig<TTop, TRight, TLeft, TBottom, TConfigSpec> const & config)
{
    SEQAN_ASSERT_EQ_MSG(source(gapsHL), source(gapsHR),
                        "Contig must be the same for left and right split alignment.");

    auto tmp = _splitAlignmentImpl(gapsHL, gapsVL, gapsHR, gapsVR, scoringScheme, config);
    return std::get<0>(tmp) + std::get<1>(tmp);
}

template <typename TSeqHL, typename TGapSpecHL, typename TSeqVL, typename TGapSpecVL,
typename TSeqHR, typename TGapSpecHR, typename TSeqVR, typename TGapSpecVR,
typename TScoreVal, typename TScoreSpec>
int splitAlignment(Gaps<TSeqHL, TGapSpecHL> & gapsHL,
                   Gaps<TSeqVL, TGapSpecVL> & gapsVL,
                   Gaps<TSeqHR, TGapSpecHR> & gapsHR,
                   Gaps<TSeqVR, TGapSpecVR> & gapsVR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    return splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, scoringScheme, AlignConfig<false, false, true, true>());
}

// Variant: banded, with Align objects.
template <typename TSequenceL, typename TAlignSpecL, typename TSequenceR, typename TAlignSpecR,
          typename TScoreVal, typename TScoreSpec,
          bool TTop, bool TRight, bool TLeft, bool TBottom, typename TConfigSpec>
int splitAlignment(Align<TSequenceL, TAlignSpecL> & alignL,
                   Align<TSequenceR, TAlignSpecR> & alignR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme,
                   AlignConfig<TTop, TRight, TLeft, TBottom, TConfigSpec> const & config,
                   int const lowerDiagonal,
                   int const upperDiagonal)
{
    SEQAN_ASSERT_EQ_MSG(source(row(alignL, 0)), source(row(alignR, 0)),
                        "Contig must be the same for left and right split alignment.");

    auto tmp = _splitAlignmentImpl(row(alignL, 0), row(alignL, 1), row(alignR, 0), row(alignR, 1),
                                   scoringScheme, config, lowerDiagonal, upperDiagonal);
    return std::get<0>(tmp) + std::get<1>(tmp);
}

template <typename TSequenceL, typename TAlignSpecL, typename TSequenceR, typename TAlignSpecR,
typename TScoreVal, typename TScoreSpec>
int splitAlignment(Align<TSequenceL, TAlignSpecL> & alignL,
                   Align<TSequenceR, TAlignSpecR> & alignR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme,
                   int const lowerDiagonal,
                   int const upperDiagonal)
{
    return splitAlignment(alignL, alignR, scoringScheme, AlignConfig<false, false, true, true>(),
                          lowerDiagonal, upperDiagonal);
}

// Variant: banded, with Gaps objects.
template <typename TSeqHL, typename TGapSpecHL, typename TSeqVL, typename TGapSpecVL,
          typename TSeqHR, typename TGapSpecHR, typename TSeqVR, typename TGapSpecVR,
          typename TScoreVal, typename TScoreSpec,
          bool TTop, bool TRight, bool TLeft, bool TBottom, typename TConfigSpec>
int splitAlignment(Gaps<TSeqHL, TGapSpecHL> & gapsHL,
                   Gaps<TSeqVL, TGapSpecVL> & gapsVL,
                   Gaps<TSeqHR, TGapSpecHR> & gapsHR,
                   Gaps<TSeqVR, TGapSpecVR> & gapsVR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme,
                   AlignConfig<TTop, TRight, TLeft, TBottom, TConfigSpec> const & config,
                   int const lowerDiagonal,
                   int const upperDiagonal)
{
    SEQAN_ASSERT_EQ_MSG(source(gapsHL), source(gapsHR),
                        "Contig must be the same for left and right split alignment.");

    auto tmp = _splitAlignmentImpl(gapsHL, gapsVL, gapsHR, gapsVR, scoringScheme, config, lowerDiagonal, upperDiagonal);
    return std::get<0>(tmp) + std::get<1>(tmp);
}

template <typename TSeqHL, typename TGapSpecHL, typename TSeqVL, typename TGapSpecVL,
          typename TSeqHR, typename TGapSpecHR, typename TSeqVR, typename TGapSpecVR,
          typename TScoreVal, typename TScoreSpec>
int splitAlignment(Gaps<TSeqHL, TGapSpecHL> & gapsHL,
                   Gaps<TSeqVL, TGapSpecVL> & gapsVL,
                   Gaps<TSeqHR, TGapSpecHR> & gapsHR,
                   Gaps<TSeqVR, TGapSpecVR> & gapsVR,
                   Score<TScoreVal, TScoreSpec> const & scoringScheme,
                   int const lowerDiagonal,
                   int const upperDiagonal)
{
    return splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, scoringScheme, AlignConfig<false, false, true, true>(),
                          lowerDiagonal, upperDiagonal);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_SPLIT_ALIGN_SPLIT_INTERFACE_H_
