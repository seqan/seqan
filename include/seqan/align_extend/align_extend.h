// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// This file contains routines to extend an existing Align object
// ==========================================================================

#ifndef INCLUDE_ALIGN_ALIGN_EXTEND_H
#define INCLUDE_ALIGN_ALIGN_EXTEND_H

namespace seqan2 {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class AliExtContext_
// ----------------------------------------------------------------------------

// Context with memory holding objects for alignment extension
// This can be reused to prevent repeated memory allocations
template <typename TGaps0, typename TGaps1, typename TDPContext>
struct AliExtContext_
{
    typedef typename Size<TGaps0>::Type TSize;
    typedef typename Position<TGaps0>::Type TPosition;

    TGaps0 leftRow0, centerRow0, rightRow0;
    TGaps1 leftRow1, centerRow1, rightRow1;

    TDPContext dpContext;

    String<TraceSegment_<TPosition, TSize> > traceSegment;
};

template <typename TGaps0, typename TGaps1, typename TDPContext>
inline void
clear(AliExtContext_<TGaps0, TGaps1, TDPContext> & prov)
{
    // gaps don't need to be cleared, because they are always
    // re-assigned before use; dpContext, too!

    // trace segment always needs to be cleared
    clear(prov.traceSegment);
    // dpContext doesn't need to be cleared
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _reverseTrace()
// ----------------------------------------------------------------------------

// Reverse a trace string and adapt internal position.
template <typename TPosition, typename TPos2, typename TSize, typename TSpec>
void _reversePartialTrace(String<TraceSegment_<TPosition,
                                               TSize>, TSpec> & trace,
                          TPos2 const lengthH,
                          TPos2 const lengthV)
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
// Function _setUpAndRunAlignImpl()
// ----------------------------------------------------------------------------

template <typename TAliExtContext_, typename TString0, typename TString1, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl(TAliExtContext_ & alignContext,
                      TString0 const & str0,
                      TString1 const & str1,
                      Score<TScoreValue, TScoreSpec> const & scoreScheme,
                      int const /*lowerDiag*/,
                      int const /*upperDiag*/,
                      TScoreValue const /*xDrop*/,
                      TTracebackConfig const & /*gapOrientation*/,
                      False const & /*TBoolBanded*/,
                      False const & /*TBoolXDrop*/)
{
    typedef FreeEndGaps_<False, False, True, True> TFreeEndGaps;
    typedef AlignConfig2<AlignExtend_<>, DPBandConfig<BandOff>, TFreeEndGaps,
                         TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TAlignConfig;

    DPScoutState_<Default> scoutState;
    return _setUpAndRunAlignment(alignContext.dpContext, alignContext.traceSegment, scoutState, str0, str1, scoreScheme,
                                 TAlignConfig());
}

template <typename TAliExtContext_, typename TString0, typename TString1, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl(TAliExtContext_ & alignContext,
                      TString0 const & str0,
                      TString1 const & str1,
                      Score<TScoreValue, TScoreSpec> const & scoreScheme,
                      int const lowerDiag,
                      int const upperDiag,
                      TScoreValue const /*xDrop*/,
                      TTracebackConfig const & /*gapOrientation*/,
                      True const & /*TBoolBanded*/,
                      False const & /*TBoolXDrop*/)
{
    typedef FreeEndGaps_<False, False, True, True> TFreeEndGaps;
    typedef AlignConfig2<AlignExtend_<>, DPBandConfig<BandOn>, TFreeEndGaps,
                         TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TAlignConfig;

    DPScoutState_<Default> scoutState;
    return _setUpAndRunAlignment(alignContext.dpContext, alignContext.traceSegment, scoutState, str0, str1, scoreScheme,
                                 TAlignConfig(lowerDiag, upperDiag));
}

template <typename TAliExtContext_, typename TString0, typename TString1, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl(TAliExtContext_ & alignContext,
                      TString0 const & str0,
                      TString1 const & str1,
                      Score<TScoreValue, TScoreSpec> const & scoreScheme,
                      int const /*lowerDiag*/,
                      int const /*upperDiag*/,
                      TScoreValue const xDrop,
                      TTracebackConfig const & /*gapOrientation*/,
                      False const & /*TBoolBanded*/,
                      True const & /*TBoolXDrop*/)
{
    typedef FreeEndGaps_<False, False, True, True> TFreeEndGaps;
    typedef AlignConfig2<AlignExtend_<XDrop_<TScoreValue> >, DPBandConfig<BandOff>, TFreeEndGaps,
                         TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TAlignConfig;

    DPScoutState_<Terminator_<XDrop_<TScoreValue> > > scoutState(xDrop);
    return _setUpAndRunAlignment(alignContext.dpContext, alignContext.traceSegment, scoutState, str0, str1, scoreScheme,
                                 TAlignConfig());
}

template <typename TAliExtContext_, typename TString0, typename TString1, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl(TAliExtContext_ & alignContext,
                      TString0 const & str0,
                      TString1 const & str1,
                      Score<TScoreValue, TScoreSpec> const & scoreScheme,
                      int const lowerDiag,
                      int const upperDiag,
                      TScoreValue const xDrop,
                      TTracebackConfig const & /*gapOrientation*/,
                      True const & /*TBoolBanded*/,
                      True const & /*TBoolXDrop*/)
{
    typedef FreeEndGaps_<False, False, True, True> TFreeEndGaps;
    typedef AlignConfig2<AlignExtend_<XDrop_<TScoreValue> >, DPBandConfig<BandOn>, TFreeEndGaps,
                         TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TAlignConfig;

    DPScoutState_<Terminator_<XDrop_<TScoreValue> > > scoutState(xDrop);
    return _setUpAndRunAlignment(alignContext.dpContext, alignContext.traceSegment, scoutState, str0, str1, scoreScheme,
                                 TAlignConfig(lowerDiag, upperDiag));
}

// ----------------------------------------------------------------------------
// Function _extendAlignmentImpl()
// ----------------------------------------------------------------------------

template <typename TSource0, typename TSource1, typename TGapsSpec0, typename TGapsSpec1,
          typename TString0, typename TString1, typename TPos, typename TScoreValue, typename TScoreSpec,
          typename TBoolBanded, typename TBoolXDrop, typename TAliExtContext_>
inline TScoreValue
_extendAlignmentImpl(Gaps<TSource0, TGapsSpec0> & row0,
                     Gaps<TSource1, TGapsSpec1> & row1,
                     TScoreValue const & origScore,
                     TString0 const & hSeq,
                     TString1 const & vSeq,
                     Tuple<TPos, 4> const & positions,
                     ExtensionDirection const & direction,
                     int const lowerDiag,
                     int const upperDiag,
                     TScoreValue const & xDrop,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     TBoolBanded const & /**/,
                     TBoolXDrop const & /**/,
                     TAliExtContext_ & alignContext)
{
    typedef typename Infix<TSource0 const>::Type TInf0;
    typedef typename Infix<TSource1 const>::Type TInf1;

    TPos const hBeginPos    = positions[0];
    TPos const vBeginPos    = positions[1];
    TPos const hEndPos      = positions[2];
    TPos const vEndPos      = positions[3];

    SEQAN_ASSERT_EQ(infix(source(row0), beginPosition(row0), endPosition(row0)),
                    infix(hSeq, hBeginPos, hEndPos));
    SEQAN_ASSERT_EQ(infix(source(row1), beginPosition(row1), endPosition(row1)),
                    infix(vSeq, vBeginPos, vEndPos));

    bool extendLeft = ((direction & EXTEND_LEFT) && (hBeginPos > 0u) && (vBeginPos > 0u));
    bool extendRight = ((direction & EXTEND_RIGHT) && (hEndPos < length(hSeq)) && (vEndPos < length(vSeq)));

    clear(alignContext);
    alignContext.centerRow0 = row0;
    alignContext.centerRow1 = row1;

    TScoreValue leftScore   = 0;
    TScoreValue centerScore = origScore;
    TScoreValue rightScore  = 0;

    TPos newAlignLen = length(row0);

    // centerScore was set to "compute yourself" by interface function without score parameter
    if (centerScore == std::numeric_limits<TScoreValue>::min())
    {
        centerScore = 0;

        for (TPos i = 0; i < length(row0); ++i)
        {
            if ( ( isGap(row0, i)) || ( isGap(row1, i)) )
            {
                if (( i==0 ) ||
                    (isGap(row0, i-1) != isGap(row0, i)) ||
                    (isGap(row1, i-1) != isGap(row1, i)) )
                {
                    centerScore += scoreGapOpen(scoreScheme);
                }
                else
                {
                    centerScore += scoreGapExtend(scoreScheme);
                }
            }
            else
            {
                centerScore += score(scoreScheme, row0[i], row1[i]);
            }
        }
    }

    // "reset" original alignment to full length on sequences and no gaps
    assignSource(row0, infix(hSeq, 0, length(hSeq)));
    assignSource(row1, infix(vSeq, 0, length(vSeq)));

    // left
    if (extendLeft)
    {
        TInf0 inf0 = infix(hSeq, 0, hBeginPos);
        TInf1 inf1 = infix(vSeq, 0, vBeginPos);

        // reverse input
        ModifiedString<TInf0, ModReverse> const r_inf0(inf0);
        ModifiedString<TInf1, ModReverse> const r_inf1(inf1);

        leftScore = _setUpAndRunAlignImpl(alignContext, r_inf0, r_inf1, scoreScheme, lowerDiag, upperDiag, xDrop,
                                          TracebackConfig_<CompleteTrace, GapsRight>(), TBoolBanded(), TBoolXDrop());
        // un-reverve
        _reversePartialTrace(alignContext.traceSegment, length(inf0), length(inf1));

        setSource(alignContext.leftRow0, inf0);
        setSource(alignContext.leftRow1, inf1);

        _adaptTraceSegmentsTo(alignContext.leftRow0, alignContext.leftRow1, alignContext.traceSegment);

        if (length(alignContext.leftRow0) > 0)
        {
            integrateGaps(row0, alignContext.leftRow0);
            integrateGaps(row1, alignContext.leftRow1);
            setClippedBeginPosition(row0, clippedBeginPosition(alignContext.leftRow0));
            setClippedBeginPosition(row1, clippedBeginPosition(alignContext.leftRow1));

            newAlignLen += length(alignContext.leftRow0);
        }
        else
        {
            extendLeft = false;
        }
    }

    // center
    if (extendLeft)
    {
        integrateGaps(row0, alignContext.centerRow0, length(alignContext.leftRow0));
        integrateGaps(row1, alignContext.centerRow1, length(alignContext.leftRow1));
    }
    else
    {
        integrateGaps(row0, alignContext.centerRow0, hBeginPos);
        integrateGaps(row1, alignContext.centerRow1, vBeginPos);
        TPos leadGaps0 = countGaps(begin(alignContext.centerRow0));
        TPos leadGaps1 = countGaps(begin(alignContext.centerRow1));

        TPos sourceBeginPos0 = toSourcePosition(alignContext.centerRow0, leadGaps0) + hBeginPos -
                               beginPosition(alignContext.centerRow0);
        TPos sourceBeginPos1 = toSourcePosition(alignContext.centerRow1, leadGaps1) + vBeginPos -
                               beginPosition(alignContext.centerRow1);

        setClippedBeginPosition(row0, toViewPosition(row0, sourceBeginPos0) - leadGaps0);
        setClippedBeginPosition(row1, toViewPosition(row1, sourceBeginPos1) - leadGaps1);
    }

    // right
    if (extendRight)
    {
        TInf0 inf0 = infix(hSeq, hEndPos, length(hSeq));
        TInf1 inf1 = infix(vSeq, vEndPos, length(vSeq));

        clear(alignContext.traceSegment);
        rightScore = _setUpAndRunAlignImpl(alignContext, inf0, inf1, scoreScheme, lowerDiag, upperDiag, xDrop,
                                           TracebackConfig_<CompleteTrace, GapsLeft>(), TBoolBanded(), TBoolXDrop());

        setSource(alignContext.rightRow0, inf0);
        setSource(alignContext.rightRow1, inf1);
        _adaptTraceSegmentsTo(alignContext.rightRow0, alignContext.rightRow1, alignContext.traceSegment);

        if (length(alignContext.rightRow0) > 0)
        {
            integrateGaps(row0, alignContext.rightRow0);
            integrateGaps(row1, alignContext.rightRow1);

            newAlignLen += length(alignContext.rightRow0);
        }
    }

    setClippedEndPosition(row0, clippedBeginPosition(row0) + newAlignLen);
    setClippedEndPosition(row1, clippedBeginPosition(row1) + newAlignLen);

    return leftScore + centerScore + rightScore;
}

// get rows from align object
template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec,
          typename TBoolBanded, typename TBoolXDrop, typename TAliExtContext_>
inline TScoreValue
_extendAlignmentImpl(Align<TStringInfix, TAlignSpec> & align,
                     TScoreValue const & origScore,
                     TString const & hSeq,
                     TString const & vSeq,
                     Tuple<TPos, 4> const & positions,
                     ExtensionDirection const & direction,
                     int const lowerDiag,
                     int const upperDiag,
                     TScoreValue const & xDrop,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     TBoolBanded const & /**/,
                     TBoolXDrop const & /**/,
                     TAliExtContext_ & alignContext)
{
    SEQAN_ASSERT_EQ_MSG(length(rows(align)), 2u, "Only works with pairwise alignments.");
    SEQAN_ASSERT_EQ_MSG(length(row(align, 0)), length(row(align, 1)), "Invalid alignment!");

    return _extendAlignmentImpl(row(align, 0), row(align, 1), origScore, hSeq, vSeq, positions, direction, lowerDiag,
                                upperDiag, xDrop, scoreScheme, TBoolBanded(), TBoolXDrop(), alignContext);
}

// create AlignContext
template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec,
          typename TBoolBanded, typename TBoolXDrop>
inline TScoreValue
_extendAlignmentImpl(Align<TStringInfix, TAlignSpec> & align,
                     TScoreValue const & origScore,
                     TString const & hSeq,
                     TString const & vSeq,
                     Tuple<TPos, 4> const & positions,
                     ExtensionDirection const & direction,
                     int const lowerDiag,
                     int const upperDiag,
                     TScoreValue const & xDrop,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     TBoolBanded const & /**/,
                     TBoolXDrop const & /**/)
{
    if (scoreGapOpen(scoreScheme) == scoreGapExtend(scoreScheme))
    {
        typedef DPContext<DPCell_<TScoreValue, LinearGaps>, typename TraceBitMap_<TScoreValue>::Type> TDPContext;
        typedef AliExtContext_<Gaps<TStringInfix, TAlignSpec>,
                               Gaps<TStringInfix, TAlignSpec>,
                               TDPContext> TAliExtContext_;
        TAliExtContext_ alignContext;
        return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions, direction, lowerDiag, upperDiag, xDrop,
                                    scoreScheme, TBoolBanded(), TBoolXDrop(), alignContext);
    }
    else
    {
        typedef DPContext<DPCell_<TScoreValue, AffineGaps>, typename TraceBitMap_<TScoreValue>::Type> TDPContext;
        typedef AliExtContext_<Gaps<TStringInfix, TAlignSpec>,
                               Gaps<TStringInfix, TAlignSpec>,
                               TDPContext> TAliExtContext_;
        TAliExtContext_ alignContext;
        return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions, direction, lowerDiag, upperDiag, xDrop,
                                    scoreScheme, TBoolBanded(), TBoolXDrop(), alignContext);
    }
}

// ----------------------------------------------------------------------------
// Function extendAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn extendAlignment
 * @headerfile <seqan/align_extend.h>
 * @brief X-Drop extension for alignment objects.
 * @signature TScoreValue extendAlignment(align, [origScore,] hSeq, vSeq, positions, extensionDirection,
 *                                        [lowerDiag, upperDiag,] [xDrop,] scoreScheme);
 *
 * @param[in,out]  align     The @link Align @endlink object to work on.  Must be an alignment over the
 *                           @link InfixSegment infix @endlink of the <i>const</i> type of <tt>hSeq</tt>
 *                           and <tt>vSeq</tt>.  Also see section "Returned Alignment".
 * @param[in]      origScore Original score value of the alignment (optional; computed if not provided).
 * @param[in]      hSeq      Full horizontal sequence.
 * @param[in]      vSeq      Full vertical sequence.
 * @param[in]      positions A @link Tuple @endlink of length 4 with the begin and end position of the
 *                           infixes in align.
 * @param[in]      extensionDirection
 *                           The extension direction (@link ExtensionDirection @endlink).
 * @param[in]      lowerDiag Lower alignment diagonal to use (<tt>int</tt>).
 * @param[in]      upperDiag Upper alignment diagonal to use (<tt>int</tt>).
 * @param[in]      xDrop     The X-drop value to use (integral value). It only limits computation of new
 * columns in the DP-Matrix and has no influence on the diagonals (but can be combined with them).
 * @param[in]      scoringScheme
 *                           The @link Score @endlink to use.
 *
 * @return          TScoreValue
 *                           The score of the new alignment.  <tt>TScoreValue</tt> is the value type of
 *                           <tt>scoringScheme</tt>.
 *
 * @section Returned Alignment
 *
 * The resulting alignment has the infixes extended to the whole underlying sequence.  The alignment
 * is clipped to give the parts of the aligned sequences.
 *
 * @section Example
 *
 * @include demos/dox/align_extend/extend_alignment.cpp
 *
 * The output is as follows:
 *
 * @include demos/dox/align_extend/extend_alignment.cpp.stdout
 *
 * @section Remarks
 *
 * It is necessary to explicitly pass hSeq, vSeq and the positions, because the
 * original hSeq and vSeq (that Align was created on), might have been infixes,
 * (especially if they are members of a ConcatDirect set) in which cases their
 * actual begin and end positions cannot be inferred from the Align object's
 * rows' source().
 */

// NO BAND, NO XDROP
template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, std::numeric_limits<TScoreValue>::min(), hSeq, vSeq, positions, direction, 0, 0, 0, scoreScheme,
                                False(), False());
}

template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TScoreValue const & origScore,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions, direction, 0, 0, 0, scoreScheme, False(),
                                False());
}

// BAND, NO XDROP
template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                int const lowerDiag,
                int const upperDiag,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, std::numeric_limits<TScoreValue>::min(), hSeq, vSeq, positions, direction, lowerDiag, upperDiag,
                                0, scoreScheme, True(), False());
}

template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TScoreValue const & origScore,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                int const lowerDiag,
                int const upperDiag,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions, direction, lowerDiag, upperDiag, 0,
                                scoreScheme, True(), False());
}

// NO BAND, XDROP
template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                TScoreValue const & xDrop,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, std::numeric_limits<TScoreValue>::min(), hSeq, vSeq, positions, direction, 0, 0, xDrop,
                                scoreScheme, False(), True());
}

template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TScoreValue const & origScore,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                TScoreValue const & xDrop,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions, direction, 0, 0, xDrop, scoreScheme, False(),
                                True());
}

// BAND, XDROP
template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                int const lowerDiag,
                int const upperDiag,
                TScoreValue const & xDrop,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, std::numeric_limits<TScoreValue>::min(), hSeq, vSeq, positions, direction, lowerDiag, upperDiag,
                                xDrop, scoreScheme, True(), True());
}

template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TScoreValue const & origScore,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                int const lowerDiag,
                int const upperDiag,
                TScoreValue const & xDrop,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions, direction, lowerDiag, upperDiag, xDrop,
                                scoreScheme, True(), True());
}

template <typename TStringInfix, typename TAlignSpec, typename TString,
          typename TPos, typename TScoreValue, typename TScoreSpec,
          typename TAliExtContext>
inline TScoreValue
extendAlignment(Align<TStringInfix, TAlignSpec> & align,
                TAliExtContext & alignContext,
                TScoreValue const & origScore,
                TString const & hSeq,
                TString const & vSeq,
                Tuple<TPos, 4> const & positions,
                ExtensionDirection const & direction,
                int const lowerDiag,
                int const upperDiag,
                TScoreValue const & xDrop,
                Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions, direction, lowerDiag, upperDiag, xDrop,
                                scoreScheme, True(), True(), alignContext);
}

}

#endif  // INCLUDE_ALIGN_ALIGN_EXTEND_H
