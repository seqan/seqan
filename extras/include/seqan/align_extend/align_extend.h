// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Knut Reinert, FU Berlin
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

#ifndef EXTRAS_INCLUDE_ALIGN_ALIGN_EXTEND_H
#define EXTRAS_INCLUDE_ALIGN_ALIGN_EXTEND_H

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

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

template <typename TTrace, typename TString, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl(TTrace & trace,
                      TString const & str0,
                      TString const & str1,
                      Score<TScoreValue, TScoreSpec> const & scoreScheme,
                      int const /*lowerDiag*/,
                      int const /*upperDiag*/,
                      TScoreValue const /*xDrop*/,
                      TTracebackConfig const & /*gapOrientation*/,
                      False const & /*TBoolBanded*/,
                      False const & /*TBoolXDrop*/)
{
    return _setUpAndRunAlignment(trace, str0, str1, scoreScheme,
                                 AlignExtend_<>(), TTracebackConfig());
}

template <typename TTrace, typename TString, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl( TTrace & trace,
                       TString const & str0,
                       TString const & str1,
                       Score<TScoreValue, TScoreSpec> const & scoreScheme,
                       int const lowerDiag,
                       int const upperDiag,
                       TScoreValue const /*xDrop*/,
                       TTracebackConfig const & /*gapOrientation*/,
                       True const & /*TBoolBanded*/,
                       False const & /*TBoolXDrop*/)
{
    return _setUpAndRunAlignment(trace, str0, str1, scoreScheme,
                                 AlignConfig<false, false, true, true>(),
                                 lowerDiag, upperDiag, AlignExtend_<>(),
                                 TTracebackConfig());
}

template <typename TTrace, typename TString, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl(TTrace & trace,
                      TString const & str0,
                      TString const & str1,
                      Score<TScoreValue, TScoreSpec> const & scoreScheme,
                      int const /*lowerDiag*/,
                      int const /*upperDiag*/,
                      TScoreValue const xDrop,
                      TTracebackConfig const & /*gapOrientation*/,
                      False const & /*TBoolBanded*/,
                      True const & /*TBoolXDrop*/)
{
    DPScoutState_<Terminator_<XDrop_<TScoreValue> > > scoutState(xDrop);
    return _setUpAndRunAlignment(trace, scoutState, str0, str1, scoreScheme,
                                 AlignExtend_<XDrop_<TScoreValue> >(),
                                 TTracebackConfig());
}

template <typename TTrace, typename TString, typename TScoreValue,
          typename TScoreSpec, typename TTracebackConfig>
inline TScoreValue
_setUpAndRunAlignImpl(TTrace & trace,
                      TString const & str0,
                      TString const & str1,
                      Score<TScoreValue, TScoreSpec> const & scoreScheme,
                      int const lowerDiag,
                      int const upperDiag,
                      TScoreValue const xDrop,
                      TTracebackConfig const & /*gapOrientation*/,
                      True const & /*TBoolBanded*/,
                      True const & /*TBoolXDrop*/)
{
    DPScoutState_<Terminator_<XDrop_<TScoreValue> > > scoutState(xDrop);
    return _setUpAndRunAlignment(trace, scoutState, str0, str1, scoreScheme,
                                 AlignConfig<false, false, true, true>(),
                                 lowerDiag, upperDiag,
                                 AlignExtend_<XDrop_<TScoreValue> >(),
                                 TTracebackConfig());
}

// ----------------------------------------------------------------------------
// Function _extendAlignmentImpl()
// ----------------------------------------------------------------------------

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
    typedef typename Infix<TString const>::Type TInf;
    typedef Align<TStringInfix, TAlignSpec> TAlign;

    TPos const hBeginPos    = positions[0];
    TPos const vBeginPos    = positions[1];
    TPos const hEndPos      = positions[2];
    TPos const vEndPos      = positions[3];

    SEQAN_ASSERT_EQ(infix(source(row(align, 0)),
                          beginPosition(row(align, 0)),
                          endPosition(row(align, 0))),
                    infix(hSeq, hBeginPos, hEndPos));
    SEQAN_ASSERT_EQ(infix(source(row(align, 1)),
                          beginPosition(row(align, 1)),
                          endPosition(row(align, 1))),
                    infix(vSeq, vBeginPos, vEndPos));

    bool extendLeft   = ((direction & EXTEND_LEFT) && (hBeginPos > 0u)
                         && (vBeginPos > 0u));
    bool extendRight  = ((direction & EXTEND_RIGHT) && (hEndPos < length(hSeq))
                         && (vEndPos < length(vSeq)));

    TAlign leftAlign;
    resize(rows(leftAlign), 2);
    TAlign centerAlign(align);
    // copy original alignment
    TAlign rightAlign;
    resize(rows(rightAlign), 2);

    TScoreValue leftScore   = 0;
    TScoreValue centerScore = origScore;
    TScoreValue rightScore  = 0;

    if (centerScore == minValue<TScoreValue>())
    {
        centerScore = 0;
        typename Row<TAlign>::Type const & row0 = row(align, 0);
        typename Row<TAlign>::Type const & row1 = row(align, 1);

        for (TPos i = 0; i < length(row(align, 0)); ++i)
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
                centerScore += score(scoreScheme, value(row0, i),
                                     value(row1, i));
            }
        }
    }

    // "reset" original alignment to full length on sequences and no gaps
    assignSource(row(align, 0), infix(hSeq, 0, length(hSeq)));
    assignSource(row(align, 1), infix(vSeq, 0, length(vSeq)));

    // left
    if (extendLeft)
    {
        TInf inf0 = infix(hSeq, 0, hBeginPos);
        TInf inf1 = infix(vSeq, 0, vBeginPos);

        // reverse input
        ModifiedString<TInf, ModReverse> const r_inf0(inf0);
        ModifiedString<TInf, ModReverse> const r_inf1(inf1);

        typedef TraceSegment_<TPos, TPos> TTraceSegment;
        String<TTraceSegment> traceRev;

        leftScore =
            _setUpAndRunAlignImpl(traceRev, r_inf0, r_inf1, scoreScheme,
                                  lowerDiag, upperDiag, xDrop,
                                  TracebackConfig_<CompleteTrace, GapsRight>(),
                                  TBoolBanded(), TBoolXDrop());
        // un-reverve
        _reversePartialTrace(traceRev, length(inf0), length(inf1));

        assignSource(row(leftAlign, 0), inf0);
        assignSource(row(leftAlign, 1), inf1);

        _adaptTraceSegmentsTo(row(leftAlign, 0), row(leftAlign, 1), traceRev);

        if (length(row(leftAlign, 0)) > 0)
        {
            integrateAlign(align, leftAlign);
            setClippedBeginPosition(row(align, 0),
                                    clippedBeginPosition(row(leftAlign, 0)));
            setClippedBeginPosition(row(align, 1),
                                    clippedBeginPosition(row(leftAlign, 1)));
        } else
        {
            extendLeft = false;
        }
    }

    // center
    integrateAlign(align, centerAlign);
    if (!extendLeft)
    {
        // These copmutations look obscure but are due to the obscure
        // nature of the Gaps module...
        setClippedBeginPosition(row(align, 0), hBeginPos
                                - beginPosition(row(centerAlign, 0))
                                + clippedBeginPosition(row(centerAlign, 0)));
        setClippedBeginPosition(row(align, 1), vBeginPos
                                - beginPosition(row(centerAlign, 1))
                                + clippedBeginPosition(row(centerAlign, 1)));
    }

    // right
    if (extendRight)
    {

        TInf inf0 = infix(hSeq, hEndPos, length(hSeq));
        TInf inf1 = infix(vSeq, vEndPos, length(vSeq));

        typedef TraceSegment_<TPos, TPos> TTraceSegment;
        String<TTraceSegment> trace;

        rightScore =
            _setUpAndRunAlignImpl(trace, inf0, inf1, scoreScheme, lowerDiag,
                                  upperDiag, xDrop,
                                  TracebackConfig_<CompleteTrace, GapsLeft>(),
                                  TBoolBanded(), TBoolXDrop());

        assignSource(row(rightAlign, 0), inf0);
        assignSource(row(rightAlign, 1), inf1);
        _adaptTraceSegmentsTo(row(rightAlign, 0), row(rightAlign, 1), trace);

        if (length(row(rightAlign, 0)) > 0)
            integrateAlign(align, rightAlign);
    }

    setClippedEndPosition(row(align, 0), clippedBeginPosition(row(align, 0))
                          + length(row(  leftAlign, 0))
                          + length(row(centerAlign, 0))
                          + length(row( rightAlign, 0)));
    setClippedEndPosition(row(align, 1), clippedBeginPosition(row(align, 1))
                          + length(row(  leftAlign, 1))
                          + length(row(centerAlign, 1))
                          + length(row( rightAlign, 1)));

    return leftScore + centerScore + rightScore;
}

// ----------------------------------------------------------------------------
// Function extendAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn extendAlignment
 * @brief extends an Align-Object
 * @signature TScoreValue extendAlignment(align, [origScore, ] hSeq, vSeq, positions, extensionDirection, [lowerDiag, upperDiag, ] [xDrop, ] scoreScheme);
 *
 * @param[in, out]   align               "Seed-Alignment", Align object on
 * infixes of hSeq and vSeq; will be returned as extended and clipped alignment
 * on full sequences
 * @param[in]       origScore           score of the seed alignment (computed if not given)
 * @param[in]       hSeq                full horizontal sequence
 * @param[in]       vSeq                full vertical sequence
 * @param[in]       positions           the begin and end positions of infixes
 * in alignment (hBegin, vBegin, hEnd, vEnd)
 * @param[in]       extensionDirection  ExtensionDirection
 * @param[in]       lowerDiag           lower diagonal for banded extension [int]
 * (optional)
 * @param[in]       upperDiag           upper diagonal for banded extension [int]
 * (optional)
 * @param[in]       xDrop               abort extension, if all values in the
 * DP-column fall xDrop under the maximum seen (optional)
 * @param[in]       scoreScheme     SeqAn Score object [Score]
 *
 * @return          TScoreValue         the score of the new alignment
 * @headerfile seqan/align_extend.h
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
    return _extendAlignmentImpl(align, minValue<TScoreValue>(), hSeq, vSeq,
                                positions, direction, 0, 0, 0, scoreScheme,
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
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions,
                                direction, 0, 0, 0, scoreScheme, False(),
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
    return _extendAlignmentImpl(align, minValue<TScoreValue>(), hSeq, vSeq,
                                positions, direction, lowerDiag, upperDiag, 0,
                                scoreScheme, True(), False());
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
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions,
                                direction, lowerDiag, upperDiag, 0, scoreScheme,
                                True(), False());
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
    return _extendAlignmentImpl(align, minValue<TScoreValue>(), hSeq, vSeq,
                                positions, direction, 0, 0, xDrop, scoreScheme,
                                False(), True());
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
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions,
                                direction, 0, 0, xDrop, scoreScheme, False(),
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
    return _extendAlignmentImpl(align, minValue<TScoreValue>(), hSeq, vSeq,
                                positions, direction, lowerDiag, upperDiag,
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
    return _extendAlignmentImpl(align, origScore, hSeq, vSeq, positions,
                                direction, lowerDiag, upperDiag, xDrop,
                                scoreScheme, True(), True());
}

}

#endif
