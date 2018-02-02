// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef APP_YARA_FIND_VERIFIER_H_
#define APP_YARA_FIND_VERIFIER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Verifier
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec>
struct Verifier
{
    typedef typename InfixOnValue<THaystack const>::Type                        THaystackInfix;
    typedef String<GapAnchor<int> >                                             TGapAnchors;
    typedef AnchorGaps<TGapAnchors>                                             TAnchorGaps;
    typedef typename Size<THaystackInfix>::Type                                 TSize;
    typedef typename Position<THaystackInfix>::Type                             TPosition;
    typedef TraceSegment_<TPosition, TSize>                                     TTraceSegment;
    typedef String<TTraceSegment>                                               TTrace;
    typedef DPScoutState_<Default>                                              TDPState;
    typedef DPContext<DPCell_<int, AffineGaps>, typename TraceBitMap_<>::Type>  TDPContext;

    // Thread-private data.
    TGapAnchors     contigAnchors;
    TGapAnchors     readAnchors;
    TTrace          traceSegments;
    TDPState        dpScoutState;
    TDPContext      dpContext;

    // Shared-memory read-only data.
    THaystack const &   haystack;

    Verifier(THaystack const & haystack) :
        haystack(haystack)
    {}
};

// ----------------------------------------------------------------------------
// Function verify<AffineGaps>()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle,
          typename THaystackPos, typename TErrors, typename TDelegate>
inline void
verify(Verifier<THaystack, TNeedle, AffineGaps> & me,
       TNeedle const & needle,
       THaystackPos haystackBegin,
       THaystackPos haystackEnd,
       TErrors maxErrors,
       TErrors maxIndels,
       TDelegate & delegate)
{
    typedef Verifier<THaystack, TNeedle, AffineGaps>    TVerifier;
    typedef typename TVerifier::THaystackInfix          THaystackInfix;
    typedef typename TVerifier::TAnchorGaps             TAnchorGaps;
    typedef Gaps<THaystackInfix, TAnchorGaps>           TContigGaps;
    typedef Gaps<TNeedle const, TAnchorGaps>            TReadGaps;
    typedef typename Size<THaystack>::Type              TCount;

    typedef AlignConfig<true, false, false, true>                       TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type         TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps> TAlignConfig2;

    THaystackInfix haystackInfix = infix(me.haystack, haystackBegin, haystackEnd);

    if (empty(haystackInfix)) return;

    clear(me.contigAnchors);
    clear(me.readAnchors);
    TContigGaps contigGaps(haystackInfix, me.contigAnchors);
    TReadGaps readGaps(needle, me.readAnchors);

    clear(me.traceSegments);
    int errors = _setUpAndRunAlignment(me.dpContext,
                                       me.traceSegments,
                                       me.dpScoutState,
                                       source(contigGaps),
                                       source(readGaps),
                                       Score<int>(0, -1000, -999, -1001),
                                       TAlignConfig2()) / -999;
    _adaptTraceSegmentsTo(contigGaps, readGaps, me.traceSegments);

    // PUBLIC INTERFACE
//    int errors = globalAlignment(contigGaps, readGaps,
//                                 Score<int>(0, -1000, -999, -1001),           // Match, mismatch, extend, open.
//                                 AlignConfig<true, false, false, true>(),     // Top, left, right, bottom.
//                                 Gotoh()) / -999;

    clipSemiGlobal(contigGaps, readGaps);

    TCount gapOpens = countGapOpens(contigGaps) + countGapOpens(readGaps);
    TCount gapExtensions = countGapExtensions(contigGaps) + countGapExtensions(readGaps);
    TCount gaps = gapOpens + gapExtensions;
    TCount events = errors + gapOpens - gapExtensions;

    if (events <= maxErrors && gaps <= maxIndels)
    {
        THaystackPos matchBegin = posAdd(haystackBegin, beginPosition(contigGaps));
        THaystackPos matchEnd = posAdd(haystackBegin, endPosition(contigGaps));
        delegate(matchBegin, matchEnd, errors);
    }
}

#endif  // #ifndef APP_YARA_FIND_VERIFIER_H_
