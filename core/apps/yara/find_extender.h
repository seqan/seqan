// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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

#ifndef APP_YARA_FIND_EXTENDER_H_
#define APP_YARA_FIND_EXTENDER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Extender
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec>
struct Extender
{
    typedef typename Infix<THaystack const>::Type       THaystackInfix;
    typedef ModifiedString<THaystackInfix, ModReverse>  THaystackInfixRev;
    typedef typename Infix<TNeedle const>::Type         TNeedleInfix;
    typedef ModifiedString<TNeedleInfix, ModReverse>    TNeedleInfixRev;

    typedef Finder<THaystackInfix>                      TFinderRight;
    typedef Finder<THaystackInfixRev>                   TFinderLeft;
    typedef PatternState_<TNeedleInfix, TSpec>          TPatternRight;
    typedef PatternState_<TNeedleInfixRev, TSpec>       TPatternLeft;

    THaystack const &   haystack;
//    TFinderRight        finderRight;
//    TFinderLeft         finderLeft;
    TPatternRight       patternRight;
    TPatternLeft        patternLeft;

    Extender(THaystack const & haystack) :
        haystack(haystack)
    {}
};

// ----------------------------------------------------------------------------
// Function _extendLeft()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec,
          typename THaystackInfix, typename TNeedleInfix, typename THaystackPos, typename TErrors, typename TMaxErrors>
inline bool _extendLeft(Extender<THaystack, TNeedle, TSpec> & extender,
                        THaystackInfix & haystackInfix,
                        TNeedleInfix & needleInfix,
                        THaystackPos & matchBegin,
                        TErrors & needleErrors,
                        TMaxErrors maxErrors)
{
    typedef Extender<THaystack, TNeedle, TSpec>         TExtender;
    typedef typename TExtender::THaystackInfixRev       THaystackInfixRev;
    typedef typename TExtender::TNeedleInfixRev         TNeedleInfixRev;
    typedef typename TExtender::TFinderLeft             TFinderLeft;

    typedef typename Value<THaystack>::Type             THaystackString;
    typedef typename Size<THaystackString>::Type        THaystackSize;

    // Lcp trick.
    THaystackSize lcp = 0;
    {  // TODO(holtgrew): Workaround to storing and returning copies in host() for nested infixes/modified strings. This is ugly and should be fixed later.
        TNeedleInfixRev needleInfixRev(needleInfix);
        THaystackInfixRev haystackInfixRev(haystackInfix);
        lcp = lcpLength(haystackInfixRev, needleInfixRev);
    }
    if (lcp == length(needleInfix))
    {
        matchBegin.i2 -= lcp;
        return true;
    }
    setEndPosition(haystackInfix, endPosition(haystackInfix) - lcp);
    setEndPosition(needleInfix, endPosition(needleInfix) - lcp);

    TErrors remainingErrors = maxErrors - needleErrors;
    TErrors minErrors = remainingErrors + 1;
    THaystackSize endPos = 0;

    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Align.
    TNeedleInfixRev needleInfixRev(needleInfix);
    THaystackInfixRev haystackInfixRev(haystackInfix);
    TFinderLeft finder(haystackInfixRev);
    extender.patternLeft.leftClip = remainingErrors;

    while (find(finder, needleInfixRev, extender.patternLeft, -static_cast<int>(remainingErrors)))
    {
        TErrors currentErrors = -getScore(extender.patternLeft);

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = position(finder) + 1;
        }
    }

    needleErrors += minErrors;
    matchBegin.i2 -= endPos + lcp;

    return needleErrors <= maxErrors;
}

// ----------------------------------------------------------------------------
// Function _extendRight()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec,
          typename THaystackInfix, typename TNeedleInfix, typename THaystackPos, typename TErrors, typename TMaxErrors>
inline bool _extendRight(Extender<THaystack, TNeedle, TSpec> & extender,
                         THaystackInfix & haystackInfix,
                         TNeedleInfix & needleInfix,
                         THaystackPos & matchEnd,
                         TErrors & needleErrors,
                         TMaxErrors maxErrors)
{
    typedef Extender<THaystack, TNeedle, TSpec>         TExtender;
    typedef typename TExtender::TFinderRight            TFinderRight;

    typedef typename Value<THaystack>::Type             THaystackString;
    typedef typename Size<THaystackString>::Type        THaystackSize;

    // Lcp trick.
    THaystackSize lcp = lcpLength(haystackInfix, needleInfix);
    if (lcp == length(needleInfix))
    {
        matchEnd.i2 += lcp;
        return true;
    }
    else if (lcp == length(haystackInfix))
    {
        needleErrors += length(needleInfix) - length(haystackInfix);
        matchEnd.i2 += length(needleInfix);
        return needleErrors <= maxErrors;
    }
    setBeginPosition(haystackInfix, beginPosition(haystackInfix) + lcp);
    setBeginPosition(needleInfix, beginPosition(needleInfix) + lcp);

    // NOTE(esiragusa): Uncomment this to disable lcp trick.
//    THaystackPos lcp = 0;

    TErrors remainingErrors = maxErrors - needleErrors;
    TErrors minErrors = remainingErrors + 1;
    THaystackSize endPos = 0;

    // NOTE(esiragusa): Comment this to disable lcp trick.
    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Remove last base.
    THaystackInfix haystackPrefix(haystackInfix);
    TNeedleInfix needlePrefix(needleInfix);
    setEndPosition(haystackPrefix, endPosition(haystackPrefix) - 1);
    setEndPosition(needlePrefix, endPosition(needlePrefix) - 1);

    // Align.
    TFinderRight finder(haystackPrefix);
    extender.patternRight.leftClip = remainingErrors;

    while (find(finder, needlePrefix, extender.patternRight, -static_cast<int>(remainingErrors)))
    {
        THaystackSize currentEnd = position(finder) + 1;
        TErrors currentErrors = -getScore(extender.patternRight);

        // Compare last base.
        if (getValue(haystackInfix, currentEnd) != back(needleInfix))
            if (++currentErrors > remainingErrors)
                continue;

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = currentEnd;
        }
    }

    needleErrors += minErrors;
    matchEnd.i2 += endPos + lcp + 1;

    return needleErrors <= maxErrors;
}

// ----------------------------------------------------------------------------
// Function extend()
// ----------------------------------------------------------------------------
// TODO(esiragusa): _setScoreThreshold(me, maxErrors)
// TODO(esiragusa): extend(me, haystackInfix, needleInfix(needle=host), needleErrors, delegate)

template <typename THaystack, typename TNeedle, typename TSpec,
          typename THaystackPos, typename TNeedlePos, typename TErrors, typename TMaxErrors, typename TDelegate>
inline void
extend(Extender<THaystack, TNeedle, TSpec> & extender,
       TNeedle const & needle,
       THaystackPos haystackBegin,
       THaystackPos haystackEnd,
       TNeedlePos needleBegin,
       TNeedlePos needleEnd,
       TErrors needleErrors,
       TMaxErrors maxErrors,
       TDelegate & delegate)
{
    typedef Extender<THaystack, TNeedle, TSpec>         TExtender;
    typedef typename TExtender::THaystackInfix          THaystackInfix;
    typedef typename TExtender::TNeedleInfix            TNeedleInfix;
    typedef typename Id<THaystack>::Type                THaystackId;

    typedef typename Value<THaystack>::Type             THaystackString;
    typedef typename Size<THaystackString>::Type        THaystackSize;

    THaystackId haystackId = getValueI1(haystackBegin);
    THaystackSize haystackLength = length(extender.haystack[haystackId]);
    TNeedlePos needleLength = length(needle);

    // Extend left.
    THaystackPos matchBegin = haystackBegin;

    if (needleBegin > 0)
    {
        THaystackPos haystackLeftBegin = THaystackPos(haystackId, 0);
        if (haystackBegin.i2 > needleBegin + maxErrors - needleErrors)
            haystackLeftBegin.i2 = haystackBegin.i2 - (needleBegin + maxErrors - needleErrors);

        THaystackInfix haystackLeft = infix(extender.haystack, haystackLeftBegin, haystackBegin);
        TNeedleInfix needleLeft = infix(needle, 0, needleBegin);

        if (!_extendLeft(extender, haystackLeft, needleLeft, matchBegin, needleErrors, maxErrors)) return;
    }

    // Extend right.
    THaystackPos matchEnd = haystackEnd;

    if (needleEnd < needleLength)
    {
        THaystackPos haystackRightEnd = THaystackPos(haystackId, haystackLength);
        if (haystackRightEnd.i2 > haystackBegin.i2 + needleLength - needleBegin + maxErrors - needleErrors)
            haystackRightEnd.i2 = haystackBegin.i2 + needleLength - needleBegin + maxErrors - needleErrors;

        THaystackInfix haystackRight = infix(extender.haystack, haystackEnd, haystackRightEnd);
        TNeedleInfix needleRight = infix(needle, needleEnd, needleLength);

        if (!_extendRight(extender, haystackRight, needleRight, matchEnd, needleErrors, maxErrors)) return;
    }

//    delegate(extender);
    delegate(matchBegin, matchEnd, needleErrors);
}

#endif  // #ifndef APP_YARA_FIND_EXTENDER_H_
