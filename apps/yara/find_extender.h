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

#ifndef APP_YARA_FIND_EXTENDER_H_
#define APP_YARA_FIND_EXTENDER_H_

//#define YARA_PRINT_ALIGN

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Extender
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec = void>
struct Extender
{
    typedef typename InfixOnValue<THaystack const>::Type THaystackInfix;
    typedef ModifiedString<THaystackInfix, ModReverse>   THaystackInfixRev;
    typedef typename InfixOnValue<TNeedle const>::Type   TNeedleInfix;
    typedef ModifiedString<TNeedleInfix, ModReverse>     TNeedleInfixRev;

    typedef AlignTextBanded<FindPrefix,
                            NMatchesNone_,
                            NMatchesNone_>               TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                TAlgorithm;

    typedef Finder<THaystackInfix>                       TFinderRight;
    typedef Finder<THaystackInfixRev>                    TFinderLeft;
    typedef PatternState_<TNeedleInfix, TAlgorithm>      TPatternRight;
    typedef PatternState_<TNeedleInfixRev, TAlgorithm>   TPatternLeft;

    THaystack const &   haystack;
//    TFinderRight        finderRight;
//    TFinderLeft         finderLeft;
    TPatternRight       patternRight;
    TPatternLeft        patternLeft;

    Extender(THaystack const & haystack) :
        haystack(haystack)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _checkSeed()
// ----------------------------------------------------------------------------

template <typename THaystackInfix, typename TNeedleInfix, typename TErrors, typename TMaxErrors>
inline bool _checkSeed(THaystackInfix & haystackInfix,
                       TNeedleInfix & needleInfix,
                       TErrors & matchErrors,
                       TMaxErrors maxErrors)
{
    typedef typename Iterator<THaystackInfix, Standard>::Type   THaystackIt;
    typedef typename Iterator<TNeedleInfix, Standard>::Type     TNeedleIt;

    if (length(haystackInfix) != length(needleInfix)) return false;

    THaystackIt hIt = begin(haystackInfix, Standard());
    THaystackIt hEnd = end(haystackInfix, Standard());
    TNeedleIt nIt = begin(needleInfix, Standard());

    for (; hIt != hEnd && matchErrors <= maxErrors; ++hIt, ++nIt)
        matchErrors += !ordEqual(value(hIt), value(nIt));

    return hIt == hEnd && matchErrors <= maxErrors;
}

// ----------------------------------------------------------------------------
// Function _extendLeft()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec,
          typename THaystackInfix, typename TNeedleInfix, typename THaystackPos, typename TErrors, typename TMaxErrors>
inline bool _extendLeft(Extender<THaystack, TNeedle, TSpec> & extender,
                        THaystackInfix & haystackInfix,
                        TNeedleInfix & needleInfix,
                        THaystackPos & matchBegin,
                        TErrors & matchErrors,
                        TMaxErrors maxErrors)
{
    typedef Extender<THaystack, TNeedle, TSpec>                 TExtender;
    typedef typename TExtender::THaystackInfixRev               THaystackInfixRev;
    typedef typename TExtender::TNeedleInfixRev                 TNeedleInfixRev;
    typedef typename TExtender::TFinderLeft                     TFinderLeft;

    typedef typename Value<THaystack>::Type                     THaystackString;
    typedef typename Size<THaystackString>::Type                THaystackSize;

    // Lcp trick.
    THaystackSize lcp = 0;
    {
        // FIXME(holtgrew): Workaround to storing and returning copies in host() for nested infixes/modified strings.
        TNeedleInfixRev needleInfixRev(needleInfix);
        THaystackInfixRev haystackInfixRev(haystackInfix);
        lcp = lcpLength(haystackInfixRev, needleInfixRev);
    }
    if (lcp == length(needleInfix))
    {
		matchBegin = posAdd(matchBegin, -(typename MakeSigned<THaystackSize>::Type)lcp);
        return true;
    }
    setEndPosition(haystackInfix, endPosition(haystackInfix) - lcp);
    setEndPosition(needleInfix, endPosition(needleInfix) - lcp);

    TErrors remainingErrors = maxErrors - matchErrors;
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

    matchErrors += minErrors;
	matchBegin = posAdd(matchBegin, -(typename MakeSigned<THaystackSize>::Type)(endPos + lcp));

    return matchErrors <= maxErrors;
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
                         TErrors & matchErrors,
                         TMaxErrors maxErrors)
{
    typedef Extender<THaystack, TNeedle, TSpec>                 TExtender;
    typedef typename TExtender::TFinderRight                    TFinderRight;

    typedef typename Value<THaystack>::Type                     THaystackString;
    typedef typename Size<THaystackString>::Type                THaystackSize;

    // Lcp trick.
    THaystackSize lcp = lcpLength(haystackInfix, needleInfix);
    if (lcp == length(needleInfix))
    {
        matchEnd = posAdd(matchEnd, lcp);
        return true;
    }
    else if (lcp == length(haystackInfix))
    {
        matchErrors += length(needleInfix) - length(haystackInfix);
        matchEnd = posAdd(matchEnd, lcp);
        return matchErrors <= maxErrors;
    }
    setBeginPosition(haystackInfix, beginPosition(haystackInfix) + lcp);
    setBeginPosition(needleInfix, beginPosition(needleInfix) + lcp);

    // NOTE(esiragusa): Uncomment this to disable lcp trick.
//    THaystackPos lcp = 0;

    TErrors remainingErrors = maxErrors - matchErrors;
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

    matchErrors += minErrors;
    matchEnd = posAdd(matchEnd, endPos + lcp + 1);

    return matchErrors <= maxErrors;
}

// ----------------------------------------------------------------------------
// Function extend()
// ----------------------------------------------------------------------------
// TODO(esiragusa): _setScoreThreshold(me, maxErrors)
// TODO(esiragusa): extend(me, haystackInfix, needleInfix(needle=host), matchErrors, delegate)

template <typename THaystack,typename TNeedle, typename TSpec,
          typename THaystackPos, typename TNeedlePos, typename TDistance, 
          typename TErrors, typename TMaxErrors, typename TDelegate>
inline void
extend(Extender<THaystack, TNeedle, TSpec> & extender,
       TNeedle const & needle,
       THaystackPos haystackBegin,
       THaystackPos haystackEnd,
       TNeedlePos needleBegin,
       TNeedlePos needleEnd,
       TDistance /* tag */,
       TErrors /* needleErrors */,
       TMaxErrors maxErrors,
       TDelegate && delegate)
{
    typedef Extender<THaystack, TNeedle, TSpec>             TExtender;
    typedef typename TExtender::THaystackInfix              THaystackInfix;
    typedef typename TExtender::TNeedleInfix                TNeedleInfix;
    typedef typename Id<THaystack>::Type                    THaystackId;

    typedef typename Value<THaystack>::Type                 THaystackString;
    typedef typename Size<THaystackString>::Type            THaystackSize;

    THaystackId haystackId = getSeqNo(haystackBegin);
    THaystackSize haystackLength = length(extender.haystack[haystackId]);
    TNeedlePos needleLength = length(needle);

    TErrors matchErrors = 0;

#ifdef YARA_PRINT_ALIGN
    std::cerr << std::endl << "======================================================================" << std::endl;
    std::cerr << "SEED: " << Pair<THaystackPos>(haystackBegin, haystackEnd) << " -- " << 
                Pair<TNeedlePos>(needleBegin, needleEnd) << std::endl;

    {
        Align<Dna5String, AnchorGaps<>> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), infix(extender.haystack, haystackBegin, haystackEnd));
        assignSource(row(align, 1), infix(needle, needleBegin, needleEnd));

        globalAlignment(row(align, 0), row(align, 1), Score<short, EditDistance>(),
                        AlignConfig<false, false, false, false>(),
                        -(int)maxErrors, (int)maxErrors);

        std::cerr << row(align, 0) << std::endl;
        std::cerr << row(align, 1) << std::endl;
//        std::cerr << align << std::endl;
    }
#endif // YARA_PRINT_ALIGN

    // Check seed due to Ns randomization in the index.
    THaystackPos haystackLeftEnd = haystackEnd;
    TNeedlePos needleLeftEnd = needleEnd;

    if (IsSameType<TDistance, HammingDistance>::VALUE)
    {
        THaystackInfix haystackSeed = infix(extender.haystack, haystackBegin, haystackEnd);
        TNeedleInfix needleSeed = infix(needle, needleBegin, needleEnd);
        if (!_checkSeed(haystackSeed, needleSeed, matchErrors, maxErrors)) 
            return;

        haystackLeftEnd = haystackBegin;
        needleLeftEnd = needleBegin;
    }

    // Extend left.
    THaystackPos matchBegin = haystackLeftEnd;

    if (needleLeftEnd > 0)
    {
        THaystackPos haystackLeftBegin = haystackBegin;
        THaystackSize haystackLeftOffset = needleBegin + (maxErrors - matchErrors);
        setSeqOffset(haystackLeftBegin, 0);
        if (getSeqOffset(haystackBegin) > haystackLeftOffset)
            setSeqOffset(haystackLeftBegin, getSeqOffset(haystackBegin) - haystackLeftOffset);

        THaystackInfix haystackLeft = infix(extender.haystack, haystackLeftBegin, haystackLeftEnd);
        TNeedleInfix needleLeft = infix(needle, 0, needleLeftEnd);

#ifdef YARA_PRINT_ALIGN
        std::cerr << "----------------------------------------------------------------------" << std::endl;
        std::cerr << "LEFT: " << Pair<THaystackPos>(haystackLeftBegin, haystackEnd) << " -- " 
                  << Pair<TNeedlePos>(0, needleEnd) << std::endl;

        typename TExtender::THaystackInfixRev haystackLeftRev(haystackLeft);
        typename TExtender::TNeedleInfixRev needleLeftRev(needleLeft);
        Align<Dna5String, AnchorGaps<>> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), haystackLeftRev);
        assignSource(row(align, 1), needleLeftRev);
#endif // YARA_PRINT_ALIGN

        if (!_extendLeft(extender, haystackLeft, needleLeft, matchBegin, matchErrors, maxErrors)) 
            return;

#ifdef YARA_PRINT_ALIGN
        globalAlignment(row(align, 0), row(align, 1), Score<short, EditDistance>(),
                        AlignConfig<false, false, false, true>(),
                        -(int)maxErrors, (int)maxErrors);
//        clipSemiGlobal(row(align, 0), row(align, 1));
        std::cerr << row(align, 0) << std::endl;
        std::cerr << row(align, 1) << std::endl;
        std::cerr << "TOTAL ERRORS: " << matchErrors << std::endl;
#endif // YARA_PRINT_ALIGN
    }

    // Extend right.
    THaystackPos matchEnd = haystackEnd;

    if (needleEnd < needleLength)
    {
        THaystackSize haystackRightOffset = needleLength - needleBegin + (maxErrors - matchErrors);
        THaystackPos haystackRightEnd = THaystackPos(haystackId, haystackLength);
        if (getSeqOffset(haystackRightEnd) > getSeqOffset(haystackBegin) + haystackRightOffset)
            setSeqOffset(haystackRightEnd, getSeqOffset(haystackBegin) + haystackRightOffset);

        THaystackInfix haystackRight = infix(extender.haystack, haystackEnd, haystackRightEnd);
        TNeedleInfix needleRight = infix(needle, needleEnd, needleLength);

#ifdef YARA_PRINT_ALIGN
        std::cerr << "----------------------------------------------------------------------" << std::endl;
        std::cerr << "RIGHT: " << Pair<THaystackPos>(haystackEnd, haystackRightEnd) << " -- " 
                  << Pair<TNeedlePos>(needleEnd, needleLength) << std::endl;

        Align<Dna5String, AnchorGaps<>> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), haystackRight);
        assignSource(row(align, 1), needleRight);
#endif // YARA_PRINT_ALIGN

        if (!_extendRight(extender, haystackRight, needleRight, matchEnd, matchErrors, maxErrors)) 
            return;

#ifdef YARA_PRINT_ALIGN
        globalAlignment(row(align, 0), row(align, 1), Score<short, EditDistance>(),
                        AlignConfig<false, false, false, true>());
//        clipSemiGlobal(row(align, 0), row(align, 1));
        std::cerr << row(align, 0) << std::endl;
        std::cerr << row(align, 1) << std::endl;
        std::cerr << "TOTAL ERRORS: " << matchErrors << std::endl;
#endif // YARA_PRINT_ALIGN
    }

#ifdef YARA_PRINT_ALIGN
    std::cerr << "MATCH: " << Pair<THaystackPos>(matchBegin, matchEnd) << " # " << matchErrors << std::endl;
#endif // YARA_PRINT_ALIGN

    delegate(matchBegin, matchEnd, matchErrors);
}

#endif  // #ifndef APP_YARA_FIND_EXTENDER_H_
