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
    typedef typename Infix<THaystack const>::Type   THaystackInfix;
    typedef Finder<THaystackInfix>                  TFinder;
    typedef Pattern<TNeedle const, TSpec>           TPattern;

    THaystack const &   haystack;
    TFinder             finder;
    TPattern            pattern;

    Verifier(THaystack const & haystack) :
        haystack(haystack)
    {
        _patternMatchNOfPattern(pattern, false);
        _patternMatchNOfFinder(pattern, false);
    }
};

// ----------------------------------------------------------------------------
// Class Verifier
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec>
struct Verifier<THaystack, TNeedle, Filter<TSpec> >
{
    typedef typename Infix<THaystack const>::Type           THaystackInfix;
    typedef Finder<THaystackInfix>                          TFinder;
    typedef StringSet<TNeedle, Segment<TNeedle> >           TSeedsSet;
    typedef Pattern<TSeedsSet, TSpec>                       TPattern;

    THaystack const &   haystack;
    TFinder             finder;
    TSeedsSet           seeds;
    TPattern            pattern;

    Verifier(THaystack const & haystack) :
        haystack(haystack)
    {
//        _patternMatchNOfPattern(pattern, false);
//        _patternMatchNOfFinder(pattern, false);
    }
};

// ----------------------------------------------------------------------------
// Function verify()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec,
          typename THaystackPos, typename TErrors, typename TDelegate>
inline void
verify(Verifier<THaystack, TNeedle, TSpec> & verifier,
       TNeedle const & needle,
       THaystackPos haystackBegin,
       THaystackPos haystackEnd,
       TErrors maxErrors,
       TDelegate & delegate)
{
    typedef Verifier<THaystack, TNeedle, TSpec>         TVerifier;
    typedef typename TVerifier::THaystackInfix          THaystackInfix;

    THaystackInfix haystackInfix = infix(verifier.haystack, haystackBegin, haystackEnd);

    clear(verifier.finder);
    setHost(verifier.finder, haystackInfix);
    setHost(verifier.pattern, needle);

    // TODO(esiragusa): Enumerate all minima.
    bool paired = false;
    while (find(verifier.finder, verifier.pattern, -static_cast<int>(maxErrors)))
        paired = true;

//    if (paired) delegate(verifier);
    if (paired) delegate(haystackBegin, haystackEnd, maxErrors);
}

// ----------------------------------------------------------------------------
// Function verify()
// ----------------------------------------------------------------------------

template <typename THaystack, typename TNeedle, typename TSpec,
          typename THaystackPos, typename TErrors, typename TDelegate>
inline void
verify(Verifier<THaystack, TNeedle, Filter<TSpec> > & verifier,
       TNeedle & needle,
       THaystackPos haystackBegin,
       THaystackPos haystackEnd,
       TErrors maxErrors,
       TDelegate & delegate)
{
    typedef Verifier<THaystack, TNeedle, TSpec>         TVerifier;
    typedef typename TVerifier::THaystackInfix          THaystackInfix;
    typedef typename Size<TNeedle>::Type                TNeedleSize;
    typedef uint64_t                                    TWord;

    TNeedleSize seedsCount = maxErrors + 1;
    TNeedleSize seedLength = BitsPerValue<TWord>::VALUE / seedsCount;
//    TNeedleSize needleLength = length(needle);
//    TNeedleSize seedLength = needleLength / seedsCount;

    clear(verifier.seeds);
    reserve(verifier.seeds, seedsCount, Exact());
    setHost(verifier.seeds, needle);

    for (TNeedleSize seedId = 0; seedId < seedsCount; ++seedId)
        appendInfixWithLength(verifier.seeds, seedId * seedLength, seedLength, Exact());

    THaystackInfix haystackInfix = infix(verifier.haystack, haystackBegin, haystackEnd);

    clear(verifier.finder);
    setHost(verifier.finder, haystackInfix);
    setHost(verifier.pattern, verifier.seeds);


    // TODO(esiragusa): Enumerate all minima.
    bool paired = false;
    while (find(verifier.finder, verifier.pattern))
        paired = true;

//    if (paired) delegate(verifier);
    if (paired) delegate(haystackBegin, haystackEnd, maxErrors);
}

#endif  // #ifndef APP_YARA_FIND_VERIFIER_H_
