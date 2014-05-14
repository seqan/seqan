// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_FIND_INDEX_LAMBDA_H
#define SEQAN_FIND_INDEX_LAMBDA_H

namespace seqan {

// ============================================================================
// Metafunction
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultFind<Index>
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TPattern>
struct DefaultFind<Index<THaystack, THaystackSpec>, TPattern>
{
    typedef Backtracking<Exact> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FindTraits_<Index, Backtracking<Exact> >
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TPattern, typename TSpec>
struct FindTraits_<Index<THaystack, THaystackSpec>, TPattern, Backtracking<Exact, TSpec> >
{
    typedef typename Iterator<Index<THaystack, THaystackSpec>, TopDown<> >::Type    TextIterator;
    typedef typename Iterator<TPattern, Rooted>::Type                               PatternIterator;
    typedef unsigned char                                                           Score;
};

// ----------------------------------------------------------------------------
// Metafunction FindTraits_<Index, Backtracking<HammingDistance> >
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TPattern, typename TSpec>
struct FindTraits_<Index<THaystack, THaystackSpec>, TPattern, Backtracking<HammingDistance, TSpec> >
{
    typedef typename Iterator<Index<THaystack, THaystackSpec>, TopDown<ParentLinks<> > >::Type  TextIterator;
    typedef typename Iterator<TPattern, Rooted>::Type                                           PatternIterator;
    typedef unsigned char                                                                       Score;
};

// ----------------------------------------------------------------------------
// Metafunction FindTraits_<Index, Index, Backtracking<Exact> >
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec, typename TSpec>
struct FindTraits_<Index<THaystack, THaystackSpec>, Index<TNeedle, TNeedleSpec>, Backtracking<Exact, TSpec> >
{
    typedef typename Iterator<Index<THaystack, THaystackSpec>, TopDown<> >::Type    TextIterator;
    typedef typename Iterator<Index<TNeedle, TNeedleSpec>, TopDown<> >::Type        PatternIterator;
    typedef unsigned char                                                           Score;
};

// ----------------------------------------------------------------------------
// Metafunction FindTraits_<Index, Index, Backtracking<HammingDistance> >
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec, typename TSpec>
struct FindTraits_<Index<THaystack, THaystackSpec>, Index<TNeedle, TNeedleSpec>, Backtracking<HammingDistance, TSpec> > :
    FindTraits_<Index<THaystack, THaystackSpec>, Index<TNeedle, TNeedleSpec>, Backtracking<Exact, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction FindTraits_<Index, Index, Backtracking<EditDistance> >
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec, typename TSpec>
struct FindTraits_<Index<THaystack, THaystackSpec>, Index<TNeedle, TNeedleSpec>, Backtracking<EditDistance, TSpec> > :
    FindTraits_<Index<THaystack, THaystackSpec>, Index<TNeedle, TNeedleSpec>, Backtracking<Exact, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction FindState_<Backtracking<HammingDistance> >
// ----------------------------------------------------------------------------

template <typename TIndex, typename TPattern, typename TSpec>
struct FindState_<TIndex, TPattern, Backtracking<HammingDistance, TSpec> >
{
    typedef typename FindTraits_<TIndex, TPattern, Backtracking<HammingDistance, TSpec> >::TextIterator  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _findImpl(..., Backtracking<Exact>)
// ----------------------------------------------------------------------------

template <typename TState, typename THaystack, typename TIndexSpec, typename TNeedle,
          typename TThreshold, typename TDelegate, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImpl(TState & /* state */,
          Index<THaystack, TIndexSpec> & index,
          TNeedle const & needle,
          TThreshold /* threshold */,
          TDelegate && delegate,
          Backtracking<Exact, TSpec>)
{
    typedef Index<THaystack, TIndexSpec>                    TText;
    typedef Backtracking<Exact, TSpec>                      TAlgorithm;
    typedef FindTraits_<TText, TNeedle const, TAlgorithm>   Traits;
    typedef typename Traits::TextIterator                   TTextIt;

    TTextIt textIt(index);

    if (goDown(textIt, needle))
    {
        delegate(textIt, TThreshold());
    }
}

// ----------------------------------------------------------------------------
// Function _findStateInit(..., Backtracking<HammingDistance>)
// ----------------------------------------------------------------------------

template <typename TTextIt, typename TIndex, typename TNeedle, typename TThreshold, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findStateInit(TTextIt & textIt,
               TIndex & index,
               TNeedle const & /* needle */,
               TThreshold /* threshold */,
               Backtracking<HammingDistance, TSpec>)
{
    textIt.index = &index;
    goRoot(textIt);
}

// ----------------------------------------------------------------------------
// Function _findImpl(..., Backtracking<HammingDistance>)
// ----------------------------------------------------------------------------

template <typename TTextIt, typename TIndex, typename TNeedle,
          typename TThreshold, typename TDelegate, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImpl(TTextIt & textIt,
          TIndex & /* index */,
          TNeedle const & needle,
          TThreshold threshold,
          TDelegate && delegate,
          Backtracking<HammingDistance, TSpec>)
{
    typedef Backtracking<HammingDistance, TSpec>            TAlgorithm;
    typedef FindTraits_<TIndex, TNeedle const, TAlgorithm>  Traits;
    typedef typename Traits::PatternIterator                TPatternIt;

    TPatternIt patternIt = begin(needle);
    TThreshold errors = 0;

    do
    {
        // Exact case.
        if (errors == threshold)
        {
            if (goDown(textIt, suffix(needle, position(patternIt))))
            {
                delegate(textIt, errors);
            }

            goUp(textIt);
        }

        // Approximate case.
        else if (errors < threshold)
        {
            // Base case.
            if (atEnd(patternIt))
            {
                delegate(textIt, errors);
            }

            // Recursive case.
            else if (goDown(textIt))
            {
                errors += !ordEqual(parentEdgeLabel(textIt), value(patternIt));
                goNext(patternIt);
                continue;
            }
        }

        // Backtrack.
        do
        {
            // Termination.
            if (isRoot(textIt)) break;

            goPrevious(patternIt);
            errors -= !ordEqual(parentEdgeLabel(textIt), value(patternIt));
        }
        while (!goRight(textIt) && goUp(textIt));

        // Termination.
        if (isRoot(textIt)) break;

        errors += !ordEqual(parentEdgeLabel(textIt), value(patternIt));
        goNext(patternIt);
    }
    while (true);
}

// ----------------------------------------------------------------------------
// Function find(index, index, errors, [](...){}, Backtracking<TDistance>());
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec,
          typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
inline void
find(Index<THaystack, THaystackSpec> & text,
     Index<TNeedle, TNeedleSpec> & pattern,
     TThreshold threshold,
     TDelegate && delegate,
     Backtracking<TDistance, TSpec>)
{
    typedef Index<THaystack, THaystackSpec>         TText;
    typedef Index<TNeedle, TNeedleSpec>             TPattern;
    typedef Backtracking<TDistance, TSpec>          TAlgorithm;
    typedef Finder_<TText, TPattern, TAlgorithm>    TFinder;

    TFinder finder;

    // This lambda is stored locally because _find() does not accept rvalue delegates.
    std::function<void(TFinder const &)> delegator = [&](...)
    {
        delegate(_textIterator(finder), _patternIterator(finder), _getScore(finder));
    };

    _find(finder, text, pattern, threshold, delegator);
}

// ----------------------------------------------------------------------------
// Function find(index, wotdIndex, errors, [](...){}, Backtracking<TDistance>());
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec,
          typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
inline void
find(Index<THaystack, THaystackSpec> & text,
     Index<TNeedle, IndexWotd<TNeedleSpec> > & pattern,
     TThreshold threshold,
     TDelegate && delegate,
     Backtracking<TDistance, TSpec>)
{
    typedef Index<THaystack, THaystackSpec>                 TText;
    typedef Index<TNeedle, IndexWotd<TNeedleSpec> >         TPattern;
    typedef Backtracking<TDistance, TSpec>                  TAlgorithm;

    typedef Finder<TText, TAlgorithm>                       TFinder_;
    typedef Pattern<TPattern, TAlgorithm>                   TPattern_;

    TFinder_ _finder(text);
    TPattern_ _pattern(pattern, length(back(host(pattern))));

    while (_resume(_finder, _pattern, threshold))
    {
        delegate(_finder.index_iterator, _pattern.index_iterator, _pattern.prefix_aligner.errors);
    }
}

}

#endif  // #ifndef SEQAN_FIND_INDEX_LAMBDA_H
