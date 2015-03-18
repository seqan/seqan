// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _findImpl(index, sequence, Backtracking<Exact>)
// ----------------------------------------------------------------------------

template <typename TState, typename TIndex, typename TNeedle,
          typename TThreshold, typename TDelegate, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImpl(TState & /* state */,
          TIndex & index,
          TNeedle const & needle,
          TThreshold /* threshold */,
          TDelegate && delegate,
          Backtracking<Exact, TSpec>)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type TIndexIt;

    TIndexIt indexIt(index);

    if (goDown(indexIt, needle))
    {
        delegate(indexIt, TThreshold());
    }
}

// ----------------------------------------------------------------------------
// Function _findImpl(trieIndex, sequence, Backtracking<Edit/HammingDistance>)
// ----------------------------------------------------------------------------

template <typename TIndexIt, typename TNeedle, typename TNeedleIt,
          typename TThreshold, typename TDelegate, typename TDistance>
inline void
_findBacktracking(TIndexIt indexIt,
                  TNeedle const & needle,
                  TNeedleIt needleIt,
                  TThreshold errors,
                  TThreshold threshold,
                  TDelegate && delegate,
                  TDistance)
{
    // Exact case.
    if (errors == threshold)
    {
        if (goDown(indexIt, suffix(needle, position(needleIt, needle))))
            delegate(indexIt, errors);
    }
    // Approximate case.
    else if (errors < threshold)
    {
        // Base case.
        if (atEnd(needleIt, needle))
        {
            delegate(indexIt, errors);
        }
        // Recursive case.
        else
        {
            // Insertion.
            if (IsSameType<TDistance, EditDistance>::VALUE)
            {
                _findBacktracking(indexIt, needle, needleIt + 1,
                                  static_cast<TThreshold>(errors + 1), threshold, delegate, TDistance());
            }

            if (goDown(indexIt))
            {
                do
                {
                    // Mismatch.
                    TThreshold delta = !ordEqual(parentEdgeLabel(indexIt), value(needleIt));
                    _findBacktracking(indexIt, needle, needleIt + 1,
                                      static_cast<TThreshold>(errors + delta), threshold, delegate, TDistance());

                    // Deletion.
                    if (IsSameType<TDistance, EditDistance>::VALUE)
                    {
                        _findBacktracking(indexIt, needle, needleIt,
                                          static_cast<TThreshold>(errors + 1), threshold, delegate, TDistance());
                    }
                }
                while (goRight(indexIt));
            }
        }
    }
}

template <typename TState, typename TIndex, typename TNeedle,
          typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
SEQAN_FUNC_ENABLE_IF(And<Is<StringTrieConcept<TIndex> >, IsSequence<TNeedle> >, void)
_findImpl(TState & /* indexIt */,
          TIndex & index,
          TNeedle const & needle,
          TThreshold threshold,
          TDelegate && delegate,
          Backtracking<TDistance, TSpec>)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type       TIndexIt;
    typedef typename Iterator<TNeedle const, Standard>::Type  TNeedleIt;

    TIndexIt indexIt(index);
    TNeedleIt needleIt = begin(needle, Standard());
    TThreshold errors = 0;

    _findBacktracking(indexIt, needle, needleIt, errors, threshold, delegate, TDistance());
}

// ----------------------------------------------------------------------------
// Function _findImpl(treeIndex, sequence, Backtracking<HammingDistance>());
// ----------------------------------------------------------------------------

template <typename TState, typename TIndex, typename TNeedle,
          typename TThreshold, typename TDelegate, typename TSpec>
SEQAN_FUNC_ENABLE_IF(And<Is<StringTreeConcept<TIndex> >, IsSequence<TNeedle> >, void)
_findImpl(TState & /* indexIt */,
          TIndex & index,
          TNeedle const & needle,
          TThreshold threshold,
          TDelegate && delegate,
          Backtracking<HammingDistance, TSpec>)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type     TIndexIt;

    TIndexIt indexIt(index);

    // This lambda is stored locally because approximateStringSearch() does not accept rvalue delegates.
    std::function<void(TNeedle const &, TIndexIt const &)> delegator = [&](TNeedle const &, TIndexIt const & it)
    {
        // NOTE(esiragusa): the approximateStringSearch() functor doesn't know the score...
        delegate(it, threshold);
    };

    approximateStringSearch(delegator, needle, indexIt, threshold);
}

// ----------------------------------------------------------------------------
// Function find(index, index, Backtracking<TDistance>());
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
    typedef Index<THaystack, THaystackSpec>                     TTextIndex;
    typedef Index<TNeedle, TNeedleSpec>                         TPatternIndex;
    typedef typename Iterator<TTextIndex, TopDown<> >::Type     TTextIndexIt;
    typedef typename Iterator<TPatternIndex, TopDown<> >::Type  TPatternIndexIt;

    TTextIndexIt textIt(text);
    TPatternIndexIt patternIt(pattern);

    // This lambda is stored locally because approximateTreeSearch() does not accept rvalue delegates.
    std::function<void(TPatternIndexIt const &, TTextIndexIt const &)> delegator =
    [&](TPatternIndexIt const & pIt, TTextIndexIt const & tIt)
    {
        // NOTE(esiragusa): the approximateTreeSearch() functor doesn't know the score...
        delegate(tIt, pIt, threshold);
    };

    approximateTreeSearch(delegator, patternIt, textIt, threshold);
}

// ----------------------------------------------------------------------------
// Function find(index, index, errors, [](...){}, Backtracking<EditDistance>());
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec,
          typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
inline void
find(Index<THaystack, THaystackSpec> & text,
     Index<TNeedle, TNeedleSpec> & pattern,
     TThreshold threshold,
     TDelegate && delegate,
     Backtracking<EditDistance, TSpec>)
{
    typedef Backtracking<EditDistance, TSpec>               TAlgorithm;
    typedef Index<THaystack, THaystackSpec>                 TText;
    typedef Index<TNeedle, TNeedleSpec>                     TPattern;
    typedef Finder_<TText, TPattern, TAlgorithm>            TFinder;

    TFinder finder;

    // This lambda is stored locally because _find() does not accept rvalue delegates.
    std::function<void(TFinder const &)> delegator = [&](TFinder const &)
    {
        delegate(_textIterator(finder), _patternIterator(finder), _getScore(finder));
    };

    _find(finder, text, pattern, threshold, delegator);
}

}

#endif  // #ifndef SEQAN_FIND_INDEX_LAMBDA_H
