// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Approximate string matching via backtracking on VSTrees
// ==========================================================================

#ifndef SANDBOX_ESIRAGUSA_INCLUDE_SEQAN_FIND_BACKTRACKING_DFS_H_
#define SANDBOX_ESIRAGUSA_INCLUDE_SEQAN_FIND_BACKTRACKING_DFS_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TTextIterator, typename TPatternIterator, typename TSize, typename TErrors, typename TDelegate>
inline void _search(TTextIterator textIt,
                    TPatternIterator patternIt,
                    TSize const patternsLength,
                    TErrors const maxErrors,
                    TSize currentDepth,
                    TErrors errors,
                    TDelegate & delegate)
{
    typedef typename EdgeLabel<TTextIterator>::Type         TTextLabel;
    typedef typename EdgeLabel<TPatternIterator>::Type      TPatternLabel;

#ifdef SEQAN_DEBUG
    std::cout << "representative: " << representative(textIt) << std::endl;
    std::cout << "repLength:      " << repLength(textIt) << std::endl;
    std::cout << "current depth:  " << currentDepth << std::endl;
#endif

    // NOTE(esiragusa): Using repLength() is fine only for tries.
    SEQAN_ASSERT_EQ(repLength(textIt), currentDepth);
    SEQAN_ASSERT_EQ(repLength(patternIt), currentDepth);
    
    // An acceptance state was reached.
    if (currentDepth == patternsLength)
    {
        onMatch(delegate, textIt, patternIt, errors);
    }
    else
    {
        // Visit all children of pattern.
        if (currentDepth < patternsLength && goDown(patternIt))
        {
            do
            {
                TTextLabel patternLabel = parentEdgeLabel(patternIt);
                TTextIterator textChildIt = textIt;

                // Search the corresponding children in the text.
                if (goDown(textChildIt, patternLabel))
                {
#ifdef SEQAN_DEBUG
                    std::cout << "pattern:        " << patternLabel << std::endl;
                    std::cout << "errors:         " << static_cast<unsigned>(errors) << std::endl;
#endif
                    _search(textChildIt, patternIt, patternsLength, maxErrors, currentDepth + 1, errors, delegate);

#ifdef SEQAN_DEBUG
                    std::cout << "back to depth:  " << currentDepth << std::endl;
#endif
                }

            } while (goRight(patternIt));
        }
    }
}

template <typename TTextIterator, typename TPatternIterator, typename TSize, typename TErrors, typename TDelegate>
inline void _dfs(TTextIterator textIt,
                 TPatternIterator patternIt,
                 TSize const patternsLength,
                 TErrors const maxErrors,
                 TSize currentDepth,
                 TErrors errors,
                 TDelegate & delegate)
{
    typedef typename EdgeLabel<TTextIterator>::Type         TTextLabel;
    typedef typename EdgeLabel<TPatternIterator>::Type      TPatternLabel;

#ifdef SEQAN_DEBUG
    std::cout << "representative: " << representative(textIt) << std::endl;
    std::cout << "repLength:      " << repLength(textIt) << std::endl;
    std::cout << "current depth:  " << currentDepth << std::endl;
#endif

    // NOTE(esiragusa): Using repLength() is fine only for tries.
    SEQAN_ASSERT_EQ(repLength(textIt), currentDepth);
    SEQAN_ASSERT_EQ(repLength(patternIt), currentDepth);

    // An acceptance state was reached.
    if (currentDepth == patternsLength)
    {
        onMatch(delegate, textIt, patternIt, errors);
    }
    else
    {
        // Visit all children of text and pattern.
        if (goDown(textIt) && (currentDepth < patternsLength && goDown(patternIt)))
        {
            // Visit all children of text.
            do
            {
                TTextLabel textLabel = parentEdgeLabel(textIt);
                TPatternIterator patternChildIt = patternIt;

                // Visit all children of pattern.
                // NOTE(esiragusa): Each children could be visited more than once.
                do
                {
                    TTextLabel patternLabel = parentEdgeLabel(patternChildIt);

                    // Align edge labels.
                    TErrors distance = (textLabel == patternLabel) ? 0 : 1;
                    TErrors newErrors = errors + distance;

#ifdef SEQAN_DEBUG
                    std::cout << "text:           " << textLabel << std::endl;
                    std::cout << "pattern:        " << patternLabel << std::endl;
                    std::cout << "distance:       " << static_cast<unsigned>(distance) << std::endl;
                    std::cout << "errors:         " << static_cast<unsigned>(newErrors) << std::endl;
#endif

                    if (newErrors < maxErrors)
                        _dfs(textIt, patternChildIt, patternsLength, maxErrors, currentDepth + 1, newErrors, delegate);
                    else
                        _search(textIt, patternChildIt, patternsLength, maxErrors, currentDepth + 1, newErrors, delegate);

#ifdef SEQAN_DEBUG
                    std::cout << "back to depth:  " << currentDepth << std::endl;
#endif
                } while (goRight(patternChildIt));

            } while (goRight(textIt));
        }
    }
}

// NOTE(esiragusa): TTextIndex and TPatternIndex must be tries
template <typename TTextIndex, typename TPatternIndex, typename TSize, typename TErrors, typename TDelegate>
void find(TTextIndex & text,
          TPatternIndex & pattern,
          TSize patternsLength,
          TErrors errors,
          TDelegate & delegate,
          DfsPreorder const & /* tag */)
{
    typedef TopDown<>                                               TIteratorSpec;
    typedef typename Iterator<TTextIndex, TIteratorSpec>::Type      TTextIterator;
    typedef typename Iterator<TPatternIndex, TIteratorSpec>::Type   TPatternIterator;
    typedef String<TPatternIterator>                                TPatternIteratorString;

    TTextIterator textIt(text);
    TPatternIterator patternIt(pattern);

    // NOTE(esiragusa): This is necessary since tail recursion was optimised in _dfs().
    if (errors > 0)
        _dfs(textIt, patternIt, patternsLength, errors, static_cast<TSize>(0), static_cast<TErrors>(0), delegate);
    else
        _search(textIt, patternIt, patternsLength, errors, static_cast<TSize>(0), static_cast<TErrors>(0), delegate);
}

}

#endif  // #ifndef SANDBOX_ESIRAGUSA_INCLUDE_SEQAN_FIND_BACKTRACKING_DFS_H_
