// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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

#ifndef TESTS_FIND_BACKTRACKING_EXP_H_
#define TESTS_FIND_BACKTRACKING_EXP_H_

#include <seqan/index.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function testFinder()                                             [FinderTester]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TFinderSpec>
inline void
testFinder(FinderTester<TText, TPattern, TSpec> & tester, Finder_<TText, TPattern, TFinderSpec> const & finder)
{
    typedef typename Fibre<TText, FibreSA>::Type const                          TTextSAFibre;
    typedef typename Infix<TTextSAFibre>::Type                                  TTextOccurrences;
    typedef typename Size<TText>::Type                                          TTextSize;

    TTextOccurrences textOccurrences = getOccurrences(_textIterator(finder));

    TTextSize textOccurrencesCount = length(textOccurrences);

    for (TTextSize i = 0; i < textOccurrencesCount; ++i)
    {
        addResult(tester, textOccurrences[i], 0, _getScore(finder));
#ifdef SEQAN_DEBUG
        std::cout << "text:           " << textOccurrences[i] << std::endl;
        std::cout << "score:          " << static_cast<unsigned>(_getScore(finder)) << std::endl;
        std::cout << std::endl;
#endif
    }
}

template <typename TText, typename TPattern, typename TPatternIndexSpec, typename TSpec, typename TFinderSpec>
inline void
testFinder(FinderTester<TText, Index<TPattern, TPatternIndexSpec>, TSpec> & tester,
           Finder_<TText, Index<TPattern, TPatternIndexSpec>, TFinderSpec> const & finder)
{
    typedef Index<TPattern, TPatternIndexSpec>                                  TPatternIndex;
    typedef typename Fibre<TText, FibreSA>::Type const                          TTextSAFibre;
    typedef typename Fibre<TPatternIndex, FibreSA>::Type const                  TPatternSAFibre;
    typedef typename Infix<TTextSAFibre>::Type                                  TTextOccurrences;
    typedef typename Infix<TPatternSAFibre>::Type                               TPatternOccurrences;
    typedef typename Size<TText>::Type                                          TTextSize;
    typedef typename Size<TPatternIndex>::Type                                  TPatternSize;

    TTextOccurrences textOccurrences = getOccurrences(_textIterator(finder));
    TPatternOccurrences patternOccurrences = getEmptyEdges(_patternIterator(finder));

    TTextSize textOccurrencesCount = length(textOccurrences);
    TPatternSize patternOccurrencesCount = length(patternOccurrences);

    for (TTextSize i = 0; i < textOccurrencesCount; ++i)
        for (TPatternSize j = 0; j < patternOccurrencesCount; ++j)
        {
            addResult(tester, textOccurrences[i], patternOccurrences[j], _getScore(finder));
#ifdef SEQAN_DEBUG
            std::cout << "text:           " << textOccurrences[i] << std::endl;
            std::cout << "pattern:        " << patternOccurrences[j] << std::endl;
            std::cout << "score:          " << static_cast<unsigned>(_getScore(finder)) << std::endl;
            std::cout << std::endl;
#endif
        }
}

// ============================================================================
// Types
// ============================================================================

typedef CharString                              TText;
typedef StringSet<CharString>                   TPattern;
typedef Index<TText, IndexSa<> >                TTextIndex;
typedef Index<TPattern, IndexSa<> >             TPatternIndex;
typedef EditDistance                            TDistance;

// ============================================================================
// Tests
// ============================================================================

// ----------------------------------------------------------------------------
// Test test_find_backtracking_multiple_banana_vs_ada_ana
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_find_backtracking_multiple_edit_banana_vs_ada_ana)
{
    typedef CharString                              TText;
    typedef StringSet<CharString>                   TPattern;
    typedef Index<TText, IndexSa<> >                TTextIndex;
    typedef Index<TPattern, IndexSa<> >             TPatternIndex;
    typedef EditDistance                            TDistance;

    typedef Backtracking<TDistance>                                         TBacktracking;
    typedef FinderTester<TTextIndex, TPatternIndex, TBacktracking>          TTester;
    typedef Finder_<TTextIndex, TPatternIndex, TBacktracking>               TFinder;

    typedef typename Fibre<TTextIndex, FibreSA>::Type                       TTextSAFibre;
    typedef typename Fibre<TPatternIndex, FibreSA>::Type                    TPatternSAFibre;
    typedef typename Value<TTextSAFibre>::Type                              TTextSAPos;
    typedef typename Value<TPatternSAFibre>::Type                           TPatternSAPos;

    TText text = "banana";

    TPattern pattern;
    appendValue(pattern, "ada");
    appendValue(pattern, "ana");

    TTextIndex textIndex(text);

    TPatternIndex patternIndex(pattern);
    indexCreate(patternIndex, FibreSA(), Trie());

    TTester tester;
    TFinder finder;

    addSolution(tester, TTextSAPos(3), TPatternSAPos(0, 0), 1);
    addSolution(tester, TTextSAPos(1), TPatternSAPos(0, 0), 1);
    addSolution(tester, TTextSAPos(3), TPatternSAPos(1, 0), 1);
    addSolution(tester, TTextSAPos(1), TPatternSAPos(1, 0), 1);
    addSolution(tester, TTextSAPos(0), TPatternSAPos(1, 0), 1);
    addSolution(tester, TTextSAPos(4), TPatternSAPos(1, 0), 1);
    addSolution(tester, TTextSAPos(2), TPatternSAPos(1, 0), 1);

    _find(finder, textIndex, patternIndex, 1, tester);
    test(tester);
}

// ----------------------------------------------------------------------------
// Test test_find_backtracking_multiple_hamming_banana_vs_ada_ana
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_find_backtracking_multiple_hamming_banana_vs_ada_ana)
{
    typedef CharString                              TText;
    typedef StringSet<CharString>                   TPattern;
    typedef Index<TText, IndexSa<> >                TTextIndex;
    typedef Index<TPattern, IndexSa<> >             TPatternIndex;
    typedef HammingDistance                         TDistance;

    typedef Backtracking<TDistance>                                     TBacktracking;
    typedef FinderTester<TTextIndex, TPatternIndex, TBacktracking>      TTester;
    typedef Finder_<TTextIndex, TPatternIndex, TBacktracking>           TFinder;

    typedef typename Fibre<TTextIndex, FibreSA>::Type                   TTextSAFibre;
    typedef typename Fibre<TPatternIndex, FibreSA>::Type                TPatternSAFibre;
    typedef typename Value<TTextSAFibre>::Type                          TTextSAPos;
    typedef typename Value<TPatternSAFibre>::Type                       TPatternSAPos;

    TText text = "banana";

    TPattern pattern;
    appendValue(pattern, "ada");
    appendValue(pattern, "ana");

    TTextIndex textIndex(text);

    TPatternIndex patternIndex(pattern);
    indexCreate(patternIndex, FibreSA(), Trie());

    TTester tester;
    TFinder finder;

    addSolution(tester, TTextSAPos(3), TPatternSAPos(0, 0), 1);
    addSolution(tester, TTextSAPos(1), TPatternSAPos(0, 0), 1);
    addSolution(tester, TTextSAPos(3), TPatternSAPos(1, 0), 0);
    addSolution(tester, TTextSAPos(1), TPatternSAPos(1, 0), 0);

    _find(finder, textIndex, patternIndex, 1, tester);
    test(tester);
}

// ----------------------------------------------------------------------------
// Test test_find_backtracking_single_hamming_banana_vs_ada
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_find_backtracking_single_hamming_banana_vs_ada)
{
    typedef CharString                              TText;
    typedef CharString                              TPattern;
    typedef Index<TText, IndexSa<> >                TTextIndex;
    typedef HammingDistance                         TDistance;

    typedef Backtracking<TDistance>                                 TBacktracking;
    typedef FinderTester<TTextIndex, TPattern, TBacktracking>       TTester;
    typedef Finder_<TTextIndex, TPattern, TBacktracking>            TFinder;

    typedef typename Fibre<TTextIndex, FibreSA>::Type               TTextSAFibre;
    typedef typename Value<TTextSAFibre>::Type                      TTextSAPos;

    TText text = "banana";

    TPattern pattern = "ada";

    TTextIndex textIndex(text);

    TTester tester;
    TFinder finder;

    addSolution(tester, TTextSAPos(3), 0, 1);
    addSolution(tester, TTextSAPos(1), 0, 1);

//    _find(finder, textIndex, pattern, 1, tester);
//    test(tester);
}

// ----------------------------------------------------------------------------
// Test test_find_backtracking_single_edit_banana_vs_ada
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_find_backtracking_single_edit_banana_vs_ada)
{
    typedef CharString                              TText;
    typedef CharString                              TPattern;
    typedef Index<TText, IndexSa<> >                TTextIndex;
    typedef EditDistance                            TDistance;

    typedef Backtracking<TDistance>                             TBacktracking;
    typedef FinderTester<TTextIndex, TPattern, TBacktracking>   TTester;
    typedef Finder_<TTextIndex, TPattern, TBacktracking>        TFinder;

    typedef typename Fibre<TTextIndex, FibreSA>::Type           TTextSAFibre;
    typedef typename Value<TTextSAFibre>::Type                  TTextSAPos;

    TText text = "banana";

    TPattern pattern = "ana";

    TTextIndex textIndex(text);

    TTester tester;
    TFinder finder;

    addSolution(tester, TTextSAPos(3), 0, 1);
    addSolution(tester, TTextSAPos(1), 0, 1);
    addSolution(tester, TTextSAPos(0), 0, 1);
    addSolution(tester, TTextSAPos(4), 0, 1);
    addSolution(tester, TTextSAPos(2), 0, 1);


//    _find(finder, textIndex, pattern, 1, tester);
//    test(tester);
}

#endif  // TESTS_FIND_BACKTRACKING_EXP_H_
