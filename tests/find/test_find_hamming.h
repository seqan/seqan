// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Helper functions for testing k-mismatch / hamming distance aproximate
// string search algorithms.
// ==========================================================================

#ifndef TESTS_FIND_TEST_FIND_HAMMING_H_
#define TESTS_FIND_TEST_FIND_HAMMING_H_

using namespace seqan;

// Helper function that iterates over all values of a given length of
// a given string and calls a functor on it.
template <typename TString, typename TFunctor>
void iterateOverStrings(size_t len, TString &string, TFunctor &functor) {
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename Iterator<TString>::Type TIterator;
    resize(string, len, 0);

    bool done = false;  // Flag for breaking loop.
    while (!done) {
        functor(string);
        // Compute next string.
        if (ordValue(string[0]) + 1 == ValueSize<TAlphabet>::VALUE) {
            // Overflow.
            TIterator it = begin(string);
            while (it != end(string)) {
                *it = 0;
                it += 1;
                if (it == end(string)) {
                    done = true;
                } else if (ordValue(*it) + 1 < ValueSize<TAlphabet>::VALUE) {
                    (*it)++;
                    break;
                }
            }
        } else {
            // No Overflow.
            string[0]++;
        }
    }
}

/*
// Called for all needles in ForAllHaystacksFunctor.
template <typename TString> struct ForAllNeedlesFunctor {
    TString &haystack;
    ForAllNeedlesFunctor(TString &haystack_) : haystack(haystack_) {}

    void operator()(TString &needle) {
        // Everything within here is called for all haystacks and all
        // needles.  We now iterate over all error values and try to
        // find the needle in the haystack, allowing the given number
        // of errors.
        std::cout << "haystack = " << haystack << ", needle = " << needle << std::endl;
        for (size_t dist = 0; dist < length(needle); ++dist) {
            std::cout << "dist = " << dist << std::endl;
            // The finder and pattern for the naive algorithm.
            Finder<TString> simpleFinder(haystack);
            Pattern<TString, HammingSimple> simplePattern(needle, -dist);
            // The finder and pattern for the Horspool algorithm.
            Finder<TString> horspoolFinder(haystack);
            Pattern<TString, HammingHorspool> horspoolPattern(needle, -dist);
            int i = 0;
            while (true) {
                bool resSimple = find(simpleFinder, simplePattern);
                bool resHorspool = find(horspoolFinder, horspoolPattern);
                SEQAN_ASSERT_EQ_MSG(resSimple, resHorspool, "i = %d", i);
                if (!resSimple || !resHorspool)
                    break;
                SEQAN_ASSERT_EQ_MSG(position(simpleFinder) + length(needle), endPosition(simpleFinder), "i = %d", i);
                SEQAN_ASSERT_EQ_MSG(position(simpleFinder), position(horspoolFinder), "i = %d", i);
                SEQAN_ASSERT_EQ_MSG(endPosition(simpleFinder), endPosition(horspoolFinder), "i = %d", i);
                SEQAN_ASSERT_EQ_MSG(getScore(simplePattern), getScore(horspoolPattern), "i = %d", i);
                i += 1;
            }
        }
    }
};


// Called for all haystacks in testFindApproximateHamming.
template <typename TString>struct ForAllHaystacksFunctor {
    size_t needle_len_min;
    size_t needle_len_max;

    ForAllHaystacksFunctor(size_t needle_len_min_, size_t needle_len_max_)
            : needle_len_min(needle_len_min_),
              needle_len_max(needle_len_max_) {}

    void operator()(TString &haystack) {
        std::cout << "Haystack == " << haystack << std::endl;
        TString needle;
        ForAllNeedlesFunctor<TString> forAllNeedles(haystack);
        for (size_t m = needle_len_min; m <= needle_len_max; ++m) {
            std::cout << "m = " << m << std::endl;
            iterateOverStrings(m, needle, forAllNeedles);
        }
    }
};


// Runs tests for approximate string search with hamming distance.
//
// The test works as follows: The function will generate all needles
// of the given alphabet of lengths [needle_len_min, needle_len_max]
// and all haystacks of length [haystack_len_min, haystack_len_max].
// It will then build a list of expected match end positions using a
// brute force finder and then build the same list for the finder
// specialization TFinderSpec.  A mismatch in position lists is
// reported as an error.
//
// The finder specialization is selected with TFinderSpec, the
// alphabet with TAlphabet.
template <typename TFinderSpec, typename TAlphabet>
void testFindApproximateHamming(
        size_t needle_len_min = 2, size_t needle_len_max = 2,
        size_t haystack_len_min = 3, size_t haystack_len_max = 4) {
    // Two containers for the expected and found positions;
    typedef typename Position<String<TAlphabet> >::Type TPos;
    String<TPos> expected;
    String<TPos> result;
    // Build all needles and haystacks, run the searches, compare the results.
    typedef String<TAlphabet> TString;
    typedef typename Iterator<TString>::Type TIterator;
    for (size_t n = haystack_len_min; n <= haystack_len_max; ++n) {
        TString haystack;
        std::cout << "n == " << n << std::endl;
        ForAllHaystacksFunctor<TString> forAllHaystacks(
                needle_len_min, needle_len_max);
        iterateOverStrings(n, haystack, forAllHaystacks);
    }
}
*/

#endif  // TESTS_FIND_TEST_FIND_HAMMING_H_

