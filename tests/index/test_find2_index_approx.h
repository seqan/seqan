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
// Author: Christopher Pockrandt <github@cpockrandt.de>
// ==========================================================================

#ifndef TESTS_FIND2_INDEX_APPROX_H_
#define TESTS_FIND2_INDEX_APPROX_H_

#include <seqan/index.h>
#include "test_index_helpers.h"

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

template <typename TChar, typename TConfig, typename TRand>
void generateText(String<TChar, TConfig> & text, TRand & rng, unsigned length)
{
    unsigned alphabetSize = ValueSize<TChar>::VALUE;
    unsigned minChar = MinValue<TChar>::VALUE;

    resize(text, length);

    for (unsigned i = 0; i < length; ++i)
    {
        text[i] = rng() % alphabetSize - minChar;
    }
}

template <typename T>
inline void printVector(std::vector<T> const & v, std::ostream & stream)
{
    stream << "[";
    for (unsigned i = 0; i < v.size(); ++i)
    {
        stream << (std::is_same<T, uint8_t>::value ? (unsigned) v[i] : v[i])
               << ((i < v.size() - 1) ? ", " : "");
    }
    stream << "]" << std::endl;
}

inline void printSearch(Search const & s, std::ostream & stream)
{
    stream << "SearchScheme (Pi): ";
    printVector(s.pi, stream);
    stream << "SearchScheme (L): ";
    printVector(s.l, stream);
    stream << "SearchScheme (U): ";
    printVector(s.u, stream);
    stream << "SearchScheme (BL): ";
    printVector(s.blocklength, stream);
}

inline void _getErrorDistributions(std::vector<uint8_t> l,
                                   std::vector<uint8_t> u,
                                   std::vector<std::vector<uint8_t> > & res,
                                   uint8_t const e)
{
    if (l.size() == 0)
    {
        std::vector<uint8_t> emptyVector;
        res.push_back(emptyVector);
        return;
    }

    uint8_t l1 = l[0];
    uint8_t u1 = u[0];
    l.erase(l.begin());
    u.erase(u.begin());

    for (uint8_t i = std::max(e, l1); i <= u1; ++i)
    {
        std::vector<std::vector<uint8_t> > _res;
        _getErrorDistributions(l, u, _res, i);
        for (auto & _resElem : _res)
        {
            _resElem.insert(_resElem.begin(), i - e);
            res.push_back(_resElem);
        }
    }
}

// Compute all possible error distributions given a search.
// The result is ordered as the search (s.pi)
inline void getErrorDistributions(Search const & s,
                                  std::vector<std::vector<uint8_t> > & res)
{
    _getErrorDistributions(s.l, s.u, res, 0u);
}

// Reorder blocks s.t. they are in a sequential order (from left to right)
template <typename T>
inline void orderVector(Search const & s, std::vector<T> & v)
{
    std::vector<T> v_tmp = v;
    for (uint8_t i = 0; i < s.pi.size(); ++i)
    {
        uint8_t index = s.pi[i] - 1;
        v[index] = v_tmp[i];
    }
}

// Reorder blocks s.t. they are in a sequential order (from left to right)
// Blocklength is stored as absolute values instead of cumulative values.
inline void getOrderedSearch(Search const & s, Search & os)
{
    for (uint8_t i = 0; i < s.pi.size(); ++i)
    {
        uint8_t index = s.pi[i] - 1;
        os.pi[index] = s.pi[i];
        os.l[index] = s.l[i];
        os.u[index] = s.u[i];
        os.blocklength[index] = s.blocklength[i] - ((i > 0) ? s.blocklength[i - 1] : 0);
    }
}

// Trivial backtracking that finds *all* matches with given distance.
template <typename TDelegate, typename TText, typename TIndex,
          typename TIndexSpec, typename TText2, typename TNeedleIter>
inline void trivialSearch(TDelegate & delegate,
                          Iter<Index<TText, TIndex>, VSTree<TopDown<TIndexSpec> > > it,
                          TText2 const & needle,
                          TNeedleIter const needleIt,
                          uint8_t const errorsLeft,
                          bool const indels)
{
    if (errorsLeft == 0)
    {
        if (goDown(it, suffix(needle, position(needleIt, needle)), Rev()))
        {
            delegate(it);
        }
    }
    else
    {
        if (atEnd(needleIt, needle))
        {
            delegate(it);
            if (!(indels && errorsLeft > 0))
                return;
        }
        if (indels && !atEnd(needleIt, needle))
        {
            // Insertion
            trivialSearch(delegate, it, needle, needleIt + 1, errorsLeft - 1, indels);
        }
        if (goDown(it, Rev()))
        {
            do
            {
                // Match / Mismatch
                if (!atEnd(needleIt, needle))
                {
                    auto delta = !ordEqual(parentEdgeLabel(it, Rev()), value(needleIt));
                    trivialSearch(delegate, it, needle, needleIt + 1, errorsLeft - delta, indels);
                }

                if (indels)
                {
                    // Deletion
                    trivialSearch(delegate, it, needle, needleIt, errorsLeft - 1, indels);
                }
            } while (goRight(it, Rev()));
        }
    }
}

// Compute random blocklengths (order: left to right)
inline void setRandomBlockLength(SearchScheme & ss, std::mt19937 & rng, uint16_t needleLength)
{
    std::vector<uint8_t> blocklength(ss[0].pi.size());

    // Set minimum length for each block considerung all searches.
    uint8_t maxErrors = 0;
    for (Search const & s : ss)
    {
        maxErrors = std::max(maxErrors, s.u.back());
    }
    // TODO: this enforces to be the needles much longer than necessary!
    for (uint8_t i = 0; i < blocklength.size(); ++i)
    {
        blocklength[i] = maxErrors;
    }
    needleLength -= maxErrors * blocklength.size();

    // Randomly distribute the rest on all blocks
    while (needleLength > 0)
    {
        ++blocklength[rng() % blocklength.size()];
        --needleLength;
    }

    _schemeSearchSetBlockLength(ss, blocklength);
    _schemeSearchInit(ss);
}

template <typename TText, typename TIndex, typename TIndexSpec>
inline void testSearch(std::mt19937 & rng,
                       Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                       Search const & s,
                       Search const & os,
                       unsigned const needleLength,
                       std::vector<uint8_t> const & errorDistribution,
                       bool const indels,
                       time_t const seed)
{
    typedef typename Value<TText>::Type TChar;
    typedef Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > >  TIt;

    TText & text = indexText(container(it.fwdIter));

    unsigned pos = rng() % (length(text) - needleLength + 1);
    TText origNeedle(infix(text, pos, pos + needleLength));

    // Modify needle s.t. it has errors matching errorDistribution.
    TText needle(origNeedle);
    uint8_t cumulativeBlocklength = 0;
    for (uint8_t block = 0; block < s.pi.size(); ++block)
    {
        uint8_t blocklength = os.blocklength[block];
        if (errorDistribution[block] > blocklength)
        {
            std::stringstream stream;
            stream << "Error in block " << (unsigned) block << "(+ 1): "
                   << (unsigned) errorDistribution[block]
                   << " errors cannot fit into a block of length "
                   << (unsigned) blocklength << "." << std::endl;
            stream << "Error Distribution: ";
            printVector(errorDistribution, stream);
            printSearch(s, stream);
            printSearch(os, stream);
            SEQAN_ASSERT_FAIL(toCString(stream.str()));
        }

        // choose random positions in needle that will be a mismatch/indel
        // repeat until all error positions are unique
        std::vector<uint8_t> errorPositions(errorDistribution[block]);
        do
        {
            clear(errorPositions);
            for (uint8_t error = 0; error < errorDistribution[block]; ++error)
            {
                errorPositions.push_back(rng() % blocklength);
            }
            sort(errorPositions.begin(), errorPositions.end());
        } while (adjacent_find(errorPositions.begin(), errorPositions.end()) != errorPositions.end());

        // construct needle with chosen error positions
        for (unsigned error = 0; error < errorPositions.size(); ++error)
        {
            unsigned pos = errorPositions[error] + cumulativeBlocklength;
            TChar newChar;
            do
            {
                newChar = TChar(rng() % ValueSize<TChar>::VALUE);
            } while(needle[pos] == newChar);
            needle[pos] = newChar;
        }
        cumulativeBlocklength += blocklength;
    }

    std::vector<unsigned> hits, expectedHitsSS, expectedHitsTrivial;
    auto delegate = [&hits](TIt const & it)
    {
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
        {
            hits.push_back(getOccurrences(it)[i]);
        }
    };

    // Find all hits using search schemes.
    _schemeSearch(delegate, it, needle, s, indels);
    for (unsigned hit : hits)
    {
        if (infix(text, hit, hit + needleLength) == origNeedle)
        {
            expectedHitsSS.push_back(hit);
        }
    }

    // Find all hits using trivial backtracking.
    clear(hits);
    uint8_t maxErrors = s.u.back();
    trivialSearch(delegate, it, needle, begin(needle, Standard()), maxErrors, indels);
    for (unsigned hit : hits)
    {
        // filter only correct error distributions
        if (origNeedle == infix(text, hit, hit + needleLength))
        {
            bool distributionOkay = true;
            unsigned leftRange = 0, rightRange = 0;
            for (unsigned block = 0; block < s.pi.size(); ++block)
            {
                rightRange += os.blocklength[block];

                uint8_t errors = 0;
                for (unsigned i = leftRange; i < rightRange; ++i)
                {
                    if (hit + i >= length(text))
                    {
                        ++errors;
                    }
                    else
                    {
                        errors += !ordEqual(needle[i], text[hit + i]);
                    }
                }
                if (errors != errorDistribution[block])
                {
                    distributionOkay = false;
                }
                leftRange += os.blocklength[block];
            }
            if (distributionOkay)
            {
                expectedHitsTrivial.push_back(hit);
            }
        }
    }

    // eliminate duplicates
    sort(expectedHitsSS.begin(), expectedHitsSS.end());
    sort(expectedHitsTrivial.begin(), expectedHitsTrivial.end());
    expectedHitsSS.erase(unique(expectedHitsSS.begin(), expectedHitsSS.end()), expectedHitsSS.end());
    expectedHitsTrivial.erase(unique(expectedHitsTrivial.begin(), expectedHitsTrivial.end()), expectedHitsTrivial.end());

    if (expectedHitsSS != expectedHitsTrivial)
    {
        std::stringstream stream;
        stream << "Seed: " << seed << std::endl
               << "Text: " << text << std::endl
               << "ErrorDistribution: ";
        printVector(errorDistribution, stream);
        stream << "Original: " << origNeedle << std::endl
               << "Modified: " << needle << std::endl
               << "ExpectedHitsSS: ";
        printVector(expectedHitsSS, stream);
        stream << "ExpectedHitsTrivial: ";
        printVector(expectedHitsTrivial, stream);
        printSearch(s, stream);
        SEQAN_ASSERT_FAIL(toCString(stream.str()));
    }
}

inline void testSearchScheme(SearchScheme & ss, bool const indels)
{
    typedef DnaString TText;
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
    typedef Index<TText, BidirectionalIndex<FMIndex<void, TMyFastConfig> > > TIndex;
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    time_t seed = ts.tv_nsec;
    std::mt19937 rng(seed);

    TText text;
    SearchScheme os(ss.size());
    std::vector<std::vector<std::vector<uint8_t> > > errorDistributions(ss.size());

    // Calculate all error distributions and sort each of them (from left to right).
    unsigned errors = 0;
    for (uint8_t schemeId = 0; schemeId < ss.size(); ++schemeId)
    {
        os[schemeId] = ss[schemeId]; // TODO: better init!
        getErrorDistributions(ss[schemeId], errorDistributions[schemeId]);
        for (std::vector<uint8_t> & resElem : errorDistributions[schemeId])
        {
            orderVector(ss[schemeId], resElem);
        }
        errors = std::max(errors, (unsigned) ss[schemeId].u.back());
    }

    for (unsigned textLength = 10; textLength < 10000; textLength *= 10)
    {
        generateText(text, rng, textLength);
        TIndex index(text);
        TIter it(index);
        for (unsigned needleLength = std::max(5lu, ss[0].pi.size() * errors); needleLength < std::min(16u, textLength); ++needleLength)
        {
            setRandomBlockLength(ss, rng, needleLength);

            for (uint8_t schemeId = 0; schemeId < ss.size(); ++schemeId)
            {
                getOrderedSearch(ss[schemeId], os[schemeId]);
            }

            // TODO: no parallelization in tests?
            #pragma omp parallel for
            for (uint8_t schemeId = 0; schemeId < ss.size(); ++schemeId)
            {
                for (auto & errorDistribution : errorDistributions[schemeId])
                {
                    testSearch(rng, it, ss[schemeId], os[schemeId], needleLength, errorDistribution, indels, seed);
                }
            }
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

// ----------------------------------------------------------------------------
// Test test_find2_index_approx_hamming
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_find2_index_approx_hamming)
{
    for (unsigned errors = 0; errors < schemes.size(); ++errors)
    {
        unsigned minErrors = 0; // TODO
        // for (unsigned minErrors = 0; minErrors <= errors; ++minErrors)
        // {
            testSearchScheme(schemes[minErrors][errors], true);
        // }
    }
}

// ----------------------------------------------------------------------------
// Test test_find2_index_approx_hamming
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_find2_index_approx_edit)
{
    for (unsigned errors = 0; errors < schemes.size(); ++errors)
    {
        unsigned minErrors = 0; // TODO
        // for (unsigned minErrors = 0; minErrors <= errors; ++minErrors)
        // {
            testSearchScheme(schemes[minErrors][errors], true);
        // }
    }
}

#endif  // TESTS_FIND2_INDEX_APPROX_H_
