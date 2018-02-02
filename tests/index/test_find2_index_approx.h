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

#include <mutex>
#include <ctime>

#include <seqan/index.h>

#include "test_index_helpers.h"

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

template <typename TChar, typename TConfig>
void generateText(String<TChar, TConfig> & text, unsigned const length)
{
    unsigned alphabetSize = ValueSize<TChar>::VALUE;
    unsigned minChar = MinValue<TChar>::VALUE;

    resize(text, length);

    for (unsigned i = 0; i < length; ++i)
    {
        text[i] = std::rand() % alphabetSize - minChar;
    }
}

template< typename TContainer>
inline void printVector(const TContainer & c, std::ostream & stream)
{
    using T = std::decay_t<decltype(c[0])>;
    stream << "[";
    for (uint8_t i = 0; i < c.size(); ++i)
    {
        stream << (std::is_same<T, uint8_t>::value ? (unsigned) c[i] : c[i])
               << ((i < c.size() - 1) ? ", " : "");
    }
    stream << "]" << std::endl;
}

template <size_t nbrBlocks>
inline void printSearch(OptimalSearch<nbrBlocks> const & s, std::ostream & stream)
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

inline void _getErrorDistributions(std::vector<std::vector<uint8_t> > & res,
                                   std::vector<uint8_t> l,
                                   std::vector<uint8_t> u,
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
        _getErrorDistributions(_res, l, u, i);
        for (auto & _resElem : _res)
        {
            _resElem.insert(_resElem.begin(), i - e);
            res.push_back(_resElem);
        }
    }
}

// Compute all possible error distributions given a search.
// The result is ordered as the search (s.pi)
template <size_t nbrBlocks>
inline void getErrorDistributions(std::vector<std::vector<uint8_t> > & res,
                                  OptimalSearch<nbrBlocks> const & s)
{
    std::vector<uint8_t> l(s.l.begin(), s.l.end());
    std::vector<uint8_t> u(s.u.begin(), s.u.end());
    _getErrorDistributions(res, l, u, 0u);
}

// Reorder blocks s.t. they are in a sequential order (from left to right)
template <size_t nbrBlocks, typename T>
inline void orderVector(OptimalSearch<nbrBlocks> const & s, std::vector<T> & v)
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
template <size_t nbrBlocks>
inline void getOrderedSearch(OptimalSearch<nbrBlocks> const & s, OptimalSearch<nbrBlocks> & os)
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
inline void trivialSearch(TDelegate && delegate,
                          Iter<Index<TText, TIndex>, VSTree<TopDown<TIndexSpec> > > it,
                          TText2 const & needle,
                          TNeedleIter const needleIt,
                          uint8_t const errorsLeft,
                          bool const indels)
{
    if (errorsLeft == 0)
    {
        if (goDown(it, suffix(needle, position(needleIt, needle)), Rev()))
            delegate(it);
        return;
    }

    if (atEnd(needleIt, needle))
    {
        delegate(it);
        if (!(indels && errorsLeft > 0))
            return;
    }

    // Insertion
    if (indels && !atEnd(needleIt, needle))
        trivialSearch(delegate, it, needle, needleIt + 1, errorsLeft - 1, indels);

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

            // Deletion
            if (indels)
                trivialSearch(delegate, it, needle, needleIt, errorsLeft - 1, indels);

        } while (goRight(it, Rev()));
    }
}

// Compute random blocklengths (order: left to right)
template <size_t nbrBlocks, size_t N>
inline void setRandomBlockLength(std::array<OptimalSearch<nbrBlocks>, N> & ss, uint32_t needleLength)
{
    std::vector<uint32_t> blocklength(ss[0].pi.size());

    // Set minimum length for each block considerung all searches.
    uint8_t maxErrors = 0;
    for (auto const & s : ss)
        maxErrors = std::max(maxErrors, s.u.back());
    // NOTE(cpockrandt): this enforces to be the needles much longer than necessary!
    for (uint8_t i = 0; i < blocklength.size(); ++i)
        blocklength[i] = maxErrors;
    needleLength -= maxErrors * blocklength.size();

    // Randomly distribute the rest on all blocks
    while (needleLength > 0)
    {
        ++blocklength[std::rand() % blocklength.size()];
        --needleLength;
    }

    _optimalSearchSchemeSetBlockLength(ss, blocklength);
    _optimalSearchSchemeInit(ss);
}

template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks, typename TDistanceTag>
inline void testOptimalSearch(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                              OptimalSearch<nbrBlocks> const & s,
                              OptimalSearch<nbrBlocks> const & os,
                              unsigned const needleLength,
                              std::vector<uint8_t> const & errorDistribution,
                              time_t const seed,
                              TDistanceTag const & /**/)
{
    typedef typename Value<TText>::Type TChar;
    typedef Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > >  TIt;

    TText & text = indexText(container(it.fwdIter));

    unsigned pos = std::rand() % (length(text) - needleLength + 1);
    TText origNeedle(infix(text, pos, pos + needleLength));

    // Modify needle s.t. it has errors matching errorDistribution.
    TText needle(origNeedle);
    uint32_t cumulativeBlocklength = 0;
    for (uint8_t block = 0; block < s.pi.size(); ++block)
    {
        uint32_t blocklength = os.blocklength[block];
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
                errorPositions.push_back(std::rand() % blocklength);
            sort(errorPositions.begin(), errorPositions.end());
        } while (adjacent_find(errorPositions.begin(), errorPositions.end()) != errorPositions.end());

        // construct needle with chosen error positions
        for (unsigned error = 0; error < errorPositions.size(); ++error)
        {
            unsigned pos = errorPositions[error] + cumulativeBlocklength;
            TChar newChar;
            do
            {
                newChar = TChar(std::rand() % ValueSize<TChar>::VALUE);
            } while(needle[pos] == newChar);
            needle[pos] = newChar;
        }
        cumulativeBlocklength += blocklength;
    }

    std::vector<unsigned> hits, expectedHitsSS, expectedHitsTrivial;
    auto delegate = [&hits](TIt const & it, TText const & /*needle*/, uint8_t /*errors*/)
    {
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            hits.push_back(getOccurrences(it)[i]);
    };
    auto delegateTrivial = [&hits](TIt const & it)
    {
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            hits.push_back(getOccurrences(it)[i]);
    };

    // Find all hits using search schemes.
    _optimalSearchScheme(delegate, it, needle, s, TDistanceTag());
    for (unsigned hit : hits)
        if (infix(text, hit, hit + needleLength) == origNeedle)
            expectedHitsSS.push_back(hit);

    // Find all hits using trivial backtracking.
    clear(hits);
    uint8_t maxErrors = s.u.back();
    trivialSearch(delegateTrivial, it, needle, begin(needle, Standard()), maxErrors, std::is_same<TDistanceTag, EditDistance>::value);
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
                    if (hit + i >= length(text))
                        ++errors;
                    else
                        errors += !ordEqual(needle[i], text[hit + i]);
                if (errors != errorDistribution[block])
                    distributionOkay = false;
                leftRange += os.blocklength[block];
            }
            if (distributionOkay || std::is_same<TDistanceTag, EditDistance>::value)
                expectedHitsTrivial.push_back(hit);
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

template <size_t nbrBlocks, size_t N, typename TDistanceTag>
inline void testOptimalSearchScheme(std::array<OptimalSearch<nbrBlocks>, N> ss, TDistanceTag const & /**/)
{
    typedef DnaString TText;
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
    typedef Index<TText, BidirectionalIndex<FMIndex<void, TMyFastConfig> > > TIndex;
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    time_t seed = std::time(nullptr);
    std::srand(seed);

    TText text;
    std::array<OptimalSearch<nbrBlocks>, N> os;
    std::vector<std::vector<std::vector<uint8_t> > > errorDistributions(ss.size());

    // Calculate all error distributions and sort each of them (from left to right).
    uint8_t errors = 0;
    for (uint8_t schemeId = 0; schemeId < ss.size(); ++schemeId)
    {
        os[schemeId] = ss[schemeId];
        getErrorDistributions(errorDistributions[schemeId], ss[schemeId]);
        for (std::vector<uint8_t> & resElem : errorDistributions[schemeId])
            orderVector(ss[schemeId], resElem);
        errors = std::max(errors, ss[schemeId].u.back());
    }

    for (uint32_t textLength = 10; textLength < 10000; textLength *= 10)
    {
        generateText(text, textLength);
        TIndex index(text);
        TIter it(index);
        for (uint16_t needleLength = std::max(static_cast<size_t>(5), ss[0].pi.size() * errors); needleLength < std::min(16u, textLength); ++needleLength)
        {
            setRandomBlockLength(ss, needleLength);

            for (uint8_t schemeId = 0; schemeId < ss.size(); ++schemeId)
                getOrderedSearch(ss[schemeId], os[schemeId]);

            for (uint8_t schemeId = 0; schemeId < ss.size(); ++schemeId)
                for (auto & errorDistribution : errorDistributions[schemeId])
                    testOptimalSearch(it, ss[schemeId], os[schemeId], needleLength, errorDistribution, seed, TDistanceTag());
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
    testOptimalSearchScheme(OptimalSearchSchemes<0, 0>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 1>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 2>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 3>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 4>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 1>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 2>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 3>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 4>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<2, 2>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<2, 3>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<2, 4>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<3, 3>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<3, 4>::VALUE, HammingDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<4, 4>::VALUE, HammingDistance());
}

// ----------------------------------------------------------------------------
// Test test_find2_index_approx_edit
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_find2_index_approx_edit)
{
    testOptimalSearchScheme(OptimalSearchSchemes<0, 0>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 1>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 2>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 3>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<0, 4>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 1>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 2>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 3>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<1, 4>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<2, 2>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<2, 3>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<2, 4>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<3, 3>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<3, 4>::VALUE, EditDistance());
    testOptimalSearchScheme(OptimalSearchSchemes<4, 4>::VALUE, EditDistance());
}

SEQAN_DEFINE_TEST(test_find2_index_approx_small_test)
{
    DnaString genome(
        "GAGAGGCCACTCGCAGGATTAAGTCAATAAGTTAATGGCGTCGGCTTCCTGGTATGTAGTACGACGCCCACAGTGACCTCATCGGTGCATTTCCTCATCGTAG"
        "GCGGAACGGTAGACACAAGGCATGATGTCAAATCGCGACTCCAATCCCAAGGTCGCAAGCCTATATAGGAACCCGCTTATGCCCTCTAATCCCGGACAGACCC"
        "CAAATATGGCATAGCTGGTTGGGGGTACCTACTAGGCACAGCCGGAAGCA");
    StringSet<DnaString> needles;
    appendValue(needles, "GGGGTTAT");
    appendValue(needles, "CTAGCTAA");

    std::set<unsigned> hits, expectedHits {186, 226, 227, 234};

    std::mutex mtx;
    auto delegate = [&hits, &mtx](auto & iter, DnaString const & /*pattern*/, uint8_t /*errors*/)
    {
        std::lock_guard<std::mutex> lck(mtx);
        for (auto occ : getOccurrences(iter))
            hits.insert(occ);
    };
    Index<DnaString, BidirectionalIndex<FMIndex<> > > index(genome);

    find<1, 2>(delegate, index, needles, EditDistance(), Serial());
    SEQAN_ASSERT(hits == expectedHits);

    hits.clear();
    find<1, 2>(delegate, index, needles, EditDistance(), Parallel());
    SEQAN_ASSERT(hits == expectedHits);
}

#endif  // TESTS_FIND2_INDEX_APPROX_H_
