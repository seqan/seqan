// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// Approximate String matching via search schemes on a substring index.
// ==========================================================================

#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_H_

namespace seqan {

template <size_t N>
struct Search
{
    std::array<uint8_t, N> pi; // oder of the blocks. permutation of [1..n]
    std::array<uint8_t, N> l; // minimum number of errors at the end of the corresponding block
    std::array<uint8_t, N> u; // maximum number of errors at the end of the corresponding block

    std::array<uint32_t, N> blocklength { }; // cumulated values / prefix sums
    uint32_t startPos = 0;
    // first character of the needle starts with position 1 (not with 0)
    // if initialDirection is true, startPos is one character to the left of the block that is searched first
    // otherwise, startPos is one character right from the block that is searched first

    bool initialDirection = false; // true <=> goToRight
};

template <size_t min, size_t max, typename TVoidType = void>
struct SearchSchemes;
// {
//     static_assert(false, "The error configuration is invalid or there exists no such search scheme.");
// };

template <typename TVoidType>
struct SearchSchemes<0, 0, TVoidType>
{
    static constexpr std::array<Search<1>, 1> VALUE { {{ {{1}}, {{0}}, {{0}} }} };
};

template <typename TVoidType>
constexpr std::array<Search<1>, 1> SearchSchemes<0, 0, TVoidType>::VALUE;

template <typename TVoidType>
struct SearchSchemes<0, 1, TVoidType>
{
    static constexpr std::array<Search<2>, 2> VALUE
    {{
        { {{1, 2}}, {{0, 0}}, {{0, 1}} },
        { {{2, 1}}, {{0, 1}}, {{0, 1}} }
    }};
};

template <typename TVoidType>
constexpr std::array<Search<2>, 2> SearchSchemes<0, 1, TVoidType>::VALUE;

template <typename TVoidType>
struct SearchSchemes<0, 2, TVoidType>
{
    static constexpr std::array<Search<4>, 3> VALUE
    {{
        { {{2, 1, 3, 4}}, {{0, 0, 1, 1}}, {{0, 0, 2, 2}} },
        { {{3, 2, 1, 4}}, {{0, 0, 0, 0}}, {{0, 1, 1, 2}} },
        { {{4, 3, 2, 1}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}} }
    }};
};

template <typename TVoidType>
constexpr std::array<Search<4>, 3> SearchSchemes<0, 2, TVoidType>::VALUE;

template <typename TVoidType>
struct SearchSchemes<0, 3, TVoidType>
{
    static constexpr std::array<Search<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 2, 2}}, {{0, 0, 3, 3, 3}} },
        { {{4, 3, 2, 1, 5}}, {{0, 0, 0, 0, 0}}, {{1, 1, 2, 2, 3}} },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3}}, {{0, 2, 2, 3, 3}} }
    }};
};

template <typename TVoidType>
constexpr std::array<Search<5>, 3> SearchSchemes<0, 3, TVoidType>::VALUE;

// TODO
// NOTE constexpr removed because of blocklength (cannot be const)
// template <typename TVoidType>
// struct SearchSchemes<0, 4>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };

template <typename TVoidType>
struct SearchSchemes<1, 1, TVoidType>
{
    static constexpr std::array<Search<2>, 2> VALUE
    {{
        { {{1, 2}}, {{0, 1}}, {{0, 1}} },
        { {{2, 1}}, {{0, 1}}, {{0, 1}} }
    }};
};

// template <typename TVoidType>
// struct SearchSchemes<1, 2>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<1, 3>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<1, 4>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<2, 2>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<2, 3>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<2, 4>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<3, 3>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<3, 4>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };
//
// template <typename TVoidType>
// struct SearchSchemes<4, 4>
// {
//     static constexpr std::array<Search<X>, Y> VALUE
//     {{
//     }};
// };

// SearchScheme scheme_0_4
// {
//     { {1, 2, 3}, {0, 0, 3}, {1, 4, 5}, {0, 0, 0}, 0, false },
//     { {2, 3, 1}, {0, 0, 0}, {2, 3, 5}, {0, 0, 0}, 0, false },
//     { {3, 2, 1}, {0, 3, 5}, {0, 5, 5}, {0, 0, 0}, 0, false }
// };

// Given the blocklengths (absolute, not cumulative values), assign it to all
// Searches in a SearchScheme. The order of blocklength has to be from left to
// right (regarding blocks)
template <size_t nbrBlocks, size_t N>
inline void _schemeSearchSetBlockLength(std::array<Search<nbrBlocks>, N> & ss, std::vector<uint32_t> const & blocklength)
{
    for (auto & s : ss)
    {
        for (uint8_t i = 0; i < s.blocklength.size(); ++i)
        {
            s.blocklength[i] = blocklength[s.pi[i]-1]
                               + ((i > 0) ? s.blocklength[i-1] : 0);
        }
    }
}

// requires blocklength to be already set!
template <size_t nbrBlocks, size_t N>
inline void _schemeSearchInit(std::array<Search<nbrBlocks>, N> & ss)
{
    // check whether 2nd block is on the left or right and choose initialDirection accordingly
    // (more efficient since we do not have to switch directions and thus have better caching performance)
    // for that we need to slightly modify search()
    for (auto & s : ss)
    {
        s.initialDirection = s.pi[0] < s.pi[1]; // ascending blockorder -> goToRight
        s.startPos = !s.initialDirection; // + 1 if we go left (right border of block), + 0 if we go right (left border of block)
        for (uint8_t i = 0; i < s.pi.size(); ++i)
        {
            if (s.pi[i] < s.pi[0] + !s.initialDirection) // x < y for goRight and x <= y (i.e. x < y + 1) for goLeft
            {
                s.startPos += s.blocklength[i] - ((i > 0) ? s.blocklength[i-1] : 0);
            }
        }
    }
}

template <size_t min, size_t max>
inline void _schemeSearchComputeFixedBlocklength(SearchSchemes<min, max> & ss, uint32_t const needleLength)
{
    uint8_t blocks = ss[0].pi.size();
    uint32_t blocklength = needleLength / blocks;
    uint8_t rest = needleLength - blocks * blocklength;
    std::vector<uint32_t> blocklengths;
    for (uint8_t i = 0; i < blocks; ++i)
    {
        blocklengths.push_back(blocklength + (i < rest));
    }

    _schemeSearchSetBlockLength(ss, blocklengths);
    _schemeSearchInit(ss);
}

template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, size_t nbrBlocks, typename TDir>
inline void _schemeSearchDeletion(TDelegate & delegate,
                                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                  TText2 const & needle,
                                  uint32_t needleLeftIt,
                                  uint32_t needleRightIt,
                                  uint8_t const errors,
                                  Search<nbrBlocks> const & s,
                                  TDir const /**/,
                                  uint8_t const blockIndex)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {
        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
        {
            _schemeSearch(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, Rev(), blockIndex2, true);
        }
        else
        {
            _schemeSearch(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, Fwd(), blockIndex2, true);
        }
    }

    if (maxErrorsLeftInBlock > 0)
    {
        if (goDown(iter, TDir()))
        {
            do
            {
                _schemeSearchDeletion(delegate, iter, needle, needleLeftIt, needleRightIt, errors + 1, s, TDir(), blockIndex);
            } while (goRight(iter, TDir()));
        }
    }
}

template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, size_t nbrBlocks, typename TDir>
inline void _schemeSearchChildren(TDelegate & delegate,
                                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                  TText2 const & needle,
                                  uint32_t const needleLeftIt,
                                  uint32_t const needleRightIt,
                                  uint8_t const errors,
                                  Search<nbrBlocks> const & s,
                                  TDir const /**/,
                                  uint8_t const blockIndex,
                                  bool const indels,
                                  uint8_t minErrorsLeftInBlock)
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightIt - needleLeftIt - 1/*(needleLeftIt != needleRightIt)*/);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()), needle[goToRight ? needleRightIt - 1 : needleLeftIt - 1]);

            // TODO: this is not optimal yet! we have more edges than in the theoretical model,
            // since we go down an edge before we check whether it can even work out!
            if (!indels && minErrorsLeftInBlock > 0 && charsLeft - 1 < minErrorsLeftInBlock - delta)
            {
                continue;
            }

            auto needleLeftIt2 = needleLeftIt - !goToRight;
            auto needleRightIt2 = needleRightIt + goToRight;

            if (needleRightIt - needleLeftIt == s.blocklength[blockIndex])
            {
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                if (indels)
                {
                    _schemeSearchDeletion(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, TDir(), blockIndex);
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                    {
                        _schemeSearch(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, Rev(), blockIndex2, indels);
                    }
                    else
                    {
                        _schemeSearch(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, Fwd(), blockIndex2, indels);
                    }
                }
            }
            else
            {
                _schemeSearch(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, TDir(), blockIndex, indels);
            }

            // Deletion
            if (indels)
            {
                _schemeSearch(delegate, iter, needle, needleLeftIt, needleRightIt, errors + 1, s, TDir(), blockIndex, indels);
            }
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, size_t nbrBlocks, typename TDir>
inline void _schemeSearchExact(TDelegate & delegate,
                               Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                               TText2 const & needle,
                               uint32_t const needleLeftIt,
                               uint32_t const needleRightIt,
                               uint8_t const errors,
                               Search<nbrBlocks> const & s,
                               TDir const /**/,
                               uint8_t const blockIndex,
                               bool const indels)
{
    bool goToRight2 = s.pi[blockIndex + 1] > s.pi[blockIndex]; // TODO: segfault? value doesn't matter for last block, but we should code it better anyway
    if (std::is_same<TDir, Rev>::value)
    {
        uint32_t infixPosLeft = needleRightIt - 1;
        uint32_t infixPosRight = needleLeftIt + s.blocklength[blockIndex] - 1;

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
        {
            return;
        }

        if (goToRight2)
        {
            _schemeSearch(delegate, iter, needle, needleLeftIt, infixPosRight + 2, errors, s, Rev(), std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), indels);
        }
        else
        {
            _schemeSearch(delegate, iter, needle, needleLeftIt, infixPosRight + 2, errors, s, Fwd(), std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), indels);
        }
    }
    else
    {
        // TODO: transform loop s.t. it does not need to use signed integers
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightIt - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftIt - 1;

        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight], TDir()))
            {
                return;
            }
            --infixPosRight;
        }
        if (goToRight2)
        {
            _schemeSearch(delegate, iter, needle, infixPosLeft, needleRightIt, errors, s, Rev(), std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), indels);
        }
        else
        {
            _schemeSearch(delegate, iter, needle, infixPosLeft, needleRightIt, errors, s, Fwd(), std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), indels);
        }
    }
}

template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, size_t nbrBlocks, typename TDir>
inline void _schemeSearch(TDelegate & delegate,
                          Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                          TText2 const & needle,
                          uint32_t const needleLeftIt,
                          uint32_t const needleRightIt,
                          uint8_t const errors,
                          Search<nbrBlocks> const & s,
                          TDir const /**/,
                          uint8_t const blockIndex,
                          bool const indels)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    // Done.
    if (minErrorsLeftInBlock == 0 && needleLeftIt == 0 && needleRightIt == length(needle) + 1) // NOTE: switch to iterator syntax
    {
        delegate(iter);
    }

    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && /*!(blockIndex > 0 && */needleRightIt - needleLeftIt - 1 /*==*/ != s.blocklength[blockIndex]/*)*/) // TODO: why blockIndex > 0? vermutlich wegen needleRightIt == nedleLeftIt
    {
        _schemeSearchExact(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, TDir(), blockIndex, indels);
    }
    // Approximate search in current block.
    else // if (s.blocklength[blockIndex] - (needleRightIt - needleLeftIt - (needleLeftIt != needleRightIt)) >= minErrorsLeftInBlock)
    {
        // Insertion
        if (indels) // && !(needleLeftIt == 0 && needleRightIt == length(needle) + 1) && charsLeft > 0 // no insertions at the end of the block when we didnt increment the blockid (for deletions)
        {
            bool const goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftIt2 = needleLeftIt - !goToRight; // TODO: can we use an unsigned int?
            uint32_t const needleRightIt2 = needleRightIt + goToRight;

            if (needleRightIt - needleLeftIt == s.blocklength[blockIndex])
            {
                // TODO: check!
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                _schemeSearchDeletion(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + 1, s, TDir(), blockIndex);
            }
            else
            {
                _schemeSearch(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + 1, s, TDir(), blockIndex, indels);
            }
        }
        _schemeSearchChildren(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, TDir(), blockIndex, indels, minErrorsLeftInBlock);
    }
}

template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, size_t nbrBlocks>
inline void _schemeSearch(TDelegate & delegate,
                          Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                          TText2 const & needle,
                          Search<nbrBlocks> const & s,
                          bool const indels)
{
    if (s.initialDirection)
    {
        _schemeSearch(delegate, it, needle, s.startPos, s.startPos + 1, 0, s, Rev(), 0, indels);
    }
    else
    {
        _schemeSearch(delegate, it, needle, s.startPos - 1, s.startPos, 0, s, Fwd(), 0, indels);
    }
}

// ----------------------------------------------------------------------------
// Function _find()
// ----------------------------------------------------------------------------

// TODO: think about design
template <size_t minErrors, size_t maxErrors, typename TText, typename TIndexSpec, typename TPattern, typename TDelegate, typename TDistanceTag>
inline void
_find(Index<TText, TIndexSpec> & index,
      TPattern const && pattern,
      TDelegate & delegate,
      Tag<TDistanceTag> /**/)
{
    auto scheme = SearchSchemes<minErrors, maxErrors>::VALUE;
    _schemeSearchComputeFixedBlocklength(scheme, length(pattern));
    Iter<Index<TText, TIndexSpec>, VSTree<TopDown<> > > it(index);
    search(delegate, it, pattern, scheme, IsSameType<Tag<TDistanceTag>, EditDistance>::VALUE, false);
}

template <size_t minErrors, size_t maxErrors, typename TText, typename TIndexSpec, typename TStringSetSpec, typename TPattern, typename TDelegate, typename TDistanceTag, typename TParallelTag>
inline void
_find(Index<TText, TIndexSpec> & index,
      StringSet<TPattern, TStringSetSpec> const && patterns,
      TDelegate & delegate,
      Tag<TDistanceTag> /**/,
      Tag<TParallelTag> /**/)
{
    SEQAN_OMP_PRAGMA(parallel for if (IsSameType<Tag<TParallelTag>, Parallel>::VALUE))
    for (size_t seqNo = 0; seqNo < length(patterns); ++seqNo)
    {
        _find<minErrors, maxErrors>(index, patterns[seqNo], delegate, TDistanceTag());
    }
}

template <size_t minErrors, size_t maxErrors, typename TText, typename TIndexSpec, typename TStringSetSpec, typename TPattern, typename TDelegate, typename TDistanceTag>
inline void
_find(Index<TText, TIndexSpec> & index,
      StringSet<TPattern, TStringSetSpec> const && patterns,
      TDelegate & delegate,
      Tag<TDistanceTag> /**/)
{
    _find<minErrors, maxErrors>(index, patterns, delegate, TDistanceTag(), Serial());
}

}

#endif  // #ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_H_
