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
// Approximate String matching via search schemes on a substring index.
// ==========================================================================

#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_H_

namespace seqan {

template <size_t N>
struct OptimalSearch
{
    std::array<uint8_t, N> pi; // order of the blocks. permutation of [1..n]
    std::array<uint8_t, N> l; // minimum number of errors at the end of the corresponding block
    std::array<uint8_t, N> u; // maximum number of errors at the end of the corresponding block

    std::array<uint32_t, N> blocklength; // cumulated length of the blocks in Search Scheme order
    //NOTE (svnbgnk) added additional information about search schemes depending on the read length
    //These values are not set to Zero during the creation of Optimal Search Schemes
    std::array<uint32_t, N> chronBL;  //cumulated length of block from left
    uint32_t startPos;
};

template <size_t min, size_t max, typename TVoidType = void>
struct OptimalSearchSchemes;

template <typename TVoidType>
struct OptimalSearchSchemes<0, 0, TVoidType>
{
    static constexpr std::array<OptimalSearch<1>, 1> VALUE { {{ {{1}}, {{0}}, {{0}}, {{0}}, 0 }} };
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<1>, 1> OptimalSearchSchemes<0, 0, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<0, 1, TVoidType>
{
    static constexpr std::array<OptimalSearch<2>, 2> VALUE
    {{
        { {{1, 2}}, {{0, 0}}, {{0, 1}}, {{0, 0}}, 0 },
        { {{2, 1}}, {{0, 1}}, {{0, 1}}, {{0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<2>, 2> OptimalSearchSchemes<0, 1, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<0, 2, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 3> VALUE
    {{
        { {{2, 1, 3, 4}}, {{0, 0, 1, 1}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{3, 2, 1, 4}}, {{0, 0, 0, 0}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 3> OptimalSearchSchemes<0, 2, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<0, 3, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 2, 2}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1, 5}}, {{0, 0, 0, 0, 0}}, {{1, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3}}, {{0, 2, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<0, 3, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<0, 4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 4}}, {{0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 0}}, {{2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 3, 3}}, {{0, 0, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<0, 4, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<1, 1, TVoidType>
{
    static constexpr std::array<OptimalSearch<2>, 2> VALUE
    {{
        { {{1, 2}}, {{0, 1}}, {{0, 1}}, {{0, 0}}, 0 },
        { {{2, 1}}, {{0, 1}}, {{0, 1}}, {{0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<2>, 2> OptimalSearchSchemes<1, 1, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<1, 2, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 3> VALUE
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 1}}, {{0, 0, 0, 1}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1}}, {{0, 0, 1, 1}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 3> OptimalSearchSchemes<1, 2, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<1, 3, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 2, 2}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1, 5}}, {{0, 0, 0, 0, 1}}, {{1, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3}}, {{0, 2, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<1, 3, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<1, 4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 4}}, {{0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 1}}, {{2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 3, 3}}, {{0, 0, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<1, 4, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<2, 2, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 3> VALUE
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 1}}, {{0, 0, 1, 2}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1}}, {{0, 0, 0, 2}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 3> OptimalSearchSchemes<2, 2, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<2, 3, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 2, 2}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1, 5}}, {{0, 0, 0, 0, 2}}, {{1, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3}}, {{0, 2, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<2, 3, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<2, 4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 4}}, {{0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 2}}, {{2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 3, 3}}, {{0, 0, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<2, 4, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<3, 3, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 2, 3}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1, 5}}, {{0, 0, 0, 0, 3}}, {{1, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3}}, {{0, 2, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<3, 3, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<3, 4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 4}}, {{0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 3}}, {{2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 3, 3}}, {{0, 0, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<3, 4, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemes<4, 4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 4}}, {{0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 4}}, {{2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 3, 4}}, {{0, 0, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> OptimalSearchSchemes<4, 4, TVoidType>::VALUE;

// Given the blocklengths (absolute, not cumulative values), assign it to all
// OptimalSearches in a OptimalSearchScheme. The order of blocklength has to be from left to
// right (regarding blocks)
template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeSetBlockLength(std::array<OptimalSearch<nbrBlocks>, N> & ss,
                                               std::vector<uint32_t> const & blocklength)
{
    for (OptimalSearch<nbrBlocks> & s : ss)
        for (uint8_t i = 0; i < s.blocklength.size(); ++i)
            s.blocklength[i] = blocklength[s.pi[i]-1] + ((i > 0) ? s.blocklength[i-1] : 0);
}

// requires blocklength to be already set!
template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeInit(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    // check whether 2nd block is on the left or right and choose initialDirection accordingly
    //NOTE (svnbngk) done
    // (more efficient since we do not have to switch directions and thus have better caching performance)
    // for that we need to slightly modify search()
    for (OptimalSearch<nbrBlocks> & s : ss)
    {
        bool initialDirectionRight = s.pi[1] > s.pi[0];
        if(initialDirectionRight){
            s.startPos = 0;
            for (uint8_t i = 1; i < s.pi.size(); ++i)
                if (s.pi[i] < s.pi[0])
                    s.startPos += s.blocklength[i] - s.blocklength[i-1];
        }else{
            s.startPos = s.blocklength[s.pi.size() - 1];
            for (uint8_t i = 0; i < s.pi.size(); ++i)
                if(s.pi[i] > s.pi[0])
                    s.startPos -= s.blocklength[i] - s.blocklength[i-1];
        }
    }
}

template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeComputeFixedBlocklength(std::array<OptimalSearch<nbrBlocks>, N> & ss, uint32_t const needleLength)
{
    uint8_t blocks = ss[0].pi.size();
    uint32_t blocklength = needleLength / blocks;
    uint8_t rest = needleLength - blocks * blocklength;
    std::vector<uint32_t> blocklengths;
    for (uint8_t i = 0; i < blocks; ++i)
        blocklengths.push_back(blocklength + (i < rest));

    _optimalSearchSchemeSetBlockLength(ss, blocklengths);
    _optimalSearchSchemeInit(ss);
}

//NOTE (svnbngk) added new function to calculate added parameters from OptimalSearch
template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeComputeChronBlocklength(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    for (OptimalSearch<nbrBlocks> & s : ss){
        s.chronBL[s.pi[0] - 1]  = s.blocklength[0];
        for (int j = 1; j < nbrBlocks; ++j)
            s.chronBL[s.pi[j] - 1] = s.blocklength[j] -  s.blocklength[j - 1];
        for (int j = 1; j < nbrBlocks; ++j)
            s.chronBL[j] += s.chronBL[j - 1];
    }
}


// Compare potential occurrences directly to genome if the range on the index is small enough.
template <typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void inTextVerification(TDelegateD & delegateDirect,
                  Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                  TNeedle const & needle,
                  uint32_t const needleLeftPos,
                  uint32_t const needleRightPos,
                  uint8_t const errors,
                  OptimalSearch<nbrBlocks> const & s,
                  uint8_t const blockIndex,
                  TDir const & /**/)
{
    auto const & genome = indexText(*iter.fwdIter.index);
    uint32_t needleL = length(needle);
    uint32_t blocks = s.pi.size();

    std::vector<uint32_t> blockStarts(blocks - blockIndex);
    std::vector<uint32_t> blockEnds(blocks - blockIndex);
    for(uint32_t j = blockIndex; j < s.pi.size(); ++j){
        uint32_t blockStart = (s.pi[j] - 1 == 0) ? 0 : s.chronBL[s.pi[j] - 2];
        blockStarts[j - blockIndex] = blockStart;
        blockEnds[j - blockIndex] = s.chronBL[s.pi[j] - 1];
    }

    //modifie blockStart or blockEnd if we are already inside a block
    if(std::is_same<TDir, Rev>::value){
        if(needleRightPos - 1 > blockStarts[0] && needleRightPos - 1 < blockEnds[0])
            blockStarts[0] = needleRightPos - 1;
    }else{
        if(needleLeftPos > blockStarts[0] && needleLeftPos < blockEnds[0])
            blockEnds[0] = needleLeftPos;
    }

    for(uint32_t i = iter.fwdIter.vDesc.range.i1; i < iter.fwdIter.vDesc.range.i2; ++i){
        bool valid = true;
        Pair<uint16_t, uint32_t> sa_info = iter.fwdIter.index->sa[i];
        //dont need look at the reverse index in this case since i dont use mappability
        uint32_t chromlength = length(genome[sa_info.i1]);
        if(!(needleLeftPos <= sa_info.i2 && chromlength - 1 >= sa_info.i2 - needleLeftPos + needleL - 1))
            continue;

        sa_info.i2 = sa_info.i2 - needleLeftPos;
        uint8_t errors2 = errors;
	//iterate over each block according to search scheme
        for(uint32_t j = 0; j < blockStarts.size(); ++j){
            // compare bases to needle
            for(uint32_t k = blockStarts[j]; k <  blockEnds[j]; ++k){
                if(needle[k] != genome[sa_info.i1][sa_info.i2 + k])
                    ++errors2;
            }
            if(errors2 < s.l[blockIndex + j] || errors2 > s.u[blockIndex + j]){
                valid = false;
                break;
            }
        }
        if(valid)
            delegateDirect(sa_info, needle, errors2);
    }
}

template <typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         TCondition & itvCondition,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {
        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, Rev(), EditDistance());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, Fwd(), EditDistance());
        }
    }

    if (maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
            _optimalSearchSchemeDeletion(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, TDir());
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeChildren(TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         TCondition & itvCondition,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & /**/,
                                         TDistanceTag const & /**/)
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()),
                                   needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);

            // NOTE (cpockrandt): this might not be optimal yet! we have more edges than in the theoretical model,
            // since we go down an edge before we check whether it can even work out!
            if (!std::is_same<TDistanceTag, EditDistance>::value && minErrorsLeftInBlock > 0 &&
                charsLeft + delta < minErrorsLeftInBlock + 1u) // charsLeft - 1 < minErrorsLeftInBlock - delta
            {
                continue;
            }

            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;

            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                if (std::is_same<TDistanceTag, EditDistance>::value)
                {
                    _optimalSearchSchemeDeletion(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                    {
                        _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Rev(), TDistanceTag());
                    }
                    else
                    {
                        _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, Fwd(), TDistanceTag());
                    }
                }
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,blockIndex, TDir(), TDistanceTag());
            }

            // Deletion
            if (std::is_same<TDistanceTag, EditDistance>::value)
            {
                _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, TDir(), TDistanceTag());
            }
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeExact(TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      TCondition & itvCondition,
                                      Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const & /**/,
                                      TDistanceTag const & /**/)
{
    bool goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] :
    s.pi[blockIndex] > s.pi[blockIndex - 1];
    if (std::is_same<TDir, Rev>::value)
    {
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            return;

        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), TDistanceTag());
        }
    }
    else
    {
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightPos - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftPos - 1;

        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight], TDir()))
                return;
            --infixPosRight;
        }
        if (goToRight2)
        {
            _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                 std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), TDistanceTag());
        }
    }
}

template <typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 TCondition & itvCondition,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 TDir const & /**/,
                                 TDistanceTag const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    // Done.
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        delegate(iter, needle, errors);
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExact(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    // Approximate search in current block.
    // (s.blocklength[blockIndex]-(needleRightPos-needleLeftPos-(needleLeftPos!=needleRightPos))>=minErrorsLeftInBlock)
    else
    {
        // Insertion
        if (std::is_same<TDistanceTag, EditDistance>::value)
        {
            bool const goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t const needleRightPos2 = needleRightPos + goToRight;

            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                _optimalSearchSchemeDeletion(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, TDir());
            }
            else
            {
                _optimalSearchScheme(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, TDir(), TDistanceTag());
            }
        }
        //use lambda function to determine if In Text Search should be used
        if(itvCondition(iter, needleLeftPos, needleRightPos, errors, s, blockIndex))
        {
            inTextVerification(delegateDirect, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir());
            return;
        }

        _optimalSearchSchemeChildren(delegate, delegateDirect, itvCondition, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

template <typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 TCondition & itvCondition,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TDistanceTag const & /**/)
{
    //NOTE (svnbngk) search as long as possible in one Direction from the beginning
    if(s.pi[1] > s.pi[0])
        _optimalSearchScheme(delegate, delegateDirect, itvCondition, it, needle, s.startPos, s.startPos + 1, 0, s, 0, Rev(), TDistanceTag());
    else
        _optimalSearchScheme(delegate, delegateDirect, itvCondition, it, needle, s.startPos, s.startPos + 1, 0, s, 0, Fwd(), TDistanceTag());
}

template <typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 TCondition & itvCondition,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const & /**/)
{
    for (auto & s : ss)
        _optimalSearchScheme(delegate, delegateDirect, itvCondition, it, needle, s, TDistanceTag());
}

//NOTE added this function so find function can be run on single needle without ITV condition
template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(TDelegate & delegate,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 TNeedle const & needle,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const & /**/)
{
    auto dummy = [](auto /*iter*/, uint32_t /*needleLeftPos*/, uint32_t /*needleRightPos*/, uint8_t /*errors*/, auto /*s*/, uint8_t const /*blockIndex*/)
    {
        return(false);
    };
    auto delegateDummy = [](Pair<uint16_t, uint32_t> const & /*pos*/, DnaString const & /*needle*/, uint8_t const /*errors*/)
    {
    };

    for (auto & s : ss)
        _optimalSearchScheme(delegate, delegateDummy, dummy, it, needle, s, TDistanceTag());
}


// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndexSpec,
          typename TChar, typename TStringSpec,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     TCondition & itvCondition,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     String<TChar, TStringSpec> const & needle,
     TDistanceTag const & /**/)
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(delegate,  delegateDirect, itvCondition, it, needle, scheme, TDistanceTag());
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag,
          typename TParallelTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     TCondition & itvCondition,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & /**/,
     TParallelTag const & /**/)
{
    typedef typename Iterator<StringSet<TNeedle, TStringSetSpec> const, Rooted>::Type TNeedleIt;
    typedef typename Reference<TNeedleIt>::Type                                       TNeedleRef;
    iterate(needles, [&](TNeedleIt const & needleIt)
    {
        TNeedleRef needle = value(needleIt);
        find<minErrors, maxErrors>(delegate, delegateDirect, itvCondition, index, needle, TDistanceTag());
    },
    Rooted(), TParallelTag());
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <size_t minErrors, size_t maxErrors,
          typename TDelegate, typename TDelegateD, typename TCondition,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     TDelegateD & delegateDirect,
     TCondition & itvCondition,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & /**/)
{
    if (std::is_same<TDistanceTag, EditDistance>::value)
        std::cerr << "In Text Verification was not implemented for EditDistance" << "\n";
    find<minErrors, maxErrors>(delegate, delegateDirect, itvCondition, index, needles, TDistanceTag(), Serial());
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <size_t minErrors, size_t maxErrors,
          typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & /**/)
{
    auto dummy = [](auto /*iter*/, uint32_t /*needleLeftPos*/, uint32_t /*needleRightPos*/, uint8_t /*errors*/, auto s, uint8_t const /*blockIndex*/)
    {
        return(false);
    };

       auto delegateDummy = [](Pair<uint16_t, uint32_t> const & /*pos*/, DnaString const & /*needle*/, uint8_t const /*errors*/)
    {
    };
    find<minErrors, maxErrors>(delegate, delegateDummy, dummy, index, needles, TDistanceTag(), Serial());
}

}

#endif  // #ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_H_
