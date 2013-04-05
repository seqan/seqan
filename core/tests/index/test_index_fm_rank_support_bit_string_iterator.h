// ==========================================================================
//                                    rsbs
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef TEST_INDEX_FM_RANK_SUPPORT_BIT_STRING_ITERATOR_H_
#define TEST_INDEX_FM_RANK_SUPPORT_BIT_STRING_ITERATOR_H_

using namespace seqan;

SEQAN_DEFINE_TEST(test_rsbs_iterator_get_value)
{
    unsigned length_ = 1000;
    String<unsigned> controlString;
    resize(controlString, length_);

    Rng<MersenneTwister> rng(SEED);
    controlString[0] = pickRandomNumber(rng) % 2;

    for (unsigned i = 1; i < length_; ++i)
    {
        controlString[i] = pickRandomNumber(rng) % 2;
    }

    typedef RankSupportBitString<void> TRankSupportBitString;
    TRankSupportBitString bitString(controlString);
    Iterator<RankSupportBitString<void> >::Type defaultIt = begin(bitString);

    for (unsigned i = 0; i < length(bitString); ++i)
    {
        SEQAN_ASSERT(isBitSet(bitString, i) == getValue(defaultIt));
        SEQAN_ASSERT(isBitSet(bitString, i) == getValue(defaultIt, Bit()));
        SEQAN_ASSERT(getRank(bitString, i) == getValue(defaultIt, Rank()));
        ++defaultIt;
    }
}

#endif  // TEST_INDEX_FM_RANK_SUPPORT_BIT_STRING_H_
