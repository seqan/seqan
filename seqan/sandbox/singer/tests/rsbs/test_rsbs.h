// ==========================================================================
//                                    rsbs
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_SINGER_TESTS_RSBS_TEST_RSBS_H_
#define SANDBOX_SINGER_TESTS_RSBS_TEST_RSBS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>
#include <seqan/fm_sequence.h>

using namespace seqan;
const int SEED = 42;

// A test for strings.
SEQAN_DEFINE_TEST(test_rsbs_defaultConstructor)
{
    typedef RankSupportBitString<> TRankSupportBitString;
    TRankSupportBitString bitString;

    SEQAN_ASSERT_EQ(length(bitString), 0u);
    SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreBits())), 0u);
    SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreBlocks())), 0u);
    SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreSuperBlocks())), 0u);
}

SEQAN_DEFINE_TEST(test_rsbs_resize)
{
    typedef RankSupportBitString<> TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type         TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Value<TFibreBits>::Type               TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type              TFibreBlocksValue;

    TRankSupportBitString bitString;
    TFibreBlocksValue const bitsPerValue_ = BitsPerValue<TFibreBitsValue>::VALUE;

    resize(bitString, 0u);
    SEQAN_ASSERT_EQ(length(bitString), 0u);
    SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreBits())), 0u);
    SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreBlocks())), 0u);
    SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreSuperBlocks())), 0u);
    for (unsigned i = 1; i < 100000; ++i)
    {
        resize(bitString, i);
        SEQAN_ASSERT_EQ(length(bitString), i);
        SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreBits())), (i + bitsPerValue_ - 1) / bitsPerValue_);
        SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreBlocks())), (i + bitsPerValue_ - 1) / bitsPerValue_);
        SEQAN_ASSERT_EQ(length(getFibre(bitString, FibreSuperBlocks())), ((i + bitsPerValue_ - 1) / bitsPerValue_ + bitsPerValue_ - 1) / bitsPerValue_);
    }

    resize(bitString, 0);
    SEQAN_ASSERT_EQ(length(bitString), 0u);

    resize(bitString, 100000);
    for (unsigned i = 0; i < length(bitString); ++i)
        SEQAN_ASSERT_EQ(getRank(bitString, i), 0u);
    
    resize(bitString, 200000, 1);
    for (unsigned i = 0; i < 100000; ++i)
        SEQAN_ASSERT_EQ(getRank(bitString, i), 0u);
    for (unsigned i = 100000; i < length(bitString); ++i)
        SEQAN_ASSERT_EQ(getRank(bitString, i), i - 99999);
}

SEQAN_DEFINE_TEST(test_rsbs_getBuPos)
{
    typedef RankSupportBitString<> TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type         TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Value<TFibreBits>::Type               TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type              TFibreBlocksValue;

    TRankSupportBitString bitString;
    TFibreBlocksValue const bitsPerValue_ = BitsPerValue<TFibreBitsValue>::VALUE;

    resize(bitString, 100000);
    for (unsigned i = 0; i < length(bitString); ++i)
    {
        SEQAN_ASSERT_EQ(getBlockPos_(bitString, i), i / bitsPerValue_);
    }
}

SEQAN_DEFINE_TEST(test_rsbs_getSBuPos)
{
    typedef RankSupportBitString<> TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type         TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Value<TFibreBits>::Type               TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type              TFibreBlocksValue;

    TRankSupportBitString bitString;
    TFibreBlocksValue const bitsPerValue_ = BitsPerValue<TFibreBitsValue>::VALUE;

    resize(bitString, 100000);
    for (unsigned i = 0; i < length(bitString); ++i)
    {
        SEQAN_ASSERT_EQ(getSuperBlockPos_(bitString, i), i / (bitsPerValue_ * bitsPerValue_));
    }

}

SEQAN_DEFINE_TEST(test_rsbs_getPosInBu)
{
    typedef RankSupportBitString<> TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type         TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Value<TFibreBits>::Type               TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type              TFibreBlocksValue;

    TRankSupportBitString bitString;
    TFibreBlocksValue const bitsPerValue_ = BitsPerValue<TFibreBitsValue>::VALUE;

    resize(bitString, 100000);
    for (unsigned i = 0; i < length(bitString); ++i)
    {
        SEQAN_ASSERT_EQ(getPosInBlock_(bitString, i), i % bitsPerValue_);
    }

}


// A test for strings.
SEQAN_DEFINE_TEST(test_rsbs_getSetBit)
{
    unsigned length_ = 100000;
    RankSupportBitString<> bitString;
    resize(bitString, length_);

    String<bool> controlBitString;
    resize(controlBitString, length_);


    Rng<MersenneTwister> rng(SEED);
    for (unsigned i = 0; i < length(bitString); ++i)
    {
        controlBitString[i] = pickRandomNumber(rng) % 2;
        setBit(bitString, i, controlBitString[i]);
        SEQAN_ASSERT_EQ(getBit(bitString, i), controlBitString[i]);
    }
}

// A test for strings.
SEQAN_DEFINE_TEST(test_rsbs_append_value)
{
    unsigned length_ = 100000;
    RankSupportBitString<> bitString;

    String<bool> controlBitString;
    resize(controlBitString, length_);

    Rng<MersenneTwister> rng(SEED);
    for (unsigned i = 0; i < length_; ++i)
    {
        controlBitString[i] = pickRandomNumber(rng) % 2;
        appendValue(bitString, controlBitString[i]);
        SEQAN_ASSERT_EQ(length(bitString), i + 1);
        SEQAN_ASSERT_EQ(getBit(bitString, i), controlBitString[i]);
    }
}

// A test for strings.
SEQAN_DEFINE_TEST(test_rsbs_rank)
{
    unsigned length_ = 100000;
    RankSupportBitString<> bitString;

    String<unsigned int> controlRankString;
    resize(controlRankString, length_);

    Rng<MersenneTwister> rng(SEED);
    bool bit = pickRandomNumber(rng) % 2;
    controlRankString[0] = bit;
    appendValue(bitString, bit);
    SEQAN_ASSERT_EQ(length(bitString), 1u);
    SEQAN_ASSERT_EQ(getRank(bitString, 0), controlRankString[0]);

    for (unsigned i = 1; i < length_; ++i)
    {
        bit = pickRandomNumber(rng) % 2;
        controlRankString[i] = controlRankString[i - 1] + bit;
        appendValue(bitString, bit);
        SEQAN_ASSERT_EQ(length(bitString), i + 1);
        SEQAN_ASSERT_EQ(getRank(bitString, i), controlRankString[i]);
    }
}

// Note: updateRank_ is tested in test_rsbs_append_value
SEQAN_DEFINE_TEST(test_rsbs_update_ranks_)
{
    unsigned length_ = 100000;
    RankSupportBitString<> bitString;
    resize(bitString, length_);

    String<bool> controlBitString;
    resize(controlBitString, length_);

    String<unsigned> controlRankString;
    resize(controlRankString, length_);

    Rng<MersenneTwister> rng(SEED);
    controlBitString[0] = pickRandomNumber(rng) % 2;
    controlRankString[0] = controlBitString[0];
    setBit(bitString, 0, controlBitString[0]);
    updateRanks_(bitString);
    SEQAN_ASSERT_EQ(getRank(bitString, 0), controlRankString[0]);

    for (unsigned i = 1; i < length(bitString); ++i)
    {
        controlBitString[i] = pickRandomNumber(rng) % 2;
        setBit(bitString, i, controlBitString[i]);
        updateRanks_(bitString);
        controlRankString[i] = controlRankString[i - 1] + controlBitString[i];
        SEQAN_ASSERT_EQ(getRank(bitString, i), controlRankString[i]);
    }

    for (unsigned i = length(bitString) / 2 ; i < length(bitString); ++i)
    {
        controlBitString[i] = pickRandomNumber(rng) % 2;
        setBit(bitString, i, controlBitString[i]);
        controlRankString[i] = controlRankString[i - 1] + controlBitString[i];
    }
    
    updateRanks_(bitString, 0);
    for (unsigned i = 0 / 2 ; i < length(bitString); ++i)
    {
        SEQAN_ASSERT_EQ(getRank(bitString, i), controlRankString[i]);
    }
}

SEQAN_DEFINE_TEST(test_rsbs_constructor)
{
    unsigned length_ = 100000;
    String<unsigned> controlString;
    resize(controlString, length_);

    String<unsigned> controlRankString;
    resize(controlRankString, length_);

    Rng<MersenneTwister> rng(SEED);
    controlString[0] = pickRandomNumber(rng) % 2;
    controlRankString[0] = controlString[0];

    for (unsigned i = 1; i < length_; ++i)
    {
        controlString[i] = pickRandomNumber(rng) % 2;
        controlRankString[i] = controlRankString[i - 1] + controlString[i];
    }

    typedef RankSupportBitString<> TRankSupportBitString;
    TRankSupportBitString bitString(controlString);

    SEQAN_ASSERT_EQ(getRank(bitString, 0), controlRankString[0]);
    for (unsigned i = 1; i < length(bitString); ++i)
    {
        SEQAN_ASSERT_EQ(getRank(bitString, i), controlRankString[i]);
    }
}

SEQAN_DEFINE_TEST(test_rsbs_equalOperator)
{
    unsigned length_ = 10000;
    RankSupportBitString<> bitString;
	RankSupportBitString<> otherBitString;

    String<bool> controlBitString;
    resize(controlBitString, length_);

    Rng<MersenneTwister> rng(SEED);
    for (unsigned i = 0; i < length_; ++i)
    {
    	controlBitString[i] = pickRandomNumber(rng) % 2;
    	appendValue(bitString, controlBitString[i]);
    	appendValue(otherBitString, controlBitString[i]);
    }

	SEQAN_ASSERT_EQ(bitString.bits, otherBitString.bits);
	SEQAN_ASSERT_EQ(bitString.blocks, otherBitString.blocks);
	SEQAN_ASSERT_EQ(bitString.superBlocks, otherBitString.superBlocks);
	SEQAN_ASSERT_EQ(bitString.length_, otherBitString.length_);
}


SEQAN_DEFINE_TEST(test_rsbs_assignOperator)
{
    unsigned length_ = 100000;
    RankSupportBitString<> bitString;

    String<bool> controlBitString;
    resize(controlBitString, length_);

    Rng<MersenneTwister> rng(SEED);
    for (unsigned i = 0; i < length_; ++i)
    {
    	controlBitString[i] = pickRandomNumber(rng) % 2;
    	appendValue(bitString, controlBitString[i]);
    }

	RankSupportBitString<> otherBitString = bitString;
	SEQAN_ASSERT(bitString == otherBitString);
}

#endif  // SANDBOX_SINGER_TESTS_RSBS_TEST_RSBS_H_
