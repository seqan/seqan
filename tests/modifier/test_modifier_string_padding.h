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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_MODIFIER_MODIFIER_STRING_PADDING_H_
#define TESTS_MODIFIER_MODIFIER_STRING_PADDING_H_

SEQAN_DEFINE_TEST(test_modified_string_padding_construction)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    ModifiedString<DnaString, ModPadding> modString(seq);

    SEQAN_ASSERT(modString._host == &seq);
    SEQAN_ASSERT(modString._cargo._paddedValue == 'A');
    SEQAN_ASSERT(modString._cargo._numPaddedChar == 0u);
    SEQAN_ASSERT(modString._cargo._remainingSteps == 0u);
}

SEQAN_DEFINE_TEST(test_modified_string_padding_expand)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    ModifiedString<DnaString, ModPadding> modString(seq);
    ModifiedString<DnaString, ModPadding> modString2(seq);

    expand(modString, 5);
    expand(modString2, 5, 'C');

    for (unsigned i = 0; i < 5; ++i)
    {
        SEQAN_ASSERT(modString[length(seq) + i] == 'A');
        SEQAN_ASSERT(modString2[length(seq) + i] == 'C');
    }
}

SEQAN_DEFINE_TEST(test_modified_string_padding_length)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    ModifiedString<DnaString, ModPadding> modString(seq);

    SEQAN_ASSERT(length(modString) == length(seq));
    expand(modString, 5);
    SEQAN_ASSERT(length(modString) == length(seq) + 5);
}

SEQAN_DEFINE_TEST(test_modified_string_padding_begin)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    ModifiedString<DnaString, ModPadding> modString(seq);

    auto it = begin(modString, Standard());
    SEQAN_ASSERT(host(it) == begin(seq, Rooted()));
    SEQAN_ASSERT_EQ(cargo(it)._remainingSteps, 0u);
    SEQAN_ASSERT_EQ(cargo(it)._numPaddedChar, 0u);

    expand(modString, 5);
    it = begin(modString, Standard());
    SEQAN_ASSERT(host(it) == begin(seq, Rooted()));
    SEQAN_ASSERT_EQ(cargo(it)._remainingSteps, 5u);
    SEQAN_ASSERT_EQ(cargo(it)._numPaddedChar, 5u);
}

SEQAN_DEFINE_TEST(test_modified_string_padding_end)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    ModifiedString<DnaString, ModPadding> modString(seq);

    auto itEnd = end(modString, Standard());
    SEQAN_ASSERT(host(itEnd) == end(seq, Rooted()));
    SEQAN_ASSERT_EQ(cargo(itEnd)._remainingSteps, 0u);
    SEQAN_ASSERT_EQ(cargo(itEnd)._numPaddedChar, 0u);

    expand(modString, 5);
    itEnd = end(modString, Standard());
    SEQAN_ASSERT(host(itEnd) == end(seq, Rooted()));
    SEQAN_ASSERT_EQ(cargo(itEnd)._remainingSteps, 0u);
    SEQAN_ASSERT_EQ(cargo(itEnd)._numPaddedChar, 5u);
}

SEQAN_DEFINE_TEST(test_modified_string_padding_difference)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    ModifiedString<DnaString, ModPadding> modString(seq);

    expand(modString, 5);

    auto itB = begin(modString, Standard());
    auto itE = end(modString, Standard());

    auto pos = 0;
    for (auto it = itB; it != itE; ++it, ++pos)
        SEQAN_ASSERT(it - itB == pos);
}

SEQAN_DEFINE_TEST(test_modified_string_padding_iterator)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    ModifiedString<DnaString, ModPadding> modString(seq);

    expand(modString, 5, static_cast<Dna>('T'));

    // Test begin, end, increment, decrement and dereference
    auto itSeq = begin(seq, Standard());
    auto it = begin(modString, Standard());
    for (; it != end(modString, Standard()); ++it)
    {
        if (itSeq < end(seq, Standard()))
            SEQAN_ASSERT_EQ(*it, *(itSeq++));
        else
            SEQAN_ASSERT_EQ(*it, static_cast<Dna>('T'));
    }

    while (it != begin(modString, Standard()))
    {
        --it;
        if (it - begin(modString, Standard()) < static_cast<typename Difference<decltype(it)>::Type>(length(seq)))
            SEQAN_ASSERT_EQ(*it, *(--itSeq));
        else
            SEQAN_ASSERT_EQ(*it, static_cast<Dna>('T'));
    }

    // Test advance.
    SEQAN_ASSERT(it == begin(modString, Standard()));
    it += 4;
    SEQAN_ASSERT_EQ(*it, seq[4]);
    it -= 2;
    SEQAN_ASSERT_EQ(*it, seq[2]);
    it += 13;
    SEQAN_ASSERT_EQ(*it, static_cast<Dna>('T'));
    it -= 15;
    SEQAN_ASSERT_EQ(*it, seq[0]);
    it = it + 17;
    SEQAN_ASSERT_EQ(*it, static_cast<Dna>('T'));
    it += 3;
    SEQAN_ASSERT(it == end(modString, Standard()));
    it -= 5;
    SEQAN_ASSERT_EQ(*it, static_cast<Dna>('T'));
    SEQAN_ASSERT_EQ(*(it - 3), seq[12]);
}

SEQAN_DEFINE_TEST(test_modified_string_padding_defect_2190)
{
    using namespace seqan;

    DnaString seq = "ACGTGGATAGCATCG";
    auto seqInf = infix(seq, 0, length(seq));

    auto test_const = [](auto const & modifier)
    {
        // using TRef = typename Reference<decltype(modifier)>::Type;
        auto x = value(modifier, 1);
        SEQAN_ASSERT_EQ(x, 'C');
    };
    { // Test working case, when reference of const modifier gives back a const reference to the value
        ModifiedString<decltype(seq), ModPadding> modString(seq);
        test_const(modString);
    }

    {  // Test defect, when reference of const modifier gives back a non-const reference to the value
        ModifiedString<decltype(seqInf), ModPadding> modString(seqInf);
        test_const(modString);
    }
}

#endif  // #ifndef TESTS_MODIFIER_MODIFIER_STRING_PADDING_H_
