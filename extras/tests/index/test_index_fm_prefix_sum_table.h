// ==========================================================================
//                               fm_index_beta
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

#ifndef TEST_PREFIX_SUM_TABLE_BETA_H_
#define TEST_PREFIX_SUM_TABLE_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
template <typename TPrefixSumTable>
void prefixSumTableConstructor(TPrefixSumTable & /*tag*/)
{
	typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef String<TChar> TText;

    int minValue = MinValue<TChar>::VALUE;

    TText text;
    generateText(text);

    TPrefixSumTable pst(text);

    Index<TText> index(text);
    typename Iterator<Index<TText>, TopDown<ParentLinks<> > >::Type it(index);

    unsigned sum = 0;
    for (int i = 1; (unsigned)i < length(pst) - 1; ++i)
    {
        goDown(it, (TChar)(i - 1 - minValue));
        sum += countOccurrences(it);
        SEQAN_ASSERT_EQ(getPrefixSum(pst, i), sum); 
        goUp(it);
    }
    SEQAN_ASSERT_EQ(getPrefixSum(pst, length(pst) - 1), length(text)); 
}

template <typename TPrefixSumTable>
void prefixSumTableGetAlphabetSize(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef String<TChar> TText;
    
    TText text = "A";
    TPrefixSumTable pst(text);
    SEQAN_ASSERT(getAlphabetSize(pst) == (unsigned)ValueSize<TChar>::VALUE);
}

template <typename TPrefixSumTable>
void prefixSumTableGetCharacterPosition(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef String<TChar> TText;

    int minValue = MinValue<TChar>::VALUE;
    int maxValue = MaxValue<TChar>::VALUE;

    TText text = "A";
    TPrefixSumTable pst(text);
    for (int i = minValue; i <= maxValue; ++i)
        SEQAN_ASSERT_EQ(getCharacterPosition(pst, (TChar)i),(unsigned)(i - minValue));
}

template <typename TPrefixSumTable>
void prefixSumTableGetCharacter(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    int minValue = MinValue<TChar>::VALUE;
    int maxValue = MaxValue<TChar>::VALUE;

    TText text = "A";
    TPrefixSumTable pst(text);
    for (int i = minValue; i <= maxValue; ++i)
    {
        SEQAN_ASSERT_EQ(getCharacter(pst, i), (TUChar)i);
    }
}

template <typename TPrefixSumTable>
void prefixSumTableGetPivotPosition(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef String<TChar> TText;

    {
        TText text;
        TPrefixSumTable pst(text);
        unsigned pivotPos = _getPivotPosition(pst, 0u, ValueSize<TChar>::VALUE - 1);
        SEQAN_ASSERT_EQ(pivotPos, 1u); 
    }
    {
        TText text = "AAA";
        TPrefixSumTable pst(text);

        unsigned pivotPos = _getPivotPosition(pst, 0u, ValueSize<TChar>::VALUE - 1);
        SEQAN_ASSERT_EQ(pivotPos, 1u); 
    }
    {
        TText text;
        generateText(text);
        TPrefixSumTable pst(text);
        for (unsigned i = 0; i < ValueSize<TChar>::VALUE - 2; ++i)
        {
            unsigned pivotPos = _getPivotPosition(pst, i, ValueSize<TChar>::VALUE - 1);
            long currentSum = std::abs((long)(pst[pivotPos] - pst[i]) - (long)(pst[ValueSize<TChar>::VALUE] - pst[pivotPos]));
            if(pivotPos > 0)
            {
                long leftSum = std::abs((long)(pst[pivotPos - 1] - pst[i]) - (long)(pst[ValueSize<TChar>::VALUE] - pst[pivotPos - 1]));
                SEQAN_ASSERT_LEQ(currentSum, leftSum);
            }
            if(pivotPos < ValueSize<TChar>::VALUE)
            {
                long rightSum = std::abs((long)(pst[pivotPos + 1] - pst[i]) - (long)(pst[ValueSize<TChar>::VALUE] - pst[pivotPos + 1]));
                SEQAN_ASSERT_LEQ(currentSum, rightSum);
            }
        }
        for (unsigned i = 0; i < ValueSize<TChar>::VALUE - 2; ++i)
        {
            unsigned pivotPos = _getPivotPosition(pst, 0u, ValueSize<TChar>::VALUE - 1 - i);
            long currentSum = std::abs((long)(pst[pivotPos] - pst[0]) - (long)(pst[ValueSize<TChar>::VALUE - i] - pst[pivotPos]));
            if(pivotPos > 0)
            {
                long leftSum = std::abs((long)(pst[pivotPos - 1] - pst[0]) - (long)(pst[ValueSize<TChar>::VALUE - i] - pst[pivotPos - 1]));
                SEQAN_ASSERT_LEQ(currentSum, leftSum);
            }
            if(pivotPos < ValueSize<TChar>::VALUE - 1 - i)
            {
                long rightSum = std::abs((long)(pst[pivotPos + 1] - pst[0]) - (long)(pst[ValueSize<TChar>::VALUE - i] - pst[pivotPos + 1]));
                SEQAN_ASSERT_LEQ(currentSum, rightSum);
            }
        }

    }
}

template <typename TPrefixSumTable>
void prefixSumTableDetermineDollarSubstitute(TPrefixSumTable & /*tag*/)
{
    typedef typename RemoveConst<TPrefixSumTable>::Type TNonConstPrefixSumTable;
	typedef typename CharacterValue<TNonConstPrefixSumTable>::Type TChar;

	{
		String<TChar> text = "ACGTNACGTNACGTNNN";
		TPrefixSumTable prefixSumTable(text);
		TChar character;
		_determineDollarSubstitute(prefixSumTable, character);
		SEQAN_ASSERT_EQ(character, TChar('A'));
	}
	{
		String<TChar> text = "AGTNAGTNAGTNNN";
		TPrefixSumTable prefixSumTable(text);
		TChar character;
		_determineDollarSubstitute(prefixSumTable, character);
		SEQAN_ASSERT_EQ(character, TChar('A'));
	}
	{
		String<TChar> text = "ACGTNAGTNAGTNNN";
		TPrefixSumTable prefixSumTable(text);
		TChar character;
		_determineDollarSubstitute(prefixSumTable, character);
		SEQAN_ASSERT_EQ(character, TChar('C'));
	}
	{
		String<TChar> text = "ACGTNACGTNACGT";
		TPrefixSumTable prefixSumTable(text);
		TChar character;
		_determineDollarSubstitute(prefixSumTable, character);
		SEQAN_ASSERT_EQ(character, TChar('N'));
	}
}

template <typename TPrefixSumTable>
void prefixSumTableGetValue(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);

    typename Reference<TPrefixSumTable>::Type temp = prefixSum(pst, 0);
    temp = 11;

    SEQAN_ASSERT_EQ(getValue(pst, 0), 11u);
}



template <typename TPrefixSumTable>
void _prefixSumTableInsertDollar(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);

    TPrefixSumTable const pstConst = pst;

    unsigned numDollar = 10;
    _insertDollar(pst, numDollar);

    for (unsigned i = 0; i < length(pst); ++i)
        SEQAN_ASSERT_EQ(getPrefixSum(pst, i), getPrefixSum(pstConst, i) + numDollar);
}


template <typename TPrefixSumTable>
void prefixSumTableLength(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);

    SEQAN_ASSERT_EQ(length(pst), length(pst.entries));
}

template <typename TPrefixSumTable>
void prefixSumTablePrefixSum(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);

    typename Reference<TPrefixSumTable>::Type temp = prefixSum(pst, 0);
    temp = 11;

    SEQAN_ASSERT_EQ(getPrefixSum(pst, 0), 11u);
}

template <typename TPrefixSumTable>
void prefixSumTableResize(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);
    resize(pst, 1u);

    SEQAN_ASSERT_EQ(length(pst), 1u);
}

template <typename TPrefixSumTable>
void prefixSumTableSetPrefixSum(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);

    SEQAN_ASSERT_NEQ(getPrefixSum(pst, 0), 100u);
    setPrefixSum(pst, 100, 0);
    SEQAN_ASSERT_EQ(getPrefixSum(pst, 0), 100u);
}

template <typename TPrefixSumTable>
void prefixSumTableValue(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);

    typename Reference<TPrefixSumTable>::Type temp = value(pst, 0);
    temp = 11;

    SEQAN_ASSERT_EQ(getPrefixSum(pst, 0), 11u);
}

template <typename TPrefixSumTable>
void prefixSumTableOpenSave(TPrefixSumTable & /*tag*/)
{
    typedef typename CharacterValue<TPrefixSumTable>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef String<TChar> TText;

    TText text = "ACGT";
    TPrefixSumTable pst(text);

    CharString tempFilename = SEQAN_TEMP_FILENAME();
    save(pst, toCString(tempFilename));

    TPrefixSumTable openPst;
    open(openPst, toCString(tempFilename));

    SEQAN_ASSERT(pst == openPst);
}

SEQAN_DEFINE_TEST(prefix_sum_table_constructor)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableConstructor(tag);

    PrefixSumTable<signed char, void> charTag;
    prefixSumTableConstructor(charTag);

    //PrefixSumTable<Dna5, void> const constTag;
    //prefixSumTableConstructor(constTag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_get_alphabet_size)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableGetAlphabetSize(tag);

    PrefixSumTable<signed char, void> charTag;
    prefixSumTableGetAlphabetSize(charTag);

    //PrefixSumTable<Dna5, void> const constTag;
    //prefixSumTableConstructor(constTag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_get_character_position)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableGetCharacterPosition(tag);

    PrefixSumTable<signed char, void> charTag;
    prefixSumTableGetCharacterPosition(charTag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_get_character)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableGetCharacter(tag);

    PrefixSumTable<signed char, void> charTag;
    prefixSumTableGetCharacter(charTag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_get_pivot_position)
{
    using namespace seqan;

    PrefixSumTable<signed char, void> tag;
    prefixSumTableGetPivotPosition(tag);

    //PrefixSumTable<Dna5, void> const constTag;
    //prefixSumTableGetPivotPosition(constTag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_determine_dollar_substitute)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableDetermineDollarSubstitute(tag);

    PrefixSumTable<Dna5, void> const constTag;
    prefixSumTableDetermineDollarSubstitute(constTag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_get_value)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableGetValue(tag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_insert_dollar_)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    _prefixSumTableInsertDollar(tag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_length)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableLength(tag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_prefix_sum)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTablePrefixSum(tag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_resize)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableResize(tag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_set_prefix_sum)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableSetPrefixSum(tag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_value)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableValue(tag);
}

SEQAN_DEFINE_TEST(prefix_sum_table_open_save)
{
    using namespace seqan;

    PrefixSumTable<Dna5, void> tag;
    prefixSumTableOpenSave(tag);
}

#endif  // TEST_PREFIX_SUM_TABLE_BETA_H_
