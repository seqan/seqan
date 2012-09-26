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

#ifndef TEST_COMPRESSED_SA_BETA_H_
#define TEST_COMPRESSED_SA_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/fm_sequence.h>
#include "test_fm_index_beta.h"

using namespace seqan;

// template <typename TCompressedSA>
// void compressedSaCompressionFactor(TCompressedSA & /*tag*/)
// { 
//     TCompressedSA compressedSA;
// 
//     assignCompressionFactor(compressedSA, 3u);
// 
//     SEQAN_ASSERT_EQ(getCompressionFactor(compressedSA), 3u);
// }

// template <typename TCompressedSA>
// void compressedSaAssignValue(TCompressedSA & /*tag*/)
// { 
//     TCompressedSA compressedSA;
//     assignCompressionFactor(compressedSA, 10u);
//     resize(compressedSA, 30);
// 
//     assignValue(compressedSA, 0, 0);
//     assignValue(compressedSA, 1, 2);
//     assignValue(compressedSA, 2, 3);
// 
//     SEQAN_ASSERT_EQ(getValue(getFibre(compressedSA, FibreSparseString()), 0), 0u);
//     SEQAN_ASSERT_EQ(getValue(getFibre(compressedSA, FibreSparseString()), 1), 2u);
//     SEQAN_ASSERT_EQ(getValue(getFibre(compressedSA, FibreSparseString()), 2), 3u);
// }

template <typename TCompressedSA>
void compressedSaClearLengthResize(TCompressedSA & /*tag*/)
{ 
    TCompressedSA compressedSA;
    SEQAN_ASSERT_EQ(length(compressedSA), 0u);

    resize(compressedSA, 30);
    SEQAN_ASSERT_EQ(length(compressedSA), 30u);

    clear(compressedSA);
    SEQAN_ASSERT_EQ(length(compressedSA), 0u);
}

template <typename TCompressedSA>
void compressedSaEmpty(TCompressedSA & /*tag*/)
{ 
    TCompressedSA compressedSA;

    SEQAN_ASSERT_EQ(empty(compressedSA), true);

    resize(compressedSA, 30);

    SEQAN_ASSERT_EQ(empty(compressedSA), false);

    clear(compressedSA);
    SEQAN_ASSERT_EQ(empty(compressedSA), true);
}

template <typename TCompressedSA>
void compressedSaCompressedSaCreate(TCompressedSA & /*tag*/)
{ 
    TCompressedSA compressedSA;


    String<unsigned> fullSA;
    appendValue(fullSA, 8);
    appendValue(fullSA, 10);
    appendValue(fullSA, 6);
    appendValue(fullSA, 5);
    appendValue(fullSA, 2);
    appendValue(fullSA, 9);
    appendValue(fullSA, 1);
    appendValue(fullSA, 0);
    appendValue(fullSA, 7);
    appendValue(fullSA, 3);

    compressedSaCreate(compressedSA, fullSA, 3u);

    SEQAN_ASSERT_EQ(length(getFibre(getFibre(compressedSA, FibreSparseString()), FibreValueString())), 4u);

    SEQAN_ASSERT_EQ(getValue(getFibre(compressedSA, FibreSparseString()), 2), 6u);
    SEQAN_ASSERT_EQ(getValue(getFibre(compressedSA, FibreSparseString()), 5), 9u);
    SEQAN_ASSERT_EQ(getValue(getFibre(compressedSA, FibreSparseString()), 7), 0u);
    SEQAN_ASSERT_EQ(getValue(getFibre(compressedSA, FibreSparseString()), 9), 3u);

    typedef typename Fibre<typename Fibre<TCompressedSA, FibreSparseString>::Type, FibreIndicatorString>::Type TIndicatorString;

    TIndicatorString & indicatorString = getFibre(getFibre(compressedSA, FibreSparseString()), FibreIndicatorString());

    SEQAN_ASSERT_EQ(getBit(indicatorString, 0), false);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 1), false);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 2), true);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 3), false);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 4), false);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 5), true);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 6), false);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 7), true);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 8), false);
    SEQAN_ASSERT_EQ(getBit(indicatorString, 9), true);
}

template <typename TCompressedSA>
void compressedSaGetFibre(TCompressedSA & /*tag*/)
{ 
    TCompressedSA compressedSA;

    resize(compressedSA, 3);
    SEQAN_ASSERT(getFibre(compressedSA, FibreSparseString()) == compressedSA.sparseString);
}

template <typename TIndex>
void compressedSaGetNextPos_(TIndex & /*tag*/)
{ 
    typedef typename Value<TIndex>::Type TText;
    typedef typename SAValue<TText>::Type TSAValue;
    typedef typename Fibre<TIndex, FibreSA>::Type TCompressedSA;
    typedef typename Fibre<TIndex, FibreLfTable>::Type TLfTable;
    typedef typename Fibre<TLfTable, FibreOccTable>::Type TOccTable;

    typedef String<TSAValue> TSAString;

    TText text;
    generateText(text);

    TSAString sa;
    resize(sa, length(text));
    createSuffixArray(sa, text, Skew7());

    TIndex index(text, 3);

    unsigned pos, pos2;
    TCompressedSA & compressedSA = getFibre(index, FibreSA());
    TOccTable occTable = getFibre(getFibre(index, FibreLfTable()), FibreOccTable());

    //static_cast<Nothing>(occTable);

    for(unsigned i = 1; i < length(text); ++i)
    {
        if (!dollarPosition(occTable, i))
        {
            pos = i;
            pos2 = pos;
            while(compressedSA[pos] % 3 != 0)
            {
                SEQAN_ASSERT(getNextPos_(compressedSA, pos) == false);
                pos2 = lfMapping(getFibre(index, FibreLfTable()), pos2);
            }
            SEQAN_ASSERT(getNextPos_(compressedSA, pos) == true);
            SEQAN_ASSERT_EQ(pos, pos2);
        }
    }
}

template <typename TCompressedSA>
void compressedSaSetLfTable(TCompressedSA & /*tag*/)
{ 
    TCompressedSA compressedSA;

    CharString dummyLfTable;

    setLfTable(compressedSA, dummyLfTable);

    SEQAN_ASSERT(compressedSA.lfTable == &dummyLfTable);
}

template <typename TIndex>
void compressedSaValueAccess(TIndex & /*tag*/)
{
    typedef typename Value<TIndex>::Type TText;
    typedef typename SAValue<TText>::Type TSAValue;

    typedef String<TSAValue> TSAString;

    TText text;
    generateText(text);

    TSAString sa;
    resize(sa, length(text));
    createSuffixArray(sa, text, Skew7());

    TIndex index(text);

    for(unsigned i = 0; i < length(text); ++i)
    {
        SEQAN_ASSERT_EQ(sa[i], getFibre(index, FibreSA())[i+1]);
    }

}

// SEQAN_DEFINE_TEST(compressed_sa_compression_factor)
// {
//     using namespace seqan;
// 
//     CompressedSA<SparseString<String<unsigned int>, void >, CharString, void> tag;
// 
//     compressedSaCompressionFactor(tag);
// }

// SEQAN_DEFINE_TEST(compressed_sa_assign_value)
// {
//     using namespace seqan;
// 
//     CompressedSA<SparseString<String<unsigned int>, void >, CharString, void> tag;
// 
//     compressedSaAssignValue(tag);
// }

SEQAN_DEFINE_TEST(compressed_sa_clear_length_resize)
{
    using namespace seqan;

    CompressedSA<SparseString<String<unsigned int>, void >, CharString, void> tag;

    compressedSaClearLengthResize(tag);
}

SEQAN_DEFINE_TEST(compressed_sa_empty)
{
    using namespace seqan;

    CompressedSA<SparseString<String<unsigned int>, void >, CharString, void> tag;

    compressedSaEmpty(tag);
}

SEQAN_DEFINE_TEST(compressed_sa_compressed_sa_create)
{
    using namespace seqan;

    CompressedSA<SparseString<String<unsigned int>, void >, CharString, void> tag;

    compressedSaCompressedSaCreate(tag);
}

SEQAN_DEFINE_TEST(compressed_sa_get_fibre)
{
    using namespace seqan;

    CompressedSA<SparseString<String<unsigned int>, void >, CharString, void> tag;

    compressedSaGetFibre(tag);
}

SEQAN_DEFINE_TEST(compressed_sa_get_next_pos_)
{
    using namespace seqan;

    typedef Dna TChar;
    typedef String<TChar> TText;

    Index<TText, FmIndex<WaveletTreeBased<FmiDollarSubstituted<> >, void > > tag;

    compressedSaGetNextPos_(tag);
}

SEQAN_DEFINE_TEST(compressed_sa_set_lf_table)
{
    using namespace seqan;

    CompressedSA<SparseString<String<unsigned int>, void >, CharString, void> tag;

    compressedSaSetLfTable(tag);
}

SEQAN_DEFINE_TEST(compressed_sa_value_access)
{
    using namespace seqan;

    typedef Dna TChar;
    typedef String<TChar> TText;

    Index<TText, FmIndex<WaveletTreeBased<FmiDollarSubstituted<> >, void > > tag;

    compressedSaValueAccess(tag);

}


#endif // TEST_COMPRESSED_SA_BETA_H_

