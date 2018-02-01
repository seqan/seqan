// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_INDEX_TEST_INDEX_FM_SPARSE_STRING_H_
#define TESTS_INDEX_TEST_INDEX_FM_SPARSE_STRING_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

// ==========================================================================
// Tests
// ========================================================================== 

// template <typename TSparseString>
// void sparseStringCompressionFactor(TSparseString & /*tag*/)
// { 
//     TSparseString sparseString;
// 
//     assignCompressionFactor(sparseString, 3u);
// 
//     SEQAN_ASSERT_EQ(getCompressionFactor(sparseString), 3u);
// }

// template <typename TSparseString>
// void sparseStringAssignGetValue(TSparseString & /*tag*/)
// { 
//     TSparseString sparseString;
//     resize(sparseString, 30);
// 
//     assignValue(sparseString, 0, 0);
//     assignValue(sparseString, 10, 2);
//     assignValue(sparseString, 20, 3);
// 
//     SEQAN_ASSERT_EQ(getValue(sparseString, 0), 0u);
//     SEQAN_ASSERT_EQ(getValue(sparseString, 10), 2u);
//     SEQAN_ASSERT_EQ(getValue(sparseString, 20), 3u);
//     SEQAN_ASSERT_EQ(getValue(sparseString, 21), DefaultValue<TSparseString>::VALUE);
// }

template <typename TSparseString>
void sparseStringGetValue(TSparseString & /*tag*/)
{ 
    TSparseString sparseString;
    resize(getFibre(sparseString, FibreValues()), 3);

    assignValue(getFibre(sparseString, FibreValues()), 0, 0);
    assignValue(getFibre(sparseString, FibreValues()), 1, 2);
    assignValue(getFibre(sparseString, FibreValues()), 2, 3);

    String<bool> ind;
    resize(ind, 30, false);
    ind[0] = true;
    ind[10] = true;
    ind[20] = true;
    createRankDictionary(getFibre(sparseString, FibreIndicators()), ind);
    
    SEQAN_ASSERT_EQ(getValue(sparseString, 0), 0u);
    SEQAN_ASSERT_EQ(getValue(sparseString, 10), 2u);
    SEQAN_ASSERT_EQ(getValue(sparseString, 20), 3u);
    SEQAN_ASSERT_EQ(getValue(sparseString, 21), +DefaultValue<TSparseString>::VALUE);
}

template <typename TSparseString>
void sparseStringClearLengthResize(TSparseString & /*tag*/)
{ 
    TSparseString sparseString;

    SEQAN_ASSERT_EQ(length(sparseString), 0u);

    resize(sparseString, 30);
    SEQAN_ASSERT_EQ(length(sparseString), 30u);

    clear(sparseString);
    SEQAN_ASSERT_EQ(length(sparseString), 0u);
}

template <typename TSparseString>
void sparseStringEmpty(TSparseString & /*tag*/)
{ 
    TSparseString sparseString;

    SEQAN_ASSERT_EQ(empty(sparseString), true);

    resize(sparseString, 30);

    SEQAN_ASSERT_EQ(empty(sparseString), false);

    clear(sparseString);
    SEQAN_ASSERT_EQ(empty(sparseString), true);
}

// template <typename TSparseString>
// void sparseStringValue(TSparseString & /*tag*/)
// { 
//     typedef typename Fibre<TSparseString, FibreValues>::Type TValueString;
//     typedef typename Value<TValueString>::Type TValue;
// 
//     TSparseString sparseString;
//     resize(sparseString, 30);
// 
//     assignValue(sparseString, 0, 0);
//     assignValue(sparseString, 1, 2);
//     assignValue(sparseString, 2, 3);
// 
//     SEQAN_ASSERT_EQ(getValue(sparseString, 0), 0u);
//     SEQAN_ASSERT_EQ(getValue(sparseString, 1), 2u);
//     SEQAN_ASSERT_EQ(getValue(sparseString, 2), 3u);
// 
//     TValue & temp = value(sparseString, 0);
//     temp = 5;
//     TValue & temp2 = value(sparseString, 2);
//     temp2 = 10;
// 
//     SEQAN_ASSERT_EQ(getValue(sparseString, 0), 5u);
//     SEQAN_ASSERT_EQ(getValue(sparseString, 2), 10u);
// 
// }

// SEQAN_DEFINE_TEST(sparse_string_compression_factor)
// {
//     using namespace seqan;
// 
//     SparseString<String<unsigned int> > tag;
// 
//     sparseStringCompressionFactor(tag);
// }

SEQAN_DEFINE_TEST(sparse_string_get_value)
{
    using namespace seqan;

    SparseString<String<unsigned int> > tag;

    sparseStringGetValue(tag);
}

SEQAN_DEFINE_TEST(sparse_string_clear_length_resize)
{
    using namespace seqan;

    SparseString<String<unsigned int> > tag;

    sparseStringClearLengthResize(tag);
}

SEQAN_DEFINE_TEST(sparse_string_empty)
{
    using namespace seqan;

    SparseString<String<unsigned int>, void > tag;

    sparseStringEmpty(tag);
}

SEQAN_DEFINE_TEST(sparse_string_get_fibre)
{
    using namespace seqan;

    SparseString<String<unsigned int>, void > tag;

    sparseStringEmpty(tag);
}

// SEQAN_DEFINE_TEST(sparse_string_value)
// {
//     using namespace seqan;
// 
//     SparseString<String<unsigned int>, void > tag;
// 
//     sparseStringValue(tag);
// }

// ==========================================================================
// Functions
// ========================================================================== 

SEQAN_BEGIN_TESTSUITE(test_fm_index_sparse_string)
{
    SEQAN_CALL_TEST(sparse_string_get_value);
    SEQAN_CALL_TEST(sparse_string_clear_length_resize);
    SEQAN_CALL_TEST(sparse_string_empty);
    SEQAN_CALL_TEST(sparse_string_get_fibre);
}
SEQAN_END_TESTSUITE

#endif // TESTS_INDEX_TEST_INDEX_FM_SPARSE_STRING_H_
