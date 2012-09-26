// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the HistogramSet data structure.
// ==========================================================================

#ifndef TEST_SYNOPSIS_TEST_SYNOPSIS_HISTOGRAM_SET_H_
#define TEST_SYNOPSIS_TEST_SYNOPSIS_HISTOGRAM_SET_H_

SEQAN_DEFINE_TEST(test_synopsis_histogram_set_construct)
{
    using namespace seqan;

    HistogramSet<unsigned> hSet(10, 2);
}

SEQAN_DEFINE_TEST(test_synopsis_histogram_set_clear)
{
    using namespace seqan;

    // Clear all histograms.
    {
        HistogramSet<unsigned> hSet(10, 2);
        setValue(hSet, 1, 1, 1);
        SEQAN_ASSERT_EQ(getValue(hSet, 1, 1), 1u);
        clear(hSet);
        SEQAN_ASSERT_EQ(getValue(hSet, 0, 0), 0u);
    }

    // Clear one histogram.
    {
        HistogramSet<unsigned> hSet(10, 2);
        setValue(hSet, 1, 1, 1);
        setValue(hSet, 2, 1, 1);
        clear(hSet, 2);
        SEQAN_ASSERT_EQ(getValue(hSet, 1, 1), 1u);
        SEQAN_ASSERT_EQ(getValue(hSet, 2, 1), 0u);
    }
}

SEQAN_DEFINE_TEST(test_synopsis_histogram_set_value)
{
    using namespace seqan;

    HistogramSet<unsigned> hSet(10, 2);
    SEQAN_ASSERT_EQ(getValue(hSet, 0, 0), 0u);
    value(hSet, 1, 1) = 1;
    SEQAN_ASSERT_EQ(getValue(hSet, 1, 1), 1u);
}

SEQAN_DEFINE_TEST(test_synopsis_histogram_set_set_value)
{
    using namespace seqan;

    HistogramSet<unsigned> hSet(10, 2);
    SEQAN_ASSERT_EQ(getValue(hSet, 0, 0), 0u);
    setValue(hSet, 1, 1, 1);
    SEQAN_ASSERT_EQ(getValue(hSet, 1, 1), 1u);
}

SEQAN_DEFINE_TEST(test_synopsis_histogram_set_increment_value)
{
    using namespace seqan;

    HistogramSet<unsigned> hSet(10, 2);
    SEQAN_ASSERT_EQ(getValue(hSet, 0, 0), 0u);
    incrementValue(hSet, 1, 1, 1);
    SEQAN_ASSERT_EQ(getValue(hSet, 1, 1), 1u);
    incrementValue(hSet, 1, 1, 3);
    SEQAN_ASSERT_EQ(getValue(hSet, 1, 1), 4u);
}

#endif  // TEST_SYNOPSIS_TEST_SYNOPSIS_HISTOGRAM_SET_H_
