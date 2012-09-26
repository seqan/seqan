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
// Tests for the synopsis module: Mining of large data sets.
// ==========================================================================

#include "config.h"

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/synopsis.h>

#include "test_synopsis_counter_buckets.h"
#include "test_synopsis_hot_list.h"
#include "test_synopsis_histogram_set.h"

SEQAN_BEGIN_TESTSUITE(test_template)
{
    // ----------------------------------------------------------------------
    // Tests for support data structures.
    // ----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_synopsis_counter_buckets_simple);

    // ----------------------------------------------------------------------
    // Tests for the frequency mining data structures.
    // ----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_synopsis_frequency_hot_list_frequent);
    SEQAN_CALL_TEST(test_synopsis_frequency_hot_list_frequent_large);

    SEQAN_CALL_TEST(test_synopsis_frequency_hot_list_lossy_counting_large);

    SEQAN_CALL_TEST(test_synopsis_frequency_hot_list_space_saving_large);

    // ----------------------------------------------------------------------
    // Tests for HistogramSet.
    // ----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_synopsis_histogram_set_construct);
    SEQAN_CALL_TEST(test_synopsis_histogram_set_clear);
    SEQAN_CALL_TEST(test_synopsis_histogram_set_value);
    SEQAN_CALL_TEST(test_synopsis_histogram_set_set_value);
    SEQAN_CALL_TEST(test_synopsis_histogram_set_increment_value);
}
SEQAN_END_TESTSUITE
