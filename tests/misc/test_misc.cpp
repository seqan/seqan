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

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <ctime>
#define SEQAN_DEBUG


#include <seqan/basic.h>

#include <seqan/map.h>

#include <seqan/misc/edit_environment.h>
#include <seqan/misc/base.h>
#include <seqan/misc/dequeue.h>
#include <seqan/misc/map.h>
#include <seqan/misc/set.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>
#include <seqan/misc/terminal.h>

#include "test_misc_interval_tree.h"
#include "test_misc_accumulators.h"
#include "test_misc_edit_environment.h"
#include "test_misc_bit_twiddling.h"

using namespace std;
using namespace seqan;

SEQAN_DEFINE_TEST(test_misc_terminal_get_terminal_size)
{
    using namespace seqan;
    
    unsigned cols = std::numeric_limits<unsigned>::max(), rows = std::numeric_limits<unsigned>::max();
    bool succ = getTerminalSize(cols, rows);

#if !defined(STDLIB_VS)
    SEQAN_ASSERT(succ);
    SEQAN_ASSERT_NEQ(cols, std::numeric_limits<unsigned>::max());
    SEQAN_ASSERT_NEQ(rows, std::numeric_limits<unsigned>::max());
#else  // #if !defined(STDLIB_VS)
    SEQAN_ASSERT_NOT(succ);
#endif  // #if !defined(STDLIB_VS)
}

SEQAN_BEGIN_TESTSUITE(test_misc) {
    SEQAN_CALL_TEST(test_misc_terminal_get_terminal_size);

    // Test Bit Twiddling.
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_char);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_signed_char);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_unsigned_char);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_short);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_signed_short);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_unsigned_short);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_int);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_signed_int);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_unsigned_int);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_long);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_signed_long);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_unsigned_long);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_long_long);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_signed_long_long);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_unsigned_long_long);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_int8);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_uint8);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_int16);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_uint16);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_int32);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_uint32);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_int64);
    SEQAN_CALL_TEST(test_misc_bit_twiddling_pop_count_uint64);

    // Test IntervalTree class
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_QueryAtBoundary);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_IntervalTree__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_IntervalTreeFromIterator__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_NonFullLength__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_AddInterval__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_TreeStructure__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_FindIntervalExcludeTouching__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_FindNoInterval__int);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_GraphMap__int_ComputeCenter_StoreIntervals);
    SEQAN_CALL_TEST(Interval_Tree__IntervalTreeTest_FindIntervalsIntervals__int_ComputeCenter);

    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_int_average);
    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_int_count);
    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_int_sum);
    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_int_clear);

    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_double_average);
    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_double_count);
    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_double_sum);
    SEQAN_CALL_TEST(test_misc_accumulators_average_accumulator_double_clear);

    SEQAN_CALL_TEST(test_misc_edit_environment_string_enumerator_hamming);
    SEQAN_CALL_TEST(test_misc_edit_environment_string_enumerator_iterator_hamming);
    SEQAN_CALL_TEST(test_misc_edit_environment_string_enumerator_edit);
    SEQAN_CALL_TEST(test_misc_edit_environment_string_enumerator_iterator_edit);
}
SEQAN_END_TESTSUITE

