// ==========================================================================
//                     test_alignment_dp_trace_segment.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_TRACE_SEGMENT_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_TRACE_SEGMENT_H_

#include <seqan/basic.h>

#include <seqan/align.h>

template <typename TPosition, typename TSize>
void
testAlignmentTracebackTraceSegmentsConstructor()
{
    using namespace seqan;

    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    { // test default ctor
        TTraceSegment traceSegment;

        SEQAN_ASSERT_EQ(traceSegment._horizontalBeginPos, (TPosition) 0);
        SEQAN_ASSERT_EQ(traceSegment._verticalBeginPos, (TPosition) 0);
        SEQAN_ASSERT_EQ(traceSegment._length, (TSize) 0);
        SEQAN_ASSERT_EQ(traceSegment._traceValue, TraceBitMap_<>::NONE);
    }

    { // test copy ctor
        TTraceSegment traceSegment;
        traceSegment._horizontalBeginPos = 10;
        traceSegment._verticalBeginPos = 3;
        traceSegment._length = 5;
        traceSegment._traceValue = TraceBitMap_<>::DIAGONAL;

        TTraceSegment traceSegment2(traceSegment);

        SEQAN_ASSERT_EQ(traceSegment2._horizontalBeginPos, (TPosition) 10);
        SEQAN_ASSERT_EQ(traceSegment2._verticalBeginPos, (TPosition) 3);
        SEQAN_ASSERT_EQ(traceSegment2._length, (TSize) 5);
        SEQAN_ASSERT_EQ(traceSegment2._traceValue, TraceBitMap_<>::DIAGONAL);

        TTraceSegment traceSegment3 = traceSegment;

        SEQAN_ASSERT_EQ(traceSegment3._horizontalBeginPos, (TPosition) 10);
        SEQAN_ASSERT_EQ(traceSegment3._verticalBeginPos, (TPosition) 3);
        SEQAN_ASSERT_EQ(traceSegment3._length, (TSize) 5);
        SEQAN_ASSERT_EQ(traceSegment3._traceValue, TraceBitMap_<>::DIAGONAL);
    }

    { // test additional ctor
        TTraceSegment traceSegment(12, 13, 8, TraceBitMap_<>::VERTICAL);

        SEQAN_ASSERT_EQ(traceSegment._horizontalBeginPos, (TPosition) 12);
        SEQAN_ASSERT_EQ(traceSegment._verticalBeginPos, (TPosition) 13);
        SEQAN_ASSERT_EQ(traceSegment._length, (TSize) 8);
        SEQAN_ASSERT_EQ(traceSegment._traceValue, TraceBitMap_<>::VERTICAL);
    }
}

template <typename TPosition, typename TSize>
void
testAlignmentTracebackTraceSegmentsAssignment()
{
    using namespace seqan;

    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    { // test assignment
        TTraceSegment traceSegment;
        traceSegment._horizontalBeginPos = 10;
        traceSegment._verticalBeginPos = 3;
        traceSegment._length = 5;
        traceSegment._traceValue = TraceBitMap_<>::DIAGONAL;

        TTraceSegment traceSegment2;

        SEQAN_ASSERT_EQ(traceSegment2._horizontalBeginPos, (TPosition) 0);
        SEQAN_ASSERT_EQ(traceSegment2._verticalBeginPos, (TPosition) 0);
        SEQAN_ASSERT_EQ(traceSegment2._length, (TSize) 0);
        SEQAN_ASSERT_EQ(traceSegment2._traceValue, TraceBitMap_<>::NONE);

        traceSegment2 = traceSegment;

        SEQAN_ASSERT_EQ(traceSegment2._horizontalBeginPos, (TPosition) 10);
        SEQAN_ASSERT_EQ(traceSegment2._verticalBeginPos, (TPosition) 3);
        SEQAN_ASSERT_EQ(traceSegment2._length, (TSize) 5);
        SEQAN_ASSERT_EQ(traceSegment2._traceValue, TraceBitMap_<>::DIAGONAL);
    }
}

template <typename TPosition, typename TSize>
void
testAlignmentTracebackTraceSegmentsCompare()
{
    using namespace seqan;

    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    TTraceSegment traceSegment;
    traceSegment._horizontalBeginPos = 10;
    traceSegment._verticalBeginPos = 3;
    traceSegment._length = 5;
    traceSegment._traceValue = TraceBitMap_<>::DIAGONAL;

    TTraceSegment traceSegment2(traceSegment);

    SEQAN_ASSERT(traceSegment2 == traceSegment);
    traceSegment._traceValue = TraceBitMap_<>::HORIZONTAL;
    SEQAN_ASSERT(traceSegment2 !=  traceSegment);
}


template <typename TPosition, typename TSize>
void
testAlignmentTracebackTraceSegmentsPosition()
{
    using namespace seqan;

    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef typename Position<TTraceSegment>::Type TPosition_;
    bool result = +IsSameType<TPosition_, TPosition>::VALUE;
    SEQAN_ASSERT(result);
}

template <typename TPosition, typename TSize>
void
testAlignmentTracebackTraceSegmentsSize()
{
    using namespace seqan;

    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef typename Size<TTraceSegment>::Type TSize_;
    bool result = +IsSameType<TSize_, TSize>::VALUE;
    SEQAN_ASSERT(result);
}

template <typename TTarget>
void testAlignmentTracebackRecordTrace(TTarget & target)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;

    TTraceValue tv1 = TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL;
    TTraceValue tv2 = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL;
    TTraceValue tv3 = TraceBitMap_<>::HORIZONTAL;
    _recordSegment(target, 0, 0, 3, tv1);
    _recordSegment(target, 0, 3, 5, tv2);
    _recordSegment(target, 5, 8, 3, tv3);
    _recordSegment(target, 8, 8, 0, TraceBitMap_<>::DIAGONAL);



    SEQAN_ASSERT_EQ(target[0]._horizontalBeginPos, 0);
    SEQAN_ASSERT_EQ(target[0]._verticalBeginPos, 0);
    SEQAN_ASSERT_EQ(target[0]._length, 3);
    SEQAN_ASSERT_EQ(target[0]._traceValue, TraceBitMap_<>::DIAGONAL);
    SEQAN_ASSERT_EQ(target[1]._horizontalBeginPos, 0);
    SEQAN_ASSERT_EQ(target[1]._verticalBeginPos, 3);
    SEQAN_ASSERT_EQ(target[1]._length, 5);
    SEQAN_ASSERT_EQ(target[1]._traceValue, TraceBitMap_<>::VERTICAL);
    SEQAN_ASSERT_EQ(target[2]._horizontalBeginPos, 5);
    SEQAN_ASSERT_EQ(target[2]._verticalBeginPos, 8);
    SEQAN_ASSERT_EQ(target[2]._length, 3);
    SEQAN_ASSERT_EQ(target[2]._traceValue, TraceBitMap_<>::HORIZONTAL);

    SEQAN_ASSERT_EQ(length(target), 3u);
}

void testAlignmentTracebackTraceSegmentGetBeginHorizontal()
{
    using namespace seqan;

    typedef TraceSegment_<int, size_t> TTraceSegment;

    TTraceSegment traceSegm;

    traceSegm._horizontalBeginPos = 5;
    traceSegm._verticalBeginPos = 12;
    traceSegm._length = 7;

    SEQAN_ASSERT_EQ(_getBeginHorizontal(traceSegm), 5);

    traceSegm._horizontalBeginPos = 7;
    SEQAN_ASSERT_EQ(_getBeginHorizontal(traceSegm), 7);
}

void testAlignmentTracebackTraceSegmentGetBeginVertical()
{
    using namespace seqan;

    typedef TraceSegment_<int, size_t> TTraceSegment;

    TTraceSegment traceSegm;

    traceSegm._horizontalBeginPos = 5;
    traceSegm._verticalBeginPos = 12;
    traceSegm._length = 7;

    SEQAN_ASSERT_EQ(_getBeginVertical(traceSegm), 12);

    traceSegm._verticalBeginPos = 7;
    SEQAN_ASSERT_EQ(_getBeginVertical(traceSegm), 7);
}

void testAlignmentTracebackTraceSegmentGetEndHorizontal()
{
    using namespace seqan;

    typedef TraceSegment_<int, size_t> TTraceSegment;

    TTraceSegment traceSegm;

    traceSegm._horizontalBeginPos = 5;
    traceSegm._verticalBeginPos = 12;
    traceSegm._length = 7;
    traceSegm._traceValue = TraceBitMap_<>::HORIZONTAL;

    SEQAN_ASSERT_EQ(_getEndHorizontal(traceSegm), 12);

    traceSegm._traceValue = TraceBitMap_<>::VERTICAL;
    SEQAN_ASSERT_EQ(_getEndHorizontal(traceSegm), 5);
}

void testAlignmentTracebackTraceSegmentGetEndVertical()
{
    using namespace seqan;

    typedef TraceSegment_<int, size_t> TTraceSegment;

    TTraceSegment traceSegm;

    traceSegm._horizontalBeginPos = 5;
    traceSegm._verticalBeginPos = 12;
    traceSegm._length = 7;
    traceSegm._traceValue = TraceBitMap_<>::VERTICAL;

    SEQAN_ASSERT_EQ(_getEndVertical(traceSegm), 19);

    traceSegm._traceValue = TraceBitMap_<>::HORIZONTAL;
    SEQAN_ASSERT_EQ(_getEndVertical(traceSegm), 12);
}

void testAlignmentTracebackTraceSegmentTranslateTraceValue()
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(_translateTraceValue(TraceBitMap_<>::DIAGONAL), "D");
    SEQAN_ASSERT_EQ(_translateTraceValue(TraceBitMap_<>::VERTICAL), "V");
    SEQAN_ASSERT_EQ(_translateTraceValue(TraceBitMap_<>::HORIZONTAL), "H");
    SEQAN_ASSERT_EQ(_translateTraceValue(TraceBitMap_<>::VERTICAL_OPEN), "v");
    SEQAN_ASSERT_EQ(_translateTraceValue(TraceBitMap_<>::HORIZONTAL_OPEN), "h");
    SEQAN_ASSERT_EQ(_translateTraceValue(TraceBitMap_<>::NONE), "0");
}

void testAlignmentTracebackTraceSegmentStreamOperator()
{
    using namespace seqan;

    typedef TraceSegment_<int, size_t> TTraceSegment;

    TTraceSegment traceSegm;

    traceSegm._horizontalBeginPos = 5;
    traceSegm._verticalBeginPos = 12;
    traceSegm._length = 7;
    traceSegm._traceValue = TraceBitMap_<>::DIAGONAL;

    std::stringstream ss;
    ss << traceSegm;

    String<char> testString = "D-(5, 12, 7)";
    String<char> res = ss.str();

    SEQAN_ASSERT_EQ(res, testString);
}


SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_constructor)
{
    testAlignmentTracebackTraceSegmentsConstructor<int, unsigned int>();
    testAlignmentTracebackTraceSegmentsConstructor<long, int>();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_assignment)
{
    testAlignmentTracebackTraceSegmentsAssignment<long, long>();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_position)
{
    testAlignmentTracebackTraceSegmentsPosition<unsigned long, unsigned int>();
    testAlignmentTracebackTraceSegmentsPosition<long, int>();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_size)
{
    testAlignmentTracebackTraceSegmentsSize<int, unsigned int>();
    testAlignmentTracebackTraceSegmentsSize<long long, int>();
}
SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_get_begin_horizontal)
{
    testAlignmentTracebackTraceSegmentGetBeginHorizontal();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_get_begin_vertical)
{
    testAlignmentTracebackTraceSegmentGetBeginVertical();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_get_end_horizontal)
{
    testAlignmentTracebackTraceSegmentGetEndHorizontal();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_get_end_vertical)
{
    testAlignmentTracebackTraceSegmentGetEndVertical();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_translate_trace_value)
{
    testAlignmentTracebackTraceSegmentTranslateTraceValue();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_operator_stream)
{
    testAlignmentTracebackTraceSegmentStreamOperator();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_operator_equal)
{
    testAlignmentTracebackTraceSegmentsCompare<long, long>();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_operator_unequal)
{
    testAlignmentTracebackTraceSegmentsCompare<long, long>();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_tracesegment_record_segment)
{
    seqan::StringSet<seqan::TraceSegment_<int, int> > traceSegments;
    testAlignmentTracebackRecordTrace(traceSegments);
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_TRACE_SEGMENT_H_
