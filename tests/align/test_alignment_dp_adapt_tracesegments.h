// ==========================================================================
//                  test_alignment_dp_adapt_tracesegments.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_ADAPT_TRACESEGMENTS_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_ADAPT_TRACESEGMENTS_H_

#include <seqan/align.h>

template <typename TPosition, typename TSize>
void
testAlign2TracebackTraceSegmentsConstructor()
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
testAlign2TracebackTraceSegmentsAssignment()
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
testAlign2TracebackTraceSegmentsCompare()
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
testAlign2TracebackTraceSegmentsPosition()
{
    using namespace seqan;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef typename Position<TTraceSegment>::Type TPosition_;
    bool result = +IsSameType<TPosition_, TPosition>::VALUE;
    SEQAN_ASSERT(result);
}

template <typename TPosition, typename TSize>
void
testAlign2TracebackTraceSegmentsSize()
{
    using namespace seqan;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef typename Size<TTraceSegment>::Type TSize_;
    bool result = +IsSameType<TSize_, TSize>::VALUE;
    SEQAN_ASSERT(result);
}

template <typename TTarget>
void testAlign2TracebackRecordTrace(TTarget & target)
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

//    _recordSegment(target, 8, 8, 10, TraceBitMap_<>::NONE); // note this should fail when uncommented
}

void testAlign2TraceAdaptorAdaptFile()
{
    using namespace seqan;
    std::stringstream stream;
    typedef TraceSegment_<size_t, size_t> TTraceSegment;
    String<TTraceSegment> traceSegments;

    appendValue(traceSegments, TTraceSegment(12, 8, 4, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(8, 4, 4, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(8, 3, 1, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(4, 3, 4, TraceBitMap_<>::HORIZONTAL));
    appendValue(traceSegments, TTraceSegment(1, 0, 3, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(0, 0, 1, TraceBitMap_<>::HORIZONTAL));

    String<char> seq0 = "AAAACCCCGGGG";
    String<char> seq1 = "AAAACCCCGGGG";

    _adaptTraceSegmentsTo(stream, seq0, seq1, traceSegments);
    String<char> result = stream.str();

    SEQAN_ASSERT_EQ(result[1], 'A');
    SEQAN_ASSERT_EQ(result[3], gapValue<char>());

    SEQAN_ASSERT_EQ(result[7], 'A');
    SEQAN_ASSERT_EQ(result[9], 'A');

    SEQAN_ASSERT_EQ(result[13], 'A');
    SEQAN_ASSERT_EQ(result[15], 'A');

    SEQAN_ASSERT_EQ(result[19], 'A');
    SEQAN_ASSERT_EQ(result[21], 'A');

    SEQAN_ASSERT_EQ(result[25], 'C');
    SEQAN_ASSERT_EQ(result[27], gapValue<char>());

    SEQAN_ASSERT_EQ(result[31], 'C');
    SEQAN_ASSERT_EQ(result[33], gapValue<char>());

    SEQAN_ASSERT_EQ(result[37], 'C');
    SEQAN_ASSERT_EQ(result[39], gapValue<char>());

    SEQAN_ASSERT_EQ(result[43], 'C');
    SEQAN_ASSERT_EQ(result[45], gapValue<char>());

    SEQAN_ASSERT_EQ(result[49], gapValue<char>());
    SEQAN_ASSERT_EQ(result[51], 'A');

    SEQAN_ASSERT_EQ(result[55], 'G');
    SEQAN_ASSERT_EQ(result[57], 'C');

    SEQAN_ASSERT_EQ(result[61], 'G');
    SEQAN_ASSERT_EQ(result[63], 'C');

    SEQAN_ASSERT_EQ(result[67], 'G');
    SEQAN_ASSERT_EQ(result[69], 'C');

    SEQAN_ASSERT_EQ(result[73], 'G');
    SEQAN_ASSERT_EQ(result[75], 'C');

    SEQAN_ASSERT_EQ(result[79], gapValue<char>());
    SEQAN_ASSERT_EQ(result[81], 'G');

    SEQAN_ASSERT_EQ(result[85], gapValue<char>());
    SEQAN_ASSERT_EQ(result[87], 'G');

    SEQAN_ASSERT_EQ(result[91], gapValue<char>());
    SEQAN_ASSERT_EQ(result[93], 'G');

    SEQAN_ASSERT_EQ(result[97], gapValue<char>());
    SEQAN_ASSERT_EQ(result[99], 'G');
}


void testAlign2TraceAdaptorAdaptAlign()
{
    using namespace seqan;
    typedef TraceSegment_<size_t, size_t> TTraceSegment;
    String<TTraceSegment> traceSegments;

    appendValue(traceSegments, TTraceSegment(12, 8, 4, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(8, 4, 4, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(8, 3, 1, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(4, 3, 4, TraceBitMap_<>::HORIZONTAL));
    appendValue(traceSegments, TTraceSegment(1, 0, 3, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(0, 0, 1, TraceBitMap_<>::HORIZONTAL));

    String<char> seq0 = "AAAACCCCGGGG";
    String<char> seq1 = "AAAACCCCGGGG";

    Align<String<char> > align;

    resize(rows(align), 2);
    assignSource(row(align, 0), seq0);
    assignSource(row(align, 1), seq1);

    // prepare align object to test the correctness of the adapt function.
    Align<String<char> > compareAlign;
    resize(rows(compareAlign), 2);
    assignSource(row(compareAlign, 0), seq0);
    assignSource(row(compareAlign, 1), seq1);

    insertGap(row(compareAlign, 0), 8);
    insertGaps(row(compareAlign, 0), 13, 4);

    insertGap(row(compareAlign, 1), 0);
    insertGaps(row(compareAlign, 1), 4, 4);

    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traceSegments);

    SEQAN_ASSERT_EQ(row(align, 0), row(compareAlign, 0));
    SEQAN_ASSERT_EQ(row(align, 1), row(compareAlign, 1));
}

void testAlign2TraceAdaptorAdaptFragments()
{
    using namespace seqan;
    typedef TraceSegment_<size_t, size_t> TTraceSegment;
    typedef Fragment<unsigned> TFragment;
    String<TFragment> fragmentString;


    String<TTraceSegment> traceSegments;

    appendValue(traceSegments, TTraceSegment(12, 8, 4, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(8, 4, 4, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(8, 3, 1, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(4, 3, 4, TraceBitMap_<>::HORIZONTAL));
    appendValue(traceSegments, TTraceSegment(1, 0, 3, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(0, 0, 1, TraceBitMap_<>::HORIZONTAL));

    _adaptTraceSegmentsTo(fragmentString, 0, 1, traceSegments);

    SEQAN_ASSERT_EQ(fragmentString[0] == TFragment(0, 8, 1, 4, 4), true);
    SEQAN_ASSERT_EQ(fragmentString[1] == TFragment(0, 1, 1, 0, 3), true);
}

void testAlign2TraceAdaptorAdaptAlignmentGraph()
{
    using namespace seqan;

    typedef TraceSegment_<size_t, size_t> TTraceSegment;
    typedef Graph<Alignment<StringSet<CharString, Dependent<> > > > TAlignGraph;

    StringSet<CharString> strings;
    appendValue(strings, "AAAACCCCGGGG");
    appendValue(strings, "AAAACCCCGGGG");

    TAlignGraph alignGraph(strings);

    String<TTraceSegment> traceSegments;

    appendValue(traceSegments, TTraceSegment(12, 8, 4, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(8, 4, 4, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(8, 3, 1, TraceBitMap_<>::VERTICAL));
    appendValue(traceSegments, TTraceSegment(4, 3, 4, TraceBitMap_<>::HORIZONTAL));
    appendValue(traceSegments, TTraceSegment(1, 0, 3, TraceBitMap_<>::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(0, 0, 1, TraceBitMap_<>::HORIZONTAL));

    _adaptTraceSegmentsTo(alignGraph, 0, 1, traceSegments);

    std::stringstream ss;
    ss << alignGraph;

    std::stringstream expectedSS;
    expectedSS << "Alignment matrix:\n"
               << "      0     .    :    .   \n"
               << "        AAAACCCC-GGGG----\n"
               << "         |||             \n"
               << "        -AAA----ACCCCGGGG\n\n\n";

    SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());

}

SEQAN_DEFINE_TEST(test_align2_trace_adaptor_trace_segment)
{
    testAlign2TracebackTraceSegmentsSize<int, unsigned int>();
    testAlign2TracebackTraceSegmentsSize<long long, int>();
    testAlign2TracebackTraceSegmentsPosition<unsigned long, unsigned int>();
    testAlign2TracebackTraceSegmentsPosition<long, int>();
    testAlign2TracebackTraceSegmentsConstructor<int, unsigned int>();
    testAlign2TracebackTraceSegmentsConstructor<long, int>();
    testAlign2TracebackTraceSegmentsCompare<long, long>();
    testAlign2TracebackTraceSegmentsAssignment<long, long>();
}

SEQAN_DEFINE_TEST(test_align2_trace_adaptor_record_trace_segment)
{
    seqan::StringSet<seqan::TraceSegment_<int, int> > traceSegments;
    testAlign2TracebackRecordTrace(traceSegments);
}


SEQAN_DEFINE_TEST(test_align2_trace_adaptor_adapt_file)
{
    testAlign2TraceAdaptorAdaptFile();
}
SEQAN_DEFINE_TEST(test_align2_trace_adaptor_adapt_align)
{
    testAlign2TraceAdaptorAdaptAlign();
}
SEQAN_DEFINE_TEST(test_align2_trace_adaptor_adapt_fragments)
{
    testAlign2TraceAdaptorAdaptFragments();
}

SEQAN_DEFINE_TEST(test_align2_trace_adaptor_adapt_alignment_graph)
{
    testAlign2TraceAdaptorAdaptAlignmentGraph();
}
//SEQAN_DEFINE_TEST(test_align2_trace_adaptor_adapt_vertex_descriptor)
//{
//    SEQAN_ASSERT_FAIL("TODO: implement test");
//}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_ADAPT_TRACESEGMENTS_H_
