// ==========================================================================
//                                SeqAn-Flexbar
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// Author: Benjamin Strauch <b.strauch@fu-berlin.de>
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// This file provides tests for the read trimming of seqan-flexbar which is
// based in the implementation of the original flexbar program in [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include "read_trimming.h"

SEQAN_DEFINE_TEST(sliding_window_test)
{
	// No error checking, we assume the file exists (it's a constant test file).
    seqan::CharString buffer = SEQAN_PATH_TO_ROOT();
    append(buffer, "/apps/seqan_flexbar/tests/testsample.fq");
	seqan::SeqFileIn inStream(toCString(buffer), seqan::OPEN_RDONLY);
	seqan::String<char> id;
	seqan::String<seqan::Dna5Q> seq;

	// The number of bases Trimmomatic cuts from the sample file with the same settings.
	unsigned expected[] = {12, 21, 7, 75, 100, 100, 4, 6, 3, 3, 3, 100, 3, 2, 57, 2, 3, 31, 97, 52,
					   2, 4, 14, 2, 2, 3, 2, 3, 4, 2, 70, 3, 2, 3, 3, 3, 3, 54, 2, 8, 3, 3, 98, 3};

	unsigned i=0;
	Mean m = Mean(5);
	while (!atEnd(inStream))
	{
        readRecord(id, seq, inStream);
        unsigned len = length(seq);
        unsigned res = trimRead(seq, 20, m);
        SEQAN_ASSERT_EQ(expected[i++], res);	 // Check if returned number of trimmed bases is correct.
        SEQAN_ASSERT_EQ(res, len - length(seq)); // Check if we really trimmed so many bases.
	}
}

SEQAN_DEFINE_TEST(cut_tail_test)
{
	unsigned expected[] = {9, 9, 3, 26, 69, 3, 2, 0, 0, 0, 0, 23, 1, 0, 22, 0, 0, 12, 35, 6,
					   0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 52, 0};

    seqan::CharString buffer = SEQAN_PATH_TO_ROOT();
    append(buffer, "/apps/seqan_flexbar/tests/testsample.fq");
	seqan::SeqFileIn inStream(toCString(buffer), seqan::OPEN_RDONLY);
	seqan::String<char> id;
	seqan::String<seqan::Dna5Q> seq;

	unsigned i=0;
	Tail t;
	while (!atEnd(inStream))
	{
        readRecord(id, seq, inStream);
        unsigned len = length(seq);
        unsigned res = trimRead(seq, 20, t);
        SEQAN_ASSERT_EQ(expected[i++], res);     // Check if returned number of trimmed bases is correct.
        SEQAN_ASSERT_EQ(res, len - length(seq)); // Check if we really trimmed so many bases.
	}
}

SEQAN_DEFINE_TEST(cut_bwa_test)
{
	unsigned expected[] = {8, 8, 2, 25, 68, 2, 1, 0, 0, 0, 0, 22, 0, 0, 21, 0, 0, 11, 34, 5,
					   0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 51, 0};

    seqan::CharString buffer = SEQAN_PATH_TO_ROOT();
    append(buffer, "/apps/seqan_flexbar/tests/testsample.fq");
	seqan::SeqFileIn inStream(toCString(buffer), seqan::OPEN_RDONLY);
	seqan::String<char> id;
	seqan::String<seqan::Dna5Q> seq;

	unsigned i=0;
	BWA b;
	while (!atEnd(inStream))
	{
        readRecord(id, seq, inStream);
        unsigned len = length(seq);
        unsigned res = trimRead(seq, 20, b);
        SEQAN_ASSERT_EQ(expected[i++], res);     // Check if returned number of trimmed bases is correct.
        SEQAN_ASSERT_EQ(res, len - length(seq)); // Check if we really trimmed so many bases.
	}
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(sliding_window_test);
    SEQAN_CALL_TEST(cut_tail_test);
    SEQAN_CALL_TEST(cut_bwa_test);
}
SEQAN_END_TESTSUITE
