// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/align/normalize.h>

SEQAN_DEFINE_TEST(test_align_remove_nullifying_indels_ins_del)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "TGAT");
    assignSource(row(align, 1), "TCAT");
    insertGap(row(align, 0), 1);
    insertGap(row(align, 1), 2);

    // std::cerr << align;

    bool b = removeNullifyingIndels(row(align, 0), row(align, 1));
    SEQAN_ASSERT(b);

    std::stringstream ss0, ss1;
    ss0 << row(align, 0);
    ss1 << row(align, 1);

    SEQAN_ASSERT_EQ("TGAT", ss0.str());
    SEQAN_ASSERT_EQ("TCAT", ss1.str());
}

SEQAN_DEFINE_TEST(test_align_remove_nullifying_indels_del_ins)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "TGAT");
    assignSource(row(align, 1), "TCAT");
    insertGap(row(align, 0), 2);
    insertGap(row(align, 1), 1);

    // std::cerr << align;

    bool b = removeNullifyingIndels(row(align, 0), row(align, 1));
    SEQAN_ASSERT(b);

    std::stringstream ss0, ss1;
    ss0 << row(align, 0);
    ss1 << row(align, 1);

    SEQAN_ASSERT_EQ("TGAT", ss0.str());
    SEQAN_ASSERT_EQ("TCAT", ss1.str());
}

SEQAN_DEFINE_TEST(test_align_remove_nullifying_indels_complex)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "TCATTCATTCATTCATTCAT");
    assignSource(row(align, 1), "TCATTCATTCATTCATTCAT");
    insertGaps(iter(row(align, 0), 4), 4);
    insertGaps(iter(row(align, 0), 11), 3);
    insertGaps(iter(row(align, 1), 8), 3);
    insertGaps(iter(row(align, 1), 14), 4);

    // std::cerr << align;

    bool b = removeNullifyingIndels(row(align, 0), row(align, 1));
    SEQAN_ASSERT(b);

    // std::cerr << align;

    std::stringstream ss0, ss1;
    ss0 << row(align, 0);
    ss1 << row(align, 1);

    SEQAN_ASSERT_EQ("TCAT-TCATTCATTCATTCAT", ss0.str());
    SEQAN_ASSERT_EQ("TCATTCATTCA-TTCATTCAT", ss1.str());
}

SEQAN_DEFINE_TEST(test_align_left_align_indels_nop)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "TGAT");
    assignSource(row(align, 1), "TCAT");
    insertGap(row(align, 0), 1);
    insertGap(row(align, 1), 3);

    // std::cerr << align;

    bool b = leftAlignIndels(row(align, 0), row(align, 1));
    SEQAN_ASSERT_NOT(b);

    std::stringstream ss0, ss1;
    ss0 << row(align, 0);
    ss1 << row(align, 1);

    SEQAN_ASSERT_EQ("T-GAT", ss0.str());
    SEQAN_ASSERT_EQ("TCA-T", ss1.str());
}

// Actual left-shifting of gaps.
SEQAN_DEFINE_TEST(test_align_left_align_indels_simple)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "TTAGAGAGTT");
    assignSource(row(align, 1), "TTAGAGTT");
    insertGaps(iter(row(align, 1), 6), 2);

    // std::cerr << align;

    bool b = leftAlignIndels(row(align, 0), row(align, 1));
    SEQAN_ASSERT(b);

    // std::cerr << align;

    std::stringstream ss0, ss1;
    ss0 << row(align, 0);
    ss1 << row(align, 1);

    SEQAN_ASSERT_EQ("TTAGAGAGTT", ss0.str());
    SEQAN_ASSERT_EQ("TT--AGAGTT", ss1.str());
}

// Actual left-shifting of multiple gaps.
SEQAN_DEFINE_TEST(test_align_left_align_indels_multiple)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "TTAGAGAGTTTCAGCAGG");
    assignSource(row(align, 1), "TTAGTTCACAG");
    insertGaps(iter(row(align, 1), 10), 1);
    insertGaps(iter(row(align, 1), 8), 1);
    insertGaps(iter(row(align, 1), 6), 1);
    insertGaps(iter(row(align, 1), 4), 2);
    insertGaps(iter(row(align, 1), 2), 2);

    // std::cerr << align;

    bool b = leftAlignIndels(row(align, 0), row(align, 1));
    SEQAN_ASSERT(b);

    // std::cerr << align;

    std::stringstream ss0, ss1;
    ss0 << row(align, 0);
    ss1 << row(align, 1);

    SEQAN_ASSERT_EQ("TTAGAGAGTTTCAGCAGG", ss0.str());
    SEQAN_ASSERT_EQ("TT----AG-TTCA-CA-G", ss1.str());
}

SEQAN_BEGIN_TESTSUITE(test_align_shift_indels_left)
{
    SEQAN_CALL_TEST(test_align_remove_nullifying_indels_ins_del);
    SEQAN_CALL_TEST(test_align_remove_nullifying_indels_del_ins);
    SEQAN_CALL_TEST(test_align_remove_nullifying_indels_complex);
    SEQAN_CALL_TEST(test_align_left_align_indels_nop);
    SEQAN_CALL_TEST(test_align_left_align_indels_simple);
    SEQAN_CALL_TEST(test_align_left_align_indels_multiple);
}
SEQAN_END_TESTSUITE
