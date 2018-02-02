// ==========================================================================
//                         Mason - A Read Simulator
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>

#include "sequencing.h"
#include "genomic_variants.h"

SEQAN_DEFINE_TEST(mason_tests_append_orientation_elementary_operations)
{
    // Below, we test for match, mismatch, insertion, deletion and insertion.
    {
        TCigarString cigar;

        std::pair<int, int> v = appendOperation(cigar, 'M');

        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'M');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
    {
        TCigarString cigar;

        std::pair<int, int> v = appendOperation(cigar, 'X');

        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'X');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
    {
        TCigarString cigar;

        std::pair<int, int> v = appendOperation(cigar, 'I');

        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 0);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'I');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
    {
        TCigarString cigar;

        std::pair<int, int> v = appendOperation(cigar, 'D');

        SEQAN_ASSERT_EQ(v.first, 0);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'D');
        SEQAN_ASSERT_EQ(cigar[0].count, 1u);
    }
}

SEQAN_DEFINE_TEST(mason_tests_append_orientation_combination)
{
    // Test with the combination of equal operations.
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('M', 1));
        std::pair<int, int> v = appendOperation(cigar, 'M');

        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'M');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('X', 1));
        std::pair<int, int> v = appendOperation(cigar, 'X');

        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'X');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('I', 1));
        std::pair<int, int> v = appendOperation(cigar, 'I');

        SEQAN_ASSERT_EQ(v.first, 1);
        SEQAN_ASSERT_EQ(v.second, 0);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'I');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('D', 1));
        std::pair<int, int> v = appendOperation(cigar, 'D');

        SEQAN_ASSERT_EQ(v.first, 0);
        SEQAN_ASSERT_EQ(v.second, 1);
        SEQAN_ASSERT_EQ(length(cigar), 1u);
        SEQAN_ASSERT_EQ(cigar[0].operation, 'D');
        SEQAN_ASSERT_EQ(cigar[0].count, 2u);
    }
}

SEQAN_DEFINE_TEST(mason_tests_append_orientation_canceling_out)
{
    // Test with the combination of operations that cancel each other out (I/D, D/I)
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('I', 1));
        std::pair<int, int> v = appendOperation(cigar, 'D');

        SEQAN_ASSERT_EQ(v.first, -1);
        SEQAN_ASSERT_EQ(v.second, 0);
        SEQAN_ASSERT_EQ(length(cigar), 0u);
    }
    {
        TCigarString cigar;

        appendValue(cigar, seqan::CigarElement<>('D', 1));
        std::pair<int, int> v = appendOperation(cigar, 'I');

        SEQAN_ASSERT_EQ(v.first, 0);
        SEQAN_ASSERT_EQ(v.second, -1);
        SEQAN_ASSERT_EQ(length(cigar), 0u);
    }
}

SEQAN_DEFINE_TEST(mason_tests_position_map_inversion)
{
    typedef PositionMap::TInterval TInterval;

    PositionMap positionMap;

    // Inversion: --1000-->|<--1000--|--1000-->
    GenomicInterval gi1(   0, 1000,    0, 1000, '+', GenomicInterval::NORMAL);
    GenomicInterval gi2(1000, 2000, 1000, 2000, '-', GenomicInterval::INVERTED);
    GenomicInterval gi3(2000, 3000, 2000, 3000, '+', GenomicInterval::NORMAL);

    // Build interval tree.
    seqan::String<TInterval> intervals;
    appendValue(intervals, TInterval(gi1.svBeginPos, gi1.svEndPos, gi1));
    appendValue(intervals, TInterval(gi2.svBeginPos, gi2.svEndPos, gi2));
    appendValue(intervals, TInterval(gi3.svBeginPos, gi3.svEndPos, gi3));
    createIntervalTree(positionMap.svIntervalTree, intervals);

    // Add breakpoints.
    positionMap.svBreakpoints.insert(std::make_pair(0, 0));
    positionMap.svBreakpoints.insert(std::make_pair(gi1.svEndPos, 1));
    positionMap.svBreakpoints.insert(std::make_pair(gi2.svEndPos, 2));
    positionMap.svBreakpoints.insert(std::make_pair(gi3.svEndPos, 3));

    // Tests for overlapsWithBreakpoint()
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(   0, 1000));
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(1000, 2000));
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(2000, 3000));

    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint( 999, 1001));
    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint(1999, 2001));
    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint(2999, 3001));

    // Tests for getGenomicInterval()
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval(   0));
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval( 999));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1000));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1999));
    SEQAN_ASSERT(gi3 == positionMap.getGenomicInterval(2000));
    SEQAN_ASSERT(gi3 == positionMap.getGenomicInterval(2999));

    // Tests for toSmallVarInterval().
    typedef std::pair<int, int> TPair;

    TPair i1 = positionMap.toSmallVarInterval(0, 100);
    SEQAN_ASSERT_EQ(i1.first, 0);
    SEQAN_ASSERT_EQ(i1.second, 100);

    TPair i2 = positionMap.toSmallVarInterval(900, 1000);
    SEQAN_ASSERT_EQ(i2.first, 900);
    SEQAN_ASSERT_EQ(i2.second, 1000);

    TPair i3 = positionMap.toSmallVarInterval(1000, 1100);
    SEQAN_ASSERT_EQ(i3.first, 2000);
    SEQAN_ASSERT_EQ(i3.second, 1900);

    TPair i4 = positionMap.toSmallVarInterval(1900, 2000);
    SEQAN_ASSERT_EQ(i4.first, 1100);
    SEQAN_ASSERT_EQ(i4.second, 1000);

    TPair i5 = positionMap.toSmallVarInterval(2000, 2100);
    SEQAN_ASSERT_EQ(i5.first, 2000);
    SEQAN_ASSERT_EQ(i5.second, 2100);

    TPair i6 = positionMap.toSmallVarInterval(2900, 3000);
    SEQAN_ASSERT_EQ(i6.first, 2900);
    SEQAN_ASSERT_EQ(i6.second, 3000);
}

SEQAN_DEFINE_TEST(mason_tests_position_map_translocation)
{
    typedef PositionMap::TInterval TInterval;

    PositionMap positionMap;

    // Translocation: --A--> --B-->
    //                --B--> --A-->
    GenomicInterval gi1(   0, 1000, 1000, 2000, '+', GenomicInterval::NORMAL);
    GenomicInterval gi2(1000, 2000,    0, 1000, '-', GenomicInterval::NORMAL);

    // Build interval tree.
    seqan::String<TInterval> intervals;
    appendValue(intervals, TInterval(gi1.svBeginPos, gi1.svEndPos, gi1));
    appendValue(intervals, TInterval(gi2.svBeginPos, gi2.svEndPos, gi2));
    createIntervalTree(positionMap.svIntervalTree, intervals);

    // Add breakpoints.
    positionMap.svBreakpoints.insert(std::make_pair(0, 0));
    positionMap.svBreakpoints.insert(std::make_pair(gi1.svEndPos, 1));
    positionMap.svBreakpoints.insert(std::make_pair(gi2.svEndPos, 2));

    // Tests for overlapsWithBreakpoint()
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(   0, 1000));
    SEQAN_ASSERT_NOT(positionMap.overlapsWithBreakpoint(1000, 2000));

    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint( 999, 1001));
    SEQAN_ASSERT(positionMap.overlapsWithBreakpoint(1999, 2001));

    // Tests for getGenomicInterval()
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval(   0));
    SEQAN_ASSERT(gi1 == positionMap.getGenomicInterval( 999));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1000));
    SEQAN_ASSERT(gi2 == positionMap.getGenomicInterval(1999));

    // Tests for toSmallVarInterval().
    typedef std::pair<int, int> TPair;

    TPair i1 = positionMap.toSmallVarInterval(0, 100);
    SEQAN_ASSERT_EQ(i1.first, 1000);
    SEQAN_ASSERT_EQ(i1.second, 1100);

    TPair i2 = positionMap.toSmallVarInterval(900, 1000);
    SEQAN_ASSERT_EQ(i2.first, 1900);
    SEQAN_ASSERT_EQ(i2.second, 2000);

    TPair i3 = positionMap.toSmallVarInterval(1000, 1100);
    SEQAN_ASSERT_EQ(i3.first, 0);
    SEQAN_ASSERT_EQ(i3.second, 100);

    TPair i4 = positionMap.toSmallVarInterval(1900, 2000);
    SEQAN_ASSERT_EQ(i4.first, 900);
    SEQAN_ASSERT_EQ(i4.second, 1000);
}

SEQAN_DEFINE_TEST(mason_tests_position_map_to_original_interval)
{
    // Deletions in the variant / insertions in the reference.
    {
        // Create the following situation in the journal.
        //
        //                    1        2
        //          0         0        0
        //          :    .    :    .   :
        //      REF XX--XXXX--XXXXX
        // SMALLVAR XXXXXXXXXXXXXXX
        TJournalEntries journal;
        reinit(journal, 100);
        recordInsertion(journal,  2, 0, 2);
        recordInsertion(journal,  8, 0, 2);

        PositionMap positionMap;
        positionMap.reinit(journal);

        // Check toOriginalInterval.
        std::pair<int, int> i1 = positionMap.toOriginalInterval(3, 5);
        SEQAN_ASSERT_EQ(i1.first, 2);
        SEQAN_ASSERT_EQ(i1.second, 3);
        std::pair<int, int> i2 = positionMap.toOriginalInterval(5, 9);
        SEQAN_ASSERT_EQ(i2.first, 3);
        SEQAN_ASSERT_EQ(i2.second, 6);
        std::pair<int, int> i3 = positionMap.toOriginalInterval(3, 9);
        SEQAN_ASSERT_EQ(i3.first, 2);
        SEQAN_ASSERT_EQ(i3.second, 6);
    }

    // Insertions in the variant / deletions in the reference.
    {
        // Create the following situation in the journal.
        //
        //                    1        2
        //          0         0        0
        //          :    .    :    .   :
        //      REF XXXXXXXXXXXXXXX
        // SMALLVAR XX--XXXX--XXXXX
        TJournalEntries journal;
        reinit(journal, 100);
        recordErase(journal,  8, 10);
        recordErase(journal,  2, 4);

        PositionMap positionMap;
        positionMap.reinit(journal);

        // Check toOriginalInterval.
        std::pair<int, int> i1 = positionMap.toOriginalInterval(2, 3);
        SEQAN_ASSERT_EQ(i1.first, 4);
        SEQAN_ASSERT_EQ(i1.second, 5);
        std::pair<int, int> i2 = positionMap.toOriginalInterval(5, 6);
        SEQAN_ASSERT_EQ(i2.first, 7);
        SEQAN_ASSERT_EQ(i2.second, 10);
        std::pair<int, int> i3 = positionMap.toOriginalInterval(2, 6);
        SEQAN_ASSERT_EQ(i3.first, 4);
        SEQAN_ASSERT_EQ(i3.second, 10);
    }
}

SEQAN_DEFINE_TEST(mason_tests_position_map_original_to_small_var)
{
    // Deletions in the variant / insertions in the reference.
    {
        // Create the following situation in the journal.
        //
        //                    1        2
        //          0         0        0
        //          :    .    :    .   :
        //      REF XX--XXXX--XXXXX
        // SMALLVAR XXXXXXXXXXXXXXX
        TJournalEntries journal;
        reinit(journal, 100);
        recordInsertion(journal,  2, 0, 2);
        recordInsertion(journal,  8, 0, 2);

        PositionMap positionMap;
        positionMap.reinit(journal);

        // Check originalToSmallVarInterval.
        std::pair<int, int> i1 = positionMap.originalToSmallVarInterval(2, 3);
        SEQAN_ASSERT_EQ(i1.first, 4);
        SEQAN_ASSERT_EQ(i1.second, 5);
        std::pair<int, int> i2 = positionMap.originalToSmallVarInterval(5, 6);
        SEQAN_ASSERT_EQ(i2.first, 7);
        SEQAN_ASSERT_EQ(i2.second, 10);
        std::pair<int, int> i3 = positionMap.originalToSmallVarInterval(2, 6);
        SEQAN_ASSERT_EQ(i3.first, 4);
        SEQAN_ASSERT_EQ(i3.second, 10);
    }

    // Insertions in the variant / deletions in the reference.
    {
        // Create the following situation in the journal.
        //
        //                    1        2
        //          0         0        0
        //          :    .    :    .   :
        //      REF XXXXXXXXXXXXXXX
        // SMALLVAR XX--XXXX--XXXXX
        TJournalEntries journal;
        reinit(journal, 100);
        recordErase(journal,  8, 10);
        recordErase(journal,  2, 4);

        PositionMap positionMap;
        positionMap.reinit(journal);

        // Check originalToSmallVarInterval.
        std::pair<int, int> i1 = positionMap.originalToSmallVarInterval(3, 5);
        SEQAN_ASSERT_EQ(i1.first, 2);
        SEQAN_ASSERT_EQ(i1.second, 3);
        std::pair<int, int> i2 = positionMap.originalToSmallVarInterval(5, 9);
        SEQAN_ASSERT_EQ(i2.first, 3);
        SEQAN_ASSERT_EQ(i2.second, 6);
        std::pair<int, int> i3 = positionMap.originalToSmallVarInterval(3, 9);
        SEQAN_ASSERT_EQ(i3.first, 2);
        SEQAN_ASSERT_EQ(i3.second, 6);
    }
}

SEQAN_BEGIN_TESTSUITE(mason_tests)
{
    SEQAN_CALL_TEST(mason_tests_append_orientation_elementary_operations);
    SEQAN_CALL_TEST(mason_tests_append_orientation_combination);
    SEQAN_CALL_TEST(mason_tests_append_orientation_canceling_out);

    SEQAN_CALL_TEST(mason_tests_position_map_inversion);
    SEQAN_CALL_TEST(mason_tests_position_map_translocation);
    SEQAN_CALL_TEST(mason_tests_position_map_to_original_interval);

    SEQAN_CALL_TEST(mason_tests_position_map_original_to_small_var);
}
SEQAN_END_TESTSUITE
