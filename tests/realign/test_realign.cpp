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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the realign module.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/realign.h>

// TODO(holtgrew): Add test with alignments that do not flush left.

static const bool PRINT_REALIGNMENTS = false;
static const bool DEBUG_REALIGNMENT = false;

// Helper function that allows one to add gaps into the i-th read alignment in a fragmentStore.
template <typename TFragmentStore>
void addGaps(TFragmentStore & store, unsigned alignID, unsigned pos)
{
    typedef typename TFragmentStore::TAlignedReadStore     TAlignedReadStore;
    typedef typename seqan::Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TGapAnchors TGapAnchors;
    typedef typename TFragmentStore::TReadSeqStore         TReadSeqStore;
    typedef typename seqan::Value<TReadSeqStore>::Type     TReadSeq;
    typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<TGapAnchors> > TGaps;

    TGaps gaps(store.readSeqStore[store.alignedReadStore[alignID].readId], store.alignedReadStore[alignID].gaps);
    insertGap(gaps, pos);
}

// One read only without any gaps.
SEQAN_DEFINE_TEST(test_realign_one_read_no_gaps)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT", seqan::CharString("r0"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0, (TSize)0, 0 + length(store.readSeqStore[0]));

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 0, 0, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 120, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT-------------------\n"
            ".....................................................................................................\n";

    // std::cerr << ">>" << EXPECTED << "<<--\n>>" << ss.str() << "<<\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// One read only with some gaps (beginning, middle, end).
SEQAN_DEFINE_TEST(test_realign_one_read_with_gaps)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT", seqan::CharString("r0"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0, (TSize)0, 0 + length(store.readSeqStore[0]) + 3);

    // Add gaps at beginning, end, and in the middle.
    addGaps(store, 0, 1);
    addGaps(store, 0, 51);
    addGaps(store, 0, 102);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 0, 0, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 120, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT-------------------\n"
            ".....................................................................................................\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// Two reads that are stacked.
SEQAN_DEFINE_TEST(test_realign_two_reads_stacked_at_beginning)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT", seqan::CharString("r0"));
    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT", seqan::CharString("r1"));
    appendRead(store,           "ATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTCTCCTCTCTC", seqan::CharString("r2"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0, (TSize)0, 0 + length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0, (TSize)0, 0 + length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0, (TSize)10, 10 + length(store.readSeqStore[2]));

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 200, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 0, 0, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 200, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 200, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTCTCCTCTCTC-----------------------------------------------------------------------------------------\n"
            ".....................................................................................................\n"
            ".....................................................................................................\n"
            "          .....................................................................................................\n";

    // std::cerr << ">>" << EXPECTED << "<<--\n>>" << ss.str() << "<<\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// Simple test case where nothing should change.
SEQAN_DEFINE_TEST(test_realign_simple_case)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT", seqan::CharString("r0"));
    appendRead(store, "TATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAA", seqan::CharString("r1"));
    appendRead(store, "AGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACAT", seqan::CharString("r2"));
    appendRead(store, "AAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTA", seqan::CharString("r3"));
    appendRead(store, "ATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTT", seqan::CharString("r4"));
    appendRead(store, "TGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTTTTG", seqan::CharString("r5"));
    appendRead(store, "TACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTTTTGAGATGGAGTCTTGCTCTGTC", seqan::CharString("r6"));
    appendRead(store, "GTTAACTTATGTACATGTTTATCTTTTTTGAGATGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAATGGCGCAATCTCGACTCTCTGCAACCTCCACCT", seqan::CharString("r7"));
    appendRead(store, "TGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAATGGCGCAATCTCGACTCTCTGCAACCTCCACCTGTCCGGTTCAAGTGATTCTCCTGTCTCAGCCTC", seqan::CharString("r8"));
    appendRead(store, "GGCTCAATCTCGCCTCTCTGCAGCCTCCGCCTGTCCGGTTCAAGTGATTCTCCTGTCTCAGCCTCCCGAGTAGCTGGGATTACAGGTGCCCGCCACCATGC", seqan::CharString("r9"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0,   (TSize)0,   0 + length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0,  (TSize)29,  29 + length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0,  (TSize)52,  52 + length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0,  (TSize)57,  57 + length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0,  (TSize)63,  63 + length(store.readSeqStore[4]));
    appendAlignedRead(store, 5, 0,  (TSize)66,  66 + length(store.readSeqStore[5]));
    appendAlignedRead(store, 6, 0,  (TSize)86,  86 + length(store.readSeqStore[6]));
    appendAlignedRead(store, 7, 0, (TSize)137, 137 + length(store.readSeqStore[7]));
    appendAlignedRead(store, 8, 0, (TSize)170, 170 + length(store.readSeqStore[8]));
    appendAlignedRead(store, 9, 0, (TSize)206, 206 + length(store.readSeqStore[9]));

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 0, 0, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 400, 0, 20);

    const char * EXPECTED =
"AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTTTTGAGATGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAATGGCGCAATCTCGACTCTCTGCAACCTCCACCTGTCCGGTTCAAGTGATTCTCCTGTCTCAGCCTCCCGAGTAGCTGGGATTACAGGTGCCCGCCACCATGC---------------------------------------------------------------------------------------------\n"
".....................................................................................................                                    .....................................................................................................\n"
"                             .....................................................................................................                                        .....................................................................................................\n"
"                                                    .....................................................................................................                                                     ...T........C.........G.....G........................................................................\n"
"                                                         .....................................................................................................\n"
"                                                               .....................................................................................................\n"
"                                                                  .....................................................................................................\n"
"                                                                                      .....................................................................................................\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// A test case with some gaps strewn in.
SEQAN_DEFINE_TEST(test_realign_simple_gaps)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATT", seqan::CharString("r0"));
    appendRead(store, "TATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAA", seqan::CharString("r1"));
    appendRead(store, "AGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACAT", seqan::CharString("r2"));
    appendRead(store, "AAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTA", seqan::CharString("r3"));
    appendRead(store, "ATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTT", seqan::CharString("r4"));
    appendRead(store, "TGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTTTTG", seqan::CharString("r5"));
    appendRead(store, "TACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTTTTGAGATGGAGTCTTGCTCTGTC", seqan::CharString("r6"));
    appendRead(store, "GTTAACTTATGTACATGTTTATCTTTTTTGAGATGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAATGGCGCAATCTCGACTCTCTGCAACCTCCACCT", seqan::CharString("r7"));
    appendRead(store, "TGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAATGGCGCAATCTCGACTCTCTGCAACCTCCACCTGTCCGGTTCAAGTGATTCTCCTGTCTCAGCCTC", seqan::CharString("r8"));
    appendRead(store, "GGCTCAATCTCGCCTCTCTGCAGCCTCCGCCTGTCCGGTTCAAGTGATTCTCCTGTCTCAGCCTCCCGAGTAGCTGGGATTACAGGTGCCCGCCACCATGC", seqan::CharString("r9"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0,   (TSize)0,   0 + length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0,  (TSize)29,  29 + length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0,  (TSize)52,  52 + length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0,  (TSize)57,  57 + length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0,  (TSize)63,  63 + length(store.readSeqStore[4]));
    appendAlignedRead(store, 5, 0,  (TSize)66,  66 + length(store.readSeqStore[5]));
    appendAlignedRead(store, 6, 0,  (TSize)86,  86 + length(store.readSeqStore[6]));
    appendAlignedRead(store, 7, 0, (TSize)137, 137 + length(store.readSeqStore[7]));
    appendAlignedRead(store, 8, 0, (TSize)170, 170 + length(store.readSeqStore[8]));
    appendAlignedRead(store, 9, 0, (TSize)206, 206 + length(store.readSeqStore[9]));

    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
    {
        addGaps(store, i, 3 * i + 1);
        store.alignedReadStore[i].endPos += 1;
    }

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 0, 0, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 400, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTTAACTTATGTACATGTTTATCTTTTTTGAGATGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAATGGCGCAATCTCGACTCTCTGCAACCTCCACCTGTCCGGTTCAAGTGATTCTCCTGTCTCAGCCTCCCGAGTAGCTGGGATTACAGGTGCCCGCCACCATGC---------------------------------------------------------------------------------------------\n"
            ".....................................................................................................                                    .....................................................................................................\n"
            "                             .....................................................................................................                                        .....................................................................................................\n"
            "                                                    .....................................................................................................                                                     ...T........C.........G.....G........................................................................\n"
            "                                                         .....................................................................................................\n"
            "                                                               .....................................................................................................\n"
            "                                                                  .....................................................................................................\n"
            "                                                                                      .....................................................................................................\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// A test case with an insertion
SEQAN_DEFINE_TEST(test_realign_simple_insert)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCAT", seqan::CharString("r0"));
    appendRead(store, "ATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAG", seqan::CharString("r1"));
    appendRead(store, "GCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCA", seqan::CharString("r2"));
    appendRead(store, "ATACTCTGGAATACTATACAGTAGTTAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAA", seqan::CharString("r3"));
    appendRead(store, "ATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT", seqan::CharString("r4"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0,   (TSize)0,   0 + length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0,  (TSize)10,  10 + length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0,  (TSize)20,  20 + length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0,  (TSize)30,  30 + length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0,  (TSize)40,  40 + length(store.readSeqStore[4]));

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 0, 0, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 200, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAA----AATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT------------------------------------------------------------\n"
    "..........................................................****......................................\n"
    "          ................................................****................................................\n"
    "                    ......................................****..........................................................\n"
    "                              ............................*****...................................................................\n"
    "                                        ..................AAAG..............................................................................\n";

    if (DEBUG_REALIGNMENT)
        std::cerr << "" << ss.str() << "---\n" << EXPECTED << "\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// A test case with an insertion and a window enclosing the part comfortably.
SEQAN_DEFINE_TEST(test_realign_simple_insert_window)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCAT", seqan::CharString("r0"));
    appendRead(store, "ATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAG", seqan::CharString("r1"));
    appendRead(store, "GCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCA", seqan::CharString("r2"));
    appendRead(store, "ATACTCTGGAATACTATACAGTAGTTAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAA", seqan::CharString("r3"));
    appendRead(store, "ATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT", seqan::CharString("r4"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0,   (TSize)0,   0 + length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0,  (TSize)10,  10 + length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0,  (TSize)20,  20 + length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0,  (TSize)30,  30 + length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0,  (TSize)40,  40 + length(store.readSeqStore[4]));

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // The window comfortably starts before the insert and ends behind it.
    reAlignment(store, /*contigID=*/0, 1, 10, false, 50, 139, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 200, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAA----AATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT------------------------------------------------------------\n"
            "..........................................................****......................................\n"
            "          ................................................****................................................\n"
            "                    ......................................****..........................................................\n"
            "                              ............................*****...................................................................\n"
            "                                        ..................AAAG..............................................................................\n";

    if (DEBUG_REALIGNMENT)
        std::cerr << "" << ss.str() << "---\n" << EXPECTED << "\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// A test case with an insertion and a window being tight on the left.
SEQAN_DEFINE_TEST(test_realign_simple_insert_window_tight_left)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCAT", seqan::CharString("r0"));
    appendRead(store, "ATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAG", seqan::CharString("r1"));
    appendRead(store, "GCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCA", seqan::CharString("r2"));
    appendRead(store, "ATACTCTGGAATACTATACAGTAGTTAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAA", seqan::CharString("r3"));
    appendRead(store, "ATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT", seqan::CharString("r4"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0,   (TSize)0,   0 + length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0,  (TSize)10,  10 + length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0,  (TSize)20,  20 + length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0,  (TSize)30,  30 + length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0,  (TSize)40,  40 + length(store.readSeqStore[4]));

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 58, 139, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 200, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAA----AATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT------------------------------------------------------------\n"
            "..........................................................****......................................\n"
            "          ................................................****................................................\n"
            "                    ......................................****..........................................................\n"
            "                              ............................*****...................................................................\n"
            "                                        ..................AAAG..............................................................................\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

// A test case with an insertion and a window being tight on the right.
SEQAN_DEFINE_TEST(test_realign_simple_insert_window_tight_right)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCAT", seqan::CharString("r0"));
    appendRead(store, "ATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAG", seqan::CharString("r1"));
    appendRead(store, "GCCACAGTGTATACTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCA", seqan::CharString("r2"));
    appendRead(store, "ATACTCTGGAATACTATACAGTAGTTAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAA", seqan::CharString("r3"));
    appendRead(store, "ATACTATACAGTAGTTAAAAAGAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT", seqan::CharString("r4"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0,   (TSize)0,   0 + length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0,  (TSize)10,  10 + length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0,  (TSize)20,  20 + length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0,  (TSize)30,  30 + length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0,  (TSize)40,  40 + length(store.readSeqStore[4]));

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, 50, 130, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 200, 0, 20);

    const char * EXPECTED =
            "AACTAAATGCATCCATGTATGCCACAGTGTATACTCTGGAATACTATACAGTAGTTAA----AATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT------------------------------------------------------------\n"
            "..........................................................****......................................\n"
            "          ................................................****................................................\n"
            "                    ......................................****..........................................................\n"
            "                              ............................*****...................................................................\n"
            "                                        ..................AAAG..............................................................................\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_DEFINE_TEST(test_realign_tricky_insert_window_cuts)
{
    seqan::FragmentStore<> store;

    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    appendRead(store, "AACTAAATGCATCCATGTATGCCACAGTGT"                      "CTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCAT", seqan::CharString("r0"));
    appendRead(store,           "ATCCATGTATGCCACAGTGT"                      "CTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAG", seqan::CharString("r1"));
    appendRead(store,                     "GCCACAGTGT"                      "CTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCA", seqan::CharString("r2"));
    appendRead(store,                             "GT""CGCCATACTCGTCCGCGGAG""CTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAA", seqan::CharString("r3"));
    appendRead(store,                                 "CGCCATACTCGTCCGCGGAG""CTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT", seqan::CharString("r4"));

    typedef typename seqan::Size<typename seqan::FragmentStore<>::TAlignedReadStore>::Type TSize;
    appendAlignedRead(store, 0, 0,   (TSize)0,   0 + length(store.readSeqStore[0]) + 20);
    appendAlignedRead(store, 1, 0,  (TSize)10,  10 + length(store.readSeqStore[1]) + 20);
    appendAlignedRead(store, 2, 0,  (TSize)20,  20 + length(store.readSeqStore[2]) + 20);
    appendAlignedRead(store, 3, 0,  (TSize)28,  28 + length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0,  (TSize)30,  30 + length(store.readSeqStore[4]));

    typedef typename seqan::FragmentStore<>::TReadSeqStore     TReadSeqStore;
    typedef typename seqan::Value<TReadSeqStore>::Type         TReadSeq;
    typedef typename seqan::FragmentStore<>::TAlignedReadStore TAlignedReadStore;
    typedef typename seqan::Value<TAlignedReadStore>::Type     TAlignedRead;
    typedef typename seqan::Gaps<TReadSeq, seqan::AnchorGaps<TAlignedRead::TGapAnchors> > TGaps;

    TGaps gaps0(store.readSeqStore[0], store.alignedReadStore[0].gaps);
    insertGaps(gaps0, 30, 20);
    TGaps gaps1(store.readSeqStore[1], store.alignedReadStore[1].gaps);
    insertGaps(gaps1, 20, 20);
    TGaps gaps2(store.readSeqStore[2], store.alignedReadStore[2].gaps);
    insertGaps(gaps2, 10, 20);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "BEFORE REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    reAlignment(store, /*contigID=*/0, 1, 10, false, /*begin=*/0, /*end=*/35, DEBUG_REALIGNMENT);

    if (PRINT_REALIGNMENTS)
    {
        std::cout << "AFTER REALIGNMENT\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        printAlignment(std::cout, layout, store, 0, 0, 400, 0, 1000);
    }

    // Get ASCII art representation of realignment and compare to expected.
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::stringstream ss;
    printAlignment(ss, layout, store, 0, 0, 200, 0, 20);

    const char * EXPECTED =
"AACTAAATGCATCCATGTATGCCACAGT-GT-----------------CTCTGGAATACTATACAGTAGTTAAAATGTGGTATAGCTGAAAGTACAGTACCGAAATGCCATTGCAGAGTAGTAAGACCCCACTTCTATTAAATGAAAAGTT-------------------------------------------------\n"
"............................*..*****************...............................................................\n"
"          ..................*..*****************.........................................................................\n"
"                    ........*..*****************...................................................................................\n"
"                          ..C.CCATACTCGTCCGCGGAG.............................................................................................\n"
"                            C.CCATACTCGTCCGCGGAG.......................................................................................................\n";

    if (DEBUG_REALIGNMENT)
        std::cerr << "" << ss.str() << "---\n" << EXPECTED << "\n";

    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_BEGIN_TESTSUITE(test_realign)
{
	SEQAN_CALL_TEST(test_realign_one_read_no_gaps);
	SEQAN_CALL_TEST(test_realign_one_read_with_gaps);

	SEQAN_CALL_TEST(test_realign_two_reads_stacked_at_beginning);

	SEQAN_CALL_TEST(test_realign_simple_case);
	SEQAN_CALL_TEST(test_realign_simple_gaps);

	SEQAN_CALL_TEST(test_realign_simple_insert);
	SEQAN_CALL_TEST(test_realign_simple_insert_window);
	SEQAN_CALL_TEST(test_realign_simple_insert_window_tight_left);
	SEQAN_CALL_TEST(test_realign_simple_insert_window_tight_right);

    SEQAN_CALL_TEST(test_realign_tricky_insert_window_cuts);
}
SEQAN_END_TESTSUITE
