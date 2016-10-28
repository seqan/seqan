// ==========================================================================
//                                   rna_io
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Lily Shellhammer <lily.shellhammer@gmail.com>
// ==========================================================================

/* TO DO: Need to get rid of warnings for test. Not sure how to get rid of all of them, but I haven't looked
into it a ton just yet.
*/

#ifndef TESTS_RNA_IO_TEST_RNA_IO_H_
#define TESTS_RNA_IO_TEST_RNA_IO_H_

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/rna_io.h>
#include <seqan/graph_types.h>

// A test for connect file reading
SEQAN_DEFINE_TEST(test_rna_io_read_connect)
{
    //Path to example.ct
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.ct");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(rnaPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::RnaIOContext rnaIOContext;
    seqan::RnaRecord rnaRecord;

    readRecord(rnaRecord, rnaIOContext, iter, seqan::Connect());

    /*CHECK CONNECT FILE VALUES */

    SEQAN_ASSERT_EQ(rnaRecord.amount, 73u);
    SEQAN_ASSERT_EQ(rnaRecord.energy, -17.50);
    SEQAN_ASSERT_EQ(rnaRecord.begPos, 1);
    SEQAN_ASSERT_EQ(rnaRecord.endPos, 73);
    SEQAN_ASSERT_EQ(rnaRecord.name,"S.cerevisiae_tRNA-PHE");
    seqan::Rna5String base = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
    SEQAN_ASSERT_EQ(rnaRecord.sequence[0], base);

    seqan::TAdjacencyIterator adj_it(rnaRecord.graph, 0);
    SEQAN_ASSERT_EQ(value(adj_it), 71u);

    /* CHECK DEFAULT VALUES */

    //SEQAN_ASSERT_EQ(rnaRecord.qual, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.offset, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.seqpos, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.annotation, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.comment, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.reactivity, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.reactivity_error, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.xsel, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.xsel_refine, /**/);
    
}

// A test for dot bracket rna file reading.
SEQAN_DEFINE_TEST(test_rna_io_read_dot_bracket)
{
    //Path to example.ct
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.dt");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(rnaPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::RnaIOContext rnaIOContext;
    seqan::RnaRecord rnaRecord;

    readRecord(rnaRecord, rnaIOContext, iter, seqan::DotBracket());

    /*CHECK DOTBRACKET FILE VALUES */

    SEQAN_ASSERT_EQ(rnaRecord.amount, 73u);
    SEQAN_ASSERT_EQ(rnaRecord.energy, -17.50);
    SEQAN_ASSERT_EQ(rnaRecord.begPos, 1);
    SEQAN_ASSERT_EQ(rnaRecord.endPos, 73);
    SEQAN_ASSERT_EQ(rnaRecord.name,"S.cerevisiae_tRNA-PHE");
    SEQAN_ASSERT_EQ(rnaRecord.sequence[0], "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA");

    seqan::TAdjacencyIterator adj_it(rnaRecord.graph, 0);
    SEQAN_ASSERT_EQ(value(adj_it), 71u);
    
    /* CHECK DEFAULT VALUES */

    //SEQAN_ASSERT_EQ(rnaRecord.qual, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.offset, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.seqpos, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.annotation, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.comment, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.reactivity, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.reactivity_error, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.xsel, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.xsel_refine, /**/);

}

///////////////////BPSEQ TEST NOT COMPLETE////////////////////////
SEQAN_DEFINE_TEST(test_rna_io_read_bpseq)
{
    //Path to example.ct
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.bpseq");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(rnaPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::RnaIOContext rnaIOContext;
    seqan::RnaRecord rnaRecord;
    //readHeader(rnaHeader, rnaIOCOntext, iter, seqan::Bpseq());
    readRecord(rnaRecord, rnaIOContext, iter, seqan::Bpseq());

    // CHECK BPSEQ FILE VALUES 
/*
    SEQAN_ASSERT_EQ(rnaRecord.amount, 73u);
    SEQAN_ASSERT_EQ(rnaRecord.begPos, 1);
    SEQAN_ASSERT_EQ(rnaRecord.endPos, 73);
    //SEQAN_ASSERT_EQ(rnaRecord.name,"S.cerevisiae_tRNA-PHE");
    SEQAN_ASSERT_EQ(rnaRecord.sequence[0], "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA");
    seqan::TAdjacencyIterator adj_it(rnaRecord.graph, 0);
    SEQAN_ASSERT_EQ(value(adj_it), 71u);
*/
    // CHECK DEFAULT VALUES 

    //SEQAN_ASSERT_EQ(rnaRecord.qual, /**/ //);
    //SEQAN_ASSERT_EQ(rnaRecord.offset, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.seqpos, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.annotation, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.comment, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.reactivity, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.reactivity_error, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.xsel, /**/);
    //SEQAN_ASSERT_EQ(rnaRecord.xsel_refine, /**/);

}

SEQAN_DEFINE_TEST(test_rna_write_connect_record)
{
    seqan::RnaRecord record;
    //set values
    record.amount = 8u;
    record.begPos = 1;
    record.endPos = 8;
    record.name = "S.cerevisiae_tRNA-PHE";
    record.energy = -17.5;
    appendValue(record.sequence, (seqan::Rna5String)"GCGGAUUU");

    for (unsigned idx = 0; idx < record.amount; ++idx)
        addVertex(record.graph);
    for (unsigned idx = 0; idx < 4; ++idx)
        addEdge(record.graph, idx, 7u - idx, 1.);


    // Write Connect records to string stream.String<char> out;
    seqan::String<char> out;
    writeRecord(out, record, seqan::Connect());

    // Compare string stream to expected value.
    seqan::String<char> expected = "8 ENERGY = \t-17.5\tS.cerevisiae_tRNA-PHE\n";
     append(expected, " 1 G\t0\t2\t8\t1\n");
     append(expected, " 2 C\t1\t3\t7\t2\n");
     append(expected, " 3 G\t2\t4\t6\t3\n");
     append(expected, " 4 G\t3\t5\t5\t4\n");
     append(expected, " 5 A\t4\t6\t4\t5\n");
     append(expected, " 6 U\t5\t7\t3\t6\n");
     append(expected, " 7 U\t6\t8\t2\t7\n");
     append(expected, " 8 U\t7\t9\t1\t8\n");
    SEQAN_ASSERT_EQ(out, expected);
}

SEQAN_DEFINE_TEST(test_rna_write_dot_bracket_record)
{
    seqan::RnaRecord record;
    //set values
    record.amount = 8u;
    record.begPos = 1;
    record.endPos = 8;
    record.name = "S.cerevisiae_tRNA-PHE";
    record.energy = -17.5;
    appendValue(record.sequence, (seqan::Rna5String)"GCGGAUUU");

    for (unsigned idx = 0; idx < record.amount; ++idx)
        addVertex(record.graph);
    for (unsigned idx = 0; idx < 4; ++idx)
        addEdge(record.graph, idx, 7u - idx, 1.);

    // Write Connect records to string stream.String<char> out;
    seqan::String<char> out;
    writeRecord(out, record, seqan::DotBracket());

    // Compare string stream to expected value.
    seqan::String<char> expected = ">S.cerevisiae_tRNA-PHE /1-8\n";
    append(expected, "GCGGAUUU\n");
    append(expected, "(((()))) (-17.5)\n");
    SEQAN_ASSERT_EQ(out, expected);
}

#endif  // TESTS_RNA_IO_TEST_RNA_IO_H_
