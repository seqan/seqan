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
// Authors: Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_RNA_IO_TEST_RNA_IO_H_
#define TESTS_RNA_IO_TEST_RNA_IO_H_

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/rna_io.h>
#include <seqan/graph_types.h>
#include <seqan/align.h>

// ----------------------------------------------------------------------------
// Connect File I/O
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rna_io_read_connect)
{
    // Path to example.ct
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.ct");
    seqan::RnaStructFileIn inputfile(seqan::toCString(rnaPath), seqan::OPEN_RDONLY);
    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, inputfile);

    SEQAN_ASSERT(rnaRecord.hasUndefinedID());
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 73u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.fixedGraphs[0].energy, -17.50f);
    SEQAN_ASSERT_EQ(rnaRecord.name, "S.cerevisiae_tRNA-PHE");
    seqan::Rna5String base = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
    SEQAN_ASSERT_EQ(rnaRecord.sequence, base);
    seqan::RnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 0);
    SEQAN_ASSERT_EQ(value(adj_it), 71u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "");
    // UNUSED: seqID, align, bppMatrGraphs, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_connect)
{
    seqan::RnaRecord record{};
    //set values
    record.name = "S.cerevisiae_tRNA-PHE";
    record.sequence = "GCGGAUUU";
    record.seqLen = length(record.sequence);
    seqan::RnaStructureGraph graph;

    for (typename seqan::Size<seqan::Rna5String>::Type idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph.inter);

    for (unsigned idx = 0u; idx < 4u; ++idx)
        addEdge(graph.inter, idx, 7u - idx, 1.0);

    graph.energy = -17.5f;
    append(record.fixedGraphs, graph);

    // Write records to string stream.String<char> out;
    seqan::CharString outstr;
    writeRecord(outstr, record, seqan::Connect());

    // Compare string stream to expected value.
    seqan::String<char> expected = "8\tENERGY = -17.5\tS.cerevisiae_tRNA-PHE\n"
            " 1\tG\t0\t2\t8\t1\n"
            " 2\tC\t1\t3\t7\t2\n"
            " 3\tG\t2\t4\t6\t3\n"
            " 4\tG\t3\t5\t5\t4\n"
            " 5\tA\t4\t6\t4\t5\n"
            " 6\tU\t5\t7\t3\t6\n"
            " 7\tU\t6\t8\t2\t7\n"
            " 8\tU\t7\t9\t1\t8\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// DotBracket File I/O
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rna_io_read_dot_bracket)
{
    //Path to example.dbn
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.dbn");
    seqan::RnaStructFileIn inputfile(seqan::toCString(rnaPath), seqan::OPEN_RDONLY);
    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, inputfile);

    SEQAN_ASSERT(rnaRecord.hasUndefinedID());
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 73u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.fixedGraphs[0].energy, -17.50f);
    SEQAN_ASSERT_EQ(rnaRecord.name, "S.cerevisiae_tRNA-PHE M10740");
    seqan::Rna5String base = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
    SEQAN_ASSERT_EQ(rnaRecord.sequence, base);
    seqan::RnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 0);
    SEQAN_ASSERT_EQ(value(adj_it), 71u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "");
    // UNUSED: seqID, align, bppMatrGraphs, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_dot_bracket)
{
    seqan::RnaRecord record{};
    //set values
    record.name = "S.cerevisiae_tRNA-PHE";
    record.sequence = "GCGGAUUU";
    record.seqLen = length(record.sequence);
    seqan::RnaStructureGraph graph;

    for (typename seqan::Size<seqan::Rna5String>::Type idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph.inter);

    for (unsigned idx = 0u; idx < 3u; ++idx)
        addEdge(graph.inter, idx, 7u - idx, 1.0);

    graph.energy = -17.5f;
    append(record.fixedGraphs, graph);

    // Write records to string stream.String<char> out;
    seqan::CharString outstr;
    writeRecord(outstr, record, seqan::DotBracket());

    // Compare string stream to expected value.
    seqan::String<char> expected = ">S.cerevisiae_tRNA-PHE/1-8\n"
            "GCGGAUUU\n"
            "(((..))) (-17.5)\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// Vienna File I/O
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rna_io_read_vienna)
{
    //Path to example.dbv
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.dbv");
    seqan::RnaStructFileIn inputfile(seqan::toCString(rnaPath), seqan::OPEN_RDONLY);
    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, inputfile);

    SEQAN_ASSERT(rnaRecord.hasUndefinedID());
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 66u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.fixedGraphs[0].energy, 0.0f);
    SEQAN_ASSERT_EQ(rnaRecord.name, "gi-12082738-AF304460");
    seqan::Rna5String base = "AGUCUUAUACACAAUGGUAAGCCAGUGGUAGUAAAGGUAUAAGAAAUUUGCUACUAUGUUACUGAA";
    SEQAN_ASSERT_EQ(rnaRecord.sequence, base);
    seqan::RnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 3);
    SEQAN_ASSERT_EQ(value(adj_it), 42u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "");
    // UNUSED: seqID, align, bppMatrGraphs, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_vienna)
{
    seqan::RnaRecord record{};
    //set values
    record.name = "S.cerevisiae_tRNA-PHE";
    record.sequence = "GCGGAUUU";
    record.seqLen = length(record.sequence);
    seqan::RnaStructureGraph graph;

    for (typename seqan::Size<seqan::Rna5String>::Type idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph.inter);

    for (unsigned idx = 0u; idx < 3u; ++idx)
        addEdge(graph.inter, idx, 7u - idx, 1.0);

    append(record.fixedGraphs, graph);

    // Write records to string stream.String<char> out;
    seqan::CharString outstr;
    writeRecord(outstr, record, seqan::Vienna());

    // Compare string stream to expected value.
    seqan::String<char> expected = ">S.cerevisiae_tRNA-PHE/1-8\n"
        "GCGGAUUU\n"
        "(((..)))\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// Stockholm File I/O
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rna_io_read_stockholm)
{
    //Path to example.sth
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.sth");
    seqan::RnaStructFileIn inputfile(seqan::toCString(rnaPath), seqan::OPEN_RDONLY);
    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, inputfile);

    SEQAN_ASSERT(rnaRecord.hasUndefinedID());
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 74u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.fixedGraphs[0].energy, 0.0f);
    SEQAN_ASSERT_EQ(rnaRecord.name, "trna");
    seqan::Rna5String base = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
    SEQAN_ASSERT_EQ(stringSet(rnaRecord.align)[0], base);
    base = "UCCGUGAUAGUUUAAUGGUCAGAAUGGGCGCUUGUCGCGUGCCAGAUCGGGGUUCAAUUCCCCGUCGCGGAG";
    SEQAN_ASSERT_EQ(stringSet(rnaRecord.align)[2], base);
    SEQAN_ASSERT_EQ(rnaRecord.seqID[0], "DF6280");
    SEQAN_ASSERT_EQ(rnaRecord.seqID[2], "DD6280");

    seqan::RnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 10);
    SEQAN_ASSERT_EQ(value(adj_it), 24u);
    // UNUSED: sequence, bppMatrGraphs, quality, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_stockholm)
{
    seqan::RnaRecord record{};
    // set values
    record.seqLen = 9u;
    record.name = "trna";
    record.comment = "alignment of 1415 tRNAs";
    seqan::Rna5String seq1 = "GCGGAUUU";
    seqan::Rna5String seq2 = "UCCGAUAUA";
    // create alignment
    seqan::resize(seqan::rows(record.align), 2);
    seqan::assignSource(seqan::row(record.align, 0), seq1);
    seqan::assignSource(seqan::row(record.align, 1), seq2);
    seqan::insertGap(seqan::row(record.align, 0), 4);
    // set sequence identifiers
    seqan::appendValue(record.seqID, seqan::CharString{"seq0"});
    seqan::appendValue(record.seqID, seqan::CharString{"seq1"});

    seqan::RnaStructureGraph graph;
    for (typename seqan::Size<seqan::Rna5String>::Type idx = 0u; idx < record.seqLen; ++idx)
        addVertex(graph.inter);

    for (unsigned idx = 0u; idx < 4u; ++idx)
        addEdge(graph.inter, idx, 7u - idx, 1.0);

    append(record.fixedGraphs, graph);

    // Write records to string stream.String<char> out;
    seqan::String<char> outstr;
    writeRecord(outstr, record, seqan::Stockholm());

    // Compare string stream to expected value.
    seqan::String<char> expected = "# STOCKHOLM 1.0\n"
            "#=GF ID      trna\n"
            "#=GF DE      alignment of 1415 tRNAs\n\n"
            "seq0        \tGCGG-AUUU\n"
            "seq1        \tUCCGAUAUA\n"
            "#=GC SS_cons\t(((()))).\n"
            "//\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// Bpseq File I/O
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rna_io_read_bpseq)
{
    // Path to example.bpseq
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.bpseq");
    seqan::RnaStructFileIn inputfile(seqan::toCString(rnaPath), seqan::OPEN_RDONLY);
    seqan::RnaRecord rnaRecord;
    readRecord(rnaRecord, inputfile);

    SEQAN_ASSERT(rnaRecord.hasUndefinedID());
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 50u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(rnaRecord.fixedGraphs[0].energy, 0.0f);
    seqan::Rna5String base = "GGGCCGGGCGCGGUGGCGCGCGCCUGUAGUCCCAGCUACUCGGGAGGCUC";
    SEQAN_ASSERT_EQ(rnaRecord.sequence, base);
    seqan::RnaAdjacencyIterator adj_it(rnaRecord.fixedGraphs[0].inter, 1);
    SEQAN_ASSERT_EQ(value(adj_it), 48u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "");
    seqan::CharString comment = " PDB ID 1E8O Signal Recognition Particle (SRP) RNA ";
    SEQAN_ASSERT_EQ(rnaRecord.comment, comment);
    seqan::CharString name = "A header line beginning with # is for comments not for actual structure information.";
    SEQAN_ASSERT_EQ(rnaRecord.name, name);
    // UNUSED: seqID, align, bppMatrGraphs, typeID, reactivity, reactError
}

SEQAN_DEFINE_TEST(test_rna_io_write_bpseq)
{
    seqan::RnaRecord record{};
    //set values
    record.name = "S.cerevisiae_tRNA-PHE";
    record.sequence = "GCGGAUUU";
    record.seqLen = length(record.sequence);
    seqan::RnaStructureGraph graph;

    for (typename seqan::Size<seqan::Rna5String>::Type idx = 0; idx < record.seqLen; ++idx)
        addVertex(graph.inter);

    for (unsigned idx = 0u; idx < 4u; ++idx)
        addEdge(graph.inter, idx, 7u - idx, 1.0);

    append(record.fixedGraphs, graph);

    // Write records to string stream.String<char> out;
    seqan::String<char> outstr;
    seqan::writeRecord(outstr, record, seqan::Bpseq());

    // Compare string stream to expected value.
    seqan::String<char> expected = "# S.cerevisiae_tRNA-PHE\n"
            "1\tG\t8\n2\tC\t7\n3\tG\t6\n4\tG\t5\n5\tA\t4\n6\tU\t3\n7\tU\t2\n8\tU\t1\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

// ----------------------------------------------------------------------------
// Ebpseq File I/O
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_rna_io_read_ebpseq)
{
    // Path to example.ebpseq
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.ebpseq");
    seqan::RnaStructFileIn inputfile(seqan::toCString(rnaPath), seqan::OPEN_RDONLY);
    seqan::RnaStructContents contents;
    seqan::readRecords(contents, inputfile, 100u);

    seqan::RnaRecord & rnaRecord = contents.records[0];
    seqan::RnaStructureGraph & fgraph = rnaRecord.fixedGraphs[0];
    seqan::RnaStructureGraph & bgraph = rnaRecord.bppMatrGraphs[1];

    SEQAN_ASSERT_EQ(rnaRecord.recordID, 0u);
    SEQAN_ASSERT_EQ(rnaRecord.seqLen, 8u);
    SEQAN_ASSERT_EQ(rnaRecord.offset, 1u);
    SEQAN_ASSERT_EQ(fgraph.energy, 0.0f);
    SEQAN_ASSERT_EQ(rnaRecord.sequence, "AGUCCGUC");
    seqan::RnaAdjacencyIterator adj_it(fgraph.inter, 1);
    SEQAN_ASSERT_EQ(value(adj_it), 6u);
    seqan::RnaAdjacencyIterator adj_it2(rnaRecord.fixedGraphs[1].inter, 1);
    SEQAN_ASSERT_EQ(value(adj_it2), 3u);
    SEQAN_ASSERT_EQ(rnaRecord.quality, "H{#!7T@g");
    SEQAN_ASSERT_EQ(rnaRecord.comment, "");
    SEQAN_ASSERT_EQ(rnaRecord.name, "Name of first sequence");
    SEQAN_ASSERT_EQ(fgraph.specs, "First fixed structure computed with a tool such as Ipknot, dotknot, RNAfold etc..");
    SEQAN_ASSERT_EQ(bgraph.specs, "Second base pair probability matrix computed with a tool such RNAfold");
    seqan::RnaAdjacencyIterator adj_it3(bgraph.inter, 1);
    SEQAN_ASSERT_EQ(value(adj_it3), 3u);
    SEQAN_ASSERT_EQ(getCargo(findEdge(bgraph.inter, 1, 3)), 0.6);
    SEQAN_ASSERT_EQ(getCargo(findEdge(bgraph.inter, 6, 1)), 0.3);
    SEQAN_ASSERT_EQ(getCargo(findEdge(bgraph.inter, 2, 5)), 0.8);
    SEQAN_ASSERT_EQ(getCargo(findEdge(bgraph.inter, 0, 3)), 0.002);
    SEQAN_ASSERT_EQ(rnaRecord.typeID[0], 1u);

    rnaRecord = contents.records[2];
    SEQAN_ASSERT_EQ(rnaRecord.reactivity[0][1], 46.9128f);
    SEQAN_ASSERT_EQ(rnaRecord.reactivity[1][1], 0.3916f);
    SEQAN_ASSERT_EQ(rnaRecord.reactError[0][1], 13.8533f);
    SEQAN_ASSERT_EQ(rnaRecord.reactError[1][1], 0.1568f);
    SEQAN_ASSERT_EQ(length(rnaRecord.typeID), 2u);
    SEQAN_ASSERT_EQ(rnaRecord.typeID[1], 1u);
    // UNUSED: seqID, align
}

SEQAN_DEFINE_TEST(test_rna_io_write_ebpseq)
{
    seqan::RnaHeader header{};
    seqan::appendValue(header.seqLabels, "a sequence");
    seqan::appendValue(header.fixLabels, "the fixed structure");
    seqan::appendValue(header.bppLabels, "the bpp structure");
    seqan::appendValue(header.typeLabels, "fictional data");

    seqan::RnaRecord record{};
    record.recordID = 0u;
    record.seqLen = 3u;
    record.offset = 1u;
    record.quality = "H7T";
    record.sequence = "AGU";
    record.comment = "these are artificial RNA structures";
    seqan::appendValue(record.typeID, 0u);

    seqan::RnaStructureGraph fgraph;
    fgraph.specs = "the fixed structure";
    for (typename seqan::Size<seqan::Rna5String>::Type idx = 0; idx < record.seqLen; ++idx)
        addVertex(fgraph.inter);

    seqan::addEdge(fgraph.inter, 0u, 2u, 1.0);
    append(record.fixedGraphs, fgraph);

    seqan::RnaStructureGraph bgraph;
    bgraph.specs = "the bpp structure";
    for (typename seqan::Size<seqan::Rna5String>::Type idx = 0u; idx < record.seqLen; ++idx)
        addVertex(bgraph.inter);
    seqan::addEdge(bgraph.inter, 0u, 2u, 0.8);
    seqan::addEdge(bgraph.inter, 1u, 2u, 0.2);
    seqan::append(record.bppMatrGraphs, bgraph);

    seqan::String<float> reactivity;
    seqan::appendValue(reactivity, 2.4f);
    seqan::appendValue(reactivity, 1.9f);
    seqan::appendValue(reactivity, 4.0f);
    seqan::appendValue(record.reactivity, reactivity);

    seqan::clear(reactivity);
    seqan::appendValue(reactivity, 0.03f);
    seqan::appendValue(reactivity, 0.10f);
    seqan::appendValue(reactivity, 0.08f);
    seqan::appendValue(record.reactError, reactivity);

    // Write records to string stream.String<char> out;
    seqan::RnaIOContext context;
    seqan::String<char> outstr;
    seqan::writeHeader(outstr, header, context, seqan::Ebpseq());
    seqan::writeRecord(outstr, record, context, seqan::Ebpseq());

    // Compare string stream to expected value.
    seqan::String<char> expected = "## S1: a sequence\n"
        "## F1: the fixed structure\n"
        "## M1: the bpp structure\n"
        "## T1: fictional data\n"
        "# S1 T1\n"
        "# I\tNT\tQU\tR1\tRE1\tF1\tM1\n"
        "1\tA\tH\t2.4\t0.03\t3\t<3/0.8>\n"
        "2\tG\t7\t1.9\t0.1\t0\t<3/0.2>\n"
        "3\tU\tT\t4\t0.08\t1\t<2/0.2 | 1/0.8>\n";
    SEQAN_ASSERT_EQ(outstr, expected);
}

SEQAN_DEFINE_TEST(test_rna_io_convert)
{
    seqan::CharString rnaPath = SEQAN_PATH_TO_ROOT();
    append(rnaPath, "/tests/rna_io/example.dbv");
    seqan::RnaStructFileIn inputfile(seqan::toCString(rnaPath), seqan::OPEN_RDONLY);

    seqan::RnaStructContents filecontents;
    seqan::readRecords(filecontents, inputfile, 100u);
    SEQAN_ASSERT_EQ(length(filecontents.records), 3u);

    std::string outpath = std::string(SEQAN_TEMP_FILENAME()) + ".ebpseq";
    seqan::RnaStructFileOut outputfile(seqan::toCString(outpath), seqan::OPEN_WRONLY);
    seqan::writeRecords(outputfile, filecontents);
}

#endif  // TESTS_RNA_IO_TEST_RNA_IO_H_
