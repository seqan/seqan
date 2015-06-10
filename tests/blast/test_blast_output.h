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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Tests for the blast module
// ==========================================================================

#ifndef SEQAN_TESTS_TEST_BLAST_OUTPUT_H_
#define SEQAN_TESTS_TEST_BLAST_OUTPUT_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/blast.h>

using namespace seqan;

template <typename TFile, typename TScore, typename TRecords, BlastProgram p, BlastTabularSpec h>
inline void
test_blast_write_do(TFile & file,
                    BlastIOContext<TScore, p, h> & context,
                    TRecords const & records,
                    int const format,
                    int const,
                    BlastReport const & /**/)
{
    if (format < 3)
    {
        writeHeader(file, context, BlastReport());
        for (auto const & r : records)
            writeRecord(file, context, r, BlastReport());
        writeFooter(file, context, BlastReport());
    } else
    {
        BlastReportFileOut<BlastIOContext<TScore, p, h>> out(context, file, BlastReport());
        writeHeader(out);
        for (auto const & r : records)
            writeRecord(out, r);
        writeFooter(out);
    }
}

template <typename TFile, typename TScore, typename TRecords, BlastProgram p, BlastTabularSpec h>
inline void
test_blast_write_do(TFile & file,
                    BlastIOContext<TScore, p, h> & context,
                    TRecords const & records,
                    int const format,
                    int const custom,
                    BlastTabular const & /**/)
{
    switch (format)
    {
        case 0: // iterate over matches
        case 1: // iterate over headers and matches
            if ((format == 1) && (custom <= 1))
                writeHeader(file, context, BlastTabular()); // noop for TABULARs
            for (auto const & r : records)
            {
                if ((format == 1) && (custom <= 1))
                    _writeRecordHeader(file, context, r, BlastTabular());

                for (auto const & m : r.matches)
                {
                    if (custom <= 1)
                        _writeMatch(file, context, m, BlastTabular());
                    else
                        writeMatch(file,
                                    BlastTabularLL(),
                                    m.qId,
                                    m.sId,
                                    m.alignStats.alignmentLength,
                                    m.alignStats.numMismatches,
                                    m.alignStats.numGaps,
                                    m.eValue,
                                    m.bitScore);
                }
            }
            if ((format == 1) && (custom <= 1))
                writeFooter(file, context, BlastTabular());
            break;
        case 2: // iteratre over records
            writeHeader(file, context, BlastTabular()); // noop for TABULARs
            for (auto const & r : records)
                writeRecord(file, context, r, BlastTabular());
            writeFooter(file, context, BlastTabular());
            break;
        case 3: // formatted file out
        {
            BlastTabularFileOut<BlastIOContext<TScore, p, h>> out(context, file, BlastTabular());
            writeHeader(out); // noop for TABULARs
            for (auto const & r : records)
                writeRecord(out, r);
            writeFooter(out);
        } break;
    }
}

template <typename TFile, typename TScore, typename TFormat, BlastProgram p, BlastTabularSpec h>
inline void
test_blast_write_record_match(TFile & file,
                              int const format,
                              int const customFields,
                              BlastIOContext<TScore, p, h> & context,
                              TFormat const &)
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;
    typedef BlastRecord<CharString, CharString, uint32_t, TAlign> TBlastRecord;

    StringSet<String<AminoAcid>, Owner<ConcatDirect<>>> queries;
    StringSet<CharString, Owner<ConcatDirect<>>> qIds;
    StringSet<String<AminoAcid>> subjects;
    StringSet<CharString> sIds;

    resize(subjects, 2);
    resize(sIds, 2);

    subjects[0] =
    "SSITEEKHIPHKEQDKDAEFLSKEALKTHMTENVLQMDRRAVQDPSTSFLQLLKAKGLLG"
    "LPDYEVNLADVNSPGFRKVAYAQTKPRRLCFPNGGTRRGSFIMDTAVVVMVSLRYVNIGK"
    "VIFPGATDVSEGEDEFWAGLPQAYGCLATEFLCIHIAIYSWIHVQSSRYDDMNASVIRAK"
    "LNLAVITSWTQLIQAEKETI";

    subjects[1] =
    "GATRDSKGNAVITSFTQARLRVYADLLGPYWIILHVIELTGVGNTGQKCTLNHMGTYAVF"
    "DLKQPPATNDLGLPKPCFIGFDIQNELAIGTVGHSEAVIAAFTQRDRLEERAESKQSLAR"
    "PVISPKLIAEVSTVLESALNQMYSSLGFYRVERAEDYAQPRKLCVVKKKSFNCLNADIWL"
    "EYRMEDQKSVPKVFKIMMDD";

    sIds[0] = "Subject_Numero_Uno";
    sIds[1] = "Subject_Numero_Dos";

    appendValue(queries, "VAYAQPRKLCYP");
    appendValue(queries, "NAYPRUTEIN");
    appendValue(queries, "AVITSFTQ");

    appendValue(qIds, "Query_Numero_Uno with args");
    appendValue(qIds, "Query_Numero_Dos with args");
    appendValue(qIds, "Query_Numero_Tres with args");

    setScoreGapOpenBlast(context.scoringScheme, -11);
    setScoreGapExtend(context.scoringScheme, -1);
    SEQAN_ASSERT(isValid(context.scoringScheme));

    String<TBlastRecord> records;
    resize(records, 3);

    context.dbName = "The Foo Database";
    context.dbTotalLength = length(subjects[0]) + length(subjects[1]);
    context.dbNumberOfSeqs = 2;

    for (int q = 0; q < 3; ++q)
    {
        records[q].qId = qIds[q];
        records[q].qLength = length(queries[q]);

        if (q == 1) // we want the second record to have no matches
            continue;

        for (int s = 0; s < 2; ++s)
        {
            records[q].matches.emplace_back(qIds[q], sIds[s]);
            TBlastMatch & m = records[q].matches.back();
            resize(rows(m.align), 2);
            assignSource(row(m.align, 0), queries[q]);
            assignSource(row(m.align, 1), subjects[s]);

            localAlignment(m.align, seqanScheme(context.scoringScheme));

            m.qStart = beginPosition(row(m.align, 0));
            m.qEnd   = endPosition(row(m.align, 0));
            m.sStart = beginPosition(row(m.align, 1));
            m.sEnd   = endPosition(row(m.align, 1));

            m.qLength = length(queries[q]);
            m.sLength = length(subjects[s]);

            m.qFrameShift = 1;
            m.sFrameShift = 1;

            computeAlignmentStats(m, context);
            computeBitScore(m, context);
            computeEValue(m, context);
        }
    }

    // sort by bit-score
    for (auto & r : records)
        r.matches.sort();

    test_blast_write_do(file, context, records, format, customFields, TFormat());
}

template <typename TScore, BlastProgram p, BlastTabularSpec h>
void test_blast_write_tabular_impl(int const format, // 0 only matches; 1 matches+headers; 2 records; 3 formattedFile
                                   int const customFields, // 0 defaults; 1 customfields; 2 lowlevel impl
                                   BlastIOContext<TScore, p, h> & context)
{
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    test_blast_write_record_match(fstream, format, customFields, context, BlastTabular());

    std::string contents;
    resize(contents, fstream.tellg());
    fstream.seekg(0, std::ios::beg);
    fstream.read(&contents[0], contents.size());
    fstream.close();

    // Generate the output string for comparison (since this includes the SeqAn version we don't want it completely
    // hardcoded (and also to save lots of redundant code for the immense combinations of options)

    std::string versionLine;
    if (!context.legacyFormat)
        versionLine.append("# BLASTP 2.2.26+ ");
    else
        versionLine.append("# BLASTP 2.2.26 ");
    versionLine.append("[I/O Module of SeqAn-");
    versionLine.append(std::to_string(SEQAN_VERSION_MAJOR));
    versionLine.append(".");
    versionLine.append(std::to_string(SEQAN_VERSION_MINOR));
    versionLine.append(".");
    versionLine.append(std::to_string(SEQAN_VERSION_PATCH));
    versionLine.append(", http://www.seqan.de]\n");

    std::string fieldsLine = "# Fields: ";
    if (context.legacyFormat)
    {
        fieldsLine.append(BlastMatchField<>::legacyColumnLabels);
    }
    else
    {
        switch(customFields)
        {
            case 0:
                fieldsLine.append(BlastMatchField<>::columnLabels[0]);
                break;
            case 1:
                fieldsLine.append(BlastMatchField<>::columnLabels[0]);
                fieldsLine.append(", score, query/sbjct frames");
                break;
            case 2:
                fieldsLine.append("Query id, Subject id, alignment length, mismatches, gaps, e-value, bit score");
                break;
        }
    }
    fieldsLine.append("\n");


    std::string compString;
    if ((context.tabularSpec == BlastTabularSpec::HEADER) && (format > 0))
    {
        compString = versionLine;
        compString.append("# Query: Query_Numero_Uno with args\n"
                          "# Database: The Foo Database\n");
        compString.append(fieldsLine);
        if (!context.legacyFormat)
            compString.append("# 2 hits found\n");
    }
    // legacy format counts gap characters to mismatches as well
    unsigned mismatches = context.legacyFormat ? 4 : 2;

    if (customFields == 0)
    {
        compString.append("Query_Numero_Uno\tSubject_Numero_Uno\t71.43\t14\t");
        compString.append(std::to_string(mismatches));
        compString.append(                                                    "\t1\t1\t12\t79\t92\t5e-04\t23.1\n"
                          "Query_Numero_Uno\tSubject_Numero_Dos\t100.00\t8\t0\t0\t3\t10\t157\t164\t0.001\t22.3\n");
    }
    else if (customFields == 1)
        compString.append(
        "Query_Numero_Uno\tSubject_Numero_Uno\t71.43\t14\t2\t1\t1\t12\t79\t92\t5e-04\t23.1\t48\t0/0\n"
        "Query_Numero_Uno\tSubject_Numero_Dos\t100.00\t8\t0\t0\t3\t10\t157\t164\t0.001\t22.3\t46\t0/0\n");
    else
        compString.append(
        "Query_Numero_Uno with args\tSubject_Numero_Uno\t14\t2\t2\t0.000534696\t23.0978\n"
        "Query_Numero_Uno with args\tSubject_Numero_Dos\t8\t0\t0\t0.000912053\t22.3274\n");

    if ((context.tabularSpec == BlastTabularSpec::HEADER) && (format > 0))
    {
        // header for empty record
        compString.append(versionLine);
        compString.append("# Query: Query_Numero_Dos with args\n"
                          "# Database: The Foo Database\n");
        if (context.legacyFormat)
            compString.append(fieldsLine);
        if (!context.legacyFormat)
            compString.append("# 0 hits found\n");

        // header for
        compString.append(versionLine);
        compString.append("# Query: Query_Numero_Tres with args\n"
                          "# Database: The Foo Database\n");
        compString.append(fieldsLine);
        if (!context.legacyFormat)
            compString.append("# 2 hits found\n");
    }

    if (customFields == 0)
        compString.append("Query_Numero_Tres\tSubject_Numero_Dos\t100.00\t8\t0\t0\t1\t8\t10\t17\t0.007\t18.9\n"
                          "Query_Numero_Tres\tSubject_Numero_Uno\t87.50\t8\t1\t0\t1\t8\t184\t191\t0.026\t16.9\n");
    else if (customFields == 1)
        compString.append("Query_Numero_Tres\tSubject_Numero_Dos\t100.00\t8\t0\t0\t1\t8\t10\t17\t0.007\t18.9\t37\t0/0\n"
                          "Query_Numero_Tres\tSubject_Numero_Uno\t87.50\t8\t1\t0\t1\t8\t184\t191\t0.026\t16.9\t32\t0/0\n");
    else
        compString.append("Query_Numero_Tres with args\tSubject_Numero_Dos\t8\t0\t0\t0.00672262\t18.8606\n"
                          "Query_Numero_Tres with args\tSubject_Numero_Uno\t8\t1\t0\t0.0255459\t16.9346\n");

    if ((context.tabularSpec == BlastTabularSpec::HEADER) && (!context.legacyFormat) && (format >= 1))
        compString.append("# BLAST processed 3 queries\n");

    if (contents != compString)
    {
        for (uint32_t i = 0; i < length(contents); ++i)
        {
            if (contents[i] != compString[i])
            {
                std::cout << contents.substr(0,i) << "\n";
                std::cout << "CONT: \"" << contents[i] << "\"\n";
                std::cout << "COMP: \"" << compString[i] << "\"\n";
                break;
            }
        }
    }
    SEQAN_ASSERT_EQ(contents, compString);
}

// TEST MATCH, WITHOUT HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_match_tabular_run_time_context_args)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;
    context.tabularSpec = BlastTabularSpec::NO_HEADER;

    test_blast_write_tabular_impl(0, 0, context);
}

SEQAN_DEFINE_TEST(test_blast_write_match_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    test_blast_write_tabular_impl(0, 0, context);
}

SEQAN_DEFINE_TEST(test_blast_write_match_tabular_legacy)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    context.legacyFormat = true;
    test_blast_write_tabular_impl(0, 0, context);
}

// TEST MATCH, WITH HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_match_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    test_blast_write_tabular_impl(0, 0, context);
}

SEQAN_DEFINE_TEST(test_blast_write_match_tabular_with_header_legacy)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    context.legacyFormat = true;
    test_blast_write_tabular_impl(0, 0, context);
}

// TEST HEADER, WITHOUT HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_header_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    test_blast_write_tabular_impl(1, 0, context);
}

SEQAN_DEFINE_TEST(test_blast_write_header_tabular_legacy)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    context.legacyFormat = true;
    test_blast_write_tabular_impl(1, 0, context);
}

// TEST HEADER, WITH HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_header_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    test_blast_write_tabular_impl(1, 0, context);
}

SEQAN_DEFINE_TEST(test_blast_write_header_tabular_with_header_legacy)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    context.legacyFormat = true;
    test_blast_write_tabular_impl(1, 0, context);
}

// TEST MATCH, WITHOUT HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_match_customfields_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(0, 1, context);
}

// TEST MATCH, WITH HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_match_customfields_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(0, 1, context);
}

// TEST HEADER, WITHOUT HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_header_customfields_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(1, 1, context);
}

// TEST HEADER, WITH HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_header_customfields_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(1, 1, context);
}

// TEST MATCH, WITHOUT HEADER, CUSTOM COLUMNS THROUGH VARIADIC TEMPLATES IMPLEMENTATION

SEQAN_DEFINE_TEST(test_blast_write_match_lowlevel_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::UNKNOWN, BlastTabularSpec::UNKNOWN> context;
    // the BlastIOContext will not be used downstream at all since this uses the low-level implementation
    // that is independent of context and other data structures
    test_blast_write_tabular_impl(0, 2, context);
}

// TEST RECORD, WITHOUT HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_record_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    test_blast_write_tabular_impl(2, 0, context);
}

SEQAN_DEFINE_TEST(test_blast_write_record_tabular_legacy)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    context.legacyFormat = true;
    test_blast_write_tabular_impl(2, 0, context);
}

// TEST RECORD, WITH HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_record_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    test_blast_write_tabular_impl(2, 0, context);
}

SEQAN_DEFINE_TEST(test_blast_write_record_tabular_with_header_legacy)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    context.legacyFormat = true;
    test_blast_write_tabular_impl(2, 0, context);
}

// TEST RECORD, WITHOUT HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_record_customfields_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(2, 1, context);
}

// TEST RECORD, WITH HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_record_customfields_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(2, 1, context);
}

// TEST FORMATTEDFILE, WITHOUT HEADER

SEQAN_DEFINE_TEST(test_blast_write_formatted_file_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    test_blast_write_tabular_impl(3, 0, context);
}

// TEST FORMATTEDFILE, WITH HEADER

SEQAN_DEFINE_TEST(test_blast_write_formatted_file_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    test_blast_write_tabular_impl(3, 0, context);
}

// TEST FORMATTEDFILE, WITHOUT HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_formatted_file_customfields_tabular)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(3, 1, context);
}

// TEST FORMATTEDFILE, WITH HEADER, CUSTOM FIELDS

SEQAN_DEFINE_TEST(test_blast_write_formatted_file_customfields_tabular_with_header)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::HEADER> context;
    append(context.fields, BlastMatchField<>::Enum::SCORE);
    append(context.fields, BlastMatchField<>::Enum::FRAMES);
    test_blast_write_tabular_impl(3, 1, context);
}

// PAIRWISE FORMAT

template <typename TString, typename TContext>
void test_blast_compare_report_outputs(TString const & output, TContext const & context)
{
    std::string compString;
    if (!context.legacyFormat)
        compString.append("BLASTP 2.2.26+ ");
    else
        compString.append("BLASTP 2.2.26 ");
    compString.append("[I/O Module of SeqAn-");
    compString.append(std::to_string(SEQAN_VERSION_MAJOR));
    compString.append(".");
    compString.append(std::to_string(SEQAN_VERSION_MINOR));
    compString.append(".");
    compString.append(std::to_string(SEQAN_VERSION_PATCH));
    compString.append(", http://www.seqan.de]\n"
    "\n\n"
    "Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,\n"
    "Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),\n"
    "\"Gapped BLAST and PSI-BLAST: a new generation of protein database search\n"
    "programs\",  Nucleic Acids Res. 25:3389-3402.\n"
    "\n\n\n"
    "Reference for SeqAn: Doering, A., D. Weese, T. Rausch, K. Reinert (2008): SeqAn --\n"
    "An efficient, generic C++ library for sequence analysis. BMC Bioinformatics,\n"
    "9(1), 11. BioMed Central Ltd. doi:10.1186/1471-2105-9-11\n"
    "\n"
    "\n"
    "\n"
    "Database: The Foo Database\n"
    "           2 sequences; 400 total letters\n"
    "\n"
    "\n"
    "Query= Query_Numero_Uno with args\n"
    "\n"
    "Length=12\n"
//     "\n"
//     "\n"
    "                                                                   Score     E\n"
    "Sequences producing significant alignments:                       (Bits)  Value\n"
    "\n"
    "Subject_Numero_Uno                                                   23  0.0005\n"
    "Subject_Numero_Dos                                                   22  0.0009\n"
    "\n"
    "ALIGNMENTS\n"
    "> Subject_Numero_Uno\n"
    "Length=200\n"
    "\n"
    " Score =  23.1 bits (48), Expect =  0.0005\n"
    " Identities = 10/14 (71%), Positives = 12/14 (86%), Gaps = 2/14 (14%)\n"
    "\n"
    "Query  1   VAYAQ--PRKLCYP  12\n"
    "           VAYAQ  PR+LC+P\n"
    "Sbjct  79  VAYAQTKPRRLCFP  92\n"
    "\n"
    "\n"
    "> Subject_Numero_Dos\n"
    "Length=200\n"
    "\n"
    " Score =  22.3 bits (46), Expect =  0.0009\n"
    " Identities = 8/8 (100%), Positives = 8/8 (100%), Gaps = 0/8 (0%)\n"
    "\n"
    "Query  3    YAQPRKLC  10 \n"
    "            YAQPRKLC\n"
    "Sbjct  157  YAQPRKLC  164\n"
    "\n"
    "\n"
    "\n"
    "Lambda     K      H\n"
    "   0.267   0.0410   0.1400\n"
    "\n"
    "Gapped\n"
    "Lambda     K      H\n"
    "   0.267   0.0410   0.1400\n"
    "\n"
    "Effective search space used: 4800\n"
    "\n"
    "\n"
    "Query= Query_Numero_Dos with args\n"
    "\n"
    "Length=10\n"
    "\n"
    "\n"
    "***** No hits found *****\n"
    "\n"
    "\n"
    "\n"
    "Lambda     K      H\n"
    "   0.267   0.0410   0.1400\n"
    "\n"
    "Gapped\n"
    "Lambda     K      H\n"
    "   0.267   0.0410   0.1400\n"
    "\n"
    "Effective search space used: 4000\n"
    "\n"
    "\n"
    "Query= Query_Numero_Tres with args\n"
    "\n"
    "Length=8\n"
    "                                                                   Score     E\n"
    "Sequences producing significant alignments:                       (Bits)  Value\n"
    "\n"
    "Subject_Numero_Dos                                                   18  0.007\n"
    "Subject_Numero_Uno                                                   16  0.03\n"
    "\n"
    "ALIGNMENTS\n"
    "> Subject_Numero_Dos\n"
    "Length=200\n"
    "\n"
    " Score =  18.9 bits (37), Expect =  0.007\n"
    " Identities = 8/8 (100%), Positives = 8/8 (100%), Gaps = 0/8 (0%)\n"
    "\n"
    "Query  1   AVITSFTQ  8 \n"
    "           AVITSFTQ\n"
    "Sbjct  10  AVITSFTQ  17\n"
    "\n"
    "\n"
    "> Subject_Numero_Uno\n"
    "Length=200\n"
    "\n"
    " Score =  16.9 bits (32), Expect =  0.03\n"
    " Identities = 7/8 (88%), Positives = 8/8 (100%), Gaps = 0/8 (0%)\n"
    "\n"
    "Query  1    AVITSFTQ  8  \n"
    "            AVITS+TQ\n"
    "Sbjct  184  AVITSWTQ  191\n"
    "\n"
    "\n"
    "\n"
    "Lambda     K      H\n"
    "   0.267   0.0410   0.1400\n"
    "\n"
    "Gapped\n"
    "Lambda     K      H\n"
    "   0.267   0.0410   0.1400\n"
    "\n"
    "Effective search space used: 3200\n"
    "\n"
    "\n"
    "  Database: The Foo Database\n"
    "  Number of letters in database: 400\n"
    "  Number of sequences in database:  2\n"
    "\n"
    "\n"
    "\n"
    "Matrix: BLOSUM62\n"
    "Gap Penalties: Existence: 11, Extension: 1\n\n");

    if (output != compString)
    {
        for (uint32_t i = 0; i < length(output); ++i)
        {
            if (output[i] != compString[i])
            {
                std::cout << output.substr(0,i) << "\n";
                std::cout << "CONT: \"" << output[i] << "\"\n";
                std::cout << "COMP: \"" << compString[i] << "\"\n";
                break;
            }
        }
    }
    SEQAN_ASSERT_EQ(output, compString);

//     std::cout << "<span style=\"font-size:80%\"><table>\n"
//                  "<tr><th>index</th>"
//                  "<th>Enum</th>"
//                  "<th>optionLabels</th>"
//                  "<th>columnLabels</th>"
//                  "<th>descriptions</th>"
//                  "<th>implemented</th></tr>\n";
//     for (int i = 0; i < 45; ++i)
//     {
//         std::cout << "<tr><td>"
//                   << i
//                   << "</td><td>"
//                   // enum have to inserted manually
//                   << "</td><td>"
//                   << BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::optionLabels[i]
//                   << "</td><td>"
//                   << BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::columnLabels[i]
//                   << "</td><td>"
//                   << BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::descriptions[i]
//                   << "</td><td>"
//                   << (BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::implemented[i]
//                     ? "&#9745;"
//                     : "&#9744;")
//                   << "</td></tr>\n";
//     }
//     std::cout << "</table></span>\n";
}

SEQAN_DEFINE_TEST(test_blast_write_pairwise)
{
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer,
                         std::ios_base::in |
                         std::ios_base::out |
                         std::ios_base::binary |
                         std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;

    test_blast_write_record_match(fstream, 0, 0, context, BlastReport());

    std::string contents;
    resize(contents, fstream.tellg());
    fstream.seekg(0, std::ios::beg);
    fstream.read(&contents[0], contents.size());
    fstream.close();

    test_blast_compare_report_outputs(contents, context);
}

SEQAN_DEFINE_TEST(test_blast_write_pairwise_formatted_file)
{
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer,
                         std::ios_base::in |
                         std::ios_base::out |
                         std::ios_base::binary |
                         std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;

    test_blast_write_record_match(fstream, 0, 3, context, BlastReport());

    std::string contents;
    resize(contents, fstream.tellg());
    fstream.seekg(0, std::ios::beg);
    fstream.read(&contents[0], contents.size());
    fstream.close();

    test_blast_compare_report_outputs(contents, context);
}

#endif  // SEQAN_TESTS_TEST_BLAST_OUTPUT_H_
