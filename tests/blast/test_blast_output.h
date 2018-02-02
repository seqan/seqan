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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Tests for the blast module
// ==========================================================================

#ifndef SEQAN_TESTS_TEST_BLAST_OUTPUT_H_
#define SEQAN_TESTS_TEST_BLAST_OUTPUT_H_

using namespace seqan;

template <typename TFile, typename TScore, typename TRecords, BlastProgram p, BlastTabularSpec h>
inline void
_testBlastOutputWriteFile(TFile & file,
                          BlastIOContext<TScore, p, h> & context,
                          TRecords const & records,
                          BlastReport const & /**/)
{
    BlastReportFileOut<BlastIOContext<TScore, p, h>> out(context, file, BlastReport());
    writeHeader(out);
    for (auto const & r : records)
        writeRecord(out, r);
    writeFooter(out);
}

template <typename TFile, typename TScore, typename TRecords, BlastProgram p, BlastTabularSpec h>
inline void
_testBlastOutputWriteFile(TFile & file,
                          BlastIOContext<TScore, p, h> & context,
                          TRecords const & records,
                          BlastTabular const & /**/)
{
    BlastTabularFileOut<BlastIOContext<TScore, p, h>> out(context, file, BlastTabular());
    writeHeader(out); // noop for TABULARs
    for (auto const & r : records)
        writeRecord(out, r);
    writeFooter(out);
}

template <typename TFile, typename TScore, typename TRecords, BlastProgram p, BlastTabularSpec h>
inline void
_testBlastOutputWriteFile(TFile & file,
                          BlastIOContext<TScore, p, h> &,
                          TRecords const & records,
                          BlastTabularLL const & /**/)
{
    for (auto const & r : records)
    {
        for (auto const & m : r.matches)
        {
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
}

template <typename TFile, typename TScore, typename TFormat, BlastProgram p, BlastTabularSpec h>
inline void
_testBlastOutputGenerateContent(TFile & file,
                                BlastIOContext<TScore, p, h> & context,
                                TFormat const &)
{
    typedef StringSet<String<AminoAcid>, Owner<ConcatDirect<>>> TQueries;
    typedef StringSet<CharString, Owner<ConcatDirect<>>> TQIds;
    typedef StringSet<String<AminoAcid>> TSubjects;
    typedef StringSet<CharString> TSIds;

    typedef Gaps<String<AminoAcid>, ArrayGaps> TGaps;
    typedef BlastMatch<TGaps, TGaps, uint32_t, typename Value<TQIds>::Type, typename Value<TSIds>::Type> TBlastMatch;
    typedef BlastRecord<TBlastMatch> TBlastRecord;

    TQueries queries;
    TQIds qIds;
    TSubjects subjects;
    TSIds sIds;

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

            assignSource(m.alignRow0, queries[q]);
            assignSource(m.alignRow1, subjects[s]);

            localAlignment(m.alignRow0, m.alignRow1, seqanScheme(context.scoringScheme));

            m.qStart = beginPosition(m.alignRow0);
            m.qEnd   = endPosition(m.alignRow0);
            m.sStart = beginPosition(m.alignRow1);
            m.sEnd   = endPosition(m.alignRow1);

            m.qLength = length(queries[q]);
            m.sLength = length(subjects[s]);

            m.qFrameShift = 1;
            m.sFrameShift = 1;

            computeAlignmentStats(m, context);
            computeBitScore(m, context);
            computeEValue(m, records[q].qLength, context);
        }
    }

    // sort by bit-score (range based for-loop broken on ICC)
    for (auto it = begin(records, Standard()), itEnd = end(records, Standard());
         it != itEnd;
         ++it)
         it->matches.sort();

    _testBlastOutputWriteFile(file, context, records, TFormat());
}

template <typename TString, typename TContext>
void _testBlastOutputCheckFileWritten(TString const & contents, TContext const & context, BlastTabular const &)
{
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
        fieldsLine.append(BlastMatchField<>::columnLabels[0]);
        if ((length(context.fields) != 1) ||
            ((length(context.fields) == 1) && (context.fields[0] != BlastMatchField<>::Enum::STD)))
            fieldsLine.append(", score, query/sbjct frames");
    }
    fieldsLine.append("\n");

    std::string compString;
    if (context.tabularSpec == BlastTabularSpec::COMMENTS)
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

    if ((length(context.fields) == 1) && (context.fields[0] == BlastMatchField<>::Enum::STD))
    {
        compString.append("Query_Numero_Uno\tSubject_Numero_Uno\t71.43\t14\t");
        compString.append(std::to_string(mismatches));
        compString.append(                                                    "\t1\t1\t12\t79\t92\t5e-04\t23.1\n"
                          "Query_Numero_Uno\tSubject_Numero_Dos\t100.00\t8\t0\t0\t3\t10\t157\t164\t0.001\t22.3\n");
    }
    else
        compString.append(
        "Query_Numero_Uno\tSubject_Numero_Uno\t71.43\t14\t2\t1\t1\t12\t79\t92\t5e-04\t23.1\t48\t0/0\n"
        "Query_Numero_Uno\tSubject_Numero_Dos\t100.00\t8\t0\t0\t3\t10\t157\t164\t0.001\t22.3\t46\t0/0\n");


    if (context.tabularSpec == BlastTabularSpec::COMMENTS)
    {
        // comments for empty record
        compString.append(versionLine);
        compString.append("# Query: Query_Numero_Dos with args\n"
                          "# Database: The Foo Database\n");
        if (context.legacyFormat)
            compString.append(fieldsLine);
        if (!context.legacyFormat)
            compString.append("# 0 hits found\n");

        // comments for
        compString.append(versionLine);
        compString.append("# Query: Query_Numero_Tres with args\n"
                          "# Database: The Foo Database\n");
        compString.append(fieldsLine);
        if (!context.legacyFormat)
            compString.append("# 2 hits found\n");
    }

    if ((length(context.fields) == 1) && (context.fields[0] == BlastMatchField<>::Enum::STD))
        compString.append("Query_Numero_Tres\tSubject_Numero_Dos\t100.00\t8\t0\t0\t1\t8\t10\t17\t0.007\t18.9\n"
                          "Query_Numero_Tres\tSubject_Numero_Uno\t87.50\t8\t1\t0\t1\t8\t184\t191\t0.026\t16.9\n");
    else
        compString.append("Query_Numero_Tres\tSubject_Numero_Dos\t100.00\t8\t0\t0\t1\t8\t10\t17\t0.007\t18.9\t37\t0/0\n"
                          "Query_Numero_Tres\tSubject_Numero_Uno\t87.50\t8\t1\t0\t1\t8\t184\t191\t0.026\t16.9\t32\t0/0\n");

    if ((context.tabularSpec == BlastTabularSpec::COMMENTS) && (!context.legacyFormat))
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

template <typename TString, typename TContext>
void _testBlastOutputCheckFileWritten(TString const & contents, TContext const &, BlastTabularLL const &)
{
    std::string compString;
    compString.append(
        "Query_Numero_Uno with args\tSubject_Numero_Uno\t14\t2\t2\t0.000534696\t23.0978\n"
        "Query_Numero_Uno with args\tSubject_Numero_Dos\t8\t0\t0\t0.000912053\t22.3274\n"
        "Query_Numero_Tres with args\tSubject_Numero_Dos\t8\t0\t0\t0.00672262\t18.8606\n"
        "Query_Numero_Tres with args\tSubject_Numero_Uno\t8\t1\t0\t0.0255459\t16.9346\n");
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

// PAIRWISE FORMAT
template <typename TString, typename TContext>
void _testBlastOutputCheckFileWritten(TString const & output, TContext const & context, BlastReport const &)
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
}

template <typename TScore, typename TFormat, BlastProgram p, BlastTabularSpec h>
void _testBlastOutput(BlastIOContext<TScore, p, h> & context, TFormat const & /**/)
{
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer,
                         std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    _testBlastOutputGenerateContent(fstream, context, TFormat());

    std::string contents;
    resize(contents, fstream.tellg());
    fstream.seekg(0, std::ios::beg);
    fstream.read(&contents[0], contents.size());
    fstream.close();

    // Generate the output string for comparison (since this includes the SeqAn version we don't want it completely
    // hardcoded (and also to save lots of redundant code for the immense combinations of options)
    _testBlastOutputCheckFileWritten(contents, context, TFormat());
}

// lowlevel
SEQAN_DEFINE_TEST(test_blast_write_lowlevel)
{
    BlastIOContext<> context;
    _testBlastOutput(context, BlastTabularLL());
}

// tabular
SEQAN_DEFINE_TEST(test_blast_write_tabular_without_comments)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;
    context.tabularSpec = BlastTabularSpec::NO_COMMENTS;

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_without_comments_customfields)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;
    context.tabularSpec = BlastTabularSpec::NO_COMMENTS;

    appendValue(context.fields, BlastMatchField<>::Enum::SCORE);
    appendValue(context.fields, BlastMatchField<>::Enum::FRAMES);

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_without_comments_legacy)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;
    context.tabularSpec = BlastTabularSpec::NO_COMMENTS;
    context.legacyFormat = true;

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_without_comments_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_COMMENTS> context;

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_without_comments_customfields_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_COMMENTS> context;
    appendValue(context.fields, BlastMatchField<>::Enum::SCORE);
    appendValue(context.fields, BlastMatchField<>::Enum::FRAMES);

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_without_comments_legacy_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_COMMENTS> context;
    context.legacyFormat = true;

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_comments)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;
    context.tabularSpec = BlastTabularSpec::COMMENTS;

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_comments_customfields)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;
    context.tabularSpec = BlastTabularSpec::COMMENTS;

    appendValue(context.fields, BlastMatchField<>::Enum::SCORE);
    appendValue(context.fields, BlastMatchField<>::Enum::FRAMES);

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_comments_legacy)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;
    context.tabularSpec = BlastTabularSpec::COMMENTS;
    context.legacyFormat = true;

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_comments_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::COMMENTS> context;

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_comments_customfields_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::COMMENTS> context;
    appendValue(context.fields, BlastMatchField<>::Enum::SCORE);
    appendValue(context.fields, BlastMatchField<>::Enum::FRAMES);

    _testBlastOutput(context, BlastTabular());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_comments_legacy_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::COMMENTS> context;
    context.legacyFormat = true;

    _testBlastOutput(context, BlastTabular());
}

// report
SEQAN_DEFINE_TEST(test_blast_write_report)
{
    BlastIOContext<> context;
    context.blastProgram = BlastProgram::BLASTP;

    _testBlastOutput(context, BlastReport());
}

SEQAN_DEFINE_TEST(test_blast_write_report_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP> context;

    _testBlastOutput(context, BlastReport());
}

SEQAN_DEFINE_TEST(test_blast_write_report_constexpr_dynmatrix)
{
    SelectableAminoAcidMatrix sel;
    SEQAN_ASSERT(getScoreMatrixId(sel) != AminoAcidScoreMatrixID::BLOSUM62);
    setScoreMatrixById(sel, AminoAcidScoreMatrixID::BLOSUM62);
    SEQAN_ASSERT(getScoreMatrixId(sel) == AminoAcidScoreMatrixID::BLOSUM62);

    BlastIOContext<SelectableAminoAcidMatrix, BlastProgram::BLASTP> context;
    context.scoringScheme._internalScheme = sel;

    _testBlastOutput(context, BlastReport());
}

#endif  // SEQAN_TESTS_TEST_BLAST_OUTPUT_H_
