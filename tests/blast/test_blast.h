// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for basic/translation.h
// ==========================================================================

#ifndef SEQAN_EXTRAS_TESTS_BASIC_TEST_BLAST_H_
#define SEQAN_EXTRAS_TESTS_BASIC_TEST_BLAST_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/blast.h>

using namespace seqan;

// ------------------------- STATISTICS ---------------------------------------

SEQAN_DEFINE_TEST(test_blast_scoring_scheme_conversion)
{
    Blosum62 scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);

    blastScoringScheme2seqanScoringScheme(scheme);

    SEQAN_ASSERT_EQ(scoreGapOpen(scheme), -12);
    SEQAN_ASSERT_EQ(scoreGapExtend(scheme), -1);

    seqanScoringScheme2blastScoringScheme(scheme);

    SEQAN_ASSERT_EQ(scoreGapOpen(scheme), -11);
    SEQAN_ASSERT_EQ(scoreGapExtend(scheme), -1);
}

SEQAN_DEFINE_TEST(test_blast_scoring_adapter)
{
    // Blosum62 and some general stuff
    {
        typedef Blosum62 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -1);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.267);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.041);
        SEQAN_ASSERT_EQ(getH(adapter),      0.140);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  1.900);
        SEQAN_ASSERT_EQ(getBeta(adapter), -30.000);

        setScoreGapOpen(scheme, -1);
        setScoreGapExtend(scheme, -11); // invalid parameters
        // TEST assignScoreScheme
        bool res = assignScoreScheme(adapter, scheme);
        SEQAN_ASSERT(!res);

        // TEST getScoreScheme
        Blosum62 const & schemeInt = getScoreScheme(adapter);
        SEQAN_ASSERT_EQ(scoreGapOpen(schemeInt), -1);
        SEQAN_ASSERT_EQ(scoreGapExtend(schemeInt), -11);

        // reset to valid
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -1);
        blastScoringScheme2seqanScoringScheme(scheme);
        assignScoreScheme(adapter, scheme);
        // now check values again
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.267);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.041);
        SEQAN_ASSERT_EQ(getH(adapter),      0.140);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  1.900);
        SEQAN_ASSERT_EQ(getBeta(adapter), -30.000);
    }

    // Blosum45
    {
        typedef Blosum45 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -3);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.190);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.031);
        SEQAN_ASSERT_EQ(getH(adapter),      0.095);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  2.000);
        SEQAN_ASSERT_EQ(getBeta(adapter), -38.000);
    }

    // Blosum80
    {
        typedef Blosum80 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -1);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.314);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.095);
        SEQAN_ASSERT_EQ(getH(adapter),      0.350);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  0.900);
        SEQAN_ASSERT_EQ(getBeta(adapter),  -9.000);
    }

    // Pam250
    {
        typedef Pam250 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -3);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.174);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.020);
        SEQAN_ASSERT_EQ(getH(adapter),      0.070);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  2.500);
        SEQAN_ASSERT_EQ(getBeta(adapter), -48.000);
    }

    // SimpleScore
    {
        typedef Score<int, Simple> TScheme;
        TScheme scheme;
        setScoreMatch(scheme, 2);
        setScoreMismatch(scheme, -3);
        setScoreGapOpen(scheme, -4);
        setScoreGapExtend(scheme, -2);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.610);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.350);
        SEQAN_ASSERT_EQ(getH(adapter),      0.680);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  0.900);
        SEQAN_ASSERT_EQ(getBeta(adapter),  -3.000);
    }
}

SEQAN_DEFINE_TEST(test_blast_blastmatch_stats_and_score)
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;

    TBlastMatch m;

    String<AminoAcid> src0 = "ARNDAYVBRNDCQFGCYVBQARNDCQEGEG";
    String<AminoAcid> src1 = "ARNAYVBRNDCCYCYVBQARNQEGEG";

    resize(rows(m.align), 2);
    assignSource(row(m.align, 0), src0);
    assignSource(row(m.align, 1), src1);

    typedef Blosum62 TScheme;
    TScheme scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);
    blastScoringScheme2seqanScoringScheme(scheme);

    int score = globalAlignment(m.align, scheme);
//     ARNDAYVBRNDCQFGCYVBQARNDCQEGEG
//     ||| ||||||||   ||||||||  |||||
//     ARN-AYVBRNDCCY-CYVBQARN--QEGEG

    SEQAN_ASSERT_EQ(score, 94);

    calcStatsAndScore(m, scheme);

    SEQAN_ASSERT_EQ(m.score, score);
    SEQAN_ASSERT_EQ(m.aliLength, 30u);
    SEQAN_ASSERT_EQ(m.identities, 24u);
    SEQAN_ASSERT_EQ(m.positives, 25u);
    SEQAN_ASSERT_EQ(m.mismatches, 2u);
    SEQAN_ASSERT_EQ(m.gaps, 4u);
    SEQAN_ASSERT_EQ(m.gapOpenings, 3u);
}

SEQAN_DEFINE_TEST(test_blast_blastmatch_bit_score_e_value)
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;
    typedef Blosum62 TScheme;

    TBlastMatch m;

    String<AminoAcid> src0 =
    "SSITEEKHIPHKEQDKDAEFLSKEALKTHMTENVLQMDRRAVQDPSTSFLQLLKAKGLLG"
    "LPDYEVNLADVNSPGFRKVAYAQTKPRRLCFPNGGTRRGSFIMDTAVVVMVSLRYVNIGK"
    "VIFPGATDVSEGEDEFWAGLPQAYGCLATEFLCIHIAIYSWIHVQSSRYDDMNASVIRAK"
    "LNLPGGLTLIQARGNEKETI";
    String<AminoAcid> src1 = "VAYAQPRKLCYP";

    resize(rows(m.align), 2);
    assignSource(row(m.align, 0), src0);
    assignSource(row(m.align, 1), src1);

    TScheme scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);
    blastScoringScheme2seqanScoringScheme(scheme);

    int score = localAlignment(m.align, scheme);
//         VAYAQTKPRRLCFP
//         |||||  || || |
//         VAYAQ--PRKLCYP

    SEQAN_ASSERT_EQ(score, 48);

    calcStatsAndScore(m, scheme);

    SEQAN_ASSERT_EQ(m.score, score);
    SEQAN_ASSERT_EQ(m.aliLength, 14u);
    SEQAN_ASSERT_EQ(m.identities, 10u);
    SEQAN_ASSERT_EQ(m.positives, 12u);
    SEQAN_ASSERT_EQ(m.mismatches, 2u);
    SEQAN_ASSERT_EQ(m.gaps, 2u);
    SEQAN_ASSERT_EQ(m.gapOpenings, 1u);

    // same as previous test until here (just local)

    BlastScoringAdapter<TScheme> adapter(scheme);
    SEQAN_ASSERT(isValid(adapter));

    calcBitScoreAndEValue(m.bitScore, m.eValue, m.score, length(src0),
                          length(src1), adapter);

    double epsilon = 1e-4;
    SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 23.0978), epsilon);
    epsilon = 1e-8; // more important on evalue
    SEQAN_ASSERT_LEQ(std::abs(m.eValue - 0.000267348), epsilon);
}

// ------------------------- OUTPUT ---------------------------------------

template <typename TFile, typename TDbSpecs, typename TRecords, typename TFormat>
inline void
_writeCustomColumnsImpl(TFile & file,
                       TRecords const & records,
                       TDbSpecs const & dbSpecs,
                       TFormat const & /**/)
{
    for (auto const & r : records)
    {
        writeHeader(file,
                    r.qId,
                    dbSpecs.dbName,
                    r.matches.size(),
                    TFormat(),
                    "Query id",
                    "Subject id",
                    "alignment length",
                    "mismatches",
                    "gaps",
                    "e-value",
                    "bit score");
        for (auto const & m : r.matches)
        {
             writeMatch(file,
                        TFormat(),
                        m.qId,
                        m.sId,
                        m.aliLength,
                        m.mismatches,
                        m.gaps,
                        m.eValue,
                        m.bitScore);
        }
    }
}

template <typename TFile, typename TDbSpecs, typename TRecords, BlastFormatProgram p, BlastFormatGeneration g>
inline void
_writeCustomColumns(TFile & file,
                   TRecords const & records,
                   TDbSpecs const & dbSpecs,
                   BlastFormat<BlastFormatFile::TABULAR, p, g> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeCustomColumnsImpl(file, records, dbSpecs, TFormat());
}

template <typename TFile, typename TDbSpecs, typename TRecords, BlastFormatProgram p, BlastFormatGeneration g>
inline void
_writeCustomColumns(TFile & file,
                   TRecords const & records,
                   TDbSpecs const & dbSpecs,
                   BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    _writeCustomColumnsImpl(file, records, dbSpecs, TFormat());
}

template <typename TFile, typename TDbSpecs, typename TRecords, typename TFormat>
inline void
_writeCustomColumns(TFile & /**/,
                   TRecords const & /**/,
                   TDbSpecs const & /**/,
                   TFormat const & /**/)
{
    SEQAN_ASSERT(false);
}

template <typename TFile, typename TDbSpecs, typename TRecords, typename TFormat>
inline void
_writeCustomFieldsImpl(TFile & file,
                       TRecords const & records,
                       TDbSpecs const & dbSpecs,
                       TFormat const & /**/)
{
    typedef BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum TField;
    String<TField> fields;
    append(fields, TField::STD);
    append(fields, TField::SCORE);
    append(fields, TField::FRAMES);

    for (int q = 0; q < 2; ++q)
        writeRecord(file, records[q], dbSpecs, fields, TFormat());
}

template <typename TFile, typename TDbSpecs, typename TRecords, BlastFormatProgram p>
inline void
_writeCustomFields(TFile & file,
                   TRecords const & records,
                   TDbSpecs const & dbSpecs,
                   BlastFormat<BlastFormatFile::TABULAR,
                               p,
                               BlastFormatGeneration::BLAST_PLUS> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    _writeCustomFieldsImpl(file, records, dbSpecs, TFormat());
}

template <typename TFile, typename TDbSpecs, typename TRecords, BlastFormatProgram p, BlastFormatGeneration g>
inline void
_writeCustomFields(TFile & file,
                   TRecords const & records,
                   TDbSpecs const & dbSpecs,
                   BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                               p,
                               BlastFormatGeneration::BLAST_PLUS> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    _writeCustomFieldsImpl(file, records, dbSpecs, TFormat());
}

template <typename TFile, typename TDbSpecs, typename TRecords, typename TFormat>
inline void
_writeCustomFields(TFile & /**/,
                   TRecords const & /**/,
                   TDbSpecs const & /**/,
                   TFormat const & /**/)
{
    SEQAN_ASSERT(false);
}

template <typename TFile,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
test_blast_write_record_match(TFile & file,
                              int const customFields,
                              BlastFormat<f, p, g> const & /**/)
{
    typedef BlastFormat<f, p, g> TFormat;
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;
    typedef BlastRecord<CharString, CharString, uint32_t, TAlign>
            TBlastRecord;
    typedef BlastDbSpecs<CharString> TDbSpecs;
    typedef Blosum62 TScheme;

    StringSet<String<AminoAcid>, Owner<ConcatDirect<>>> queries;
    StringSet<CharString, Owner<ConcatDirect<>>> qIds;
    StringSet<String<AminoAcid>> subjects;
    StringSet<CharString> sIds;
    CharString dbName = "The Foo Database";

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
    appendValue(queries, "AVITSFTQ");

    appendValue(qIds, "Query_Numero_Uno");
    appendValue(qIds, "Query_Numero_Dos");


    TScheme scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);
    blastScoringScheme2seqanScoringScheme(scheme);
    BlastScoringAdapter<TScheme> adapter(scheme);
    SEQAN_ASSERT(isValid(adapter));

    String<TBlastRecord> records;
    resize(records, 2);

    TDbSpecs dbSpecs(dbName);
    dbSpecs.dbTotalLength = length(subjects[0]) + length(subjects[1]);
    dbSpecs.dbNumberOfSeqs = 2;

    for (int q = 0; q < 2; ++q)
    {
        records[q].qId = qIds[q];
        records[q].qLength = length(queries[q]);

        for (int s = 0; s < 2; ++s)
        {
            records[q].matches.emplace_back(qIds[q], sIds[s]);
            TBlastMatch & m = records[q].matches.back();
            resize(rows(m.align), 2);
            assignSource(row(m.align, 0), subjects[s]);
            assignSource(row(m.align, 1), queries[q]);

            localAlignment(m.align, scheme);

            m.sStart = beginPosition(row(m.align, 0));
            m.sEnd   = endPosition(row(m.align, 0));
            m.qStart = beginPosition(row(m.align, 1));
            m.qEnd   = endPosition(row(m.align, 1));

            m.qLength = length(queries[q]);

            m.qFrameShift = 1;
            m.sFrameShift = 1;

            calcStatsAndScore(m, scheme);

            calcBitScoreAndEValue(m, dbSpecs, adapter);
        }
    }

    writeTop(file, dbSpecs, TFormat());

    if (customFields == 0)
    {
        for (int q = 0; q < 2; ++q)
            writeRecord(file, records[q], dbSpecs, TFormat());
    } else if (customFields == 1) // only TABULAR, BLAST_PLUS
    {
        _writeCustomFields(file, records, dbSpecs, TFormat());
    } else // only TABULAR
    {
        _writeCustomColumns(file, records, dbSpecs, TFormat());
    }

    writeBottom(file, dbSpecs, adapter, TFormat());

}

template <BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
void test_blast_write_tabular_impl(BlastFormat<f, p, g> const & /**/,
                                   int customFields = 0)
{
    typedef BlastFormat<f,p,g> TFormat;
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    test_blast_write_record_match(fstream, customFields, TFormat());

    std::string contents;
    resize(contents, fstream.tellg());
    fstream.seekg(0, std::ios::beg);
    fstream.read(&contents[0], contents.size());
    fstream.close();

    std::string compString;
    if (f == BlastFormatFile::TABULAR_WITH_HEADER)
    {
        std::string header = "# BLASTP I/O Module of SeqAn-";
        header.append(std::to_string(SEQAN_VERSION_MAJOR));
        header.append(".");
        header.append(std::to_string(SEQAN_VERSION_MINOR));
        header.append(".");
        header.append(std::to_string(SEQAN_VERSION_PATCH));
        header.append(" (http://www.seqan.de)\n"
        "# Query: Query_Numero_Uno\n"
        "# Database: The Foo Database\n"
        "# Fields: ");
        if (customFields == 0)
        {
            header.append(BlastMatchField<g>::columnLabels[0]);
        } else if (customFields == 1)
        {
            header.append(BlastMatchField<g>::columnLabels[0]);
            header.append(", score, query/sbjct frames");
        } else
        {
            header.append("Query id, Subject id, alignment length, mismatches, gaps, e-value, bit score");
        }
        header.append("\n");
        if (g == BlastFormatGeneration::BLAST_PLUS)
            header.append("# 2 hits found\n");
        compString.append(header);
    }

    if (customFields == 0)
        compString.append(
        "Query_Numero_Uno\tSubject_Numero_Uno\t71.43\t14\t2\t1\t1\t12\t79\t92\t5e-04\t23.1\n"
        "Query_Numero_Uno\tSubject_Numero_Dos\t100.00\t8\t0\t0\t3\t10\t157\t164\t0.001\t22.3\n");
    else if (customFields == 1)
        compString.append(
        "Query_Numero_Uno\tSubject_Numero_Uno\t71.43\t14\t2\t1\t1\t12\t79\t92\t5e-04\t23.1\t48\t+1/+1\n"
        "Query_Numero_Uno\tSubject_Numero_Dos\t100.00\t8\t0\t0\t3\t10\t157\t164\t0.001\t22.3\t46\t+1/+1\n");
    else
        compString.append(
        "Query_Numero_Uno\tSubject_Numero_Uno\t14\t2\t2\t0.000534696\t23.0978\n"
        "Query_Numero_Uno\tSubject_Numero_Dos\t8\t0\t0\t0.000912053\t22.3274\n");


    if (f == BlastFormatFile::TABULAR_WITH_HEADER)
    {
        std::string header = "# BLASTP I/O Module of SeqAn-";
        header.append(std::to_string(SEQAN_VERSION_MAJOR));
        header.append(".");
        header.append(std::to_string(SEQAN_VERSION_MINOR));
        header.append(".");
        header.append(std::to_string(SEQAN_VERSION_PATCH));
        header.append(" (http://www.seqan.de)\n"
        "# Query: Query_Numero_Dos\n" // <--- the only line different to above
        "# Database: The Foo Database\n"
        "# Fields: ");
        if (customFields == 0)
        {
            header.append(BlastMatchField<g>::columnLabels[0]);
        } else if (customFields == 1)
        {
            header.append(BlastMatchField<g>::columnLabels[0]);
            header.append(", score, query/sbjct frames");
        } else
        {
            header.append("Query id, Subject id, alignment length, mismatches, gaps, e-value, bit score");
        }
        header.append("\n");
        if (g == BlastFormatGeneration::BLAST_PLUS)
            header.append("# 2 hits found\n");
        compString.append(header);

    }
    if (customFields == 0)
        compString.append(
        "Query_Numero_Dos\tSubject_Numero_Uno\t87.50\t8\t1\t0\t1\t8\t184\t191\t0.026\t16.9\n"
        "Query_Numero_Dos\tSubject_Numero_Dos\t100.00\t8\t0\t0\t1\t8\t10\t17\t0.007\t18.9\n");
    else if (customFields == 1)
        compString.append(
        "Query_Numero_Dos\tSubject_Numero_Uno\t87.50\t8\t1\t0\t1\t8\t184\t191\t0.026\t16.9\t32\t+1/+1\n"
        "Query_Numero_Dos\tSubject_Numero_Dos\t100.00\t8\t0\t0\t1\t8\t10\t17\t0.007\t18.9\t36\t+1/+1\n");
    else
        compString.append(
        "Query_Numero_Dos\tSubject_Numero_Uno\t8\t1\t0\t0.0255459\t16.9346\n"
        "Query_Numero_Dos\tSubject_Numero_Dos\t8\t0\t0\t0.00672262\t18.8606\n");


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

// WITHOUT HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_tabular)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    test_blast_write_tabular_impl(TFormat());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_legacy)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;
    test_blast_write_tabular_impl(TFormat());
}

// WITH HEADER, DEFAULTS FIELDS

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_header)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    test_blast_write_tabular_impl(TFormat());
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_header_legacy)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;
    test_blast_write_tabular_impl(TFormat());
}

// WITHOUT HEADER, CUSTOM FIELDS (through match object)

SEQAN_DEFINE_TEST(test_blast_write_tabular_customfields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;
    test_blast_write_tabular_impl(TFormat(), 1);
}

// WITH HEADER, CUSTOM FIELDS (through match object)

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_header_customfields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    test_blast_write_tabular_impl(TFormat(), 1);
}

// WITHOUT HEADER, CUSTOM COLUMNS (arbitrary column support)

SEQAN_DEFINE_TEST(test_blast_write_tabular_customcolumns)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    test_blast_write_tabular_impl(TFormat(), 2);
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_customcolumns_legacy)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;
    test_blast_write_tabular_impl(TFormat(), 2);
}

// WITHOUT HEADER, CUSTOM COLUMNS (arbitrary column support)

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_header_customcolumns)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    test_blast_write_tabular_impl(TFormat(), 2);
}

SEQAN_DEFINE_TEST(test_blast_write_tabular_with_header_customcolumns_legacy)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;
    test_blast_write_tabular_impl(TFormat(), 2);
}

// PAIRWISE FORMAT

SEQAN_DEFINE_TEST(test_blast_write_pairwise)
{
    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                        BlastFormatProgram::BLASTP,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    test_blast_write_record_match(fstream, false, TFormat());

    std::string contents;
    resize(contents, fstream.tellg());
    fstream.seekg(0, std::ios::beg);
    fstream.read(&contents[0], contents.size());
    fstream.close();

    std::string compString = "BLASTP I/O Module of SeqAn-";
    compString.append(std::to_string(SEQAN_VERSION_MAJOR));
    compString.append(".");
    compString.append(std::to_string(SEQAN_VERSION_MINOR));
    compString.append(".");
    compString.append(std::to_string(SEQAN_VERSION_PATCH));
    compString.append(" (http://www.seqan.de)\n"
    "\n"
    "Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,\n"
    "Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),\n"
    "\"Gapped BLAST and PSI-BLAST: a new generation of protein database search\n"
    "programs\",  Nucleic Acids Res. 25:3389-3402.\n"
    "\n"
    "Reference for SeqAn: Döring, A., D. Weese, T. Rausch, K. Reinert (2008): SeqAn --\n"
    "An efficient, generic C++ library for sequence analysis. BMC Bioinformatics,\n"
    "9(1), 11. BioMed Central Ltd. doi:10.1186/1471-2105-9-11\n"
    "\n"
    "\n"
    "\n"
    "Database: The Foo Database\n"
    "           2 sequences; 400 total letters\n"
    "\n"
    "\n"
    "Query= Query_Numero_Uno\n"
    "\n"
    "Length=12\n"
    "\n"
    "\n"
    "                                                                   Score     E\n"
    "Sequences producing significant alignments:                       (Bits)  Value\n"
    "\n"
    "Subject_Numero_Uno                                                   23  0.0005\n"
    "Subject_Numero_Dos                                                   22  0.0009\n"
    "\n"
    "ALIGNMENTS\n"
    "> Subject_Numero_Uno\n"
    "Length=0\n"
    "\n"
    " Score =  23.1 bits (48), Expect =  0.0005\n"
    " Identities = 10/14 (71%), Positives = 12/14 (86%), Gaps = 2/14 (14%)\n"
    "\n"
    "Query  1   VAYAQTKPRRLCFP  14\n"
    "           VAYAQ  PR+LC+P\n"
    "Sbjct  79  VAYAQ--PRKLCYP  90\n"
    "\n"
    "> Subject_Numero_Dos\n"
    "Length=0\n"
    "\n"
    " Score =  22.3 bits (46), Expect =  0.0009\n"
    " Identities = 8/8 (100%), Positives = 8/8 (100%), Gaps = 0/8 (0%)\n"
    "\n"
    "Query  3    YAQPRKLC  10 \n"
    "            YAQPRKLC\n"
    "Sbjct  157  YAQPRKLC  164\n"
    "\n"
    "\n"
    "Query= Query_Numero_Dos\n"
    "\n"
    "Length=8\n"
    "\n"
    "\n"
    "                                                                   Score     E\n"
    "Sequences producing significant alignments:                       (Bits)  Value\n"
    "\n"
    "Subject_Numero_Uno                                                   16  0.03\n"
    "Subject_Numero_Dos                                                   18  0.007\n"
    "\n"
    "ALIGNMENTS\n"
    "> Subject_Numero_Uno\n"
    "Length=0\n"
    "\n"
    " Score =  16.9 bits (32), Expect =  0.03\n"
    " Identities = 7/8 (88%), Positives = 8/8 (100%), Gaps = 0/8 (0%)\n"
    "\n"
    "Query  1    AVITSWTQ  8  \n"
    "            AVITS+TQ\n"
    "Sbjct  184  AVITSFTQ  191\n"
    "\n"
    "> Subject_Numero_Dos\n"
    "Length=0\n"
    "\n"
    " Score =  18.9 bits (37), Expect =  0.007\n"
    " Identities = 8/8 (100%), Positives = 8/8 (100%), Gaps = 0/8 (0%)\n"
    "\n"
    "Query  1   AVITSFTQ  8 \n"
    "           AVITSFTQ\n"
    "Sbjct  10  AVITSFTQ  17\n"
    "\n"
    "\n"
    "Matrix:BLOSUM62\n"
    "Gap Penalties: Existence: -11, Extension: -1\n\n");

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

#endif  // SEQAN_EXTRAS_TESTS_BASIC_TEST_BLAST_H_
