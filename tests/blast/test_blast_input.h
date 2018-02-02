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

// Files that are being read by this implementation
#define PLUS_COMMENTS_DEFAULTS   "/tests/blast/defaultfields.m9"
#define LEGACY_COMMENTS_DEFAULTS "/tests/blast/defaultfields_legacy.m9"
#define NOCOMMENTS_DEFAULTS      "/tests/blast/defaultfields.m8"
#define PLUS_COMMENTS_CUSTOM     "/tests/blast/customfields.m9"
#define NOCOMMENTS_CUSTOM        "/tests/blast/customfields.m8"

using namespace seqan;

inline void
_test_blast_read_tabular_match_lowlevel(std::string const & path)
{
    std::string inPath = seqan::getAbsolutePath(path.c_str());

    std::ifstream ifstream(toCString(inPath),
                           std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(ifstream.is_open());

    auto it = directionIterator(ifstream, Input());

    // first line is comments
    SEQAN_ASSERT(!onMatch(it, BlastTabularLL()));

    // skip commentss and comment lines
    skipUntilMatch(it, BlastTabularLL());

    // now we should be onMatch
    SEQAN_ASSERT(onMatch(it, BlastTabularLL()));

    CharString  field1;
    std::string field2;
    double      field3;
    double      field4; // is actually int but cast to double must work, too
    unsigned    field5;
    unsigned    field6;
    unsigned    field7;
    unsigned    field8;
    unsigned    field9;
    unsigned    field10;
    double      field11;
    double      field12;

    readMatch(it, BlastTabularLL(), field1, field2, field3, field4, field5, field6,
               field7, field8, field9, field10, field11, field12);

    SEQAN_ASSERT_EQ(field1,          "SHAA004TF");
    SEQAN_ASSERT_EQ(field2,          "sp|P0A916|OMPW_SHIFL");
    SEQAN_ASSERT_EQ(field3,          50.43);
    SEQAN_ASSERT_EQ(field4,          115.0);
    if (path == LEGACY_COMMENTS_DEFAULTS)
        SEQAN_ASSERT_EQ(field5,      57u);
    else // legacy blast includes gaps in mismatches
        SEQAN_ASSERT_EQ(field5,      49u);
    SEQAN_ASSERT_EQ(field6,          2u);
    SEQAN_ASSERT_EQ(field7,          389u);
    SEQAN_ASSERT_EQ(field8,          733u);
    SEQAN_ASSERT_EQ(field9,          1u);
    SEQAN_ASSERT_EQ(field10,         107u);
    SEQAN_ASSERT_LEQ(std::abs(field11  - 1e-26), 1e-10);
    SEQAN_ASSERT_LEQ(std::abs(field12 - 108), 1e-3);

    // SEQAN_TRY reading less coluumns than are present
    readMatch(it, BlastTabularLL(), field1, field2, field3, field4, field5, field6,
               field7, field8, field9, field10);

    SEQAN_ASSERT_EQ(field1,          "SHAA004TF");
    SEQAN_ASSERT_EQ(field2,          "sp|P0A915|OMPW_ECOLI");
    SEQAN_ASSERT_EQ(field3,          50.43);
    SEQAN_ASSERT_EQ(field4,          115.0);
    if (path == LEGACY_COMMENTS_DEFAULTS) // legacy blast includes gaps in mismatches
        SEQAN_ASSERT_EQ(field5,      57u);
    else
        SEQAN_ASSERT_EQ(field5,      49u);
    SEQAN_ASSERT_EQ(field6,          2u);
    SEQAN_ASSERT_EQ(field7,          389u);
    SEQAN_ASSERT_EQ(field8,          733u);
    SEQAN_ASSERT_EQ(field9,          1u);
    SEQAN_ASSERT_EQ(field10,         107u);

    // read remaining matches
    while (onMatch(it, BlastTabularLL()))
        readMatch(it, BlastTabularLL(), field1); // only one field

    // go to last record with matches
    skipUntilMatch(it, BlastTabularLL());

    // goto second (last) match
    skipLine(it);

    // check if exceptions are properly thrown
    bool exceptThrown = false;
    SEQAN_TRY
    {
        // no strings here to take the strings
        readMatch(it, BlastTabularLL(), field3, field4, field5, field6, field7, field8, field9, field10, field11,
                  field12);
    } SEQAN_CATCH (BadLexicalCast const & e)
    {
        exceptThrown = true;
    }

    SEQAN_ASSERT(exceptThrown);
    // skip rest of line
    skipLine(it);

    ifstream.close();
}

SEQAN_DEFINE_TEST(test_blast_read_lowlevel)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match_lowlevel(PLUS_COMMENTS_DEFAULTS);
    // basic match reading should work between plus and legacy
    _test_blast_read_tabular_match_lowlevel(LEGACY_COMMENTS_DEFAULTS);
}

template <typename TContext>
void _testReadTabularWithoutComments(TContext &,
                                   std::string path,
                                   bool const defaults,
                                   bool const islegacy)
{
    typedef BlastMatchField<> TField;

    std::string inPath = seqan::getAbsolutePath(path.c_str());

    BlastRecord<> r;

    BlastTabularFileIn<TContext> fileIn(toCString(inPath));
    TContext & context = seqan::context(fileIn);

    /* Variables we will be comparing with */
    std::vector<typename TField::Enum> const fieldsDefault =
    {
        TField::Enum::STD
    };

    std::vector<typename TField::Enum> const fieldsCustom =
    {
        TField::Enum::Q_SEQ_ID,
        TField::Enum::S_SEQ_ID,
        TField::Enum::LENGTH,
        TField::Enum::N_IDENT,
        TField::Enum::MISMATCH,
        TField::Enum::POSITIVE,
        TField::Enum::GAPS,
        TField::Enum::Q_START,
        TField::Enum::Q_END,
        TField::Enum::S_START,
        TField::Enum::S_END,
        TField::Enum::FRAMES,
        TField::Enum::BIT_SCORE
    };

    /* Begin of TESTS */

    // read comments of file
    readHeader(fileIn);
    SEQAN_ASSERT(context.tabularSpec == BlastTabularSpec::NO_COMMENTS);

    // fieldsList as in-parameter
    if (defaults)
        context.fields = fieldsDefault;
    else
        context.fields = fieldsCustom;

    context.legacyFormat = islegacy;

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,                              "SHAA004TF");
    SEQAN_ASSERT_EQ(length(r.matches),                  17u);

    // TEST first match of this record extensively:
    {
        auto & m = r.matches.front();

        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TF");
        SEQAN_ASSERT_EQ(m.sId,                          "sp|P0A916|OMPW_SHIFL");
        SEQAN_ASSERT_EQ(m.alignStats.numMatches,        58u /*(115u * 50.43) / 100*/);
        SEQAN_ASSERT_EQ(m.alignStats.alignmentLength,   115u);
        if (path == LEGACY_COMMENTS_DEFAULTS) // legacy blast includes gaps in mismatches
            SEQAN_ASSERT_EQ(m.alignStats.numMismatches, 57u);
        else
            SEQAN_ASSERT_EQ(m.alignStats.numMismatches, 49u);
        SEQAN_ASSERT_EQ(m.qStart,                       389u);
        SEQAN_ASSERT_EQ(m.qEnd,                         733u);
        SEQAN_ASSERT_EQ(m.sStart,                       1u);
        SEQAN_ASSERT_EQ(m.sEnd,                         107u);

        if (defaults)
        {
            SEQAN_ASSERT_EQ(m.alignStats.numGapOpens,   2u);
            SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 108), 1e-3);
            SEQAN_ASSERT_LEQ(std::abs(m.eValue  - 1e-26), 1e-10);
        } else
        {
            SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores,    72u);
            SEQAN_ASSERT_EQ(m.qFrameShift,              2);
            SEQAN_ASSERT_EQ(m.sFrameShift,              0);
        }
        // legacy blast includes gaps in mismatches, so gaps are only computed
        // for BLAST_PLUS
        if (!context.legacyFormat)
        {
            if  (path == LEGACY_COMMENTS_DEFAULTS) // legacy output with plus reader
                SEQAN_ASSERT_EQ(m.alignStats.numGaps,   0u /*(115u - 58 - 49)*/);
            else // plus output and plus reader
                SEQAN_ASSERT_EQ(m.alignStats.numGaps,   8u /*(115u - 58 - 49)*/);
        } else
        {
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,       std::numeric_limits<decltype(m.alignStats.numGaps)>::max());
        }
    }

    // TEST basic stuff on other matches
    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TF"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_GEQ(m.bitScore,                    40.8);
        SEQAN_ASSERT_LEQ(m.bitScore,                    108.0);
    }

    if (defaults)
    {
        SEQAN_ASSERT_EQ(length(context.fields),         1u);
        SEQAN_ASSERT_EQ((uint8_t)context.fields[0],     (uint8_t)TField::Enum::STD);
    } else
    {
        SEQAN_ASSERT_EQ(length(context.fields),         13u);
        for (unsigned i = 0; i < 13u; ++i)
            SEQAN_ASSERT_EQ((uint8_t)context.fields[i], (uint8_t)fieldsCustom[i]);
    }

    readRecord(r, fileIn);

    SEQAN_ASSERT_EQ(r.qId,                              "SHAA004TR");
    SEQAN_ASSERT_EQ(length(r.matches),                  2u);

    // TEST last match of this record extensively:
    {
        auto & m = r.matches.back();

        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TR");
        SEQAN_ASSERT_EQ(m.sId,                          "sp|Q0HGZ8|META_SHESM");
        SEQAN_ASSERT_EQ(m.alignStats.numMatches,        77u /* (77u * 100) / 100 */);
        SEQAN_ASSERT_EQ(m.alignStats.alignmentLength,   77u);
        SEQAN_ASSERT_EQ(m.alignStats.numMismatches,     0u);
        SEQAN_ASSERT_EQ(m.qStart,                       232u);
        SEQAN_ASSERT_EQ(m.qEnd,                         2u);
        SEQAN_ASSERT_EQ(m.sStart,                       1u);
        SEQAN_ASSERT_EQ(m.sEnd,                         77u);

        // legacy blast includes gaps in mismatches, so gaps are only computed
        // for BLAST_PLUS
        if (!islegacy)
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,       0u); // always zero here
        else
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,       std::numeric_limits<decltype(m.alignStats.numGaps)>::max());

        if (defaults)
        {
            SEQAN_ASSERT_EQ(m.alignStats.numGapOpens,  0u);
            SEQAN_ASSERT_LEQ(std::abs(m.eValue  - 3e-42), 1e-10);
            SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 152), 1e-3);
        } else
        {
            SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores, 77u);
            SEQAN_ASSERT_EQ(m.qFrameShift,              -2);
            SEQAN_ASSERT_EQ(m.sFrameShift,              0);
        }
    }

    // TEST basic stuff on other matches
    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TR"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_EQ(m.bitScore,                     152.0);
    }

    // bottom of file
    readFooter(fileIn);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),  0u);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_without_comments)
{
    BlastIOContext<> context;
    _testReadTabularWithoutComments(context, NOCOMMENTS_DEFAULTS, true, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_without_comments_customfields)
{
    BlastIOContext<> context;
    _testReadTabularWithoutComments(context, NOCOMMENTS_CUSTOM, false, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_without_comments_legacy)
{
    BlastIOContext<> context;
    _testReadTabularWithoutComments(context, NOCOMMENTS_DEFAULTS, true, true);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_without_comments_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_COMMENTS> context;
    _testReadTabularWithoutComments(context, NOCOMMENTS_DEFAULTS, true, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_without_comments_customfields_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_COMMENTS> context;
    _testReadTabularWithoutComments(context, NOCOMMENTS_CUSTOM, false, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_without_comments_legacy_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTP, BlastTabularSpec::NO_COMMENTS> context;
    _testReadTabularWithoutComments(context, NOCOMMENTS_DEFAULTS, true, true);
}

template <typename TContext>
void _testReadTabularWithComments(TContext &,
                                std::string path,
                                bool const defaults,
                                bool const islegacy)
{
    typedef BlastMatchField<> TField;

    std::string inPath = seqan::getAbsolutePath(path.c_str());

    BlastRecord<> r;

    BlastTabularFileIn<TContext> fileIn(toCString(inPath));
    TContext & context = seqan::context(fileIn);

    /* Variables we will be comparing with */
    std::vector<typename TField::Enum> const fieldsDefault =
    {
        TField::Enum::STD
    };

    std::vector<typename TField::Enum> const fieldsCustom =
    {
        TField::Enum::Q_SEQ_ID,
        TField::Enum::S_SEQ_ID,
        TField::Enum::LENGTH,
        TField::Enum::N_IDENT,
        TField::Enum::MISMATCH,
        TField::Enum::POSITIVE,
        TField::Enum::GAPS,
        TField::Enum::Q_START,
        TField::Enum::Q_END,
        TField::Enum::S_START,
        TField::Enum::S_END,
        TField::Enum::FRAMES,
        TField::Enum::BIT_SCORE
    };

    //TODO add fieldCustom2 without BIT_SCORE

    std::string const dbString = "/tmp/uniprot_sprot.fasta";
    std::string const vString  = islegacy ? "BLASTX 2.2.26 [Sep-21-2011]" : "BLASTX 2.2.26+";

    /* Begin of TESTS */

    // read comments of file
    readHeader(fileIn);
    SEQAN_ASSERT(context.tabularSpec == BlastTabularSpec::COMMENTS);

    readRecord(r, fileIn);

    SEQAN_ASSERT_EQ(context.legacyFormat, islegacy);

    SEQAN_ASSERT_EQ(r.qId,                              "SHAA003TF  Sample 1 Mate SHAA003TR trimmed_to 27 965");
    SEQAN_ASSERT_EQ(context.dbName,                     dbString);
    SEQAN_ASSERT_EQ(context.versionString,              vString);
    SEQAN_ASSERT_EQ(length(r.matches),                  0u);
    SEQAN_ASSERT_EQ(length(context.fields),             1u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),    12u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.otherLines),         0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),  0u);

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,                              "SHAA003TR  Sample 1 Mate SHAA003TF trimmed_to 17 935");
    SEQAN_ASSERT_EQ(context.dbName,                     dbString);
    SEQAN_ASSERT_EQ(context.versionString,              vString);
    SEQAN_ASSERT_EQ(length(r.matches),                  0u);
    SEQAN_ASSERT_EQ(length(context.fields),             1u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),    12u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.otherLines),         0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),  0u);

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,                              "SHAA004TF  Sample 1 Mate SHAA004TR trimmed_to 25 828");
    SEQAN_ASSERT_EQ(context.dbName,                     dbString);
    SEQAN_ASSERT_EQ(context.versionString,              vString);
    SEQAN_ASSERT_EQ(length(r.matches),                  17u);
    SEQAN_ASSERT_EQ(length(context.otherLines),         0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),  0u);

    if (!defaults)
    {
        SEQAN_ASSERT_EQ(length(context.fields),         13u);
        SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),13u);
        for (uint8_t i = 0; i < 13; ++i)
            SEQAN_ASSERT_EQ((uint32_t)context.fields[i],(uint32_t)fieldsCustom[i]);
    } else
    {
        SEQAN_ASSERT_EQ(length(context.fields),         1u);
        SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),12u); // defaults
        SEQAN_ASSERT_EQ((uint32_t)context.fields[0],    (uint32_t)TField::Enum::STD);
    }

    // TEST first match of this record extensively:
    {
        auto & m = r.matches.front();

        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TF");
        SEQAN_ASSERT_EQ(m.sId,                          "sp|P0A916|OMPW_SHIFL");
        SEQAN_ASSERT_EQ(m.alignStats.numMatches,        58u /*(115u * 50.43) / 100*/);
        SEQAN_ASSERT_EQ(m.alignStats.alignmentLength,   115u);
        if (path == LEGACY_COMMENTS_DEFAULTS) // legacy blast includes gaps in mismatches
            SEQAN_ASSERT_EQ(m.alignStats.numMismatches, 57u);
        else
            SEQAN_ASSERT_EQ(m.alignStats.numMismatches, 49u);
        SEQAN_ASSERT_EQ(m.qStart,                       389u);
        SEQAN_ASSERT_EQ(m.qEnd,                         733u);
        SEQAN_ASSERT_EQ(m.sStart,                       1u);
        SEQAN_ASSERT_EQ(m.sEnd,                         107u);

        if (defaults)
        {
            SEQAN_ASSERT_EQ(m.alignStats.numGapOpens,   2u);
            SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 108), 1e-3);
            SEQAN_ASSERT_LEQ(std::abs(m.eValue  - 1e-26), 1e-10);
        } else
        {
            SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores,    72u);
            SEQAN_ASSERT_EQ(m.qFrameShift,              2);
            SEQAN_ASSERT_EQ(m.sFrameShift,              0);
        }
        // legacy blast includes gaps in mismatches, so gaps are only computed
        // for BLAST_PLUS
        if (!context.legacyFormat)
        {
            if  (path == LEGACY_COMMENTS_DEFAULTS) // legacy output with plus reader
                SEQAN_ASSERT_EQ(m.alignStats.numGaps,   0u /*(115u - 58 - 49)*/);
            else // plus output and plus reader
                SEQAN_ASSERT_EQ(m.alignStats.numGaps,   8u /*(115u - 58 - 49)*/);
        } else
        {
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,       std::numeric_limits<decltype(m.alignStats.numGaps)>::max());
        }
    }

    // TEST basic stuff on other matches
    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TF"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_GEQ(m.bitScore,                    40.8);
        SEQAN_ASSERT_LEQ(m.bitScore,                    108.0);
    }

    if (defaults)
    {
        SEQAN_ASSERT_EQ(length(context.fields),         1u);
        SEQAN_ASSERT_EQ((uint8_t)context.fields[0],     (uint8_t)TField::Enum::STD);
    } else
    {
        SEQAN_ASSERT_EQ(length(context.fields),         13u);
        for (unsigned i = 0; i < 13u; ++i)
            SEQAN_ASSERT_EQ((uint8_t)context.fields[i], (uint8_t)fieldsCustom[i]);
    }

    // fieldsList as in-parameter
    if (defaults)
        context.fields = fieldsDefault;
    else
        context.fields = fieldsCustom;

    context.ignoreFieldsInComments = true;

    readRecord(r, fileIn);

    SEQAN_ASSERT_EQ(r.qId,                              "SHAA004TR  Sample 1 Mate SHAA004TF trimmed_to 20 853");
    SEQAN_ASSERT_EQ(context.dbName,                     dbString);
    SEQAN_ASSERT_EQ(context.versionString,              vString);
    SEQAN_ASSERT_EQ(length(r.matches),                  2u);
    SEQAN_ASSERT_EQ(length(context.otherLines),         0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),  0u);

    // TEST last match of this record extensively:
    {
        auto & m = r.matches.back();

        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TR");
        SEQAN_ASSERT_EQ(m.sId,                          "sp|Q0HGZ8|META_SHESM");
        SEQAN_ASSERT_EQ(m.alignStats.numMatches,        77u /* (77u * 100) / 100 */);
        SEQAN_ASSERT_EQ(m.alignStats.alignmentLength,   77u);
        SEQAN_ASSERT_EQ(m.alignStats.numMismatches,     0u);
        SEQAN_ASSERT_EQ(m.qStart,                       232u);
        SEQAN_ASSERT_EQ(m.qEnd,                         2u);
        SEQAN_ASSERT_EQ(m.sStart,                       1u);
        SEQAN_ASSERT_EQ(m.sEnd,                         77u);

        // legacy blast includes gaps in mismatches, so gaps are only computed
        // for BLAST_PLUS
        if (!islegacy)
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,       0u); // always zero here
        else
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,       std::numeric_limits<decltype(m.alignStats.numGaps)>::max());

        if (defaults)
        {
            SEQAN_ASSERT_EQ(m.alignStats.numGapOpens,  0u);
            SEQAN_ASSERT_LEQ(std::abs(m.eValue  - 3e-42), 1e-10);
            SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 152), 1e-3);
        } else
        {
            SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores, 77u);
            SEQAN_ASSERT_EQ(m.qFrameShift,              -2);
            SEQAN_ASSERT_EQ(m.sFrameShift,              0);
        }
    }

    // TEST basic stuff on other matches
    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,                          "SHAA004TR"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_EQ(m.bitScore,                     152.0);
    }

    context.ignoreFieldsInComments = false;

    readRecord(r, fileIn);

    SEQAN_ASSERT_EQ(context.legacyFormat,               islegacy);

    SEQAN_ASSERT_EQ(r.qId,                              "SHAA005TR  Sample 1 Mate SHAA005TF trimmed_to 22 960");
    SEQAN_ASSERT_EQ(context.dbName,                     "");
    SEQAN_ASSERT_EQ(context.versionString,              vString);
    SEQAN_ASSERT_EQ(length(r.matches),                  0u);
    SEQAN_ASSERT_EQ(length(context.otherLines),         1u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),  2u);

    SEQAN_ASSERT_EQ(context.otherLines[0],              " Datacase: /tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.conformancyErrors[0],       "No or multiple database lines present.");

    SEQAN_ASSERT(!onRecord(fileIn));
    // bottom of file
    readFooter(fileIn);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),  0u);

}

SEQAN_DEFINE_TEST(test_blast_read_tabular_with_comments)
{
    BlastIOContext<> context;
    _testReadTabularWithComments(context, PLUS_COMMENTS_DEFAULTS, true, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_with_comments_customfields)
{
    BlastIOContext<> context;
    _testReadTabularWithComments(context, PLUS_COMMENTS_CUSTOM, false, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_with_comments_legacy)
{
    BlastIOContext<> context;
    _testReadTabularWithComments(context, LEGACY_COMMENTS_DEFAULTS, true, true);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_with_comments_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTX, BlastTabularSpec::COMMENTS> context;
    _testReadTabularWithComments(context, PLUS_COMMENTS_DEFAULTS, true, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_with_comments_customfields_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTX, BlastTabularSpec::COMMENTS> context;
    _testReadTabularWithComments(context, PLUS_COMMENTS_CUSTOM, false, false);
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_with_comments_legacy_constexpr)
{
    BlastIOContext<Blosum62, BlastProgram::BLASTX, BlastTabularSpec::COMMENTS> context;
    _testReadTabularWithComments(context, LEGACY_COMMENTS_DEFAULTS, true, true);
}
