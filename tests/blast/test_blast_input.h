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

#include <fstream>

// Files that are being read by this implementation
#define PLUS_HEADER_DEFAULTS   "/tests/blast/plus_header_defaults.blast"
#define LEGACY_HEADER_DEFAULTS "/tests/blast/legacy_header_defaults.blast"
#define PLUS_HEADER_CUSTOM     "/tests/blast/plus_header_custom.m9"
// same for HEADER and NOHEADERS:
#define NOHEADER_DEFAULTS      "/tests/blast/noheader_defaults.m8"
// same for HEADER and NOHEADERS:
#define NOHEADER_CUSTOM        "/tests/blast/noheader_custom.m8"

using namespace seqan;

template <typename TNum,
          typename std::enable_if<Is<NumberConcept<TNum>>::VALUE, int>::type = 0>
inline void
clear(TNum & num)
{
    num = 0;
}

template <typename TArg,
          typename... TArgs>
inline void
clear(TArg & arg, TArgs & ... args)
{
    clear(arg);
    clear(args...);
}

//legacy doesnt have interface with fieldList so we need to wrap around the call
// template <typename TMatch,
//           typename TIt,
//           typename TFields,
//           BlastFormatFile f,
//           BlastFormatProgram p>
// inline void
// _readMatchWrap(TMatch & match, TIt & it, BlastIOContext & context,
//                TFields const &, bool,
//                BlastFormat<f, p, BlastFormatGeneration::BLAST_LEGACY> const &)
// {
//     typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_LEGACY> BlastTabular;
//     readMatch(match, it, context, BlastTabular());
// }
//
// template <typename TMatch,
//           typename TIt,
//           typename TFields,
//           BlastFormatFile f,
//           BlastFormatProgram p>
// inline void
// _readMatchWrap(TMatch & match, TIt & it, BlastIOContext & context,
//                TFields const & fields, bool useFields,
//                BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> const &)
// {
//     typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> BlastTabular;
//     if (useFields)
//         readMatch(match, it, context, BlastTabular());
//     else
//         readMatch(match, it, context, fields, BlastTabular());
// }

void _test_blast_read_tabular_match(std::string const & path,
                                    bool const defaults,
                                    bool const withheader,
                                    bool const islegacy)
{
    // only used when defaults == false
    typedef BlastMatchField<> TField;

    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + std::string(path);

    std::ifstream ifstream(toCString(inPath),
                           std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(ifstream.is_open());

    auto it = directionIterator(ifstream, Input());
    BlastIOContext<> context;
    if (!defaults)
    {
        std::vector<typename TField::Enum> customFields
        {
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
                // actually BITSCORE is also in the file, but we are checking
                // that the implementation copes with less requested fields
                // than available ones which we defined as legal.
            }
        };
        context.fields = customFields;
    }

    if (withheader)
        setBlastTabularSpec(context, BlastTabularSpec::HEADER);
    else
        setBlastTabularSpec(context, BlastTabularSpec::NO_HEADER);

    context.legacyFormat = islegacy;

//     Iterator<std::ifstream, Rooted>::Type it = begin(ifstream);
//     auto it = directionIterator(ifstream, Bidirectional());
//     FormattedFile<> it(ifstream);

    // skip headers and comment lines
    skipUntilMatch(it, BlastTabular());

    //---- FIRST MATCH ---- //

    // now we should be onMatch
    SEQAN_ASSERT(onMatch(it, BlastTabular()));

    BlastMatch<> m;
    readMatch(m, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(m.qId,          "SHAA004TF");
    SEQAN_ASSERT_EQ(m.sId,          "sp|P0A916|OMPW_SHIFL");
    SEQAN_ASSERT_EQ(m.alignStats.numMatches,   58u /*(115u * 50.43) / 100*/);
    SEQAN_ASSERT_EQ(m.alignStats.alignmentLength,    115u);
    if (path == LEGACY_HEADER_DEFAULTS) // legacy blast includes gaps in mismatches
        SEQAN_ASSERT_EQ(m.alignStats.numMismatches,   57u);
    else
        SEQAN_ASSERT_EQ(m.alignStats.numMismatches,   49u);
    SEQAN_ASSERT_EQ(m.qStart,       389u);
    SEQAN_ASSERT_EQ(m.qEnd,         733u);
    SEQAN_ASSERT_EQ(m.sStart,       1u);
    SEQAN_ASSERT_EQ(m.sEnd,         107u);

    if (defaults)
    {
        SEQAN_ASSERT_EQ(m.alignStats.numGapOpens,  2u);
        SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 108), 1e-3);
        SEQAN_ASSERT_LEQ(std::abs(m.eValue  - 1e-26), 1e-10);
    } else
    {
        SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores,    72u);
        SEQAN_ASSERT_EQ(m.qFrameShift,  2);
        SEQAN_ASSERT_EQ(m.sFrameShift,  0);
    }
    // legacy blast includes gaps in mismatches, so gaps are only computed
    // for BLAST_PLUS
    if (!context.legacyFormat)
    {
        if  (path == LEGACY_HEADER_DEFAULTS) // legacy output with plus reader
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,     0u /*(115u - 58 - 49)*/);
        else // plus output and plus reader
            SEQAN_ASSERT_EQ(m.alignStats.numGaps,     8u /*(115u - 58 - 49)*/);
    } else
    {
        SEQAN_ASSERT_EQ(m.alignStats.numGaps, std::numeric_limits<decltype(m.alignStats.numGaps)>::max());
    }
    //---- FIRST MATCH END ---- //

    // now we should be onMatch
    SEQAN_ASSERT(onMatch(it, BlastTabular()));

    unsigned count = 0;
    if ((path == PLUS_HEADER_DEFAULTS) || (path == LEGACY_HEADER_DEFAULTS))
    {
        while (onMatch(it, BlastTabular()))
        {
//             if (defaults) // skipMatch with verification
//                 skipMatch(it, context, BlastTabular());
//             else
                skipMatch(it, context, BlastTabular());

            ++count;
        }

        // 16 further matches where skipped
        SEQAN_ASSERT_EQ(count, 16u);

        // back on header
        SEQAN_ASSERT(!onMatch(it, BlastTabular()));
    }
    else
    {
        for (int i = 0; i < 16; ++i)
        {
            SEQAN_ASSERT(onMatch(it, BlastTabular()));
//             if (defaults) // skipMatch with verification
//                 skipMatch(it, context, BlastTabular());
//             else
                skipMatch(it, context, BlastTabular());
        }
    }

    skipUntilMatch(it, BlastTabular());

    // skipMatch without verification
    skipMatch(it, context, BlastTabular());

    //---- LAST MATCH ---- //
    // read another match
//     _readMatchWrap(m, it, context, customFields, defaults, BlastTabular());
    readMatch(m, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(m.qId,          "SHAA004TR");
    SEQAN_ASSERT_EQ(m.sId,          "sp|Q0HGZ8|META_SHESM");
    SEQAN_ASSERT_EQ(m.alignStats.numMatches,   77u /* (77u * 100) / 100 */);
    SEQAN_ASSERT_EQ(m.alignStats.alignmentLength,    77u);
    SEQAN_ASSERT_EQ(m.alignStats.numMismatches,   0u);
    SEQAN_ASSERT_EQ(m.qStart,       232u);
    SEQAN_ASSERT_EQ(m.qEnd,         2u);
    SEQAN_ASSERT_EQ(m.sStart,       1u);
    SEQAN_ASSERT_EQ(m.sEnd,         77u);

    // legacy blast includes gaps in mismatches, so gaps are only computed
    // for BLAST_PLUS
    if (!islegacy)
        SEQAN_ASSERT_EQ(m.alignStats.numGaps,     0u); // always zero here
    else
        SEQAN_ASSERT_EQ(m.alignStats.numGaps, std::numeric_limits<decltype(m.alignStats.numGaps)>::max());

    if (defaults)
    {
        SEQAN_ASSERT_EQ(m.alignStats.numGapOpens,  0u);
        SEQAN_ASSERT_LEQ(std::abs(m.eValue  - 3e-42), 1e-10);
        SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 152), 1e-3);
    } else
    {
        SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores,    77u);
        SEQAN_ASSERT_EQ(m.qFrameShift, -2);
        SEQAN_ASSERT_EQ(m.sFrameShift,  0);
    }
    //---- LAST MATCH END ---- //

    // check if exceptions are properly thrown
    CharString exceptComment;
    // we are on header again
    if ((path == PLUS_HEADER_DEFAULTS) || (path == LEGACY_HEADER_DEFAULTS))
    {
        SEQAN_TRY
        {
            readMatch(m, it, context, BlastTabular());
        } SEQAN_CATCH (ParseError const & e)
        {
            exceptComment = e.what();
        }

        SEQAN_ASSERT_EQ(exceptComment, "Not on beginning of Match (you should have skipped comments).");
    }

    exceptComment = "";
    SEQAN_TRY
    {
        skipUntilMatch(it, BlastTabular());
    } SEQAN_CATCH (ParseError const & e)
    {
        exceptComment = e.what();
    }

    SEQAN_ASSERT_EQ(exceptComment, "EOF reached without finding Match.");

    ifstream.close();
}

SEQAN_DEFINE_TEST(test_blast_read_match_tabular)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match(NOHEADER_DEFAULTS, true, false, false);
    // works when header is there, too
    _test_blast_read_tabular_match(PLUS_HEADER_DEFAULTS, true, false, false);
    // basic match reading should work between plus and legacy
    _test_blast_read_tabular_match(LEGACY_HEADER_DEFAULTS, true, false, false );
}

SEQAN_DEFINE_TEST(test_blast_read_match_tabular_legacy)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match(NOHEADER_DEFAULTS, true, false, true);
    // works when header is there, too
    _test_blast_read_tabular_match(PLUS_HEADER_DEFAULTS, true, false, true);
    // basic match reading should work between plus and legacy
    _test_blast_read_tabular_match(LEGACY_HEADER_DEFAULTS, true, false, true );
}

SEQAN_DEFINE_TEST(test_blast_read_match_customfields_tabular)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match(NOHEADER_CUSTOM, false, false, false);
    // works when header is there, too
    _test_blast_read_tabular_match(PLUS_HEADER_CUSTOM, false, false, false);
    // custom columns not supported in legacy format
}

SEQAN_DEFINE_TEST(test_blast_read_match_tabular_with_header)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match(PLUS_HEADER_DEFAULTS, true, true, false);
    // basic match reading should work between plus and legacy
    _test_blast_read_tabular_match(LEGACY_HEADER_DEFAULTS, true, true, false);
}

SEQAN_DEFINE_TEST(test_blast_read_match_tabular_with_header_legacy)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match(LEGACY_HEADER_DEFAULTS, true, true, true);
    // basic match reading should work between plus and legacy
    _test_blast_read_tabular_match(PLUS_HEADER_DEFAULTS, true, true, true);
}

SEQAN_DEFINE_TEST(test_blast_read_match_customfields_tabular_with_header)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match(PLUS_HEADER_CUSTOM, false, true, false);
    // custom columns not supported in legacy format
}

inline void
_test_blast_read_tabular_match_columns(std::string const & path)
{
    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + std::string(path);

    std::ifstream ifstream(toCString(inPath),
                           std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(ifstream.is_open());

    auto it = directionIterator(ifstream, Input());
//     BlastIOContext context;

    // first line is header
    SEQAN_ASSERT(!onMatch(it, BlastTabular()));

    // skip headers and comment lines
    skipUntilMatch(it, BlastTabular());

    // now we should be onMatch
    SEQAN_ASSERT(onMatch(it, BlastTabular()));

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

    readMatch(it, BlastTabular(), field1, field2, field3, field4, field5, field6,
              field7, field8, field9, field10, field11, field12);

    SEQAN_ASSERT_EQ(field1,          "SHAA004TF");
    SEQAN_ASSERT_EQ(field2,          "sp|P0A916|OMPW_SHIFL");
    SEQAN_ASSERT_EQ(field3,          50.43);
    SEQAN_ASSERT_EQ(field4,          115.0);
    if (path == LEGACY_HEADER_DEFAULTS)
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
    readMatch(it, BlastTabular(), field1, field2, field3, field4, field5, field6,
              field7, field8, field9, field10);

    SEQAN_ASSERT_EQ(field1,          "SHAA004TF");
    SEQAN_ASSERT_EQ(field2,          "sp|P0A915|OMPW_ECOLI");
    SEQAN_ASSERT_EQ(field3,          50.43);
    SEQAN_ASSERT_EQ(field4,          115.0);
    if (path == LEGACY_HEADER_DEFAULTS) // legacy blast includes gaps in mismatches
        SEQAN_ASSERT_EQ(field5,      57u);
    else
        SEQAN_ASSERT_EQ(field5,      49u);
    SEQAN_ASSERT_EQ(field6,          2u);
    SEQAN_ASSERT_EQ(field7,          389u);
    SEQAN_ASSERT_EQ(field8,          733u);
    SEQAN_ASSERT_EQ(field9,          1u);
    SEQAN_ASSERT_EQ(field10,         107u);

    // read remaining matches
    while (onMatch(it, BlastTabular()))
        readMatch(it, BlastTabular(), field1); // only one field

    // go to last record with matches
    skipUntilMatch(it, BlastTabular());

    // goto second (last) match
    skipLine(it);

    // check if exceptions are properly thrown
    CharString exceptComment;
    SEQAN_TRY
    {
        // no strings here to take the strings
        readMatch(it, BlastTabular(), field3, field4, field5, field6, field7, field8, field9, field10, field11,
                  field12);
    } SEQAN_CATCH (BadLexicalCast const & e)
    {
        exceptComment = e.what();
    }

    SEQAN_ASSERT(startsWith(exceptComment, "Unable to convert"));
    // skip rest of line
    skipLine(it);

    exceptComment = "";
    SEQAN_TRY
    {
        skipUntilMatch(it, BlastTabular());
    } SEQAN_CATCH (ParseError const & e)
    {
        exceptComment = e.what();
    }

    SEQAN_ASSERT_EQ(exceptComment, "EOF reached without finding Match.");

    ifstream.close();
}

SEQAN_DEFINE_TEST(test_blast_read_match_customcolumns_tabular)
{
    // this testsBlastTabular#readMatch(),BlastTabular#skipMatch() and
    //BlastTabular#skipUntilMatch
    _test_blast_read_tabular_match_columns(PLUS_HEADER_DEFAULTS);
    // basic match reading should work between plus and legacy
    _test_blast_read_tabular_match_columns(LEGACY_HEADER_DEFAULTS);
}

inline void
_test_blast_read_tabular_with_header(bool custom = false)
{
    // only used when defaults == false
    typedef BlastMatchField<> TField;
    std::array<typename TField::Enum, 13> customFields
    {
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
        }
    };

    std::string inPath = std::string(SEQAN_PATH_TO_ROOT());
    inPath += (custom ? PLUS_HEADER_CUSTOM : PLUS_HEADER_DEFAULTS);

    std::ifstream ifstream(toCString(inPath),
                           std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(ifstream.is_open());

    auto it = directionIterator(ifstream, Input());

    // first line is header
    SEQAN_ASSERT(!onMatch(it, BlastTabular()));

    BlastRecord<> r;
    BlastIOContext<> context;

    // back on header
    SEQAN_ASSERT(!onMatch(it, BlastTabular()));

    readRecordHeader(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(context.legacyFormat, false);

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA003TF  Sample 1 Mate SHAA003TR trimmed_to 27 965");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26+");
    SEQAN_ASSERT_EQ(length(r.matches),
                    0u);
    SEQAN_ASSERT_EQ(length(context.fields),
                    1u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    readRecordHeader(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA003TR  Sample 1 Mate SHAA003TF trimmed_to 17 935");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26+");
    SEQAN_ASSERT_EQ(length(r.matches),
                    0u);
    SEQAN_ASSERT_EQ(length(context.fields),
                    1u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    readRecordHeader(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA004TF  Sample 1 Mate SHAA004TR trimmed_to 25 828");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26+");
    SEQAN_ASSERT_EQ(length(r.matches),
                    17u);
    if (custom)
    {
        SEQAN_ASSERT_EQ(length(context.fields),
                        13u);
        SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                        13u);
        for (uint8_t i = 0; i < 13; ++i)
            SEQAN_ASSERT_EQ((uint32_t)context.fields[i],
                            (uint32_t)customFields[i]);
    } else
    {
        SEQAN_ASSERT_EQ(length(context.fields),
                        1u);
        SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                        12u); // defaults
        SEQAN_ASSERT_EQ((uint32_t)context.fields[0],
                        (uint32_t)TField::Enum::STD);
    }
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    while (onMatch(it, BlastTabular()))
        skipMatch(it, context, BlastTabular());

    readRecordHeader(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA004TR  Sample 1 Mate SHAA004TF trimmed_to 20 853");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26+");
    SEQAN_ASSERT_EQ(length(r.matches),
                    2u);
    if (custom)
    {
        SEQAN_ASSERT_EQ(length(context.fields),
                        13u);
        SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                        13u); // not overwritten because not available
        for (uint8_t i = 0; i < 13; ++i)
            SEQAN_ASSERT_EQ((uint32_t)context.fields[i],
                            (uint32_t)customFields[i]);
    } else
    {
        SEQAN_ASSERT_EQ(length(context.fields),
                        1u);
        SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                        12u); // defaults
        SEQAN_ASSERT_EQ((uint32_t)context.fields[0],
                        (uint32_t)TField::Enum::STD);
    }
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    while (onMatch(it, BlastTabular()))
        skipMatch(it, context, BlastTabular());

    readRecordHeader(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA005TR  Sample 1 Mate SHAA005TF trimmed_to 22 960");
    SEQAN_ASSERT_EQ(context.dbName,
                    "");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26+");
    SEQAN_ASSERT_EQ(length(r.matches),
                    0u);
    SEQAN_ASSERT_EQ(length(context.fields),
                    1u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u); // not overwritten because not available
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    1u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    2u);

    SEQAN_ASSERT_EQ(context.otherLines[0],
                    "Datacase: /tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.conformancyErrors[0],
                    "No or multiple database lines present.");

    // test skipHeader
    ifstream.close();
    ifstream.open(toCString(inPath));
    SEQAN_ASSERT(ifstream.is_open());

    it = directionIterator(ifstream, Input());

    skipHeader(it, context, BlastTabular());
    skipHeader(it, context, BlastTabular());
    skipHeader(it, context, BlastTabular());
    SEQAN_ASSERT(onMatch(it, BlastTabular()));

    ifstream.close();
}

SEQAN_DEFINE_TEST(test_blast_read_header_tabular_with_header)
{
    _test_blast_read_tabular_with_header(false);
}

SEQAN_DEFINE_TEST(test_blast_read_header_customfields_tabular_with_header)
{
    _test_blast_read_tabular_with_header(true);
}

SEQAN_DEFINE_TEST(test_blast_read_header_tabular_with_header_legacy)
{
    typedef BlastMatchField<> TField;

    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + LEGACY_HEADER_DEFAULTS;

    std::ifstream ifstream(toCString(inPath),
                           std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(ifstream.is_open());

    auto it = directionIterator(ifstream, Input());

    // first line is header
    SEQAN_ASSERT(!onMatch(it, BlastTabular()));

    BlastRecord<> r;
    CharString fieldStringsConcat;
    BlastIOContext<> context;

    // back on header
    SEQAN_ASSERT(!onMatch(it, BlastTabular()));

    // all parameters, strict == false
    readRecordHeader(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(context.legacyFormat, true);

    fieldStringsConcat = concat(context.fieldsAsStrings, ", ");

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA003TF  Sample 1 Mate SHAA003TR trimmed_to 27 965");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26 [Sep-21-2011]");
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u);
    SEQAN_ASSERT_EQ(fieldStringsConcat,
                    TField::columnLabels[0]);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);

    readRecordHeader(r, it, context, BlastTabular());

    fieldStringsConcat = concat(context.fieldsAsStrings, ", ");

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA003TR  Sample 1 Mate SHAA003TF trimmed_to 17 935");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26 [Sep-21-2011]");
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u);
    SEQAN_ASSERT_EQ(fieldStringsConcat,
                    TField::columnLabels[0]);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);

    readRecordHeader(r, it, context, BlastTabular());
    fieldStringsConcat = concat(context.fieldsAsStrings, ", ");

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA004TF  Sample 1 Mate SHAA004TR trimmed_to 25 828");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26 [Sep-21-2011]");
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u);
    SEQAN_ASSERT_EQ(fieldStringsConcat,
                    TField::columnLabels[0]);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);

    while (onMatch(it, BlastTabular()))
        skipMatch(it, context, BlastTabular());

    readRecordHeader(r, it, context, BlastTabular());
    fieldStringsConcat = concat(context.fieldsAsStrings, ", ");

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA004TR  Sample 1 Mate SHAA004TF trimmed_to 20 853");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26 [Sep-21-2011]");
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u);
    SEQAN_ASSERT_EQ(fieldStringsConcat,
                    TField::columnLabels[0]);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);

    while (onMatch(it, BlastTabular()))
        skipMatch(it, context, BlastTabular());

    readRecordHeader(r, it, context, BlastTabular());

    fieldStringsConcat = concat(context.fieldsAsStrings, ", ");

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA005TR  Sample 1 Mate SHAA005TF trimmed_to 22 960");
    SEQAN_ASSERT_EQ(context.dbName,
                    ""); // couldnt be retrieved
    SEQAN_ASSERT_EQ(context.versionString,
                    "BLASTX 2.2.26 [Sep-21-2011]");
    SEQAN_ASSERT_EQ(length(context.fieldsAsStrings),
                    12u);
    SEQAN_ASSERT_EQ(fieldStringsConcat,
                    TField::columnLabels[0]);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    1u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    2u);

    SEQAN_ASSERT_EQ(context.otherLines[0],
                    "Datacase: /tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.conformancyErrors[0],
                    "No or multiple database lines present.");

    // test skipHeader
    ifstream.close();
    ifstream.open(toCString(inPath));
    SEQAN_ASSERT(ifstream.is_open());

    it = directionIterator(ifstream, Input());

    skipHeader(it, context, BlastTabular());
    skipHeader(it, context, BlastTabular());
    skipHeader(it, context, BlastTabular());
    SEQAN_ASSERT(onMatch(it, BlastTabular()));

    ifstream.close();
}

void _test_blast_read_tabular_record_noheader(bool const defaults, bool const islegacy)
{
    typedef BlastMatchField<> TField;

    std::string inPath = std::string(SEQAN_PATH_TO_ROOT());
    inPath += (defaults ? NOHEADER_DEFAULTS : NOHEADER_CUSTOM);

    std::ifstream ifstream(toCString(inPath),
                           std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(ifstream.is_open());

    auto it = directionIterator(ifstream, Input());

    BlastRecord<> r;
    BlastIOContext<> context;

    if (defaults)
    {
        context.fields =
        {
            TField::Enum::STD
        };
    }
    else
    {
        context.fields =
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
    }

    // first line is onMatch
    SEQAN_ASSERT(onMatch(it, BlastTabular()));

    readRecord(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(context.legacyFormat, islegacy);

    SEQAN_ASSERT_EQ(r.qId,              "SHAA004TF");
    SEQAN_ASSERT_EQ(length(r.matches),  17u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TF"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_GEQ(m.bitScore,    40.8);
        SEQAN_ASSERT_LEQ(m.bitScore,    108.0);
        if (!defaults)
        {
            if (m.bitScore > 100)
                SEQAN_ASSERT_EQ(m.qFrameShift, 2);
            else
                SEQAN_ASSERT_EQ(m.qFrameShift, -1);
        }
    }

    readRecord(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(context.legacyFormat, islegacy);

    SEQAN_ASSERT_EQ(r.qId,              "SHAA004TR");
    SEQAN_ASSERT_EQ(length(r.matches),  2u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TR"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_EQ(m.bitScore,    152.0);
        if (!defaults)
            SEQAN_ASSERT_EQ(m.qFrameShift, -2);
    }

    ifstream.close();
}

SEQAN_DEFINE_TEST(test_blast_read_record_tabular)
{
    _test_blast_read_tabular_record_noheader(true, false);
}

SEQAN_DEFINE_TEST(test_blast_read_record_customfields_tabular)
{
    _test_blast_read_tabular_record_noheader(false, false);
}

SEQAN_DEFINE_TEST(test_blast_read_record_tabular_legacy)
{
    _test_blast_read_tabular_record_noheader(true, true);
}

//     typedef BlastFormat<BlastFormatFile::TABULAR,
//                         BlastFormatProgram::BLASTX,
//                         BlastFormatGeneration::BLAST_LEGACY> BlastTabular;
//
//     std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) +
//                          std::string(NOHEADER_DEFAULTS);
//
//     std::ifstream ifstream(toCString(inPath),
//                            std::ios_base::in | std::ios_base::binary);
//     SEQAN_ASSERT(ifstream.is_open());
//
//     auto it = directionIterator(ifstream, Input());
//
//     BlastRecord<> r;
//     BlastIOContext context;
//
//     // first line is match
//     SEQAN_ASSERT(onMatch(it, BlastTabular()));
//
//     readRecord(r, it, context, BlastTabular());
//     SEQAN_ASSERT_EQ(r.qId,              "SHAA004TF");
//     SEQAN_ASSERT_EQ(length(r.matches),  17u);
//
//     for (auto const & m : r.matches)
//     {
//         SEQAN_ASSERT_EQ(m.qId,          "SHAA004TF"); // truncated at first space
//         // check bitscore as an example field that is both in default and custom
//         SEQAN_ASSERT_GEQ(m.bitScore,    40.8);
//         SEQAN_ASSERT_LEQ(m.bitScore,    108.0);
//     }
//
//     readRecord(r, it, context, BlastTabular());
//     SEQAN_ASSERT_EQ(r.qId,              "SHAA004TR");
//     SEQAN_ASSERT_EQ(length(r.matches),  2u);
//
//     for (auto const & m : r.matches)
//     {
//         SEQAN_ASSERT_EQ(m.qId,          "SHAA004TR"); // truncated at first space
//         // check bitscore as an example field that is both in default and custom
//         SEQAN_ASSERT_EQ(m.bitScore,    152.0);
//     }
//
//     ifstream.close();
// }

void _test_blast_read_record_withheader(const char * path, bool const defaults, bool const islegacy)
{
    typedef BlastMatchField<> TField;

    std::string inPath = std::string(SEQAN_PATH_TO_ROOT());
    inPath += path;

    std::ifstream ifstream(toCString(inPath),
                           std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(ifstream.is_open());

    auto it = directionIterator(ifstream, Input());

    BlastRecord<> r;
    BlastIOContext<> context;

    std::vector<typename TField::Enum> fieldsDefault =
    {
        {
            TField::Enum::STD
        }
    };
    std::vector<typename TField::Enum> fieldsCustom =
    {
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
        }
    };
    // first line is header
    SEQAN_ASSERT(!onMatch(it, BlastTabular()));

    // fieldsList as out-parameter, but no fields, because no matches
    readRecord(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(context.legacyFormat, islegacy);

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA003TF  Sample 1 Mate SHAA003TR trimmed_to 27 965");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(length(r.matches),
                    0u);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    // fieldsList as out-parameter, but no fields, because no matches
    readRecord(r, it, context, BlastTabular());
    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA003TR  Sample 1 Mate SHAA003TF trimmed_to 17 935");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(length(r.matches),
                    0u);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    // fieldsList as out-parameter
    readRecord(r, it, context, BlastTabular());
    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA004TF  Sample 1 Mate SHAA004TR trimmed_to 25 828");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(length(r.matches),
                    17u);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TF"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_GEQ(m.bitScore,    40.8);
        SEQAN_ASSERT_LEQ(m.bitScore,    108.0);
    }

    if (defaults)
    {
        SEQAN_ASSERT_EQ(length(context.fields),  1u);
        SEQAN_ASSERT_EQ((uint8_t)context.fields[0], (uint8_t)TField::Enum::STD);
    } else
    {
        SEQAN_ASSERT_EQ(length(context.fields),  13u);
        for (unsigned i = 0; i < 13u; ++i)
            SEQAN_ASSERT_EQ((uint8_t)context.fields[i], (uint8_t)fieldsCustom[i]);
    }

    // fieldsList as in-parameter
    if (defaults)
        context.fields = fieldsDefault;
    else
        context.fields = fieldsCustom;

    context.ignoreFieldsInHeader = true;

    readRecord(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA004TR  Sample 1 Mate SHAA004TF trimmed_to 20 853");
    SEQAN_ASSERT_EQ(context.dbName,
                    "/tmp/uniprot_sprot.fasta");

    SEQAN_ASSERT_EQ(length(r.matches),
                    2u);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    0u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    0u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TR"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_EQ(m.bitScore,    152.0);
    }

    context.ignoreFieldsInHeader = false;

    readRecord(r, it, context, BlastTabular());

    SEQAN_ASSERT_EQ(context.legacyFormat, islegacy);

    SEQAN_ASSERT_EQ(r.qId,
                    "SHAA005TR  Sample 1 Mate SHAA005TF trimmed_to 22 960");
    SEQAN_ASSERT_EQ(context.dbName,
                    "");
    SEQAN_ASSERT_EQ(length(r.matches),
                    0u);
    SEQAN_ASSERT_EQ(length(context.otherLines),
                    1u);
    SEQAN_ASSERT_EQ(length(context.conformancyErrors),
                    2u);

    SEQAN_ASSERT_EQ(context.otherLines[0],
                    "Datacase: /tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(context.conformancyErrors[0],
                    "No or multiple database lines present.");

    ifstream.close();
}

SEQAN_DEFINE_TEST(test_blast_read_record_tabular_with_header)
{
    _test_blast_read_record_withheader(PLUS_HEADER_DEFAULTS, true, false);
}

SEQAN_DEFINE_TEST(test_blast_read_record_customfields_tabular_with_header)
{
    _test_blast_read_record_withheader(PLUS_HEADER_CUSTOM, false, false);
}

SEQAN_DEFINE_TEST(test_blast_read_record_tabular_with_header_legacy)
{
    _test_blast_read_record_withheader(LEGACY_HEADER_DEFAULTS, true, true);

}
//  {
//     std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + LEGACY_HEADER_DEFAULTS;
//
//     std::ifstream ifstream(toCString(inPath),
//                            std::ios_base::in | std::ios_base::binary);
//     SEQAN_ASSERT(ifstream.is_open());
//
//     auto it = directionIterator(ifstream, Input());
//
//     BlastRecord<> r;
//     std::string dbName;
//     BlastIOContext context;
//
//     // first line is header
//     SEQAN_ASSERT(!onMatch(it, BlastTabular()));
//
//     readRecord(r, dbName, it, context, BlastTabular());
//     SEQAN_ASSERT_EQ(r.qId,              "SHAA003TF  Sample 1 Mate SHAA003TR trimmed_to 27 965");
//     SEQAN_ASSERT_EQ(dbName,           "/tmp/uniprot_sprot.fasta");
//     SEQAN_ASSERT_EQ(length(r.matches),  0u);
//
//     readRecord(r, dbName, it, context, BlastTabular());
//     SEQAN_ASSERT_EQ(r.qId,              "SHAA003TR  Sample 1 Mate SHAA003TF trimmed_to 17 935");
//     SEQAN_ASSERT_EQ(dbName,           "/tmp/uniprot_sprot.fasta");
//     SEQAN_ASSERT_EQ(length(r.matches),  0u);
//
//     readRecord(r, dbName, it, context, BlastTabular());
//     SEQAN_ASSERT_EQ(r.qId,              "SHAA004TF  Sample 1 Mate SHAA004TR trimmed_to 25 828");
//     SEQAN_ASSERT_EQ(dbName,           "/tmp/uniprot_sprot.fasta");
//     SEQAN_ASSERT_EQ(length(r.matches),  17u);
//
//     for (auto const & m : r.matches)
//     {
//         SEQAN_ASSERT_EQ(m.qId,          "SHAA004TF"); // truncated at first space
//         // check bitscore as an example field that is both in default and custom
//         SEQAN_ASSERT_GEQ(m.bitScore,    40.8);
//         SEQAN_ASSERT_LEQ(m.bitScore,    108.0);
//     }
//
//     readRecord(r, dbName, it, context, BlastTabular());
//     SEQAN_ASSERT_EQ(r.qId,              "SHAA004TR  Sample 1 Mate SHAA004TF trimmed_to 20 853");
//     SEQAN_ASSERT_EQ(dbName,           "/tmp/uniprot_sprot.fasta");
//     SEQAN_ASSERT_EQ(length(r.matches),  2u);
//
//     for (auto const & m : r.matches)
//     {
//         SEQAN_ASSERT_EQ(m.qId,          "SHAA004TR"); // truncated at first space
//         // check bitscore as an example field that is both in default and custom
//         SEQAN_ASSERT_EQ(m.bitScore,    152.0);
//     }
//
//     readRecord(r, dbName, it, context, BlastTabular());
//     SEQAN_ASSERT_EQ(r.qId,              "SHAA005TR  Sample 1 Mate SHAA005TF trimmed_to 22 960");
//     SEQAN_ASSERT_EQ(dbName,           ""); // couldn't be detected (typo)
//     SEQAN_ASSERT_EQ(length(r.matches),  0u);
// }


SEQAN_DEFINE_TEST(test_blast_read_formatted_file_tabular)
{
    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) +
                         std::string(NOHEADER_DEFAULTS);

    BlastTabularIn<> fileIn(toCString(inPath));

    guessFormat(fileIn);
    SEQAN_ASSERT(getBlastTabularSpec(context(fileIn)) == BlastTabularSpec::NO_HEADER);
    // legacy format or not and BlastProgram are detected upon reading records

    BlastRecord<> r;

    readRecord(r, fileIn);

    SEQAN_ASSERT_EQ(r.qId,              "SHAA004TF");
    SEQAN_ASSERT_EQ(length(r.matches),  17u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TF"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_GEQ(m.bitScore,    40.8);
        SEQAN_ASSERT_LEQ(m.bitScore,    108.0);
    }

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,              "SHAA004TR");
    SEQAN_ASSERT_EQ(length(r.matches),  2u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TR"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_EQ(m.bitScore,    152.0);
    }
}

SEQAN_DEFINE_TEST(test_blast_read_formatted_file_tabular_with_header)
{
    std::string inPath = std::string(SEQAN_PATH_TO_ROOT()) + PLUS_HEADER_DEFAULTS;


    BlastTabularIn<> fileIn(toCString(inPath));

    guessFormat(fileIn);
    SEQAN_ASSERT(getBlastTabularSpec(context(fileIn)) == BlastTabularSpec::HEADER);
    // legacy format or not and BlastProgram are detected upon reading records

    BlastRecord<> r;

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,                 "SHAA003TF  Sample 1 Mate SHAA003TR trimmed_to 27 965");
    SEQAN_ASSERT_EQ(context(fileIn).dbName, "/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(length(r.matches),  0u);

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,              "SHAA003TR  Sample 1 Mate SHAA003TF trimmed_to 17 935");
    SEQAN_ASSERT_EQ(context(fileIn).dbName,"/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(length(r.matches),  0u);

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,              "SHAA004TF  Sample 1 Mate SHAA004TR trimmed_to 25 828");
    SEQAN_ASSERT_EQ(context(fileIn).dbName,"/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(length(r.matches),  17u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TF"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_GEQ(m.bitScore,    40.8);
        SEQAN_ASSERT_LEQ(m.bitScore,    108.0);
    }

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,              "SHAA004TR  Sample 1 Mate SHAA004TF trimmed_to 20 853");
    SEQAN_ASSERT_EQ(context(fileIn).dbName,"/tmp/uniprot_sprot.fasta");
    SEQAN_ASSERT_EQ(length(r.matches),  2u);

    for (auto const & m : r.matches)
    {
        SEQAN_ASSERT_EQ(m.qId,          "SHAA004TR"); // truncated at first space
        // check bitscore as an example field that is both in default and custom
        SEQAN_ASSERT_EQ(m.bitScore,    152.0);
    }

    readRecord(r, fileIn);
    SEQAN_ASSERT_EQ(r.qId,              "SHAA005TR  Sample 1 Mate SHAA005TF trimmed_to 22 960");
    SEQAN_ASSERT_EQ(context(fileIn).dbName,""); // couldn't be detected (typo)
    SEQAN_ASSERT_EQ(length(r.matches),  0u);
}
