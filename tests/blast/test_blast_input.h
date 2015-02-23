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

using namespace seqan;


inline void
_writeExampleFile(std::fstream & stream)
{
    char const * string =
"# BLASTX 2.2.26+\n"
"# Query: SHAA001TF  Sample 1 Mate SHAA001TR trimmed_to 85 969\n"
"# Subject: sp|B3TN98|CCSA_BRADI Cytochrome c biogenesis protein ccsA OS=Brachypodium distachyon GN=ccsA PE=3 SV=1\n"
"# 0 hits found\n"
"# BLASTX 2.2.26+\n"
"# Query: SHAA001TF  Sample 1 Mate SHAA001TR trimmed_to 85 969\n"
"# Subject: sp|A6MM86|CCSA_BUXMI Cytochrome c biogenesis protein ccsA OS=Buxus microphylla GN=ccsA PE=3 SV=1\n"
"# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n"
"# 2 hits found\n"
"SHAA001TF\tsp|A6MM86|CCSA_BUXMI\t29.66\t118\t65\t7\t469\t783\t145\t257\t3e-05\t30.4\n"
"SHAA001TF\tsp|A6MM86|CCSA_BUXMI\t29.63\t27\t19\t0\t7\t87\t236\t262\t0.017\t21.9\n"
"# BLASTX 2.2.26+\n"
"# Query: SHAA001TF  Sample 1 Mate SHAA001TR trimmed_to 85 969\n"
"# Subject: sp|Q7YJT4|CCSA_CALFG Cytochrome c biogenesis protein ccsA OS=Calycanthus floridus var. glaucus GN=ccsA PE=3 SV=1\n"
"# 1 hits found\n"
"SHAA001TF\tsp|Q7YJT4|CCSA_CALFG\t29.63\t27\t19\t0\t7\t87\t236\t262\t0.017\t21.9\n"
"# BLASTX 2.2.26+\n"
"# Query: SHAA001TF  Sample 1 Mate SHAA001TR trimmed_to 85 969\n"
"# Subject: sp|A4QKP4|CCSA_CAPBU Cytochrome c biogenesis protein ccsA OS=Capsella bursa-pastoris GN=ccsA PE=3 SV=1\n"
"# 0 hits found\n";

    stream << string;
}

template <typename TFormat>
void _test_blast_read_tabular_basics_match(TFormat const &)
{
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    _writeExampleFile(fstream);
    fstream.seekg(0, std::ios::beg); // go back to beginning
    fstream.seekp(0, std::ios::beg); // go back to beginning
    auto it = directionIterator(fstream, Input());

//     Iterator<std::fstream, Rooted>::Type it = begin(fstream);
//     auto it = directionIterator(fstream, Bidirectional());
//     FormattedFile<> it(fstream);

    // first line is header
    SEQAN_ASSERT(!onMatch(it, TFormat()));

    // skip headers and comment lines
    skipUntilMatch(it, TFormat());

    //---- FIRST MATCH ---- //

    // now we should be onMatch
    SEQAN_ASSERT(onMatch(it, TFormat()));

    BlastMatch<> m;
    readMatch(m, it, TFormat());

    SEQAN_ASSERT_EQ(m.qId,          "SHAA001TF");
    SEQAN_ASSERT_EQ(m.sId,          "sp|A6MM86|CCSA_BUXMI");
    SEQAN_ASSERT_EQ(m.identities,   35u /*(118u * 29.66) / 100*/);
    SEQAN_ASSERT_EQ(m.aliLength,    118u);
    SEQAN_ASSERT_EQ(m.mismatches,   65u);
    SEQAN_ASSERT_EQ(m.gapOpenings,  7u);
    SEQAN_ASSERT_EQ(m.qStart,       469u);
    SEQAN_ASSERT_EQ(m.qEnd,         783u);
    SEQAN_ASSERT_EQ(m.sStart,       145u);
    SEQAN_ASSERT_EQ(m.sEnd,         257u);
    SEQAN_ASSERT_LEQ(std::abs(m.eValue  - 3e-05), 1e-10);
    SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 30.4), 1e-3);

    SEQAN_ASSERT_EQ(m.gaps,         18u /*(118u - 35 - 65)*/);

    //---- SECOND MATCH ---- //

    // now we should be onMatch
    SEQAN_ASSERT(onMatch(it, TFormat()));

    // skipMatch with verification
    skipMatch(it, false, TFormat());

    // back on header
    SEQAN_ASSERT(!onMatch(it, TFormat()));

    skipUntilMatch(it, TFormat());

    // skipMatch without verification
    skipMatch(it, TFormat());

    // check if exceptions are properly thrown
    CharString exceptComment;
    try
    {
        readMatch(m, it, TFormat());
    } catch (ParseError const & e)
    {
        exceptComment = e.what();
    }

    SEQAN_ASSERT_EQ(exceptComment, "Not on beginning of Match (you should have "
                                   "skipped comments).");

    exceptComment = "";
    try
    {
        skipUntilMatch(it, TFormat());
    } catch (ParseError const & e)
    {
        exceptComment = e.what();
    }

    SEQAN_ASSERT_EQ(exceptComment, "EOF reached without finding Match.");

    fstream.close();
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_basics_match)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        BlastFormatProgram::BLASTX,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    _test_blast_read_tabular_basics_match(TFormat());
}

SEQAN_DEFINE_TEST(test_blast_read_tabular_basics_match_legacy)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        BlastFormatProgram::BLASTX,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;
    _test_blast_read_tabular_basics_match(TFormat());
}
