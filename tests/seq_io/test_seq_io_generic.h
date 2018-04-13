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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Generic test code for seq_io module.
// ==========================================================================

#ifndef TEST_STREAM_TEST_SEQ_IO_GENERIC_H_
#define TEST_STREAM_TEST_SEQ_IO_GENERIC_H_

#include <fstream>

std::fstream* createFastAFile(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR =
"> sequenceID_with special chars an irregular linebreaks\n\
AAAACGTGCGGTTGGGCAAAAAACTTTCTTATATTCTATCTATCTTGTAGCTAGCTGTAGCTAGCTAGCATCGTAGC\n\
CCCAGAGTGTCATGCATGTCGA\n\
ACGTGTTTTTGGGGCGGTTATATATATATATATT\n\
\n\
>sequence2... with no linebreaks and no newline at end\n\
ACGTNNNNNNNCGTACTTGCTAGCTAGCTAGCTAGCTAGCATCGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTACTATCATCTACTATCTACTATCATCTACTATCATTCATCGATCGATCGATCGTACGTACGATCGATCGATCGATCGTACGATCGATGCTACGTACGTACG";


    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

std::fstream* createFastAFileProtein(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR =
"> sequenceID_with special chars an irregular linebreaks\n\
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGAT\n\
VITNLFSAIPYIGTNLV\n\
EWIWGGFSVDKATLNRFFAF\n\
HFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG\n\
>sequence2... with no linebreaks and no newline at end\n\
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX";


    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

std::fstream* createFastAFileAnnotatedProtein(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR =
"> sequenceID_with special chars an irregular linebreaks\n\
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGAT\n\
VITNLFSAIP[annotation+-*/]YIGTNLV\n\
EWIWGGFSVDKATLNRFFAF\n\
HFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG\n\
>sequence2... with no linebreaks and no newline at end\n\
GLMPFLHTSKHRSMMLRPLSQALFW[CTD(5)->gogogo]TLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX";


    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

std::fstream* createFastQFile(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR =
"@SEQ_ID\n\
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n\
+\n\
!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65\n\
@ 2ndSequence with formatting obscurities\n\
GATTTGGGGTTCAAAGC\n\
AGTATCGATCAAATAGTAAATCCATTT\n\
GTTCAACTCACAGTTT\n\
+ 2ndSequence with formatting obscurities\n\
!''*((((***+))%%%\n\
++)(%%%%).\n\
@***-+*''))**55CCF>>>>>>CCCCCCC65";
    // the second quality-string is line-broken at places to produce lines
    // beginning with '+' and '@'. The new code handles this absolutely fine!

    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

// something that is definitely not a sequence format
std::fstream* createBogusFile(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strncpy(filenameBuffer, toCString(tempFilename), 999);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR =
"asd45623\n\
-.-+´adfasddsad2342asfd8ss8g9s8g9s8fags9d8ga9s\n\
67utghjgnm\n\
!-.,()//&%$$%§§%\n\
e5bv56 u67 gdh g x#++++´\n\
..\n\
.....asdsadfasd0809809vnvbnojokjojokjokj\n\
";
    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

#endif  // TEST_STREAM_TEST_SEQ_IO_GENERIC_H_
