// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#define SEQAN_DEBUG
//#define SEQAN_TEST

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>  // Header under test.

#include "test_embl.h"
#include "test_file.h"

SEQAN_BEGIN_TESTSUITE(test_file)
{
    SEQAN_CALL_TEST(test_file_stream);
    SEQAN_CALL_TEST(test_file_cstream);

    SEQAN_CALL_TEST(test_file_raw);

    SEQAN_CALL_TEST(test_file_fasta_crlf);
    SEQAN_CALL_TEST(test_file_fasta_lf);
    SEQAN_CALL_TEST(test_file_fasta_cr);
    SEQAN_CALL_TEST(test_file_fasta_write);

    SEQAN_CALL_TEST(test_file_cgviz);
    SEQAN_CALL_TEST(test_file_fasta_align);

    SEQAN_CALL_TEST(test_file_embl);
    SEQAN_CALL_TEST(test_file_genbank);

//    SEQAN_CALL_TEST(test_file_reader_iterator);
    SEQAN_CALL_TEST(test_file_reader_string);
//    SEQAN_CALL_TEST(test_file_reader_string2_fasta);
//    SEQAN_CALL_TEST(test_file_reader_string2_embl);
//    SEQAN_CALL_TEST(test_file_reader_string2_genbank);
    SEQAN_CALL_TEST(test_file_reader_string3);

	SEQAN_CALL_TEST(test_file_embl_file);
	SEQAN_CALL_TEST(test_file_embl_meta);
}
SEQAN_END_TESTSUITE
