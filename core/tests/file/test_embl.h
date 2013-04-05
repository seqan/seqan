// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_FILE_TEST_FILE_EMBL_H_
#define TESTS_FILE_TEST_FILE_EMBL_H_

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
SEQAN_DEFINE_TEST(test_file_embl_file)
{
    char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_ROOT());
    strcat(buffer, "/core/tests/file/takifugu_scl_embl.txt");

	std::fstream strm; 
	strm.open(buffer, ios_base::in | ios_base::binary);

	String<char> line;
	String<char> feature_line;

	readLineType(strm, feature_line, "FT", Embl());
	//cout << feature_line << "\n";

	int count = 0;
	int next_pos = readFeature(feature_line, 0, line, "exon", Embl());
	while(next_pos != 0)
	{
		++count;
	//  cout << line << "\n";
		next_pos = readFeature(feature_line, next_pos, line, "exon", Embl());
	}

	SEQAN_ASSERT_EQ(count, 3);
}


SEQAN_DEFINE_TEST(test_file_embl_meta)
{
    char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_ROOT());
    strcat(buffer, "/core/tests/file/takifugu_scl_embl.txt");
    
	std::fstream strm; 
	strm.open(buffer, ios_base::in | ios_base::binary);

	String<char> line;

	String<char> feature_line;
	String<char> meta;
	readMeta(strm,meta,Embl());
	
	readLineType(meta, line, "KW", Embl());
	SEQAN_ASSERT_EQ(line, "SCL gene.");

	readLineType(meta, line, "RX", Embl());
	SEQAN_ASSERT_EQ(infix(line,0,28), "DOI; 10.1073/pnas.101532998.");
	SEQAN_ASSERT(length(line) == 46u || length(line) == 47u);

	clear(line);
	readLineType(meta, feature_line, "FT", Embl());
	//cout << feature_line << "\n";

	int count = 0;
	int next_pos = readFeature(feature_line, 0, line, "CDS", Embl());
	while(next_pos != 0)
	{
		++count;
	//  cout << line << "\n";
		next_pos = readFeature(feature_line, next_pos, line, "CDS", Embl());
	}

	SEQAN_ASSERT_EQ(count, 1);
}

#endif  // TESTS_FILE_TEST_FILE_EMBL_H_
