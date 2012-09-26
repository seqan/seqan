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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the SeqAn model store, I/O functionality.
// ==========================================================================

#include <seqan/basic.h>  // For test functionality.
#include <seqan/store.h>  // Header under test.

using namespace seqan;

SEQAN_DEFINE_TEST(test_store_io_sam) 
{
	FragmentStore<> store;
    char buffer[1023];

	// 1. LOAD CONTIGS
    strcpy(buffer, SEQAN_PATH_TO_ROOT());
    strcat(buffer, "/core/tests/store/ex1.fa");
    
	loadContigs(store, buffer);

	// 2. LOAD SAM ALIGNMENTS
    strcpy(buffer, SEQAN_PATH_TO_ROOT());
    strcat(buffer, "/core/tests/store/ex1.sam.copy");
	MultiSeqFile sam1;
	open(sam1.concat, buffer);
	split(sam1, Raw());
    
	{
		// read reference Sam from file
		std::ifstream samFile(buffer);
		SEQAN_ASSERT(samFile);
		read(samFile, store, Sam());
	}
	
	// 3. WRITE SAM ALIGNMENTS
    strcpy(buffer, SEQAN_TEMP_FILENAME());
	{
		// write Sam to temp file
		std::ofstream samFileOut(buffer);
		SEQAN_ASSERT(samFileOut);
		write(samFileOut, store, Sam());
	}
	
	// 4. COMPARE BOTH SAM FILES
	MultiSeqFile sam2;
	open(sam2.concat, buffer);
	split(sam2, Raw());
	
	SEQAN_ASSERT(!empty(sam1));
	SEQAN_ASSERT(!empty(sam2));
	for (unsigned i = 0; i < length(sam1); ++i)
	{
		if (sam1[i] != sam2[i])
		{
			std::cout << "    \t" << sam1[i] << std::endl;
			std::cout << " != \t" << sam2[i] << std::endl;
			SEQAN_ASSERT_FAIL("Files differ in line %d.", i);
		}
	}
}
