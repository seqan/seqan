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
// Author: Anne-Katrin Emde <emde@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the SeqAn moduel "refinement".
// ==========================================================================

#define SEQAN_DEBUG

#define SEQAN_VERBOSE

// Test path
#define TEST_PATH "projects/tests/refinement/"
#define LIB_PATH "core/include/seqan/refinement/"

// SeqAn Includes
#include <seqan/refinement.h>
#include <seqan/basic.h>

// Test files
#include "test_graph_impl_align.h"
#include "test_graph_match_refinement.h"

using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_refinement)
{
    // Test AlignmentGraph.
	 SEQAN_CALL_TEST(Test_Refinement_AlignmentGraphNoEdgeWeights);
	 SEQAN_CALL_TEST(Test_Refinement_AlignmentGraphEdgeWeights);
	 SEQAN_CALL_TEST(Test_Refinement_AlignmentGraphIterators);
	 SEQAN_CALL_TEST(Test_Refinement_AlignmentGraphOutput);
     SEQAN_CALL_TEST(Test_Refinement_HeaviestCommonSubsequence);
     SEQAN_CALL_TEST(Test_Refinement_OutEdgeIteratorAlignment);

    // Test Match Refinement.
    SEQAN_CALL_TEST(RefineMatchesSelfEdges);

    //SEQAN_CALL_TEST(GraphMatchRefine);
    SEQAN_CALL_TEST(RefineAlign);
    SEQAN_CALL_TEST(RefineInexactFragment);
}
SEQAN_END_TESTSUITE
