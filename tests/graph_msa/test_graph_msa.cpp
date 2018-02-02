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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: Tobias Rausch <rausch@embl.de>
// ===========================================================================
// Tests for the Graph MSA module.
// ===========================================================================

// External / STL
#include <iostream>
#include <fstream>
#include <string>

// SeqAn
#include <seqan/basic.h>
#include <seqan/stream.h>

#include <seqan/graph_msa.h>

// Test files
#include "test_graph_msa_guide_tree.h"
#include "test_graph_tcoffee.h"

SEQAN_BEGIN_TESTSUITE(test_graph_msa)
{
    // Call Tests.
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_neighbour_joining);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_weight_avg);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_avg);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_min);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_max);

    SEQAN_CALL_TEST(test_distances);
	SEQAN_CALL_TEST(test_libraries);
	SEQAN_CALL_TEST(test_external_libraries);
	SEQAN_CALL_TEST(test_triplet_extension);
	SEQAN_CALL_TEST(test_sop);
	SEQAN_CALL_TEST(test_progressive);
	SEQAN_CALL_TEST(test_reversable_fragments);	
}
SEQAN_END_TESTSUITE

