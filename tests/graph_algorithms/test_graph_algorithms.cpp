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

#define SEQAN_DEBUG
//#define SEQAN_TEST
#define SEQAN_VERBOSE

// SeqAn
#include <seqan/graph_algorithms.h>
#include "test_graph_algorithms.h"

// SeqAn Namespace
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_graph_algorithms)
{
	SEQAN_CALL_TEST(test_heap_tree);
	SEQAN_CALL_TEST(test_breadth_first_search);
	SEQAN_CALL_TEST(test_depth_first_search);
	SEQAN_CALL_TEST(test_topological_sort);
	SEQAN_CALL_TEST(test_strongly_connected_components);
	SEQAN_CALL_TEST(test_connected_components);
	SEQAN_CALL_TEST(test_prims_algorithm);
	SEQAN_CALL_TEST(test_union_find);	
	SEQAN_CALL_TEST(test_kruskals_algorithm);
	SEQAN_CALL_TEST(test_mst_all);
	SEQAN_CALL_TEST(test_dag_shortest_path);
	SEQAN_CALL_TEST(test_bellmann_ford);
	SEQAN_CALL_TEST(test_dijkstra);
	SEQAN_CALL_TEST(test_all_pairs_shortest_path);
	SEQAN_CALL_TEST(test_floyd_warshall);
	SEQAN_CALL_TEST(test_transitive_closure);
	SEQAN_CALL_TEST(test_ford_fulkerson);
	SEQAN_CALL_TEST(test_path_growing_algorithm);
	SEQAN_CALL_TEST(test_longest_increasing_subsequence);
	SEQAN_CALL_TEST(test_longest_common_subsequence);
	SEQAN_CALL_TEST(test_heaviest_increasing_subsequence);
	SEQAN_CALL_TEST(test_hmm_algorithm);	
}
SEQAN_END_TESTSUITE

