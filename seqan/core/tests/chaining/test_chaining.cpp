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

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>

//#define SEQAN_DEBUG
#define SEQAN_TEST
//#define SEQAN_NOSRAN //suppress srand

#include <seqan/sequence.h>
#include <seqan/chaining.h>

using namespace seqan;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

template< typename TContainer >
void
_generateRandomFrags(TContainer & dest,
					 int num,
					 int min,
					 int max,
					 int minwidth,
					 int maxwidth,
					 int dim )
{
	typedef typename Value< TContainer >::Type FragType;
	typename Key< FragType >::Type * left_pos;
	typename Key< FragType >::Type * right_pos;
	reserve( dest, num );
	double d_dim = static_cast< double >( dim );
	for( int i = 0; i < num; ++i )
	{
		
		double width_sum = 0;
		left_pos = new typename Key< FragType >::Type[ dim ];
		right_pos = new typename Key< FragType >::Type[ dim ];
		for( int d = 0; d < dim; ++d )
		{
			left_pos[ d ] = ( rand() % ( max - min ) ) + min;
			int width = ( rand() % ( maxwidth - minwidth ) )+ minwidth;
			width_sum += width;
			right_pos[ d ] = left_pos[ d ] + width;			
			
		}
		FragType frag( left_pos, right_pos, dim, static_cast< typename Weight< FragType >::Type >(exp(log(width_sum)/d_dim)) * 100 );

		delete[] left_pos;
		delete[] right_pos;

		appendValue( dest, frag );
	}
}

//////////////////////////////////////////////////////////////////////////////
//helper function

template <typename TChain, typename TScoring>
void _showChain(TChain & ch,
				TScoring scoring)
{
	typedef typename Iterator<TChain>::Type TIterator;
	typedef typename Value<TChain>::Type TFragment;

	TIterator it = begin(ch);
	TIterator it_end = end(ch);
	unsigned int dim = dimension(*it); //chain must not be empty!
	while (it < it_end)
	{
		TFragment & frag = *it;
		printf("%4i: ", weight(frag));
		cout << "(";
		for (unsigned int i = 0; i < dim; ++i)
		{
			printf("%2i", leftPosition(frag, i));
			if (i < (dim-1)) cout << ", ";
		}
		cout << ")(";
		for (unsigned int i = 0; i < dim; ++i)
		{
			printf("%2i", rightPosition(frag, i));
			if (i < (dim-1)) cout << ", ";
		}
		cout << ")\n";

		++it;
		if (it == it_end) break;

		printf("%4i\n", scoreChainGap(scoring, frag, *it));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoring>
void testChainer(int count, 
				 int dim,
				 TScoring scoring)
{
	String< Seed<int, MultiSeed> > fragments;
	reserve(fragments, count);
	
	_generateRandomFrags(fragments, count, 1, 3 * count, 10, 20, dim);
//_showChain(fragments, scoring);

	String< Seed<int, MultiSeed> > ch;
	reserve(ch, count+2);


	//build chain
	int chain_score = globalChaining(fragments, ch, scoring);

    //std::cout << chain_score << "\n";
//_showChain(ch, scoring);

	//verify validity of chain
    SEQAN_ASSERT_GT(length(ch), 0u);
	int sum = weight(ch[0]);
	for (unsigned int i = 1; i < length(ch); ++i)
	{
		SEQAN_ASSERT(_chainGenericChainable(ch[i-1], ch[i]));
		sum += scoreChainGap(scoring, ch[i-1], ch[i]) + weight(ch[i]);
	}
	//verify score of chain
	SEQAN_ASSERT_EQ(sum, chain_score);

	//build generic chain
	int chain_score2 = globalChaining(fragments, ch, scoring, GenericChaining());
    //std::cout << chain_score2 << "\n";
//_showChain(ch, scoring);

	//verify validity of generic chain
    SEQAN_ASSERT_GT(length(ch), 0u);
	sum = weight(ch[0]);
	for (unsigned int i = 1; i < length(ch); ++i) {
            SEQAN_ASSERT(_chainGenericChainable(ch[i-1], ch[i]));
            sum += scoreChainGap(scoring, ch[i-1], ch[i]) + weight(ch[i]);
	}
	//verify score of generic chain
	SEQAN_ASSERT_EQ(sum, chain_score2);

	//compare results of two chaining algorithms
	SEQAN_ASSERT_EQ(chain_score2, chain_score);
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_chaining_test_chainer_zero_score) {
    testChainer(1000, 2, Score<int, Zero>());
}


SEQAN_DEFINE_TEST(test_chaining_test_chainer_manhattan_score) {
    testChainer(1000, 2, Score<int, Manhattan>());
}


SEQAN_DEFINE_TEST(test_chaining_test_chainer_chain_sop_score) {
    testChainer(1000, 2, Score<int, ChainSoP>());
}

SEQAN_BEGIN_TESTSUITE(test_chaining) {
    SEQAN_CALL_TEST(test_chaining_test_chainer_zero_score);
    SEQAN_CALL_TEST(test_chaining_test_chainer_manhattan_score);
//    SEQAN_CALL_TEST(test_chaining_test_chainer_chain_sop_score);

    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/chain_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/chain_generic.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/chain_meta_fragment.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/chain_point.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/chain_wrapper_point.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/fragment.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/geom_distribution.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/range_max_tree.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/range_tree.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rmt_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rmt_common_algos.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rmt_compl_algos.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rmt_def_algos.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rmt_skip_base_element.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rmt_skip_element.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_common_algos.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_impl.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_skip_base_element.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_skip_element.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_sl_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_sl_compl_algos.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_sl_def_algos.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/rt_sl_impl.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/score_chain.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/score_chain_sop.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/score_manhattan.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/score_zero.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_base_element.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_element.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_list.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_list_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_list_dynamic.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_list_impl.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_list_iterator.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_list_type.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/skip_pool_alloc.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/tree_chain.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/tree_chain_sop.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/chaining/tree_chain_utils.h");
}
SEQAN_END_TESTSUITE
