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

#define SEQAN_DEBUG


#include<seqan/seeds.h>





using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////


void test_global_seed_chain()
{
	typedef Seed<int, SimpleSeed> TSeed;
	String<TSeed> chain;
	String<TSeed> chain2;

	SeedSet<int, SimpleSeed, DefaultScore, void> seedContainer2(1,1);
	addSeed(seedContainer2,1,1,2,3,Single());
	addSeed(seedContainer2,2,5,2,5,Single());
	addSeed(seedContainer2,3,9,2,3,Single());
	addSeed(seedContainer2,5,2,3,4,Single());
	addSeed(seedContainer2,5,8,3,4,Single());
	addSeed(seedContainer2,6,5,1,2,Single());
	addSeed(seedContainer2,8,6,1,2,Single());
	addSeed(seedContainer2,9,1,3,3,Single());
	addSeed(seedContainer2,10,4,2,3,Single());
	addSeed(seedContainer2,10,6,1,2,Single());
	addSeed(seedContainer2,10,8,2,3,Single());

	SEQAN_ASSERT_EQ(globalChaining(seedContainer2, chain), 10);
	SEQAN_ASSERT_EQ(length(chain), 4u); 

	SEQAN_ASSERT_EQ(globalChaining(seedContainer2, chain2, -1 ,13, 13), -4);

	SEQAN_ASSERT_EQ(length(chain2), 4u); 

	String<TSeed> chain3;
	String<TSeed> chain4;

	SeedSet<int, SimpleSeed, DefaultScore, void> seedContainer3(1,1);
	addSeed(seedContainer3,1,1,2,3,Single());

	SEQAN_ASSERT_EQ(globalChaining(seedContainer3, chain3), 3);
	SEQAN_ASSERT_EQ(length(chain3), 1u); 
}


void Main_GlobalSeedChain(){
	// test_global_seed_chain();
}
