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
#define SEQAN_TEST

#include<seqan/seeds.h>


using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////


void test_banded_alignment()
{	
	//Needleman Wunsch
	String<char> query = "cgtacgtgagtga";
	String<char> database = "cgatta";
	Seed<int, SimpleSeed> seed(0,0,8,5);
	Score<int, Simple> scoreMatrix(3,-3,-2);
	//setRightDiagonal(seed, -4);
	Align<String<char>,ArrayGaps> alignment;
	Segment<String<char>, InfixSegment> seg1(query, leftDim0(seed), rightDim0(seed)+1);
	Segment<String<char>, InfixSegment> seg2(database, leftDim1(seed), rightDim1(seed)+1);
	resize(rows(alignment), 2);
	assignSource(row(alignment, 0), seg1);
	assignSource(row(alignment, 1), seg2);
	//bandedAlignment(alignment, seed, 1, scoreMatrix);
	SEQAN_ASSERT_EQ(bandedAlignment(alignment, seed, 1, scoreMatrix), 6);
	SEQAN_ASSERT(row(alignment,0) == "cgtacgtga" );
	SEQAN_ASSERT(row(alignment,1) == "cg-at-t-a");

	setLeftDiagonal(seed, 1);
	Align<String<char>,ArrayGaps> alignment1b;
	Segment<String<char>, InfixSegment> seg1x(query, leftDim0(seed), rightDim0(seed)+1);
	Segment<String<char>, InfixSegment> seg2x(database, leftDim1(seed), rightDim1(seed)+1);
	resize(rows(alignment1b), 2);
	assignSource(row(alignment1b, 0), seg1x);
	assignSource(row(alignment1b, 1), seg2x);
	SEQAN_ASSERT_EQ(bandedAlignment(alignment1b, seed, 1, scoreMatrix), 6);
	SEQAN_ASSERT(row(alignment1b,0) == "cgtacgtga" );
	SEQAN_ASSERT(row(alignment1b,1) == "cg-at-t-a");
	
	String<char> query2 = "ACTTTCATTTT";
	String<char> database2 = "ACTGTTCAGGG";
	Seed<int, SimpleSeed> seed2(0,0,3,4); 
	Score<int, Simple> scoreMatrix2(2,-1,-1);
	Align<String<char>,ArrayGaps> alignment2;
	Segment<String<char>, InfixSegment> seg1b(query2, leftDim0(seed2), rightDim0(seed2)+1);
	Segment<String<char>, InfixSegment> seg2b(database2, leftDim1(seed2), rightDim1(seed2)+1);
	resize(rows(alignment2), 2);
	assignSource(row(alignment2, 0), seg1b);
	assignSource(row(alignment2, 1), seg2b);
	SEQAN_ASSERT_EQ(bandedAlignment(alignment2, seed2, 2, scoreMatrix2), 7);
	SEQAN_ASSERT(row(alignment2,0) == "ACT-T" );
	SEQAN_ASSERT(row(alignment2,1) == "ACTGT");


	//Gotoh
	setRightDiagonal(seed, -3);
	query = "cgtacgtgagtga";
	database = "cgatta";
	Score<int, Simple> scoreMatrix3(3,-2,-1,-3);

	Align<String<char>,ArrayGaps> alignment3;
	Segment<String<char>, InfixSegment> seg1c(query, leftDim0(seed), rightDim0(seed)+1);
	Segment<String<char>, InfixSegment> seg2c(database, leftDim1(seed), rightDim1(seed)+1);
	resize(rows(alignment3), 2);
	assignSource(row(alignment3, 0), seg1c);
	assignSource(row(alignment3, 1), seg2c);
	SEQAN_ASSERT_EQ(bandedAlignment(alignment3, seed, 1, scoreMatrix3), 6);
	SEQAN_ASSERT(row(alignment3,0) == "cgtacgtga" );
	SEQAN_ASSERT(row(alignment3,1) == "cg-a--tta");

	setRightDiagonal(seed, -4);
	Align<String<char>,ArrayGaps> alignment3b;
	Segment<String<char>, InfixSegment> seg1y(query, leftDim0(seed), rightDim0(seed)+1);
	Segment<String<char>, InfixSegment> seg2y(database, leftDim1(seed), rightDim1(seed)+1);
	resize(rows(alignment3b), 2);
	assignSource(row(alignment3b, 0), seg1y);
	assignSource(row(alignment3b, 1), seg2y);
	SEQAN_ASSERT_EQ(bandedAlignment(alignment3b, seed, 1, scoreMatrix3), 6);
	SEQAN_ASSERT(row(alignment3b,0) == "cgtacgtga" );
	SEQAN_ASSERT(row(alignment3b,1) == "cg-a--tta");

	Align<String<char>,ArrayGaps> alignment4;
	Segment<String<char>, InfixSegment> seg1d(query2, leftDim0(seed2), rightDim0(seed2)+1);
	Segment<String<char>, InfixSegment> seg2d(database2, leftDim1(seed2), rightDim1(seed2)+1);
	resize(rows(alignment4), 2);
	assignSource(row(alignment4, 0), seg1d);
	assignSource(row(alignment4, 1), seg2d);
	SEQAN_ASSERT_EQ(bandedAlignment(alignment4, seed2,2, scoreMatrix3), 9);
	SEQAN_ASSERT(row(alignment4,0) == "ACT-T" );
	SEQAN_ASSERT(row(alignment4,1) == "ACTGT");

}

void test_banded_chain_align()
{
	String<char> query = "ACGTCCTCGTACACCGTCTTAA";
	String<char> database = "TACGATCCACACCGCGTCT";

	Score<int, Simple> scoreMatrix(2,-1,-2);

	String<Seed<int, SimpleSeed> > seedChain1;

	appendValue(seedChain1, Seed<int, SimpleSeed>(12,10,18,18));
	appendValue(seedChain1, Seed<int, SimpleSeed>(1,2,5,7));
	
	Align<String<char>,ArrayGaps> alignment1;
	resize(rows(alignment1), 2);
	assignSource(row(alignment1, 0), query);
	assignSource(row(alignment1, 1), database);

    //cout << "Score: " << bandedChainAlignment(seedChain1, 2, alignment1, scoreMatrix) << endl;
	SEQAN_ASSERT_EQ(bandedChainAlignment(seedChain1, 2, alignment1, scoreMatrix), 11);

	//cout << alignment1 << endl;
	SEQAN_ASSERT(row(alignment1,0) == "ACGTCCTCGTACACCGTCTTAA" );
	SEQAN_ASSERT(row(alignment1,1) == "TACGATC-C--ACACCG-CGTCT");


    //------ on infixes -> only parts of the source object are used ------------

	Align<String<char>,ArrayGaps> alignment2;
	resize(rows(alignment2), 2);
	assignSource(row(alignment2, 0), query);
	assignSource(row(alignment2, 1), database);
	setClippedBeginPosition(row(alignment2, 0), 1);
	setClippedBeginPosition(row(alignment2, 1), 2);
	setBeginPosition(row(alignment2, 0), 0);
	setBeginPosition(row(alignment2, 1), 0);
	setClippedEndPosition(row(alignment2, 0), length(query));
	setClippedEndPosition(row(alignment2, 1), length(database));

    //cout << "Score: " << bandedChainAlignment(seedChain1, 2, alignment2, scoreMatrix) << endl;
	SEQAN_ASSERT_EQ(bandedChainAlignment(seedChain1, 2, alignment2, scoreMatrix), 11);

	//cout << alignment2 << endl;
	SEQAN_ASSERT(row(alignment2,0) == "CGTCCTCGTACACCGTCTTAA" );
	SEQAN_ASSERT(row(alignment2,1) == "CGATC-C--ACACCG-CGTCT");


    // ------------------- affine gap costs --------------------------------

    Score<int, Simple> scoreMatrix2(3,-2,-1, -3);

	String<Seed<int, SimpleSeed> > seedChain2;

	appendValue(seedChain2, Seed<int, SimpleSeed>(12,10,18,18));
	appendValue(seedChain2, Seed<int, SimpleSeed>(1,2,5,7));
	
	Align<String<char>,ArrayGaps> alignment3;
	resize(rows(alignment3), 2);
	assignSource(row(alignment3, 0), query);
	assignSource(row(alignment3, 1), database);

	//cout << "Score: " << bandedChainAlignment(seedChain2, 2, alignment3, scoreMatrix2) << endl;
	SEQAN_ASSERT_EQ(bandedChainAlignment(seedChain2, 2, alignment3, scoreMatrix2), 24);

	//cout << alignment3 << endl;
	SEQAN_ASSERT(row(alignment3,0) == "ACG-TCCTCGTACAC--CGTCTTAA");
	SEQAN_ASSERT(row(alignment3,1) == "TACGATCC----ACACCGCGTCT");


	// -------------------- affine gap costs and on infixes ----------------
	
    Align<String<char>,ArrayGaps> alignment4;
	resize(rows(alignment4), 2);
	assignSource(row(alignment4, 0), query);
	assignSource(row(alignment4, 1), database);
	setClippedBeginPosition(row(alignment2, 0), 1);
	setClippedBeginPosition(row(alignment2, 1), 2);
	setBeginPosition(row(alignment2, 0), 0);
	setBeginPosition(row(alignment2, 1), 0);
	setClippedEndPosition(row(alignment2, 0), length(query));
	setClippedEndPosition(row(alignment2, 1), length(database));

	//cout << "Score: " << bandedChainAlignment(seedChain2, 2, alignment4, scoreMatrix2) << endl;
	SEQAN_ASSERT_EQ(bandedChainAlignment(seedChain2, 2, alignment4, scoreMatrix2), 21);

	//cout << alignment4 << endl;
	SEQAN_ASSERT(row(alignment4,0) == "CG-TCCTCGTACAC--CGTCTTAA");
	SEQAN_ASSERT(row(alignment4,1) == "CGATCC----ACACCGCGTCT");
	
}

void Main_BandedAlign(){
	//Fragment<unsigned int, Default()> x();

	test_banded_alignment();
	test_banded_chain_align();
}
