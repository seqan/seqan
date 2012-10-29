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


#include <seqan/basic.h>
#include <seqan/seeds.h>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

//template <typename TAlgorithmSpec>
void Test_SimpleSeeds()
{

//____________________________________________________________________________
// Standard Functions
	Seed<int,SimpleSeed> seed;
	Seed<int, SimpleSeed> seed1(0,0,7);
	Seed<int, SimpleSeed> seed2(0,1,7,5);

	SEQAN_ASSERT_EQ(startDiagonal(seed1), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed1), 0);
	SEQAN_ASSERT_EQ(startDiagonal(seed2), 1);
	SEQAN_ASSERT_EQ(endDiagonal(seed2), -2);
	SEQAN_ASSERT_EQ(leftDim0(seed2), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed2), 7);
	SEQAN_ASSERT_EQ(leftDim1(seed2), 1);
	SEQAN_ASSERT_EQ(rightDim1(seed2), 5);
	SEQAN_ASSERT_EQ(leftDiagonal(seed2), 1);
	SEQAN_ASSERT_EQ(rightDiagonal(seed2), -2);
	SEQAN_ASSERT_EQ(length(seed2), 8);

	SEQAN_ASSERT_EQ(leftPosition(seed2,0), 0);
	SEQAN_ASSERT_EQ(leftPosition(seed2,1), 1);
	SEQAN_ASSERT_EQ(rightPosition(seed2,0), 7);
	SEQAN_ASSERT_EQ(rightPosition(seed2,1), 5);

	setLeftDim0(seed2,3);
	setRightDim0(seed2,9);
	setLeftDim1(seed2,5);
	setRightDim1(seed2,12);
	setLeftDiagonal(seed2,29);
	setRightDiagonal(seed2,7);
	SEQAN_ASSERT_EQ(leftDim0(seed2), 3);
	SEQAN_ASSERT_EQ(rightDim0(seed2), 9);
	SEQAN_ASSERT_EQ(leftDim1(seed2), 5);
	SEQAN_ASSERT_EQ(rightDim1(seed2), 12);
	SEQAN_ASSERT_EQ(leftDiagonal(seed2), 29);
	SEQAN_ASSERT_EQ(rightDiagonal(seed2), 7);

//____________________________________________________________________________
// Merge Algorithms
	Seed<int, SimpleSeed> seed3(0,0,7);
	Seed<int, SimpleSeed> seed4(4,5,7);
	_mergeTwoSeeds(seed3,seed4,Merge());
	SEQAN_ASSERT_EQ(startDiagonal(seed3), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed3), 1);
	SEQAN_ASSERT_EQ(startDiagonal(seed3), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed3), 1);
	SEQAN_ASSERT_EQ(leftDim0(seed3), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed3), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed3), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed3), 11);
	SEQAN_ASSERT_EQ(leftDiagonal(seed3), 1);
	SEQAN_ASSERT_EQ(rightDiagonal(seed3), 0);
	SEQAN_ASSERT_EQ(length(seed3), 11);

	
	Seed<int, SimpleSeed> seed5(0,0,7);
	_mergeTwoSeeds(seed5,4,5,7,Merge());
	SEQAN_ASSERT_EQ(startDiagonal(seed5), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed5), 1);
	SEQAN_ASSERT_EQ(startDiagonal(seed5), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed5), 1);
	SEQAN_ASSERT_EQ(leftDim0(seed5), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed5), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed5), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed5), 11);
	SEQAN_ASSERT_EQ(leftDiagonal(seed5), 1);
	SEQAN_ASSERT_EQ(rightDiagonal(seed5), 0);
	SEQAN_ASSERT_EQ(length(seed5), 11);

	Seed<int, SimpleSeed> seed6(0,0,7);
	_mergeTwoSeeds(seed6,4,5,10,11,Merge());
	SEQAN_ASSERT_EQ(startDiagonal(seed6), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed6), 1);
	SEQAN_ASSERT_EQ(startDiagonal(seed6), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed6), 1);
	SEQAN_ASSERT_EQ(leftDim0(seed6), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed6), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed6), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed6), 11);
	SEQAN_ASSERT_EQ(leftDiagonal(seed6), 1);
	SEQAN_ASSERT_EQ(rightDiagonal(seed6), 0);
	SEQAN_ASSERT_EQ(length(seed6), 11);

	Score<int,Simple> matrix(2,-1,-1);
	Seed<int, SimpleSeed> seed10(0,0,7);
	int score10 = 14;
	Seed<int, SimpleSeed> seed11(0,0,7);
	int score11 = 14;
	Seed<int, SimpleSeed> seed12(0,0,7);
	int score12 = 14;
	Seed<int, SimpleSeed> seed13(4,5,7);
	int score13 = 14;


	_mergeTwoSeedsScore(seed10, score10, seed13, score13, matrix, Manhattan(), Merge());
	SEQAN_ASSERT_EQ(leftDim0(seed10), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed10), 10);
	SEQAN_ASSERT_EQ(score10, 20);
	_mergeTwoSeedsScore(seed11, score11, 4, 5, 7, score13, matrix, Manhattan(), Merge());
	SEQAN_ASSERT_EQ(leftDim0(seed11), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed11), 10);
	SEQAN_ASSERT_EQ(score11, 20);
	_mergeTwoSeedsScore(seed12, score12, 4, 5, 10, 11, score13, matrix, Manhattan(), Merge());
	SEQAN_ASSERT_EQ(leftDim0(seed12), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed12), 10);
	SEQAN_ASSERT_EQ(score12, 20);

//____________________________________________________________________________
// Extension Algorithms



	String<Dna> query =	   "AAACCCTTTGGGTTTTT";
	String<Dna> database = "AACCCCTTTGGTGAAAAA";
	Seed<int, SimpleSeed> seed7(4,4,3);
	extendSeed(seed7, query, database, 2, MatchExtend());
	SEQAN_ASSERT_EQ(leftDim0(seed7), 3);
	SEQAN_ASSERT_EQ(rightDim0(seed7), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed7), 3);
	SEQAN_ASSERT_EQ(rightDim1(seed7), 10);

	Seed<int, SimpleSeed> seed8(4,4,3);
	extendSeed(seed8, 2, matrix, query, database, 2, UngappedXDrop());
	SEQAN_ASSERT_EQ(leftDim0(seed8), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed8), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed8), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed8), 10);

	Seed<int, SimpleSeed> seed9(4,4,3);
	extendSeed(seed9, 1, matrix, query, database, 2, GappedXDrop());
    //std::cout << infix(query, leftDim0(seed9), rightDim0(seed9)+1) << std::endl;
    //std::cout << infix(database, leftDim1(seed9), rightDim1(seed9)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed9), 12);
	SEQAN_ASSERT_EQ(leftDim1(seed9), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed9), 13);
    
	String<Dna> query1 =	"aaaacgatcgatgc";
	String<Dna> database1 = "ttttcgatcgatgcttttt";
	Seed<int, SimpleSeed> seed9a(4,4,9);
	extendSeed(seed9a, 1, matrix, query1, database1, 2, GappedXDrop());
    //std::cout << infix(query1, leftDim0(seed9a), rightDim0(seed9a)+1) << std::endl;
    //std::cout << infix(database1, leftDim1(seed9a), rightDim1(seed9a)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9a), 3);
	SEQAN_ASSERT_EQ(rightDim0(seed9a), 13);
	SEQAN_ASSERT_EQ(leftDim1(seed9a), 3);
	SEQAN_ASSERT_EQ(rightDim1(seed9a), 14);
    SEQAN_ASSERT_EQ(rightDiagonal(seed9a), 0);
    SEQAN_ASSERT_EQ(leftDiagonal(seed9a), 2);

	Seed<int, SimpleSeed> seed9b(5,5,7);
	extendSeed(seed9b, 1, matrix, query1, database1, 2, GappedXDrop());
    //std::cout << infix(query1, leftDim0(seed9b), rightDim0(seed9b)+1) << std::endl;
    //std::cout << infix(database1, leftDim1(seed9b), rightDim1(seed9b)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9b), 3);
	SEQAN_ASSERT_EQ(rightDim0(seed9b), 13);
	SEQAN_ASSERT_EQ(leftDim1(seed9b), 3);
	SEQAN_ASSERT_EQ(rightDim1(seed9b), 14);

	Seed<int, SimpleSeed> seed9c(6,6,6);
	extendSeed(seed9c, 0, matrix, query1, database1, 2, GappedXDrop());
    //std::cout << infix(query1, leftDim0(seed9c), rightDim0(seed9c)+1) << std::endl;
    //std::cout << infix(database1, leftDim1(seed9c), rightDim1(seed9c)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9c), 4);
	SEQAN_ASSERT_EQ(rightDim0(seed9c), 13);
	SEQAN_ASSERT_EQ(leftDim1(seed9c), 4);
	SEQAN_ASSERT_EQ(rightDim1(seed9c), 13);
    
	String<Dna> query2 =        "cgatcgatgcaaaaaaaaa";
	String<Dna> database2 = "ttttcgatcgatgc";
	Seed<int, SimpleSeed> seed9d(1,5,5);
	extendSeed(seed9d, 1, matrix, query2, database2, 2, GappedXDrop());
    //std::cout << infix(query2, leftDim0(seed9d), rightDim0(seed9d)+1) << std::endl;
    //std::cout << infix(database2, leftDim1(seed9d), rightDim1(seed9d)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9d), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed9d), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed9d), 3);
	SEQAN_ASSERT_EQ(rightDim1(seed9d), 13);
    SEQAN_ASSERT_EQ(rightDiagonal(seed9d), 2);
    SEQAN_ASSERT_EQ(leftDiagonal(seed9d), 5);
    
	String<Dna> query3 = "aaaaaaaaacgatcgatgcaaaaaaaaa";
	String<Dna> database3 =       "cgatcgatgccaact";
	Seed<int, SimpleSeed> seed9e(11,2,4);
	extendSeed(seed9e, 1, matrix, query3, database3, 2, GappedXDrop());
    //std::cout << infix(query3, leftDim0(seed9e), rightDim0(seed9e)+1) << std::endl;
    //std::cout << infix(database3, leftDim1(seed9e), rightDim1(seed9e)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9e), 8);
	SEQAN_ASSERT_EQ(rightDim0(seed9e), 22);
	SEQAN_ASSERT_EQ(leftDim1(seed9e), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed9e), 13);
    SEQAN_ASSERT_EQ(rightDiagonal(seed9e), -10);
    SEQAN_ASSERT_EQ(leftDiagonal(seed9e), -7);
        
	String<Dna> query4 =      "ttttagtgacgttttaaaaaa";
	String<Dna> database4 = "ccccagctgatcgtttgcccccc";
	Seed<int, SimpleSeed> seed9f(9,11,5);
	Score<int,Simple> matrix1(1,-5,-5);
	extendSeed(seed9f, 7, matrix1, query4, database4, 2, GappedXDrop());
    //std::cout << infix(query4, leftDim0(seed9f), rightDim0(seed9f)+1) << std::endl;
    //std::cout << infix(database4, leftDim1(seed9f), rightDim1(seed9f)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9f), 4);
	SEQAN_ASSERT_EQ(rightDim0(seed9f), 14);
	SEQAN_ASSERT_EQ(leftDim1(seed9f), 4);
	SEQAN_ASSERT_EQ(rightDim1(seed9f), 16);
    SEQAN_ASSERT_EQ(rightDiagonal(seed9f), 0);
    SEQAN_ASSERT_EQ(leftDiagonal(seed9f), 2);
        
	String<Dna> query5 =    "aaaaaattttgcagtgatttt";
	String<Dna> database5 = "ccccccgtttgctagtcgacccc";
	Seed<int, SimpleSeed> seed9g(7,7,5);
	extendSeed(seed9g, 7, matrix1, query5, database5, 2, GappedXDrop());
    //std::cout << infix(query5, leftDim0(seed9g), rightDim0(seed9g)+1) << std::endl;
    //std::cout << infix(database5, leftDim1(seed9g), rightDim1(seed9g)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9g), 6);
	SEQAN_ASSERT_EQ(rightDim0(seed9g), 16);
	SEQAN_ASSERT_EQ(leftDim1(seed9g), 6);
	SEQAN_ASSERT_EQ(rightDim1(seed9g), 18);
    SEQAN_ASSERT_EQ(rightDiagonal(seed9g), 0);
    SEQAN_ASSERT_EQ(leftDiagonal(seed9g), 2);
        
	String<Dna> query6 =  "ccccagctgatcgtttgcccccc";
	String<Dna> database6 = "ttttagtgacgttttaaaaaa";
	Seed<int, SimpleSeed> seed9h(11,9,5);
	extendSeed(seed9h, 7, matrix1, query6, database6, 2, GappedXDrop());
    //std::cout << infix(query6, leftDim0(seed9h), rightDim0(seed9h)+1) << std::endl;
    //std::cout << infix(database6, leftDim1(seed9h), rightDim1(seed9h)+1) << std::endl << std::endl;
	SEQAN_ASSERT_EQ(leftDim0(seed9h), 4);
	SEQAN_ASSERT_EQ(rightDim0(seed9h), 16);
	SEQAN_ASSERT_EQ(leftDim1(seed9h), 4);
	SEQAN_ASSERT_EQ(rightDim1(seed9h), 14);
    SEQAN_ASSERT_EQ(rightDiagonal(seed9h), -2);
    SEQAN_ASSERT_EQ(leftDiagonal(seed9h), 0);

	String<Dna> query7 =  "CTGACAGCGAGAAACAGTAACCAGCTAGCCT";
	String<Dna> database7 = "AGCAGCAAACAGTAAGCCAGCAGCCT";
	Seed<int, SimpleSeed> seed15(26, 21, 5);
	extendSeed(seed15, 45, Score<int, Simple>(1, -9, -9, -9), query7, database7, 2, GappedXDrop());
	SEQAN_ASSERT_EQ(leftDim0(seed15), 2);
	SEQAN_ASSERT_EQ(leftDim1(seed15), 0);
}

void Test_MultiSeeds(){

//____________________________________________________________________________
// Standard Functions
	Seed<int,ChainedSeed> seed;
	Seed<int, ChainedSeed> seed1(4,5,7);
	SEQAN_ASSERT_EQ(startDiagonal(seed1), 1);
	SEQAN_ASSERT_EQ(endDiagonal(seed1), 1);
	SEQAN_ASSERT_EQ(leftDim0(seed1), 4);
	SEQAN_ASSERT_EQ(rightDim0(seed1), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed1), 5);
	SEQAN_ASSERT_EQ(rightDim1(seed1), 11);
	SEQAN_ASSERT_EQ(leftDiagonal(seed1), 1);
	SEQAN_ASSERT_EQ(rightDiagonal(seed1), 1);
	SEQAN_ASSERT_EQ(length(seed1), 7);
	SEQAN_ASSERT_EQ(_getFirstDiag(seed1).i2, 5);

	Seed<int, ChainedSeed> const seed132(4,5,7);
	SEQAN_ASSERT_EQ(_getLastDiag(seed132).i1, 4);


	
	Seed<int,ChainedSeed> seed2(0,0,0);
	setLeftDim0(seed2,3);
	setRightDim0(seed2,9);
	setLeftDim1(seed2,2);
	setRightDim1(seed2,12);
	setLeftDiagonal(seed2,29);
	setRightDiagonal(seed2,7);
	SEQAN_ASSERT_EQ(leftDim0(seed2), 2);
	SEQAN_ASSERT_EQ(rightDim0(seed2), 12);
	SEQAN_ASSERT_EQ(leftDim1(seed2), 2);
	SEQAN_ASSERT_EQ(rightDim1(seed2), 12);
	SEQAN_ASSERT_EQ(leftDiagonal(seed2), 29);
	SEQAN_ASSERT_EQ(rightDiagonal(seed2), 7);

//____________________________________________________________________________
// Merge Algorithms
	Seed<int, ChainedSeed> seed3(0,0,7);
	Seed<int, ChainedSeed> seed4(4,5,7);
	_mergeTwoSeeds(seed3,seed4,Merge());
	SEQAN_ASSERT_EQ(startDiagonal(seed3), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed3), 1);
	SEQAN_ASSERT_EQ(startDiagonal(seed3), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed3), 1);
	SEQAN_ASSERT_EQ(leftDim0(seed3), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed3), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed3), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed3), 11);
	SEQAN_ASSERT_EQ(leftDiagonal(seed3), 1);
	SEQAN_ASSERT_EQ(rightDiagonal(seed3), 0);
	SEQAN_ASSERT_EQ(length(seed3), 11);

	
	Seed<int, ChainedSeed> seed5(0,0,7);
	_mergeTwoSeeds(seed5,4,5,7,Merge());
	SEQAN_ASSERT_EQ(startDiagonal(seed5), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed5), 1);
	SEQAN_ASSERT_EQ(startDiagonal(seed5), 0);
	SEQAN_ASSERT_EQ(endDiagonal(seed5), 1);
	SEQAN_ASSERT_EQ(leftDim0(seed5), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed5), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed5), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed5), 11);
	SEQAN_ASSERT_EQ(leftDiagonal(seed5), 1);
	SEQAN_ASSERT_EQ(rightDiagonal(seed5), 0);
	SEQAN_ASSERT_EQ(length(seed5), 11);


	Score<int,Simple> matrix(2,-1,-1);

    Seed<int, ChainedSeed> seed10(0,0,7);
	int score10 = 14;
	Seed<int, ChainedSeed> seed11(0,0,7);
	int score11 = 14;
	Seed<int, ChainedSeed> seed12(4,5,7);
	int score12 = 14;


	_mergeTwoSeedsScore(seed10, score10, 4, 5, 7, 14, matrix, Manhattan(), Merge());
	_mergeTwoSeedsScore(seed10, score10, 3, 4, 10, 20, matrix, Manhattan(), Merge());
    SEQAN_ASSERT_EQ(leftDim0(seed10), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed10), 12);
	SEQAN_ASSERT_EQ(score10, 26);
	
	_mergeTwoSeedsScore(seed11, score11, seed12, score12, matrix, Manhattan(), Merge());
	SEQAN_ASSERT_EQ(leftDim0(seed11), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed11), 10);
	SEQAN_ASSERT_EQ(score11, 20);

//____________________________________________________________________________
// Extension Algorithms

	String<Dna> query =	   "AAACCCTTTGGGTTTTT";
	String<Dna> database = "AACCCCTTTGGTGAAAAA";
	Seed<int, ChainedSeed> seed7(4,4,3);
	extendSeed(seed7, query, database, 2, MatchExtend());
	SEQAN_ASSERT_EQ(leftDim0(seed7), 3);
	SEQAN_ASSERT_EQ(rightDim0(seed7), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed7), 3);
	SEQAN_ASSERT_EQ(rightDim1(seed7), 10);

	Seed<int, ChainedSeed> seed8(4,4,3);
	extendSeed(seed8, 2, matrix, query, database, 2, UngappedXDrop());
	SEQAN_ASSERT_EQ(leftDim0(seed8), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed8), 10);
	SEQAN_ASSERT_EQ(leftDim1(seed8), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed8), 10);

	Seed<int, ChainedSeed> seed9(4,4,3);
	extendSeed(seed9, 1, matrix, query, database, 2, GappedXDrop());
	SEQAN_ASSERT_EQ(leftDim0(seed9), 0);
	SEQAN_ASSERT_EQ(rightDim0(seed9), 12);
	SEQAN_ASSERT_EQ(leftDim1(seed9), 0);
	SEQAN_ASSERT_EQ(rightDim1(seed9), 13);

//____________________________________________________________________________
// Alignment Calculation
    // TODO(holtgrew): Fix? Directly goes into ArrayGaps, would be rewrite but worth it, given that we have seeds2 soon?
    /*DISABLED since it breaks with new align module
	Align<String<Dna>, ArrayGaps> aligned;
	getAlignment(seed8, aligned, query, database, matrix);
	SEQAN_ASSERT(row(aligned, 0) == "AAACCCTTTGG");
	SEQAN_ASSERT(row(aligned, 1) == "AACCCCTTTGG");
	*/
	
	
//____________________________________________________________________________
// Score Calculation
	SEQAN_ASSERT_EQ(scoreSeed(seed8, query, database, matrix), 19);
}

void Main_Seeds(){
	Test_SimpleSeeds();
	Test_MultiSeeds();
}
