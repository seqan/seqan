#include <iostream>
#include <seqan/seeds.h>
#include <seqan/file.h>

using namespace seqan;

//define some constants
int const gaps_max = 200; //minimal sequence length for chaining
int const q_max = 13;     //start value for q
int const q_min = 9;      //minimal q
int const limit = 20;     //local seed chaining limit
int const bandwidth = 5;  //local seed chaining bandwidth
int const score_min = 30; //minimal score for local seed chaining
SimpleScore const scoring_scheme(3, -2, -1, -3); //scoring scheme
int const B = 7;          //width for banded alignment


//function laganChaining
template <typename TSeed, typename TSegment, typename TSize>
void laganChaining(std::list<TSeed> & chain, 
				   TSegment const & a,
				   TSegment const & b,
				   TSize q)
{
	if (((TSize)length(a) <= gaps_max) && ((TSize)length(b) <= gaps_max)) return;

	//Step 1: find seeds
	typedef typename Value<TSeed>::Type TPosition;
	typedef SeedSet<TPosition, SimpleSeed, DefaultScore> TSeedSet;
	TSeedSet seedset(limit, score_min, scoring_scheme);

	typedef Index< TSegment, IndexQGram<SimpleShape > > TQGramIndex;
	TQGramIndex index_qgram(b);

	typedef Finder<TQGramIndex> TFinder;
	TFinder finder(index_qgram);

	while (length(seedset) == 0)
	{
		if (q < q_min) return;

		resize(indexShape(index_qgram), q);
		for (unsigned int i = 0; i < length(a)-q+1; ++i)
		{
			while (find(finder, infix(a, i, i+q)))
			{
				typedef typename Position<TFinder>::Type TPosition;
				TPosition a_pos = beginPosition(a)+i;
				TPosition b_pos = beginPosition(b)+position(finder);
				if (!addSeed(seedset, a_pos, b_pos, q, 0, Merge()))
				if (!addSeed(seedset, a_pos, b_pos, q, host(a), host(b), bandwidth, Chaos()))
					 addSeed(seedset, a_pos, b_pos, q, Single());
			}
			clear(finder);
		}
		--q;
	}

	//Step 2: global chaining
	globalChaining(seedset, chain);
	clear(seedset);

	//Step 3: recursively fill gaps
	if (q > q_min)
	{
		std::list<TSeed> subchain;
		typedef typename std::list<TSeed>::iterator TIterator;

		TIterator it = chain.begin();
		TIterator it2 = it; 
		++it2;

		laganChaining(subchain, 
			infix(host(a), beginPosition(a), leftDim0(*it)), 
			infix(host(b), beginPosition(b), leftDim1(*it)), q);
		chain.splice(it, subchain);

		while(it2 != chain.end())
		{
			laganChaining(subchain, 
				infix(host(a), rightDim0(*it), leftDim0(*it2)), 
				infix(host(b), rightDim1(*it), leftDim1(*it2)), q);
			chain.splice(it2, subchain);

			it = it2;
			++it2;
		}

		laganChaining(subchain, 
			infix(host(a), rightDim0(*it), endPosition(a)), 
			infix(host(b), rightDim1(*it), endPosition(b)), q);
		chain.splice(it2, subchain);
	}
}

int main(int argc, const char *argv[])
{
	if (argc < 3) {	
		std::cerr << "Usage: ./lagan lagan1.fasta lagan2.fasta" << std::endl; 
		return -1; 
	}

	//load sequences
	typedef String<Dna> TString;
	TString a;
	TString b;
	std::fstream fstrm1;
	fstrm1.open(argv[1], std::ios_base::in | std::ios_base::binary);
	read(fstrm1, a, Fasta());
	fstrm1.close();
	std::fstream fstrm2;
	fstrm2.open(argv[2], std::ios_base::in | std::ios_base::binary);
	read(fstrm2, b, Fasta());
	fstrm2.close();

	////Doesn't work on Mac
	//TString a = String<Dna, FileReader<Fasta> >(file_a); // Crashes if file does not exist
	//TString b = String<Dna, FileReader<Fasta> >(file_b);


	//do LAGAN
	typedef Seed<int, SimpleSeed> TSeed;
	std::list<TSeed> chain;

	//Step 1 to 3
	laganChaining(chain, infix(a, 0, length(a)), infix(b, 0, length(b)), q_max);

	//Step 4: banded alignment
	Align<TString, ArrayGaps> alignment;
	resize(rows(alignment), 2);
	setSource(row(alignment, 0), a);
	setSource(row(alignment, 1), b);
	int score = bandedChainAlignment(chain, B, alignment, scoring_scheme);

	std::cout << "Score: " << score << std::endl;
	std::cout << alignment << std::endl;

	return 0;
}
