///A tutorial about the heaviest increasing subsequence algorithm.
#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;


int main() {
///Creation of a simple sequence and corresponding weights
	String<char> seq("zeitgeist");
	String<unsigned int> weights;
	resize(weights, length(seq), 1);
	assignProperty(weights, 2, 10);
///Out-parameter: A string of positions belonging to the heaviest increasing subsequence
	typedef Position<String<unsigned int> >::Type TPosition;
	String<TPosition> pos;
///Heaviest increasing subsequence and the corresponding weight
	unsigned int w = heaviestIncreasingSubsequence(seq, weights, pos);
///Console Output
	for(int i = 0; i< (int) length(seq); ++i) {
		std::cout << seq[i] << "(Weight=" << getProperty(weights, i) << "),";
	}
	std::cout << std::endl;
	std::cout << "His: " << std::endl;
	for(int i = length(pos)-1; i>=0; --i) {
		std::cout << seq[pos[i]] <<  ',';
	}
	std::cout << "(Weight=" << w << ')' << std::endl;
	return 0;
}
