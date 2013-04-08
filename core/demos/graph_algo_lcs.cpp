///A tutorial about the longest common subsequence algorithm.
#include <iostream>
#include <seqan/graph_algorithms.h>
#include <seqan/graph_align.h>

using namespace seqan;

int main() {
///Create two sequences
	String<char> seq1("abacx");
	String<char> seq2("baabca");
	typedef StringSet<String<char>, Dependent<> > TStringSet;
	TStringSet string_set;
	appendValue(string_set, seq1);
	appendValue(string_set, seq2);
	Graph<Alignment<TStringSet> > alignment_graph(string_set);
///Compute the longest common subsequence
	std::cout << "Score = " << globalAlignment(alignment_graph, stringSet(alignment_graph), Lcs()) << std::endl;
	std::cout << alignment_graph << std::endl;
	return 0;
}
