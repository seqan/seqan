// FRAGMENT(includes)
#include <iostream>
#include <seqan/graph_algorithms.h>
#include <seqan/basic/basic_math.h>

using namespace seqan;

// FRAGMENT(typedefs)
int main() {
	typedef LogProb<> TProbability;
	typedef Dna TAlphabet;
	typedef Size<TAlphabet>::Type TSize;
	typedef Graph<Hmm<TAlphabet, TProbability, Default> > THmm;
	typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
	typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;
	
// FRAGMENT(variables)
	Dna dnaA = Dna('A');
	Dna dnaC = Dna('C');
	Dna dnaG = Dna('G');
	Dna dnaT = Dna('T');

	THmm hmm;

// FRAGMENT(begin-state)
	TVertexDescriptor begState = addVertex(hmm);
	assignBeginState(hmm, begState);

// FRAGMENT(main-states-emissions)
	TVertexDescriptor exonState = addVertex(hmm);
	emissionProbability(hmm, exonState, dnaA) = 0.25;
	emissionProbability(hmm, exonState, dnaC) = 0.25;
	emissionProbability(hmm, exonState, dnaG) = 0.25;
	emissionProbability(hmm, exonState, dnaT) = 0.25;

	TVertexDescriptor spliceState = addVertex(hmm);
	emissionProbability(hmm, spliceState, dnaA) = 0.05;
	emissionProbability(hmm, spliceState, dnaC) = 0.0;
	emissionProbability(hmm, spliceState, dnaG) = 0.95;
	emissionProbability(hmm, spliceState, dnaT) = 0.0;

	TVertexDescriptor intronState = addVertex(hmm);
	emissionProbability(hmm, intronState, dnaA) = 0.4;
	emissionProbability(hmm, intronState, dnaC) = 0.1;
	emissionProbability(hmm, intronState, dnaG) = 0.1;
	emissionProbability(hmm, intronState, dnaT) = 0.4;

// FRAGMENT(end-state)
	TVertexDescriptor endState = addVertex(hmm);
	assignEndState(hmm, endState);

// FRAGMENT(transitions)
	addEdge(hmm, begState, exonState, 1.0);
	addEdge(hmm, exonState, exonState, 0.9);
	addEdge(hmm, exonState, spliceState, 0.1);
	addEdge(hmm, spliceState, intronState, 1.0);
	addEdge(hmm, intronState, intronState, 0.9);
	addEdge(hmm, intronState, endState, 0.1);

// FRAGMENT(print-model)
	std::cout << hmm << std::endl;

// FRAGMENT(viterbi)
	String<Dna> sequence = "CTTCATGTGAAAGCAGACGTAAGTCA";
	String<TVertexDescriptor> path;
	TProbability p = viterbiAlgorithm(hmm, sequence, path);
	std::cout << "Viterbi algorithm" << std::endl;
	std::cout << "Probability of best path: " << p << std::endl;
	std::cout << "Sequence: " << std::endl;
	for(TSize i = 0; i < length(sequence); ++i) std::cout << sequence[i] << ',';
	std::cout << std::endl;
	std::cout << "State path: " << std::endl;
	for(TSize i = 0; i < length(path); ++i) {
		std::cout << path[i];
		if (isSilent(hmm, path[i])) std::cout << " (Silent)";
		if (i < length(path) - 1) std::cout << ',';
	}
	std::cout << std::endl;

// FRAGMENT(forward-algorithm)
	std::cout << "Forward algorithm" << std::endl;
	p = forwardAlgorithm(hmm, sequence);
	std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;

// FRAGMENT(backward-algorithm)
	std::cout << "Backward algorithm" << std::endl;
	p = backwardAlgorithm(hmm, sequence);
	std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;

	return 0;
}
