///A tutorial about HMMs with silent states
#include <iostream>
#include <fstream>
#include <seqan/graph_algorithms.h>
#include <seqan/basic/basic_math.h>

using namespace seqan;

int main()
{
///HMM creation
    typedef LogProb<> TProbability;
    typedef Dna TAlphabet;
    typedef Size<TAlphabet>::Type TSize;
    typedef Graph<Hmm<TAlphabet, TProbability> > THmm;
    typedef VertexDescriptor<THmm>::Type TVertexDescriptor;

    Dna dnaA = Dna('A');
    Dna dnaC = Dna('C');
    Dna dnaG = Dna('G');
    Dna dnaT = Dna('T');

    THmm hmm;

///Begin state
    TVertexDescriptor begState = addVertex(hmm);
    assignBeginState(hmm, begState);

///Emission states
    TVertexDescriptor emitState1 = addVertex(hmm);
    emissionProbability(hmm, emitState1, dnaA) = 0.0;
    emissionProbability(hmm, emitState1, dnaC) = 0.8;
    emissionProbability(hmm, emitState1, dnaG) = 0.1;
    emissionProbability(hmm, emitState1, dnaT) = 0.1;

    TVertexDescriptor emitState2 = addVertex(hmm);
    emissionProbability(hmm, emitState2, dnaA) = 0.0;
    emissionProbability(hmm, emitState2, dnaC) = 0.2;
    emissionProbability(hmm, emitState2, dnaG) = 0.2;
    emissionProbability(hmm, emitState2, dnaT) = 0.6;

    TVertexDescriptor emitState3 = addVertex(hmm);
    emissionProbability(hmm, emitState3, dnaA) = 0.7;
    emissionProbability(hmm, emitState3, dnaC) = 0.1;
    emissionProbability(hmm, emitState3, dnaG) = 0.1;
    emissionProbability(hmm, emitState3, dnaT) = 0.1;

    TVertexDescriptor emitState4 = addVertex(hmm);
    emissionProbability(hmm, emitState4, dnaA) = 0.25;
    emissionProbability(hmm, emitState4, dnaC) = 0.25;
    emissionProbability(hmm, emitState4, dnaG) = 0.25;
    emissionProbability(hmm, emitState4, dnaT) = 0.25;

///Silent states (deletion states)
    TVertexDescriptor delState1 = addVertex(hmm, true);
    TVertexDescriptor delState2 = addVertex(hmm, true);
    TVertexDescriptor delState3 = addVertex(hmm, true);
    TVertexDescriptor delState4 = addVertex(hmm, true);

///End state
    TVertexDescriptor eState = addVertex(hmm);
    assignEndState(hmm, eState);

///Transitions
    addEdge(hmm, begState, emitState1, 0.5);
    addEdge(hmm, begState, delState1, 0.5);
    addEdge(hmm, emitState1, emitState2, 0.5);
    addEdge(hmm, emitState1, delState2, 0.5);
    addEdge(hmm, delState1, emitState2, 0.5);
    addEdge(hmm, delState1, delState2, 0.5);
    addEdge(hmm, emitState2, emitState3, 0.5);
    addEdge(hmm, emitState2, delState3, 0.5);
    addEdge(hmm, delState2, emitState3, 0.5);
    addEdge(hmm, delState2, delState3, 0.5);
    addEdge(hmm, emitState3, emitState4, 0.5);
    addEdge(hmm, emitState3, delState4, 0.5);
    addEdge(hmm, delState3, emitState4, 0.5);
    addEdge(hmm, delState3, delState4, 0.5);
    addEdge(hmm, emitState4, eState, 1.0);
    addEdge(hmm, delState4, eState, 1.0);

///Print the whole model
    std::cout << hmm << std::endl;

///Viterbi algorithm
    String<Dna> sequence = "CA";
    String<TVertexDescriptor> path;
    TProbability p = viterbiAlgorithm(path, hmm, sequence);
    std::cout << "Viterbi algorithm" << std::endl;
    std::cout << "Probability of the best path: " << p << std::endl;
    std::cout << "Sequence: " << std::endl;
    for (TSize i = 0; i < length(sequence); ++i)
        std::cout << sequence[i] << ',';
    std::cout << std::endl;
    std::cout << "State path: " << std::endl;
    for (TSize i = 0; i < length(path); ++i)
    {
        std::cout << path[i];
        if (isSilent(hmm, path[i]))
            std::cout << " (Silent)";
        if (i < length(path) - 1)
            std::cout << ',';
    }
    std::cout << std::endl;

///Forward algorithm
    std::cout << "Forward algorithm" << std::endl;
    p = forwardAlgorithm(hmm, sequence);
    std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;

///Backward algorithm
    std::cout << "Backward algorithm" << std::endl;
    p = backwardAlgorithm(hmm, sequence);
    std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;

    return 0;
}
