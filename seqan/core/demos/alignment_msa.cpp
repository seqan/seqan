///A tutorial about heuristic multiple sequence alignment.
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
/// A set of sequences to align
	typedef String<AminoAcid> TSequence;
	TSequence seq1 = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE";
	TSequence seq2 = "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK";
	TSequence seq3 = "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK";
	TSequence seq4 = "MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK";

/// Scoring: Blosum62 where gex = -1, gop = -11
	Blosum62 sc(-1, -11);

/// Create an Align object or alternatively an Alignment graph object to store the MSA
	Align<TSequence, ArrayGaps> align;
	resize(rows(align), 4);
	assignSource(row(align, 0), seq1);
	assignSource(row(align, 1), seq2);
	assignSource(row(align, 2), seq3);
	assignSource(row(align, 3), seq4);

/// Heuristic MSA
	globalMsaAlignment(align, sc);

/// Output of the MSA
	std::cout << align << std::endl;
	return 0;
}
