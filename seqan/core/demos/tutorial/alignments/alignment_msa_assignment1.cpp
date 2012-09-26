//FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
	typedef String<AminoAcid> TSequence;
	Align<TSequence> align;
	resize(rows(align), 4);
	assignSource(row(align, 0),"DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE");
	assignSource(row(align, 1),"RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK");
	assignSource(row(align, 2),"FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK");
	assignSource(row(align, 3),"HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK");

//FRAGMENT(alignment)
	globalMsaAlignment(align, Blosum80(-1, -11));
	std::cout << align << std::endl;
	
	return 0;
}
