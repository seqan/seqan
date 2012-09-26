//FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
//FRAGMENT(init)
	typedef String<AminoAcid> TSequence;
	StringSet<TSequence> seq;
	appendValue(seq,"DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE");
	appendValue(seq,"RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK");
	appendValue(seq,"FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK");
	appendValue(seq,"HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK");

//FRAGMENT(alignment)
	Graph<Alignment<StringSet<TSequence, Dependent<> > > > aliG(seq);
	globalMsaAlignment(aliG, Blosum62(-1, -11));
	std::cout << aliG << std::endl;
	
	return 0;
}
