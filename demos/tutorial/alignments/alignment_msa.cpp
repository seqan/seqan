//![main]
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
    char const * strings[4] =
    {
        "DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMA"
        "KADKARYEREMKTYIPPKGE",
        "RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQA"
        "MHREKYPNYKYRPRRKAKMLPK",
        "FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQ"
        "EFERNLARFREDHPDLIQNAKK",
        "HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQ"
        "LHMQLYPGWSARDNYGKKKKRKREK"
    };
//![main]

//![init]
    Align<String<AminoAcid> > align;
    resize(rows(align), 4);
    for (int i = 0; i < 4; ++i)
        assignSource(row(align, i), strings[i]);
//![init]

//![alignment]
    globalMsaAlignment(align, Blosum62(-1, -11));
    std::cout << align << "\n";

    return 0;
}
//![alignment]
