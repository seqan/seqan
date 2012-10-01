///A tutorial about motif search.
#include <iostream>
#include "seqan/find_motif.h"

using namespace seqan;

///Function to output found motifs.
template <typename TMotifFinder>
void printMotifs(TMotifFinder & finder)
{
	for (int i = 0; i < (int) motifCount(finder); ++i)
	{
		std::cout << i << ": " << getMotif(finder, i) << std::endl;
	}
}

int main() 
{
	std::srand((unsigned) time(NULL));

///Motif search on a small set of nucleotide sequences.
	unsigned int t = 3;		//number of input sequences
	unsigned int n = 6;		//length of sequence
	unsigned int l = 4;		//length of motif
	unsigned int d = 1;		//number of substitutions
	bool is_exact = true;	//size of Hamming distance
	unsigned int h = 0;		//size of the neighborhood considering at first

	String<DnaString> dataset;
	appendValue(dataset,DnaString("ACAGCA"));
	appendValue(dataset,DnaString("AGGCAG"));
	appendValue(dataset,DnaString("TCAGTC"));
	
///Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb1(t,l,d,is_exact,h);
	findMotif(finder_epb1,dataset,Omops());
	std::cout << getMotif(finder_epb1) << std::endl;

///Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb2(t,l,d,is_exact,h);
	findMotif(finder_epb2,dataset,Oops());
	std::cout << getMotif(finder_epb2) << std::endl;

///Application of Pms1-Zoops 
	MotifFinder<Dna, Pms1> finder_pms1(l,d,is_exact);
	findMotif(finder_pms1,dataset,Zoops());
	printMotifs(finder_pms1); 

///Application of Pmsp-Tcm
	MotifFinder<Dna, Pmsp> finder_pmsp(l,d,is_exact);
	findMotif(finder_pmsp,dataset,Tcm());
	printMotifs(finder_pmsp); 
	
///Application of PROJECTION-Oops
	unsigned int m = t*(n-l+1);
    MotifFinder<Dna, Projection> finder_proj(t,l,m,d,is_exact);
	findMotif(finder_proj, dataset, Oops());
	printMotifs(finder_proj);

///Application of PROJECTION-Omops
    MotifFinder<Dna, Projection> finder_proj_omops(t,l,m,d,is_exact);
	findMotif(finder_proj_omops, dataset, Omops());
	printMotifs(finder_proj_omops);

///Application of PROJECTION-Zoops
	MotifFinder<Dna, Projection> finder_proj_zoops(t,l,m,d,is_exact);
	findMotif(finder_proj_zoops, dataset, Zoops());
	printMotifs(finder_proj_zoops);
	
///Application of PROJECTION-Tcm
    MotifFinder<Dna, Projection> finder_proj_tcm(t,l,m,d,is_exact);
	findMotif(finder_proj_tcm, dataset, Tcm());
	printMotifs(finder_proj_tcm);

	return 0;
}

