// FRAGMENT(includes)
#include <iostream>
#include <seqan/find_motif.h>

using namespace seqan;

// FRAGMENT(typedefs)
int main() 
{
    typedef MotifFinder<Dna, Pms1> TMotifFinder;
    typedef String<DnaString> TString;
    typedef Size<TString>::Type TSize;

// FRAGMENT(sequences)
	TString dataset;
	appendValue(dataset, DnaString("ACAGCA"));
	appendValue(dataset, DnaString("AGGCAG"));
	appendValue(dataset, DnaString("TCAGTC"));

// FRAGMENT(initialization)
	TSize motifLength = 4;		//length of motif
	TSize mm = 1;	        	//number of mismatches
	bool is_exact = false;	    //occurences of motif need to have exactly mm mismatches

	TMotifFinder finder_pms1(motifLength, mm, is_exact);

// FRAGMENT(search)
	findMotif(finder_pms1, dataset, Omops());

	for (int i = 0; i < (int) motifCount(finder_pms1); ++i)
		std::cout << i << ": " << getMotif(finder_pms1, i) << std::endl;

	return 0;
}

