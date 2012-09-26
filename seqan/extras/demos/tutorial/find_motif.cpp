// FRAGMENT(includes)
#include <iostream>
#include <seqan/find_motif.h>

using namespace seqan;

// FRAGMENT(typedefs)
int main() 
{
    typedef MotifFinder<char, Projection> TMotifFinder;
    typedef String<CharString > TString;
    typedef Size<TString>::Type TSize;

// FRAGMENT(sequences)
	TString dataset;
	appendValue(dataset, CharString("hatpins"));
	appendValue(dataset, CharString("low-fat"));
	appendValue(dataset, CharString("habitat"));
    TSize seqCount = length(dataset);

// FRAGMENT(initialization)
	std::srand((unsigned) time(NULL));

    TSize seqLength = length(dataset[0]); // length of sequences
	TSize motifLength = 3;		          // length of motif
	TSize mm = 1;	                	  // number of mismatches
	bool is_exact = true;	              // occurences of motif need to have exactly mm mismatches
    TSize numPos = seqCount * (seqLength - motifLength + 1);

	TMotifFinder finder_proj(seqCount, motifLength, numPos, mm, is_exact);

// FRAGMENT(search)
	findMotif(finder_proj, dataset, Oops());

	for (int i = 0; i < (int) motifCount(finder_proj); ++i)
		std::cout << i << ": " << getMotif(finder_proj, i) << std::endl;

	return 0;
}

