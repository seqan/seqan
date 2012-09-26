// FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
// FRAGMENT(init)
	typedef String<char>			TSequence;	// sequence type
	typedef Align<TSequence,ArrayGaps>	TAlign;		// align type

	TSequence seq1 = "CDFGHC";
	TSequence seq2 = "CDEFGAHC";

	TAlign align;
	resize(rows(align), 2); 
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

// FRAGMENT(alignment)
	int score = globalAlignment(align,Score<int>(1,-1,-1,-1), Hirschberg());
	std::cout << "Score = " << score << std::endl;
	std::cout << align;
	
//FRAGMENT(typedef2)
	StringSet<TSequence> sequences;
	appendValue(sequences,seq1);
	appendValue(sequences,seq2);
	typedef StringSet<TSequence, Dependent<> > TDepStringSet;
	typedef Graph<Alignment<TDepStringSet> > TAlignGraph;
	TAlignGraph alignG(sequences);
	
//FRAGMENT(alignment2)
	score = globalAlignment(alignG,Score<int>(1,-1,-1,-1), NeedlemanWunsch());
	std::cout << "Score = " << score << std::endl;
	std::cout << alignG;
	
	return 0;
}
