// FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
	typedef String<Dna>			TSequence;	// sequence type
	typedef Align<TSequence,ArrayGaps>	TAlign;		// align type
	typedef Row<TAlign>::Type		TRow;
	typedef Iterator<TRow>::Type		TIterator;

// FRAGMENT(init)
	TSequence seq1 = "acgtacgtact";
	TSequence seq2 = "actactacgt";

	TAlign align;
	resize(rows(align), 2); 
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

// FRAGMENT(alignment)
	int score = globalAlignment(align,Score<int>(1,-1,-1,-1), Hirschberg());
	std::cout << "Score = " << score << std::endl;
	std::cout << align;
	
//FRAGMENT(iterate)
	unsigned aliLength = _max(length(row(align, 0)), length(row(align, 1)));
	for(unsigned i = 0; i < length(rows(align)); ++i)
	{
		TIterator it = iter(row(align,i), 0);
		TIterator itEnd = iter(row(align,i), aliLength);
		unsigned pos = 0;
		std::cout << "Row " << i << " contains gaps at positions: ";
		std::cout << std::endl;
		while (it != itEnd)
		{
			if(isGap(it))
 				std::cout << pos << std::endl;
			++it;
			++pos;
		}
	}

	return 0;
}
