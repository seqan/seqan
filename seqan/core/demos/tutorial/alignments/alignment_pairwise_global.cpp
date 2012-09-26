// FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
	typedef String<char>				TSequence;	// sequence type
	typedef Align<TSequence,ArrayGaps>	TAlign;		// align type
	typedef Row<TAlign>::Type			TRow;
	typedef Iterator<TRow>::Type		TIterator;

// FRAGMENT(init)
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
	
//FRAGMENT(iterate)
	TIterator it1 = begin(row(align,0));
	TIterator it1End = end(row(align,0));
	TIterator it2 = begin(row(align,1));
	TIterator it2End = end(row(align,1));

	while (it1 != it1End && it2 != it2End)
	{
		if(isGap(it1)) std::cout << "#";
		else std::cout << *it1;
		std::cout << ' ';
		if(isGap(it2)) std::cout << "#";
		else std::cout << *it2;
		std::cout << std::endl;
		++it1;
		++it2;
	}


	return 0;
}
