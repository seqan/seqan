//FRAGMENT(includes)
#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

//FRAGMENT(count-one-mers-define-function)
template <typename TString>
void countOneMers(TString const sequence) {
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TAlphabet;
	
	//FRAGMENT(count-one-mers-initialize-table)
	typedef String<TSize> TCounterString;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TCounterString counter;
	resize(counter, alphSize, 0);
	
	//FRAGMENT(count-one-mers-count-chars)	
	typedef typename Iterator<TString const>::Type TIter;
	TIter itEnd = end(sequence);
	for (TIter it = begin(sequence); it != itEnd; goNext(it))
		value(counter, ordValue(value(it))) += 1;
		
	//FRAGMENT(count-one-mers-print-chars)
	typedef typename Iterator<TCounterString>::Type TCounterIter;
	TCounterIter countIt = begin(counter);
	TCounterIter countItEnd = end(counter);
	for (TSize pos = 0; countIt != countItEnd; ++countIt, ++pos)
		if (value(countIt) > 0)
			std::cout << TAlphabet(pos) << ':' << value(countIt) << std::endl;
}

//FRAGMENT(main)
int main() {
	std::cout << "ASCII String: Hello world!" << std::endl;
	countOneMers<CharString>("Hello world!");
	
	std::cout << "DNA String: TATACGCTA" << std::endl;
	countOneMers<DnaString>("TATACGCTA");
	
	std::cout << "Peptide String: MQDRVKRPMNAFIVWSRDQRRKMALEN" << std::endl;
	countOneMers<Peptide>("MQDRVKRPMNAFIVWSRDQRRKMALEN");
	
	return 0;
}
