// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	Index<DnaString, IndexQGram<OneGappedShape> > index("CATGATTACATA");
	stringToShape(indexShape(index), "1101");

// FRAGMENT(output)
	hash(indexShape(index), "ATCA");
	for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
		std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

	return 0;
}
