// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	typedef Index<DnaString, IndexQGram< UngappedShape<3> > > TIndex;
	TIndex index("CATGATTACATA");

// FRAGMENT(output)
	hash(indexShape(index), "CAT");
	for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
		std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

	return 0;
}
