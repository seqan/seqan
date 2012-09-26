// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	String<char> myString = "abracadabra";

	typedef Index< String<char> > TMyIndex;
    
// FRAGMENT(iterator)
	TMyIndex myIndex(myString);

// FRAGMENT(iteration)
	Iterator< TMyIndex, TopDown< ParentLinks<Preorder> > >::Type myIterator(myIndex);

	while (!atEnd(myIterator))
	{
		std::cout << representative(myIterator) << std::endl;
		++myIterator;
	}

	return 0;
}
