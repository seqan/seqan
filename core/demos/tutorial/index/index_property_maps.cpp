// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	String<char> myString = "abracadabra";

	typedef Index< String<char>, IndexWotd<> > TMyIndex;
	TMyIndex myIndex(myString);
	String<int> propMap;

// FRAGMENT(iteration)
	Iterator< TMyIndex, TopDown< ParentLinks<Preorder> > >::Type myIterator(myIndex);

	int depth;
	while (!atEnd(myIterator))
	{
		if (isRoot(myIterator))
			depth = 0;
		else
			depth = getProperty(propMap, nodeUp(myIterator)) + 1;

		resizeVertexMap(myIndex, propMap);
		assignProperty(propMap, value(myIterator), depth);

		++myIterator;
	}

// FRAGMENT(output)
	goBegin(myIterator);
	while (!atEnd(myIterator))
	{
		std::cout << getProperty(propMap, value(myIterator)) << '\t' << representative(myIterator) << std::endl;
		++myIterator;
	}
	return 0;
}
