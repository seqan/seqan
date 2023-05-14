#include <iostream>
#include <seqan/index.h>

using namespace seqan2;

int main()
{
	typedef Index<CharString> TIndex;
	TIndex index("mississippi");
	Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
/*
	do {
		//...
	} while (isRoot(it));
*/
	return 0;
}
