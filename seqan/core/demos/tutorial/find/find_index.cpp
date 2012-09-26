// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main()
{
	Index<CharString> index("tobeornottobe");
	CharString needle = "be";
	Finder<Index<CharString> > finder(index);
// FRAGMENT(output)
	Pattern<CharString> pattern(needle);
	while (find(finder, pattern))
		std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

// FRAGMENT(output_short)
	clear(finder);
	while (find(finder, "be"))
		std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

	return 0;
}
