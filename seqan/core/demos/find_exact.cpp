///A tutorial about the use of exact find algorithms.
#include <iostream>
#include <seqan/find.h>

using namespace seqan;

///This function prints the positions of all occurrences of $needle$ within $haystack$.
///It uses the algorithm specified by $TAlgorithm$ that is passed as a @glos:Specialization@ to @Class.Pattern@.
template <typename TAlgorithm>
void printAllOccs(String<char>& haystack, 
				  String<char>& needle)
{
	Finder<String<char> > finder(haystack);
	Pattern<String<char>, TAlgorithm> pattern(needle);
	while (find(finder, pattern)) 
	{
		std::cout << position(finder) << ", ";
	}
	std::cout << std::endl;
}

///The main function calls $printAllOccs$ for different exact string matching algorithms.
int main() 
{
	String<char> haystack = "send more money!";
	String<char> needle = "mo";

	printAllOccs<Horspool>(haystack, needle);
	printAllOccs<BomAlgo> (haystack, needle);
	printAllOccs<BndmAlgo>(haystack, needle);
	printAllOccs<ShiftAnd>(haystack, needle);
	printAllOccs<ShiftOr> (haystack, needle);

	return 0;
}

