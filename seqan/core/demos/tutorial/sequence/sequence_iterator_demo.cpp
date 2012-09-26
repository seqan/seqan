// The comment lines containing FRAGMENT(fragment-line) are there for the
// documentation system.  You can ignore them when reading this file.
// FRAGMENT(includes)
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main() {
// FRAGMENT(metafunctions)
	String<char> str = "admn";
	Iterator<String<char> >::Type it = begin(str);
	Iterator<String<char> >::Type itEnd = end(str);
// FRAGMENT(iterators)
	while (it != itEnd) {
		std::cout << *it;
		++it;
	}
	std::cout << std::endl;
// FRAGMENT(rooted-iterators)
	Iterator<String<char>, Rooted >::Type it2 = begin(str);
	for (goBegin(it2); !atEnd(it2); goNext(it2)) {
		++value(it2);
	}
// FRAGMENT(iterator-reverse)
	goEnd(it2);
	while (!atBegin(it2)) {
		goPrevious(it2);
		std::cout << getValue(it2);
	}
	std::cout << std::endl;
// FRAGMENT(assign-value)
	assignValue(begin(str), 'X');
	std::cout << str << std::endl;
	
	return 0;
}
