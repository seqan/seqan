//FRAGMENT(includes)
#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

void countOneMers(CharString const & str) {
	//FRAGMENT(count-one-mers-initialize-table)
	String<int> table;
	resize(table, 'z' - 'a' + 1, 0);
	
	//FRAGMENT(count-one-mers-count-chars)
	for (unsigned int i = 0; i < length(str); ++i) {
		table[str[i] - 'a'] += 1;
	}
	
	//FRAGMENT(count-one-mers-print-chars)
	for (unsigned int i = 0; i < 'z' - 'a' + 1; ++i) {
		if (table[i] == 0)
			continue;
		std::cout << static_cast<char>('a' + i) << " " << table[i] << std::endl;
	}
}

//FRAGMENT(main)
int main() {
	std::cout << "String: helloworld" << std::endl;
	countOneMers("helloworld");
	
	std::cout << "String: mississippi" << std::endl;
	countOneMers("mississippi");
	
	std::cout << "String: banana" << std::endl;
	countOneMers("banana");
	
	return 0;
}
