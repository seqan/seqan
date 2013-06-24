///A tutorial about the use of the reverse modifier.
#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;


int main ()
{
// FRAGMENT(string)
	String<char> myString = "A man, a plan, a canal-Panama";
// FRAGMENT(modifier)
	ModifiedString< String<char>, ModReverse > myModifier(myString);

// FRAGMENT(output1)
	std::cout << myString << std::endl;
	std::cout << myModifier << std::endl;
// FRAGMENT(output2)
	replace(myString, 9, 9, "master ");
	std::cout << myString << std::endl;
	std::cout << myModifier << std::endl;
	return 0;
}
