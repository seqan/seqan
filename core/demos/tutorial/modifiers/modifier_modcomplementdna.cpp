// FRAGMENT(main)
#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;


int main ()
{
// FRAGMENT(mod)
	String<Dna> mySeq = "ACCGTT";

	ModifiedString< String<Dna>, ModComplementDna > myCompl(mySeq);

	std::cout << mySeq << std::endl;
	std::cout << myCompl << std::endl;
// FRAGMENT(end)
	return 0;
}
