//![main]
#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

int main()
{
//![main]
//![mod]
    String<Dna> mySeq = "ACCGTT";

    ModifiedString<String<Dna>, ModComplementDna> myCompl(mySeq);

    std::cout << mySeq << std::endl;
    std::cout << myCompl << std::endl;
//![mod]
//![end]
    return 0;
}
//![end]
