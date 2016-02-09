#include <iostream>
#include <seqan/basic.h>

using namespace seqan;

int main()
{
    typedef Dna TAlphabet;

//![main]
    unsigned bits = BitsPerValue<TAlphabet>::VALUE;
    std::cout << "Number of bits needed to store a value of type Dna: " << bits << '\n';
//![main]

    return 0;
}
