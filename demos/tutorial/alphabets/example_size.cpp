#include <iostream>
#include <seqan/basic.h>

using namespace seqan;

int main()
{
//![main]
    typedef Dna TAlphabet;

    unsigned alphSize = ValueSize<TAlphabet>::VALUE;
    std::cout << "Alphabet size of Dna: " << alphSize << '\n';
//![main]

    return 0;
}
