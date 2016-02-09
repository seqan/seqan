#include <iostream>
#include <seqan/basic.h>

using namespace seqan;

int main()
{
//![main]
    Dna a = 'A';
    Dna c = 'C';
    Dna g = 'G';
    Dna t = 'T';

    std::cout <<"A: " << (unsigned)ordValue(a) << '\n';
    std::cout <<"C: " << (unsigned)ordValue(c) << '\n';
    std::cout <<"G: " << (unsigned)ordValue(g) << '\n';
    std::cout <<"T: " << (unsigned)ordValue(t) << '\n';
//![main]

    return 0;
}
