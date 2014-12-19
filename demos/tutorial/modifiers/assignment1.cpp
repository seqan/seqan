#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

int main()
{
    typedef String<Dna> TSequence;

    TSequence seq1 = "CCCGGCATCATCC";
    TSequence seq2 = "CTTGGCATTATTC";

    std::cout << seq1 << std::endl;
    std::cout << seq2 << std::endl;
    std::cout << std::endl;

    return 0;
}
