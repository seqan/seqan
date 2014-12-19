#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

struct ConvertCT :
    public std::unary_function<Dna, Dna>
{
    inline Dna operator()(Dna x) const
    {
        if (x == 'C') return 'T';

        return x;
    }

};


int main()
{
    typedef String<Dna> TSequence;

    TSequence seq1 = "CCCGGCATCATCC";
    TSequence seq2 = "CTTGGCATTATTC";

    std::cout << seq1 << std::endl;
    std::cout << seq2 << std::endl;
    std::cout << std::endl;

    typedef ModifiedString<TSequence, ModView<ConvertCT> > TModCT;
    TModCT modCT1(seq1);
    TModCT modCT2(seq2);

    std::cout << modCT1 << std::endl;
    std::cout << modCT2 << std::endl;

    return 0;
}
