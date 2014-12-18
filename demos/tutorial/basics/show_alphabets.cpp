//![includes]
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;
//![includes]

//![showAllLettersOfMyAlphabet]
template <typename TAlphabet>
void showAllLettersOfMyAlphabet(TAlphabet const &)
{
    typedef typename Size<TAlphabet>::Type TSize;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;
    for (TSize i = 0; i < alphSize; ++i)
        std::cout << i << ',' << TAlphabet(i) << "  ";
    std::cout << std::endl;
}
//![showAllLettersOfMyAlphabet]

//![main]
int main()
{
    showAllLettersOfMyAlphabet(AminoAcid());
    showAllLettersOfMyAlphabet(Dna());
    showAllLettersOfMyAlphabet(Dna5());
    return 0;
}
//![main]
