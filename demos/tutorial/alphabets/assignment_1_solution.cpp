#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;

// We define a function which takes
// the alphabet type as an argument
template <typename TAlphabet>
void showAllLettersOfMyAlphabet(TAlphabet const &)
{
    typedef typename ValueSize<TAlphabet>::Type TSize;
    // We need to determine the alphabet size
    // using the metafunction ValueSize
    TSize alphSize = ValueSize<TAlphabet>::VALUE;
    // We iterate over all characters of the alphabet
    // and output them
    for (TSize i = 0; i < alphSize; ++i)
        std::cout << static_cast<unsigned>(i) << ',' << TAlphabet(i) << "  ";
    std::cout << std::endl;

}

int main()
{
    showAllLettersOfMyAlphabet(AminoAcid());
    showAllLettersOfMyAlphabet(Dna());
    showAllLettersOfMyAlphabet(Dna5());
    return 0;
}
