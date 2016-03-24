#include <seqan/sequence.h>
#include <iostream>

using namespace seqan;

int main()
{
//![construction]
    Iterator<DnaString>::Type           it1;  // A standard iterator
    Iterator<DnaString, Standard>::Type it2;  // Same as above
    Iterator<DnaString, Rooted>::Type   it3;  // A rooted iterator
//![construction]
    ignoreUnusedVariableWarning(it1);
    ignoreUnusedVariableWarning(it2);
    ignoreUnusedVariableWarning(it3);

    std::cout << "\n//![use-case]\n";
//![use-case]
    DnaString genome = "ACGTACGTACGT";
    typedef Iterator<DnaString>::Type TIterator;
    for (TIterator it = begin(genome); it != end(genome); ++it)
    {
        std::cout << *it;
    }
//![use-case]
    std::cout << "\n//![use-case]\n";

    return 0;
}
