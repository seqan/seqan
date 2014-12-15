///A tutorial showing how to nest modifiers.
#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

int main()
{
    String<Dna> myString = "attacgg";
///A nested modifier.
    typedef ModifiedString<String<Dna>, ModComplementDna>   TMyComplement;
    typedef ModifiedString<TMyComplement, ModReverse>       TMyReverseComplement;

///A reverse complemented string.
    TMyReverseComplement myReverseComplement(myString);
    std::cout << myString << std::endl;
    std::cout << myReverseComplement << std::endl;
    replace(myString, 1, 1, "cgt");
    std::cout << myString << std::endl;
    std::cout << myReverseComplement << std::endl;
    std::cout << DnaStringReverseComplement(myString) << std::endl;
    return 0;
}
