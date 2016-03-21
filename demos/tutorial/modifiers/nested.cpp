///A tutorial showing how to nest modifiers.
#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

int main()
{
//![string]
    String<Dna> myString = "attacgg";
//![string]
///A nested modifier.
//![complement]
    typedef ModifiedString<String<Dna>, ModComplementDna>   TMyComplement;
//![complement]
//![reverse]
    typedef ModifiedString<TMyComplement, ModReverse> TMyReverseComplement;
//![reverse]

///A reverse complemented string.
//![constructor]
    TMyReverseComplement myReverseComplement(myString);
//![constructor]
    std::cout << "//![output]" << '\n';
//![output]
    std::cout << myString << '\n';
    std::cout << myReverseComplement << '\n';

    replace(myString, 1, 1, "cgt");

    std::cout << myString << '\n';
    std::cout << myReverseComplement << '\n';
//![output]
    std::cout << "//![output]" << '\n';
    std::cout << myString << '\n';
//![alternative]
    std::cout << DnaStringReverseComplement(myString) << std::endl;
//![alternative]
    return 0;
}
