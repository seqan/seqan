///A tutorial about the use of alphabets.
#include <iostream>
#include <seqan/basic.h>
using namespace seqan;

int main()
{
///The typical alphabet is convertible to $char$.
///A conversion of a $char$ back and forth into another alphabet can, however, change the value of the $char$.
    Dna a = 'a';
    std::cout << a << std::endl;

///'f' is not in Dna5 and hence $b$ is set to 'N'.
    Dna5 b = 'f';
    std::cout << b << std::endl;

///Many SeqAn alphabet classes can be converted into each other.
    b = a;
    std::cout << b << std::endl;

    Iupac c = b;
    std::cout << c << std::endl;

    return 0;
}
