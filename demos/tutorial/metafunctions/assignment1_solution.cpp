#include <iostream>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace seqan;

template <typename T>
void checkContainerForDna(T & container)
{
    // Type Metafunction Value<>
    typedef typename Value<T>::Type TAlphType;

    // Value Metafunction IsSameType<> evaluated at compile time
    if (IsSameType<TAlphType, Dna>::VALUE)
        std::cout << "I have found a container with Dna!" <<  std::endl;
    else
        std::cout << "No Dna anywhere." <<  std::endl;
}

int main()
{
    typedef String<Dna> TDnaString;
    TDnaString dna = "AAAATTTT";

    typedef String<int> TIntString;

    TIntString numbers;
    appendValue(numbers, 1);
    appendValue(numbers, 3);

    checkContainerForDna(dna);
    checkContainerForDna(numbers);

    return 0;
}
