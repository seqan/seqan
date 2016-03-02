#include <seqan/sequence.h>

using namespace seqan;

//![func_exchange1]
void amino_exchangeFirstValues(String<AminoAcid> & str)
{
    if (length(str) < 2)
        return;
    AminoAcid temp = str[0];
    str[0] = str[1];
    str[1] = temp;
}
//![func_exchange1]

//![func_exchange2]
template <typename T>
void general_exchangeFirstValues(T & str)
{
    if (length(str) < 2)
        return;
    AminoAcid temp = str[0];
    str[0] = str[1];
    str[1] = temp;
}
//![func_exchange2]

//![func_exchange3]
template <typename T>
void exchangeFirstValues(T & str)
{
    if (length(str) < 2)
        return;
    typename Value<T>::Type temp = str[0];
    str[0] = str[1];
    str[1] = temp;
}
//![func_exchange3]

//![length]
template <typename T>
void printLenOfFixedSizeString(T const &)
{
    std::cout << LENGTH<T>::VALUE << std::endl;
}

//![length]
//![length]
int main()
{
//![length]
//![amino]
    String<AminoAcid> amino_str = "ARN";
//![amino]
    std::cout << "//![iterator]" << std::endl;
//![iterator]
    String<char> str = "I am a string";
    Iterator<String<char> >::Type it = begin(str);
    while (! atEnd(it, str))
    {
        std::cout << *it;
        ++it;
    }
    std::cout << std::endl;
//![iterator]
    std::cout << "//![iterator]" << std::endl;
    std::cout << "//![length]" << std::endl;
//![length]
    String<char, Array<100> > my_str;
    printLenOfFixedSizeString(my_str);
//![length]
    std::cout << "//![length]" << std::endl;
//![length]
    return 0;
}
//![length]
