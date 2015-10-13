#include <iostream>
#include <seqan/stream.h>  // For String output to std::cout.
#include <seqan/sequence.h>

using namespace seqan;

void replaceAs(CharString & str)
{
    typedef Iterator<CharString, Standard>::Type TIterator;

    for (TIterator it = begin(str, Standard()); it != end(str, Standard()); ++it)
    {
        if (value(it) == 'a')
            value(it) = 'X';
    }
}

int main()
{
    CharString str1 = "abcdefghijklmnopqrstuvwxyz";
    replaceAs(str1);
    std::cout << str1 << std::endl;

    CharString str2 = "Hello SeqAn!";
    replaceAs(str2);
    std::cout << str2 << std::endl;

    CharString str3 = "Hello Seqan!";
    replaceAs(str3);
    std::cout << str3 << std::endl;
    return 0;
}
