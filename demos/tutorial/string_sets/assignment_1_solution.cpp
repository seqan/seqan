#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    // Build strings
    DnaString str0 = "AAA";
    DnaString str1 = "CCC";
    DnaString str2 = "GGG";
    DnaString str3 = "TTT";
    // Build string set and append strings
    StringSet<DnaString> stringSet;
    appendValue(stringSet, str0);
    appendValue(stringSet, str1);
    appendValue(stringSet, str2);
    appendValue(stringSet, str3);
    // Print the length of the string set
    std::cout << length(stringSet) << std::endl;
    // Print all elements
    for (unsigned i = 0; i < length(stringSet); ++i)
    {
        std::cout << stringSet[i] << std::endl;
    }
    return 0;
}
