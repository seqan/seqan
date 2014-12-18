//![includes]
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;
//![includes]

//![print-strings-rec]
// Helper function for printPermutations().
void printStringsRec(String<char> & current, unsigned int pos)
{
    if (pos < length(current))
    {
        for (char c = 'a'; c <= 'z'; ++c)
        {
            current[pos] = c;
            printStringsRec(current, pos + 1);
        }
    }
    else
    {
        std::cout << current << std::endl;
    }
}
//![print-strings-rec]

//![print-strings]
// Print all strings of the alphabet {a, ..., z} of length len.
void printStrings(int len)
{
    String<char> current;
    // Build first value;
    for (int i = 0; i < len; ++i)
        current += "a";
    // Now, current == "a" * len.

    // To iterate: human -- to recurse: DIVINE.
    printStringsRec(current, 0);
}
//![print-strings]

//![main]
int main()
{
    printStrings(3);
    return 0;
}
//![main]
