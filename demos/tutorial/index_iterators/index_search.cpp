//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;
//![includes]

//![initialization]
int main()
{
    typedef Index<CharString> TIndex;
    TIndex index("How much wood would a woodchuck chuck?");
//![initialization]
//![initialization]

//![iterator]
    Iterator<TIndex, TopDown<> >::Type it(index);
//![iterator]

//![iteration]
    CharString pattern = "wood";
    while (repLength(it) < length(pattern))
    {
        // go down edge starting with the next pattern character
        if (!goDown(it, pattern[repLength(it)]))
            return 0;

        unsigned endPos = std::min((unsigned)repLength(it), (unsigned)length(pattern));
        // compare remaining edge characters with pattern
        std::cout << representative(it) << std::endl;
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos))
            return 0;
    }
//![iteration]

//![output]
    // if we get here the pattern was found
    // output match positions
    for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
        std::cout << getOccurrences(it)[i] << std::endl;

    return 0;
}
//![output]
