//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;
//![includes]

//![initialization]
int main()
{
    typedef Index<CharString> TIndex;
    TIndex index("How many wood would a woodchuck chuck.");
    Iterator<TIndex, TopDown<> >::Type it(index);
//![initialization]

//![output]
    if (goDown(it, "wood"))
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            std::cout << getOccurrences(it)[i] << std::endl;

    return 0;
}
//![output]
