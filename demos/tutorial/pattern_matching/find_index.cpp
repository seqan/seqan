//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;
//![includes]

//![initialization]
int main()
{
    Index<CharString> index("tobeornottobe");
    CharString needle = "be";
    Finder<Index<CharString> > finder(index);
//![initialization]
//![output]
    Pattern<CharString> pattern(needle);
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
//![output]

//![output_short]
    clear(finder);
    while (find(finder, "be"))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

    return 0;
}
//![output_short]
