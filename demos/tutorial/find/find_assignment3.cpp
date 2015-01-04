//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;
//![includes]

//![initialization]
int main()
{
    typedef Index<CharString, IndexQGram<UngappedShape<4>, OpenAddressing> > TIndex;
    TIndex index("tobeornottobe");
    Finder<TIndex> finder(index);
//![initialization]

//![output]
    while (find(finder, "tobe"))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

    return 0;
}
//![output]
