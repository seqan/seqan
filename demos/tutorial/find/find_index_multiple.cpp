//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;
//![includes]

//![initialization]
int main()
{
    typedef StringSet<CharString> THaystacks;
    THaystacks haystacks;
    appendValue(haystacks, "tobeornottobe");
    appendValue(haystacks, "thebeeonthecomb");
    appendValue(haystacks, "beingjohnmalkovich");

    Index<THaystacks> index(haystacks);
    Finder<Index<THaystacks> > finder(haystacks);
//![initialization]

//![output]
    clear(finder);
    while (find(finder, "be"))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

    return 0;
}
//![output]
