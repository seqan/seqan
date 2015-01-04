// The comment lines containing ![fragment-line] are there for the
// documentation system.  You can ignore them when reading this file.
//![includes]
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
//![includes]
//![metafunctions]
    String<char> str = "admn";
    Iterator<String<char> >::Type it = begin(str);
    Iterator<String<char> >::Type itEnd = end(str);
//![metafunctions]
//![iterators]
    while (it != itEnd)
    {
        std::cout << *it;
        ++it;
    }
    std::cout << std::endl;
//![iterators]
//![rooted-iterators]
    Iterator<String<char>, Rooted>::Type it2 = begin(str);
    for (goBegin(it2); !atEnd(it2); goNext(it2))
    {
        ++value(it2);
    }
//![rooted-iterators]
//![iterator-reverse]
    goEnd(it2);
    while (!atBegin(it2))
    {
        goPrevious(it2);
        std::cout << getValue(it2);
    }
    std::cout << std::endl;
//![iterator-reverse]
//![assign-value]
    assignValue(begin(str), 'X');
    std::cout << str << std::endl;

    return 0;
}
//![assign-value]
