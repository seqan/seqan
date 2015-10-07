//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

//![includes]
//![initialization]
int main()
{
    String<char> myString = "abracadabra";

    typedef Index<String<char> > TMyIndex;
    TMyIndex myIndex(myString);

//![iteration]
    Iterator<TMyIndex, TopDown<ParentLinks<Preorder> > >::Type myIterator(myIndex);

    while (!atEnd(myIterator))
    {
        std::cout << representative(myIterator) << std::endl;
        ++myIterator;
    }

    return 0;
}
//![iteration]
