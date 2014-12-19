//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

//![includes]
//![initialization]
int main()
{
    String<char> myString = "abracadabra";

    typedef Index<String<char>, IndexWotd<> > TMyIndex;
    TMyIndex myIndex(myString);
    String<int> propMap;

//![initialization]
//![iteration]
    Iterator<TMyIndex, TopDown<ParentLinks<Preorder> > >::Type myIterator(myIndex);

    int depth;
    while (!atEnd(myIterator))
    {
        if (isRoot(myIterator))
            depth = 0;
        else
            depth = getProperty(propMap, nodeUp(myIterator)) + 1;

        resizeVertexMap(propMap, myIndex);
        assignProperty(propMap, value(myIterator), depth);

        ++myIterator;
    }

//![iteration]
//![output]
    goBegin(myIterator);
    while (!atEnd(myIterator))
    {
        std::cout << getProperty(propMap, value(myIterator)) << '\t' << representative(myIterator) << std::endl;
        ++myIterator;
    }
    return 0;
}
//![output]
