//![initialization]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    StringSet<CharString> myStringSet;
    appendValue(myStringSet, "CDFGHC");
    appendValue(myStringSet, "CDEFGAHC");

    typedef Index<StringSet<CharString> > TMyIndex;
    TMyIndex myIndex(myStringSet);

//![initialization]
//![iteration]
    Iterator<TMyIndex, Mums>::Type myIterator(myIndex);

    while (!atEnd(myIterator))
    {
        std::cout << representative(myIterator) << std::endl;
        ++myIterator;
    }

    return 0;
}
//![iteration]
