//![initialization]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    StringSet<CharString> myStringSet;
    appendValue(myStringSet, "tobeornottobe");
    appendValue(myStringSet, "thebeeonthecomb");
    appendValue(myStringSet, "beingjohnmalkovich");

    typedef Index<StringSet<CharString> > TMyIndex;
    TMyIndex myIndex(myStringSet);

//![initialization]
//![iteration1]
    Iterator<TMyIndex, TopDown<ParentLinks<Postorder> > >::Type myIterator(myIndex);

    // Top-down iterators start in the root node which is not the first node of a
    // postorder DFS. Thus we have to manually go to the DFS start with goBegin
    goBegin(myIterator);
    while (!atEnd(myIterator))
    {
        std::cout << representative(myIterator) << std::endl;
        ++myIterator;
    }

//![iteration1]
//![iteration2]
    Iterator<TMyIndex, BottomUp<> >::Type myIterator2(myIndex);

    while (!atEnd(myIterator2))
    {
        std::cout << representative(myIterator2) << std::endl;
        ++myIterator2;
    }

    return 0;
}
//![iteration2]
