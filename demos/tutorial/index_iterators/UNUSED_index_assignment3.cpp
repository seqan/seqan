//![initialization]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef Index<CharString> TIndex;
    TIndex index("tobeornottobe");
    Iterator<TIndex, TopDown<ParentLinks<> > >::Type it(index);
//![initialization]
//![iteration]

    do
    {
        std::cout << representative(it) << std::endl;
        if (!goDown(it) && !goRight(it))
            while (goUp(it) && !goRight(it))
                ;
    }
    while (!isRoot(it));

    return 0;
}
//![iteration]
