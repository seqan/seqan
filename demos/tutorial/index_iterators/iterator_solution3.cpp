//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

//![includes]
//![iteration]
template <typename TIndexSpec>
void constrainedDFS()
{
    typedef Index<CharString, TIndexSpec> TIndex;
    TIndex index("mississippi");
    typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type it(index);

    do
    {
        std::cout << representative(it) << std::endl;
        if (!goDown(it) || repLength(it) > 3)
            do
            {
                if (!goRight(it))
                    while (goUp(it) && !goRight(it))
                        ;
            }
            while (repLength(it) > 3);
    }
    while (!isRoot(it));
    std::cout << std::endl;
}

int main()
{
    constrainedDFS<IndexEsa<> >();
    constrainedDFS<IndexWotd<> >();
    return 0;
}
//![iteration]
