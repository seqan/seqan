#include <seqan/index.h>

using namespace seqan;

int main ()
{
    typedef Index<CharString> TIndex;

    TIndex index("tobeornottobe");
    Iterator<TIndex, TopDown<ParentLinks<> > >::Type it(index);

    // ----------------------------------------------------------
    std::cout << "DFS pre-order traversal" << std::endl;
    do {
        // Print theletters from the root to the current node
        std::cout << representative(it) << std::endl;

        if (!goDown(it) && !goRight(it))
            while (goUp(it) && !goRight(it)) ;
    } while (!isRoot(it));

    // ----------------------------------------------------------
    std::cout << "DFS post-order traversal" << std::endl;
    goBegin(it);

    do {
        if (!goDown(it))
        {
            std::cout << representative(it) << std::endl;
            if(!goRight(it))
                while (goUp(it))
                {
                    std::cout << representative(it) << std::endl;
                    if (goRight(it))
                        break;
                }
        }
    } while (!isRoot(it));

    return 0;
}
