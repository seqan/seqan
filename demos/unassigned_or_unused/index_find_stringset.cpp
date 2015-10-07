///A tutorial about an index based finder for a StringSet.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
///First, we create a @Class.StringSet@ of 3 @Class.String.Strings@.
    StringSet<String<char> > mySet;
    resize(mySet, 3);
    mySet[0] = "tobeornottobe";
    mySet[1] = "thebeeonthecomb";
    mySet[2] = "beingjohnmalkovich";

///Then we create an @Class.Index@ of our @Class.StringSet@ and
///a @Class.Finder@ of the @Class.Index@.
    Index<StringSet<String<char> > > myIndex(mySet);
    Finder<Index<StringSet<String<char> > > > myFinder(myIndex);

///Finally we search for the string $"be"$ and output all occurrences
    std::cout << "hit at ";
    while (find(myFinder, "be"))
        std::cout << position(myFinder) << "  ";
    std::cout << std::endl;

    return 0;
}
