///A tutorial about finding supermaximal repeats.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
///We begin with a @Class.String@ to store our sequence.
    String<char> myString = "How many wood would a woodchuck chuck.";

///Then we create an @Class.Index@ of this @Class.StringSet@.
    typedef Index<String<char> > TMyIndex;
    TMyIndex myIndex(myString);

///To find supermaximal repeats, we use SeqAn's @Spec.SuperMaxRepeats Iterator@
///and set the minimum repeat length to 3.
    Iterator<TMyIndex, SuperMaxRepeats>::Type myRepeatIterator(myIndex, 3);

    while (!atEnd(myRepeatIterator))
    {
///A repeat can be represented by its length and positions it occurs at.
///@Function.getOccurrences@ returns an unordered sequence of these positions
///The length of this sequence, i.e. the repeat abundance can be obtained
///from @Function.countOccurrences@.
        for (unsigned i = 0; i < countOccurrences(myRepeatIterator); ++i)
            std::cout << getOccurrences(myRepeatIterator)[i] << ", ";

///@Function.repLength@ returns the length of the repeat string.
        std::cout << repLength(myRepeatIterator) << "   ";

///The repeat string itself can be determined with @Function.representative@.
        std::cout << "\t\"" << representative(myRepeatIterator) << '\"' << std::endl;

        ++myRepeatIterator;
    }

    return 0;
}
