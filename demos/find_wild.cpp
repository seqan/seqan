///A tutorial about the use of wildcard find algorithms.
#include <iostream>
#include <seqan/find.h>

using namespace seqan;

///This program uses the algorithm @Spec.WildShiftAnd@ to perform a wildcard search.
int main()
{
    String<char> hayst = "If you must cross a course cross cow across a crowded cow crossing, "
                         "cross the cross coarse cow across the crowded cow crossing carefully.";
    String<char> ndl = "cr?o[uw]";
///The pattern matches e.g. "cow", "crow", and "cou", but not "cros".
    Finder<String<char> > finder(hayst);
    Pattern<String<char>, WildShiftAnd> pattern(ndl);

    while (find(finder, pattern))
    {
        std::cout << position(finder) << "\n";
    }
    return 0;
}
