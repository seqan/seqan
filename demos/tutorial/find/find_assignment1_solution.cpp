#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString haystack = "Simon, send more money!";
    String<CharString> needles;
    appendValue(needles, "mo");
    appendValue(needles, "send");
    appendValue(needles, "more");

    Finder<CharString> finder(haystack);
    Pattern<String<CharString>, WuManber> pattern(needles);
    while (find(finder, pattern))
    {
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t";
        std::cout << position(pattern) << '\t' << infix(finder) << std::endl;
    }
    return 0;
}
