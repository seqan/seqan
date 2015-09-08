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

    return 0;
}
