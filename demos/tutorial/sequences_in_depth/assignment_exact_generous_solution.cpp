#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    unsigned num = 100000;
    double start;

    String<Dna> str;
    clear(str);
    start = sysTime();
    for (unsigned i = 0; i < num; ++i)
        appendValue(str, 'A', Exact());
    std::cout << "Strategy Exact() took: " << sysTime() - start << " s\n\n";

    clear(str);
    shrinkToFit(str);
    start = sysTime();
    for (unsigned i = 0; i < num; ++i)
        appendValue(str, 'A', Generous());
    std::cout << "Strategy Generous() took: " << sysTime() - start << " s\n\n";

    return 0;
}
